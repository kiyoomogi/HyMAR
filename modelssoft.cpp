#include "modelssoft.h"
#include "base/src/version.h"

namespace models {

    ModelStrainSoftening::ModelStrainSoftening() : ModelMohr() {
    }

    string ModelStrainSoftening::getProperties(void) const {
        return ModelMohr::getProperties() +
            ",table-cohesion,table-friction,table-dilation,table-tension"
            ",strain-shear-plastic,strain-tension-plastic";
    }

    uint32 ModelStrainSoftening::getMinorVersion() const {
        return version.revision_;
    }

    static std::atomic<uint64> nSS{0};
    base::Property ModelStrainSoftening::getProperty(uint32 index) const {
        if (nSS.load() == 0)
            nSS = (split(ModelMohr::getProperties(), ",")).size();
        if (index <= nSS)
            return ModelMohr::getProperty(index);
        else if (index - nSS <= 6) {
            switch (index - nSS) {
            case 1: return cTable_;
            case 2: return fTable_;
            case 3: return dTable_;
            case 4: return tTable_;
            case 5: return sHP_;
            case 6: return tHP_;
            }
        }

        return 0.0;
    }

    void ModelStrainSoftening::setProperty(uint32 index,const base::Property &p,uint32 restoreVersion) {
        ConstitutiveModel::setProperty(index, p, restoreVersion);
        if (nSS.load()==0)
            nSS = (split(ModelMohr::getProperties(), ",")).size();
        if (index <= nSS)
            ModelMohr::setProperty(index, p, restoreVersion);
        else if (index - nSS <= 6) {
            switch (index - nSS) {
            case 1: cTable_ = p.to<string>();  break;
            case 2: fTable_ = p.to<string>();  break;
            case 3: dTable_ = p.to<string>();  break;
            case 4: tTable_ = p.to<string>();  break;
            case 5: sHP_    = p.to<double>();  break;
            case 6: tHP_    = p.to<double>();  break;
            }
        }
    }

    bool ModelStrainSoftening::isPropertyReadOnly(uint32 index) const {
        if (nSS.load() == 0)
            nSS = (split(ModelMohr::getProperties(), ",")).size();
        if (index <= nSS)
            return ModelMohr::isPropertyReadOnly(index);
        else {
            auto i = index - nSS;
            if (i == 5 || i == 6)
                return true;
        }
        return false;
    }

    bool ModelStrainSoftening::isPropertyAdvanced(uint32 index) const {
        if (nSS.load() == 0)
            nSS = (split(ModelMohr::getProperties(), ",")).size();
        if (index <= nSS)
            return ModelMohr::isPropertyAdvanced(index);

        return false;
    }

    void ModelStrainSoftening::copy(const ConstitutiveModel *m) {
        const ModelStrainSoftening *mm = dynamic_cast<const ModelStrainSoftening *>(m);
        if (!mm) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        ModelMohr::copy(m);
        cTable_ = mm->cTable_;
        fTable_ = mm->fTable_;
        dTable_ = mm->dTable_;
        tTable_ = mm->tTable_;
        sHP_    = mm->sHP_;
        tHP_    = mm->tHP_;
    }

    void ModelStrainSoftening::initialize(uint32 d, State *s) {
        ConstitutiveModel::initialize(d,s);

        iCohesion_ = iFriction_ = iDilation_ = iTension_ = nullptr;

        if (cTable_.length()) iCohesion_ = s->getTableIndexFromID(cTable_);
        if (fTable_.length()) iFriction_ = s->getTableIndexFromID(fTable_);
        if (dTable_.length()) iDilation_ = s->getTableIndexFromID(dTable_);
        if (tTable_.length()) iTension_  = s->getTableIndexFromID(tTable_);
        if(iTension_ && brittle_)
            throw std::runtime_error("Internal error: flag-brittle not allowed on when assigning table-tension.");

        if (iTension_)  tension_  = s->getYFromX(iTension_, tHP_);
        if (iCohesion_) cohesion_ = s->getYFromX(iCohesion_,sHP_);
        if (iFriction_) friction_ = s->getYFromX(iFriction_,sHP_);
        if (iDilation_) dilation_ = s->getYFromX(iDilation_,sHP_);

        //if (friction_) {
        //    double dApex = cohesion_ / tan(friction_ * degrad);
        //    tension_ = tension_ < dApex ? tension_ : dApex;
        //}
        //else if (cohesion_ == 0.0)
        //    tension_ = 0.0;
        ModelMohr::updateParameters(s);
    }

    static const uint32 Dqs = 0;
    static const uint32 Dqt = 1;
    //
    void ModelStrainSoftening::run(uint32 d, State *s) {
        FP_S;
        ConstitutiveModel::run(d,s);
        FP_S;

        //bool noUpdateE = false;
        if (s->modulus_reduction_factor_ > 0.0) {
            ModelElastic::moduliReduction(s->modulus_reduction_factor_);
            //noUpdateE = true;
        }
        //ModelMohr::updateParameters(noUpdateE,&sf1_,&sf3_,s);

        if (s->state_ & shear_now) s->state_ |= shear_past;
        s->state_ &= ~shear_now;
        if (s->state_ & tension_now) s->state_ |= tension_past;
        s->state_ &= ~tension_now;
        uint32 iPlas = 0;

        FP_S;
        if (!s->sub_zone_) {
            s->working_[Dqs] = 0.0;
            s->working_[Dqt] = 0.0;
            if (s->plasticStrainActive_)
                s->stnPSI_ = SymTensor(0.0);
        }
        double dSubZoneVolume = s->getSubZoneVolume();
        FP_S;

        ModelElastic::elasticTrial(s);
        FP_S;
        s->viscous_ = true;
        if (!canFail()) return;

        DVect3 pps(0.0); // incremental principal plastic strain vector

        FP_S;
        if (supportsUniaxial() && s && s->condition_==State::CondType::D1) {
            double lams(0.0), lamt(0.0);
            double s11 = s->stnS_.s11();         
            if (s11 > 0.0) {
                double fs = nph_*s11  - csn_;
                double fsd = fs / nph_;
                double ft = s11 - tension_;
                if (fsd > 0.0 && fsd > ft) {
                    // yielding in shear only
                    iPlas = 1;
                    s->state_ |= shear_now;
                    s->stnS_.rs11() = csn_/nph_; // s11 - fs/nph_;
                    lams = fs / (e1_*nph_*nps_);
                }
                else if (ft > 0.0 && ft >= fsd) {
                    // yielding in tension only
                    iPlas = 2;
                    s->state_ |= tension_now;
                    if (brittle_) tension_ = 0.0;
                    s->stnS_.rs11() = tension_;
                    lamt = ft / e1_;
                }
            }
            else {
                double fs = -s11 - csn_;
                if (fs > 0.0) {
                    // yielding in shear only
                    iPlas = 1;
                    s->state_ |= shear_now;
                    s->stnS_.rs11() = -csn_; // s11 + fs
                    lams = fs / e1_;
                }
            }

            s->working_[Dqs] += (lams/std::sqrt(3.0)) * dSubZoneVolume;
            s->working_[Dqt] += lamt * dSubZoneVolume;

            if (iPlas) s->viscous_ = false;
        } else if (supportsPlaneStress() && s && s->condition_==State::CondType::D2) {
            double sp = (s->stnS_.s11() + s->stnS_.s22())*0.5;
            double sn = (s->stnS_.s11() - s->stnS_.s22())*0.5;
            double st = s->stnS_.s12();
            double sz = s->stnS_.s33();
            double sw = std::sqrt(sn*sn + st*st);
            double s1 = sp - sw;
            double s2 = sp + sw;
            double theta2 = (isD0(st) && isD0(sn)) ? 0.0 : std::atan2(st, sn); // in rad
            double s1n=0.0, s2n=0.0;
            double lamt1(0.0), lamt2(0.0), lams1(0.0), lams2(0.0);
            if (s2 >= sz && s1 <= sz) {
                // case A: s1 <= sz <= s2
                double ft1 = s1 - tension_;
                double ft2 = s2 - tension_;
                double fsa = -s1 + nph_*s2 - csn_; // in-plane
                if (ft1 > 0.0) {
                    // yielding in tension in two principal directions
                    iPlas = 4;
                    s->state_ |= tension_now;
                    if (brittle_) tension_ = 0.0;
                    s1n = tension_;
                    s2n = tension_;
                    double delta = e1_*e1_ - e2_*e2_;
                    assert(delta > 0.0);
                    lamt1 = (e1_*ft1 - e2_*ft2) / delta;
                    lamt2 = (e1_*ft2 - e2_*ft1) / delta;
                } else {
                    if (fsa > 0.0 && ft2 > 0.0) {
                        // yielding in tension && shear
                        iPlas = 3;
                        s->state_ |= shear_now;
                        s->state_ |= tension_now;
                        if (brittle_) tension_ = 0.0;
                        s1n = nph_*tension_ - csn_;
                        s2n = tension_;
                        double delta = e1_*e1_ - e2_*e2_;
                        double fs1 = s1 - nph_*s2n + csn_;
                        assert(delta > 0.0);
                        double eps1 = (e1_*fs1 - e2_*ft2) / delta;
                        double eps2 = (e1_*ft2 - e2_*fs1) / delta;
                        lams1 = -eps1;
                        lamt2 = eps2 - lams1*nps_;
                    } else if (fsa <= 0.0 && ft2 > 0.0) {
                        // yielding in tension only
                        iPlas = 2;
                        s->state_ |= tension_now;
                        if (brittle_) tension_ = 0.0;
                        s1n = s1 - (e2_/e1_)*ft2;
                        s2n = tension_;
                        lamt2 = ft2/e1_;
                    } else if (fsa > 0.0 && ft2 <= 0.0) {
                        // yielding in shear only
                        iPlas = 1;
                        s->state_ |= shear_now;
                        s1n = s1 + fsa*sc1_;
                        s2n = s2 + fsa*sc3_;
                        lams1 = fsa*sf1_;
                        lams2 = fsa*sf3_;
                    }
                }
            } else if (s2 < sz) {
                // Case B: s1 <= s2 < sz
                double fsb = -s1 + nph_*sz - csn_;
                if (fsb > 0.0) {
                    // yielding in shear only
                    iPlas = 1;
                    s->state_ |= shear_now;
                    s1n = s1 + fsb;  // or: nph_* sz - csn_;
                    s2n = s2 + (e2_/e1_)*fsb;
                    lams1 = fsb/e1_;
                }
            } else if (s1 > sz) {
                // Case C: sz < s1 <= s2
                double ft1 = s1 - tension_;
                double ft2 = s2 - tension_;
                double fsc = -sz + nph_*s2 - csn_;
                if (ft1 > 0.0) {
                    // yielding in tension in two principal directions
                    iPlas = 4;
                    s->state_ |= tension_now;
                    if (brittle_) tension_ = 0.0;
                    s1n = tension_;
                    s2n = tension_;
                    double delta = e1_*e1_ - e2_*e2_;
                    assert(delta > 0.0);
                    lamt1 = (e1_*ft1 - e2_*ft2) / delta;
                    lamt2 = (e1_*ft2 - e2_*ft1) / delta;
                } else {
                    if (ft2 > 0.0) {
                        // yielding in tension
                        iPlas = 2;
                        s->state_ |= tension_now;
                        if (brittle_) tension_ = 0.0;
                        s1n = s1 - (e2_/e1_)*ft2;
                        s2n = tension_;
                        lamt2 = ft2/e1_;
                    }
                    else if (fsc > 0.0) {
                        // yielding in shear
                        iPlas = 1;
                        s->state_ |= shear_now;
                        s1n = s1 - (e2_/e1_) * fsc/nph_;
                        s2n = s2 - fsc/nph_;
                        lams2 = fsc/(e1_*nph_);
                    }
                }
            }

            if (iPlas) {
                double pp = (s1n + s2n)*0.5;
                double pn = (s2n - s1n)*0.5;
                double cos2t = std::cos(theta2);
                double sin2t = std::sin(theta2);
                s->stnS_.rs11() = pp + pn*cos2t;
                s->stnS_.rs22() = pp - pn*cos2t;
                s->stnS_.rs12() = pn*sin2t;
                s->stnS_.rs33() = sz;
                s->stnS_.rs13() = 0.0;
                s->stnS_.rs23() = 0.0;
                s->viscous_ = false;

                if (iPlas == 1 || iPlas == 3) {
                    double dDepa = d1d3 * (lams1 + lams2);
                    lams1 -= dDepa;
                    lams2 -= dDepa;
                    s->working_[Dqs] += sqrt(0.5 * (lams1*lams1 + dDepa*dDepa + lams2*lams2)) * dSubZoneVolume;
                }

                if (iPlas == 2 || iPlas == 4)
                    s->working_[Dqt] += (lamt1 + lamt2) * dSubZoneVolume;
            }
        } else {
            FP_S;
            SymTensorInfo info;
            DVect3 prin = s->stnS_.getEigenInfo(&info);
            FP_S;

            double fs = -prin.x() + nph_*prin.z() - csn_;
            double fsd = fs / rc_;
            FP_S;
            double ftz = prin.z() - tension_;
            double fty = prin.y() - tension_;
            double ftx = prin.x() - tension_;
            FP_S;
            if (fsd > 0.0 && fsd >= ftz)
                ModelMohr::shearCorrection(s, &prin, &iPlas, fs, &pps);
            else if (ftz > 0.0 && ftz >= fsd)
                ModelMohr::tensionCorrection(s, &prin, &iPlas, ftz, brittle_, &pps);
            FP_S;
            ModelMohr::apexCorrection(friction_, s, &prin, &iPlas, brittle_);
            FP_S;

            if (iPlas) {
                FP_S;
                s->viscous_ = false;
                s->stnS_ = info.resolve(prin);
                //
                if (s->plasticStrainActive_) {
                    SymTensor pst = info.resolve(pps);
                    s->stnPSI_ += pst * s->getSubZoneVolume();
                }
                //
                if (iPlas == 1) {
                    double dDe1p = fs*sf1_;
                    double dDe3p = fs*sf3_;
                    double dDepa = d1d3 * (dDe1p + dDe3p);
                    dDe1p -= dDepa;
                    dDe3p -= dDepa;
                    s->working_[Dqs] += sqrt(0.5*(dDe1p*dDe1p + dDepa*dDepa + dDe3p*dDe3p)) * dSubZoneVolume;
                }
                else if (iPlas == 2) {
                    double dAux = ftz/e1_;
                    s->working_[Dqt] += dAux * dSubZoneVolume;
                }
                else if (iPlas == 3) {
                    double dAux = (ftz + fty) / (e1_ + e2_);
                    s->working_[Dqt] += dAux * dSubZoneVolume;
                }
                else if (iPlas == 4) {
                    double dAux = (ftz + fty + ftx) / (e1_ + 2.0*e2_);
                    s->working_[Dqt] += dAux * dSubZoneVolume;
                }
                FP_S;
            }
            FP_S;
        }

        if (s->sub_zone_==s->total_sub_zones_-1) {
            FP_S;
            double dAux = 1.0 / s->getZoneVolume();
            FP_S;
            if (s->overlay_==2) dAux *= 0.5;

            // accumulate plastic strain
            if (s->plasticStrainActive_)
                plasticStrain_ += s->stnPSI_ * dAux;

            sHP_ += s->working_[Dqs] * dAux;
            tHP_ += s->working_[Dqt] * dAux;

            if (s->working_[Dqs] > 0.0) {
                if (iCohesion_) cohesion_ = s->getYFromX(iCohesion_,sHP_);
                if (iFriction_) friction_ = s->getYFromX(iFriction_,sHP_);
                if (iDilation_) dilation_ = s->getYFromX(iDilation_,sHP_);
            }
            if (s->working_[Dqt] > 0.0 && iTension_)
                tension_ = s->getYFromX(iTension_,tHP_);

            FP_S;
            ModelMohr::updateParameters(s, true);
            FP_S;
        }
    }
}

// EOF
