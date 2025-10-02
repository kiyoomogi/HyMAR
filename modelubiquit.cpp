#include "modelubiquit.h"
#include "base/src/version.h"

namespace models {

    ModelUbiquitous::ModelUbiquitous(unsigned short option) : ModelMohr(option), ModelWeakplaneType() {
        FP_S;
    }

    uint32 ModelUbiquitous::getMinorVersion() const {
        return version.revision_;
    }

    string ModelUbiquitous::getProperties(void) const {
        return ModelMohr::getProperties() 
            + ",joint-cohesion,joint-friction,joint-dilation,joint-tension,"
            + ModelWeakplaneType::getProperties();       
    }

    string ModelUbiquitous::getStates(void) const {
        return ModelMohr::getStates() + ",j-shear-n,j-tension-n,j-shear-p,j-tension-p";
    }

    static std::atomic<uint64> nMC{0};
    static std::atomic<uint64> nWP{0};
    base::Property ModelUbiquitous::getProperty(uint32 index) const {
        if (not nMC)
            nMC = (split(ModelMohr::getProperties(), ",")).size();
        if (not nWP)
            nWP = (split(ModelWeakplaneType::getProperties(), ",")).size();
        if (index <= nMC)
            return ModelMohr::getProperty(index);
        else if (index <= nMC + 4) {
            switch (index - nMC) {
            case  1: return (jCohesion_);
            case  2: return (jFriction_);
            case  3: return (jDilation_);
            case  4: return (jTension_);
            }
        } else if (index - nMC - 4 <= nWP) 
            return ModelWeakplaneType::getProperty(index - nMC - 4);

        return 0.0;
    }

    void ModelUbiquitous::setProperty(uint32 index,const base::Property &p,uint32 restoreVersion) {
        ConstitutiveModel::setProperty(index, p, restoreVersion);
        if (not nMC)
            nMC = (split(ModelMohr::getProperties(), ",")).size();
        if (not nWP)
            nWP = (split(ModelWeakplaneType::getProperties(), ",")).size();
        if (index <= nMC)
            ModelMohr::setProperty(index,p,restoreVersion);
        else if (index <= nMC + 4) {
            switch (index - nMC) {
            case 1: jCohesion_ = p.to<double>(); break;
            case 2: jFriction_ = p.to<double>(); break;
            case 3: jDilation_ = p.to<double>(); break;
            case 4: jTension_ = p.to<double>();  break;
            }
        } else if (index - nMC - 4 <= nWP)
            ModelWeakplaneType::setProperty(index - nMC - 4, p, restoreVersion);
    }

    void ModelUbiquitous::copy(const ConstitutiveModel *m) {
        const ModelUbiquitous *mm = dynamic_cast<const ModelUbiquitous *>(m);
        if (!mm) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        ModelMohr::copy(m);
        jCohesion_= mm->jCohesion_;
        jFriction_= mm->jFriction_;
        jDilation_= mm->jDilation_;
        jTension_ = mm->jTension_;
        dip_      = mm->dip_;
        dd_       = mm->dd_;
        norma_    = mm->norma_;
        angle_    = mm->angle_;
    }


    void ModelUbiquitous::initialize(uint32 d,State *s) {
        ConstitutiveModel::initialize(d,s);
        if ((d!=2)&&(d!=3)) throw std::runtime_error("Illegal dimension in the ubiquitous-joint model");

        ModelMohr::updateParameters(s);
        updateJointParameters(true);
    }

    void ModelUbiquitous::run(uint32 d, State *s) {
        ConstitutiveModel::run(d,s);

        if (s->modulus_reduction_factor_ > 0.0)
            ModelMohr::moduliReduction(s->modulus_reduction_factor_);

        if (s->state_ & shear_now) s->state_ |= shear_past;
        s->state_ &= ~shear_now;
        if (s->state_ & tension_now) s->state_ |= tension_past;
        s->state_ &= ~tension_now;
        uint32 iPlas = 0;

        if (s->state_ & joint_shear_now) s->state_ |= joint_shear_past;
        s->state_ &= ~joint_shear_now;
        if (s->state_ & joint_tension_now) s->state_ |= joint_tension_past;
        s->state_ &= ~joint_tension_now;
        uint32 jPlas = 0;

        if (!s->sub_zone_ and s->plasticStrainActive_)
            s->stnPSI_ = SymTensor(0.0);

        ModelElastic::elasticTrial(s);
        s->viscous_ = true;        
        if (!canFail()) return;

        DVect3 pps(0.0); // incremental principal plastic strain vector

        SymTensorInfo info;
        DVect3 prin = s->stnS_.getEigenInfo(&info);

        double fs  = - prin.x() + nph_ * prin.z() - csn_;
        double fsd = fs / rc_;
        double ftz = prin.z() - tension_;
        if (fsd > 0.0 && fsd >= ftz)
            ModelMohr::shearCorrection(s, &prin, &iPlas, fs, &pps);
        else if (ftz > 0.0 && ftz >= fsd)
            ModelMohr::tensionCorrection(s, &prin, &iPlas, ftz, brittle_, &pps);
        ModelMohr::apexCorrection(friction_, s, &prin, &iPlas, brittle_);

        if (iPlas) {
            s->stnS_ = info.resolve(prin);
            if (s->plasticStrainActive_) {
                SymTensor pst = info.resolve(pps);
                s->stnPSI_ += pst * s->getSubZoneVolume();
            }
        }

        stressToLocal(s, norma_.axes());

        dsp11_ = 0.0; dsp33_=0.0; dsp13_=0.0; dsp23_=0.0;
        DVect3 ppsj(0.0);

        double fjt = sig_ - jTension_;
        double fjs = tau_ + sig_*qphi_ - jCohesion_;
        double fjd = fjs / jrc_;
        if (fjs > 0.0 && fjd >= fjt)
            shearCorrectionJoint(s, &jPlas, fjs, &ppsj);
        else if (fjt > 0.0 && fjt >= fjd)
            tensionCorrectionJoint(s, &jPlas, fjt, &ppsj);
        apexCorrectionJoint(s,&jPlas);

        if (jPlas > 0) {
            stressToGlobal(s, norma_.axes(), &ppsj);
            if (s->plasticStrainActive_) {
                SymTensor pst = info.resolve(ppsj);
                s->stnPSI_ += pst * s->getSubZoneVolume();
            }
        }

        if (s->isLarge()) 
            largeStrainCorrection(s, norma_.axes());

        if (s->sub_zone_ == s->total_sub_zones_ - 1) {
            if (s->plasticStrainActive_) {
                double dVal = 1.0 / s->getZoneVolume();
                if (s->overlay_ == 2) dVal *= 0.5;
                plasticStrain_ += s->stnPSI_ * dVal;
            }
        }

        if (iPlas || jPlas)
            s->viscous_ = false;
    }

    void ModelUbiquitous::scaleProperties(const double &scale,const std::vector<uint32> &props) {
        for (uint32 u=0;u<props.size();++u) {
            switch (props[u]) {
            case  5: cohesion_ *= scale;  break;
            case  6: friction_ = std::max(0.0,std::min(89.0,std::atan(std::tan(friction_*degrad)*scale)/degrad));  break;
            case  7: dilation_ = std::max(0.0,std::min(89.0,std::atan(std::tan(dilation_*degrad)*scale)/degrad));  break;
            case  8: tension_  *= scale;  break;
            case 10: jCohesion_ *= scale;  break;
            case 11: jFriction_ = std::max(0.0,std::min(89.0,std::atan(std::tan(jFriction_*degrad)*scale)/degrad));  break;
            case 12: jDilation_ = std::max(0.0,std::min(89.0,std::atan(std::tan(jDilation_*degrad)*scale)/degrad));  break;
            case 13: jTension_  *= scale;  break;
            }
            setValid(0);
        }
    }

    bool ModelUbiquitous::updateJointParameters(bool bEUpdated) {
        if (!bEUpdated) ModelElastic::updateParameters();

        if (jCohesion_ < 0.0)
            throw std::runtime_error("Ubiquitous-Joint type model: joint cohesion is not allowed less than 0.");

        if (jFriction_ >= 90.0 || jFriction_ < 0.0)
            throw std::runtime_error("Ubiquitous-Joint type model: joint friction angle is not in the valid range of 0 to 90.");

        if (jDilation_ > jFriction_)
            throw std::runtime_error("Ubiquitous-Joint type model: joint dilationn angle is not allowed greater than joint friction angle.");

        qphi_ = std::tan(jFriction_ * degrad);
        qpsi_ = std::tan(jDilation_ * degrad);
        jc1_ = 1. / (g2_ + e1_ * qphi_ * qpsi_);
        jc2_ = jc1_ * qpsi_;
        jrc_ = std::sqrt(1. + qphi_ * qphi_);

        if (jFriction_ > 0.0) {
            double dapex = jCohesion_ / qphi_;
            jTension_ = jTension_ < dapex ? jTension_ : dapex;
        }
        else if (isD0(jCohesion_))
            jTension_ = 0.0;

        double dLenTemp = norma_.norm().mag(); 
        if (isD0(dLenTemp)) {
            DVect3 norm = Orientation3::getNormFromDipDD(dip_*degrad,dd_*degrad);
            norma_.set_norm(norm);
            angle_ = Orientation3::getJointAngleFromNorm(norm) / degrad;
            dLenTemp = norma_.norm().mag();
        }
        normX_ = norma_.norm().x();
        normY_ = norma_.norm().y();
        normZ_ = norma_.norm().z();
        if (dLenTemp) {
            normX_ /= dLenTemp;
            normY_ /= dLenTemp;
            normZ_ /= dLenTemp;
        }
        return !bEUpdated;
    }

    void ModelUbiquitous::stressToLocal(State* s, const Axes3D& aAxes) {
        double e1x = aAxes.e1().x();
        double e1y = aAxes.e1().y();
        double e1z = aAxes.e1().z();
        double e2x = aAxes.e2().x();
        double e2y = aAxes.e2().y();
        double e2z = aAxes.e2().z();
        double e3x = aAxes.e3().x();
        double e3y = aAxes.e3().y();
        double e3z = aAxes.e3().z();
        double dB3x  = s->stnS_.s11()*e3x + s->stnS_.s12()*e3y + s->stnS_.s13()*e3z;
        double dB3y  = s->stnS_.s12()*e3x + s->stnS_.s22()*e3y + s->stnS_.s23()*e3z;
        double dB3z  = s->stnS_.s13()*e3x + s->stnS_.s23()*e3y + s->stnS_.s33()*e3z;
        sp13_ = e1x*dB3x + e1y*dB3y + e1z*dB3z;
        sp23_ = e2x*dB3x + e2y*dB3y + e2z*dB3z;
        sig_  = e3x*dB3x + e3y*dB3y + e3z*dB3z;
        tau_  = std::sqrt(sp13_*sp13_ + sp23_*sp23_);
    }

    void ModelUbiquitous::stressToGlobal(State* s, const Axes3D& aAxes, DVect3* psj) {
        double e1x = aAxes.e1().x();
        double e1y = aAxes.e1().y();
        double e1z = aAxes.e1().z();
        double e2x = aAxes.e2().x();
        double e2y = aAxes.e2().y();
        double e2z = aAxes.e2().z();
        double e3x = aAxes.e3().x();
        double e3y = aAxes.e3().y();
        double e3z = aAxes.e3().z();
        double dB1x = dsp11_*e1x + dsp13_*e3x;
        double dB2x = dsp11_*e1y + dsp13_*e3y;
        double dB3x = dsp11_*e1z + dsp13_*e3z;
        double dB1y = dsp11_*e2x + dsp23_*e3x;
        double dB2y = dsp11_*e2y + dsp23_*e3y;
        double dB3y = dsp11_*e2z + dsp23_*e3z;
        double dB1z = dsp13_*e1x + dsp23_*e2x + dsp33_*e3x;
        double dB2z = dsp13_*e1y + dsp23_*e2y + dsp33_*e3y;
        double dB3z = dsp13_*e1z + dsp23_*e2z + dsp33_*e3z;
        s->stnS_.rs11() += e1x*dB1x + e2x*dB1y + e3x*dB1z;
        s->stnS_.rs12() += e1x*dB2x + e2x*dB2y + e3x*dB2z;
        s->stnS_.rs13() += e1x*dB3x + e2x*dB3y + e3x*dB3z;
        s->stnS_.rs22() += e1y*dB2x + e2y*dB2y + e3y*dB2z;
        s->stnS_.rs23() += e1y*dB3x + e2y*dB3y + e3y*dB3z;
        s->stnS_.rs33() += e1z*dB3x + e2z*dB3y + e3z*dB3z;

        if (s->plasticStrainActive_) {
            dB1x = psj->x() * e3x;
            dB2x = psj->x() * e3y;
            dB3x = psj->x() * e3z;
            dB1y = psj->y() * e3x;
            dB2y = psj->y() * e3y;
            dB3y = psj->y() * e3z;
            dB1z = psj->x() * e1x + psj->y() * e2x + psj->z() * e3x;
            dB2z = psj->x() * e1y + psj->y() * e2y + psj->z() * e3y;
            dB3z = psj->x() * e1z + psj->y() * e2z + psj->z() * e3z;
            double VSubZone = s->getSubZoneVolume();
            s->stnPSI_.rs11() += (e1x * dB1x + e2x * dB1y + e3x * dB1z) * VSubZone;
            s->stnPSI_.rs12() += (e1x * dB2x + e2x * dB2y + e3x * dB2z) * VSubZone;
            s->stnPSI_.rs13() += (e1x * dB3x + e2x * dB3y + e3x * dB3z) * VSubZone;
            s->stnPSI_.rs22() += (e1y * dB2x + e2y * dB2y + e3y * dB2z) * VSubZone;
            s->stnPSI_.rs23() += (e1y * dB3x + e2y * dB3y + e3y * dB3z) * VSubZone;
            s->stnPSI_.rs33() += (e1z * dB3x + e2z * dB3y + e3z * dB3z) * VSubZone;
        }
    }

    void ModelUbiquitous::apexCorrectionJoint(State* s, uint32* jPlasticity) {
        if (jFriction_ > 0.0) {
            double japex = jCohesion_ / tan(jFriction_*degrad);
            if ((sig_ + dsp33_)  > japex) { 
                s->state_ |= joint_tension_now;
                *jPlasticity = 2;
                dsp33_ = japex - sig_;
                dsp13_ = - sp13_;
                dsp23_ = - sp23_;
                dsp11_ = dsp33_ * e2_/e1_;
            }
        }
    }

    void ModelUbiquitous::tensionCorrectionJoint(State *s,uint32 *jPlasticity,const double &fjt, DVect3* psj) {
        s->state_ |= joint_tension_now;
        if(jPlasticity) *jPlasticity = 2;
        dsp33_ = -fjt;
        dsp11_ = dsp33_ * e2_/e1_;
        if (s->plasticStrainActive_ && psj) {
            psj->rz() += fjt / e1_;
        }
    }

    void ModelUbiquitous::shearCorrectionJoint(State *s,uint32 *jPlasticity,const double &fjs, DVect3* psj) {
        s->state_ |= joint_shear_now;
        if(jPlasticity) *jPlasticity = 1;
        double dRat = 0.0;
        if (tau_ > 1.0e-6 * jCohesion_)
            dRat = fjs * jc1_ * g2_ / tau_;
        dsp13_ = - dRat * sp13_;
        dsp23_ = - dRat * sp23_;
        dsp33_ = - fjs * jc2_ * e1_;
        dsp11_ = - fjs * jc2_ * e2_;
        if (s->plasticStrainActive_ && psj) {
            psj->rx() += fjs * jc1_ * (sp13_ / tau_);
            psj->ry() += fjs * jc1_ * (sp23_ / tau_);
            psj->rz() += fjs * jc2_;
        }
    }

    static const uint32 StackX = 10;
    static const uint32 StackY = 11;
    static const uint32 StackZ = 12;
    //
    void ModelUbiquitous::largeStrainCorrection(State *s, const Axes3D &aAxes) {
        if (!s->sub_zone_) {
            s->working_[StackX] = 0.0;
            s->working_[StackY] = 0.0;
            s->working_[StackZ] = 0.0;
        }
        DVect3 rot = s->getRotation();
        double e1x = aAxes.e1().x();
        double e1y = aAxes.e1().y();
        double e1z = aAxes.e1().z();
        double e2x = aAxes.e2().x();
        double e2y = aAxes.e2().y();
        double e2z = aAxes.e2().z();
        double e3x = aAxes.e3().x();
        double e3y = aAxes.e3().y();
        double e3z = aAxes.e3().z();
        s->working_[StackX] +=   (rot.x())*e3y + (rot.y())*e3z;
        s->working_[StackY] += - (rot.x())*e3x + (rot.z())*e3z;
        s->working_[StackZ] += - (rot.y())*e3x - (rot.z())*e3y;
        double dB3x = s->stnE_.s11()*e3x + s->stnE_.s12()*e3y + s->stnE_.s13()*e3z;
        double dB3y = s->stnE_.s12()*e3x + s->stnE_.s22()*e3y + s->stnE_.s23()*e3z;
        double dB3z = s->stnE_.s13()*e3x + s->stnE_.s23()*e3y + s->stnE_.s33()*e3z;
        double dSp13 = dB3x*e1x + dB3y*e1y + dB3z*e1z;
        double dSp23 = dB3x*e2x + dB3y*e2y + dB3z*e2z;
        s->working_[StackX] -= (dSp13*e1x + dSp23*e2x);
        s->working_[StackY] -= (dSp13*e1y + dSp23*e2y);
        s->working_[StackZ] -= (dSp13*e1z + dSp23*e2z);
        if (s->sub_zone_==s->total_sub_zones_-1) {
            double dAux = 1.0 / s->total_sub_zones_;
            normX_ += s->working_[StackX] * dAux;
            normY_ += s->working_[StackY] * dAux;
            normZ_ += s->working_[StackZ] * dAux;
            norma_.set_norm(DVect3(normX_, normY_, normZ_));
            DVect2 dv2 = Orientation3::getDipDDFromNorm(DVect3(normX_,normY_,normZ_));
            dip_ = dv2.x() / degrad;
            dd_  = dv2.y() / degrad;
            angle_ = Orientation3::getJointAngleFromNorm(DVect3(normX_,normY_,normZ_)) / degrad;
        }
    }

}//namespace models
// EOF
