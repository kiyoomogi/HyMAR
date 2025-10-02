#include "modelelastic.h"
#include "models/src/convert.h"
#include "base/src/version.h"

namespace models {

    ModelElastic::ModelElastic(unsigned short option) : ConstitutiveModel(option) {
        FP_S;
    }

    uint32 ModelElastic::getMinorVersion() const {
        return version.revision_;
    }

    base::Property ModelElastic::getProperty(uint32 index) const {
        switch (index) {
        case 1: return bulk_;
        case 2: return shear_;
        case 3: { double young;
            getYPfromBS(bulk_,shear_,&young,nullptr);
            return young; }
        case 4: { double poisson;
            getYPfromBS(bulk_,shear_,nullptr,&poisson);
            return poisson; }
        }
        return(0.0);
    }

    void ModelElastic::setProperty(uint32 index,const base::Property &p,uint32 restoreVersion) {
        ConstitutiveModel::setProperty(index,p,restoreVersion);

        switch (index)  {
        case 1: bulk_  = p.to<double>(); break;
        case 2: shear_ = p.to<double>(); break;
        case 3: {
            double young(0),poisson(0);
            getYPfromBS(bulk_,shear_,&young,&poisson);
            young = p.to<double>();
            if (!restoreVersion) getBSfromYP(young,poisson,&bulk_,&shear_); } break;
        case 4: {
            double young(0),poisson(0);
            getYPfromBS(bulk_,shear_,&young,&poisson);
            poisson = p.to<double>();
            if (young <= 0.0 && !restoreVersion)
                throw std::runtime_error("Young's Modulus must be non-zero before you can specify Poisson's Ratio.");
            if (!restoreVersion) getBSfromYP(young,poisson,&bulk_,&shear_); } break;
        }
    }

    void ModelElastic::copy(const ConstitutiveModel *m) {
        ConstitutiveModel::copy(m);
        const ModelElastic *em = dynamic_cast<const ModelElastic *>(m);
        if (!em) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        bulk_ = em->bulk_;
        shear_ = em->shear_;
    }

    void ModelElastic::initialize(uint32 d,State *s) {
        ConstitutiveModel::initialize(d,s);
        updateParameters(s);
    }

    void ModelElastic::run(uint32 d,State *s) {
        ConstitutiveModel::run(d,s);

        if (s->modulus_reduction_factor_ > 0.0)
            moduliReduction(s->modulus_reduction_factor_,s);

        elasticTrial(s);
        s->viscous_ = true;
    }

    void ModelElastic::elasticTrial(State *s) {
        const double e11 = s->stnE_.s11();
        if (supportsUniaxial() && s && s->condition_==State::CondType::D1) {
            assert(isD0(s->stnE_.s12()));
            assert(isD0(s->stnE_.s13()));
            assert(isD0(s->stnE_.s23()));
            s->stnS_.rs11() += e1_ * e11;
        } else if (supportsPlaneStress() && s && s->condition_==State::CondType::D2) {
            assert(isD0(s->stnE_.s13()));
            assert(isD0(s->stnE_.s23()));
            const double e22 = s->stnE_.s22();
            s->stnS_.rs11() += e1_*e11 + e2_*e22;
            s->stnS_.rs22() += e2_*e11 + e1_*e22;
            s->stnS_.rs12() += s->stnE_.s12() * g2_;
        } else {
            const double e22 = s->stnE_.s22();
            const double e33 = s->stnE_.s33();
            s->stnS_.rs11() += e11*e1_ + (e22 + e33)*e2_;
            s->stnS_.rs22() += (e11 + e33)*e2_ + e22*e1_;
            s->stnS_.rs33() += (e11 + e22)*e2_ + e33*e1_;
            s->stnS_.rs12() += s->stnE_.s12() * g2_;
            s->stnS_.rs13() += s->stnE_.s13() * g2_;
            s->stnS_.rs23() += s->stnE_.s23() * g2_;
        }
    }

    bool ModelElastic::updateParameters(State *s, bool bEUpdated) {
        if (supportsUniaxial() && s && s->condition_==State::CondType::D1) {
            e1_ = 9.0*bulk_*shear_ / (3.0*bulk_ + shear_);
        }
        else if (supportsPlaneStress() && s && s->condition_==State::CondType::D2) {
            const double e10 = bulk_ + shear_*d4d3;
            const double e20 = bulk_ - shear_*d2d3;
            g2_ = 2.0*shear_;
            e1_ = e10 - e20*e20/e10;
            e2_ = e1_ - g2_;
        }
        else {
            e1_ = bulk_ + shear_*d4d3;
            e2_ = bulk_ - shear_*d2d3;
            g2_ = 2.0*shear_;
        }
        return !bEUpdated;
    }

    double ModelElastic::moduliReduction(const double &factor, State *s) {
        double shear_new = shear_ * factor;
        if (supportsUniaxial() && s && s->condition_==State::CondType::D1) {
            e1_ = 9.0*bulk_*shear_new / (3.0*bulk_ + shear_new);
        }
        else if (supportsPlaneStress() && s && s->condition_==State::CondType::D2) {
            double e10 = bulk_ + shear_new*d4d3;
            double e20 = bulk_ - shear_new*d2d3;
            g2_ = 2.0*shear_new;
            e1_ = e10 - e20*e20/e10;
            e2_ = e1_ - g2_;
        }
        else {
            e1_ = bulk_ + shear_new*d4d3;
            e2_ = bulk_ - shear_new*d2d3;
            g2_ = 2.0 * shear_new;
        }
        return shear_new;
    }

} // namespace models
// EOF
