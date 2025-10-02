#include "modelaniso.h"
#include "base/src/version.h"
#include <atomic>

namespace models {    
    ModelAnisotropic::ModelAnisotropic() : ModelWeakplaneType(){
    }

    uint32 ModelAnisotropic::getMinorVersion() const {
        return version.revision_;
    }

    string ModelAnisotropic::getProperties() const {
        return "young-plane,young-normal,poisson-plane,poisson-normal,shear-normal,"
               + ModelWeakplaneType::getProperties();
    }

    string ModelAnisotropic::getStates() const {
        return "";
    }

    double ModelAnisotropic::getConfinedModulus(void) const {
        return getBulkModulus() + d4d3*getShearModulus();
    }

    double ModelAnisotropic::getShearModulus(void) const {
        double G12 = isD0(Nu12_+1.0) ? 0.0 : 0.5*E1_/(1.0+Nu12_);
        return (G12 > G13_ ? G12 : G13_);
    }

    double ModelAnisotropic::getBulkModulus(void) const {
        double M = E1_ > E3_ ? E1_ : E3_;
        double N = Nu12_ > Nu13_ ? Nu12_ : Nu13_;
        return M/(3.0*(1.0-2.0*N));
    }

    static std::atomic<uint64> n{0};
    base::Property ModelAnisotropic::getProperty(uint32 index) const {
        if (not n)
            n = (split(ModelWeakplaneType::getProperties(), ",")).size();
        if (index <= 5) {
            switch (index) {
            case 1: return (E1_);
            case 2: return (E3_);
            case 3: return (Nu12_);
            case 4: return (Nu13_);
            case 5: return (G13_);
            }
        }
        else if (index - 5 <= n)
            return ModelWeakplaneType::getProperty(index - 5);

        return 0.0;
    }

    void ModelAnisotropic::setProperty(
        uint32 index, const base::Property& p, uint32 restoreVersion) {
        ConstitutiveModel::setProperty(index, p, restoreVersion);
        if (not n)
            n = (split(ModelWeakplaneType::getProperties(), ",")).size();
        if (index <= 5) {
            switch (index) {
            case 1: E1_ = p.to<double>(); break;
            case 2: E3_ = p.to<double>(); break;
            case 3: Nu12_ = p.to<double>(); break;
            case 4: Nu13_ = p.to<double>(); break;
            case 5: G13_ = p.to<double>(); break;
            }
        } else if (index - 5 <= n) {
            ModelWeakplaneType::setProperty(index - 5, p, restoreVersion);
        }
    }

    void ModelAnisotropic::copy(const ConstitutiveModel *m) {
        ConstitutiveModel::copy(m);
        const ModelAnisotropic *mm = dynamic_cast<const ModelAnisotropic *>(m);
        if (!mm) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        //
        E1_   = mm->E1_;
        E3_   = mm->E3_;
        Nu12_ = mm->Nu12_;
        Nu13_ = mm->Nu13_;
        G13_  = mm->G13_;
        //
        dip_  = mm->dip_;
        dd_   = mm->dd_;
        norma_ = mm->norma_;
        angle_= mm->angle_;
    }

    void ModelAnisotropic::initialize(uint32 d,State *s)  {
        ConstitutiveModel::initialize(d,s);

        if (E1_ <= 0.0 || E3_ <= 0.0 )
            throw std::runtime_error("Anisotropic Model - Invalid Young's modulus.");
        if (G13_ <= 0.0 )
            throw std::runtime_error("Anisotropic Model - Invalid shear modulus.");    
        if (Nu12_ >= 0.5 || Nu13_ >= 0.5 || Nu12_ < -1.0 || Nu13_ < -1.0)
            throw std::runtime_error("Anisotropic Model - Invalid Poisson ratio.");

        /* --- transversely isotropic case --- */
        double E2   = E1_;
        double Nu23 = Nu13_;
        double G23  = G13_;
        double G12 = isD0(Nu12_+1.0) ? 0.0 : 0.5*E1_/(1.0+Nu12_);
        if (E3_ <= 0.0) E3_ = E1_;

        buildMatrix(E1_,E2,E3_,Nu23,Nu13_,Nu12_,G23,G13_,G12,dd_,dip_);
    }

    void ModelAnisotropic::run(uint32 d,State *s) {
        ConstitutiveModel::run(d,s);
        AnisotropicElasticTrial(s);
        s->viscous_ = false;
    }

    void ModelAnisotropic::AnisotropicElasticTrial(State *s) {
        double DG12 = 2.0*s->stnE_.s12();
        double DG13 = 2.0*s->stnE_.s13();
        double DG23 = 2.0*s->stnE_.s23();
        double DS11 = s->stnE_.s11();
        double DS22 = s->stnE_.s22();
        double DS33 = s->stnE_.s33();
        s->stnS_.rs11() += A11_*DS11 + A12_*DS22 + A13_*DS33 + A14_*DG12 + A15_*DG13 + A16_*DG23;
        s->stnS_.rs22() += A12_*DS11 + A22_*DS22 + A23_*DS33 + A24_*DG12 + A25_*DG13 + A26_*DG23;
        s->stnS_.rs33() += A13_*DS11 + A23_*DS22 + A33_*DS33 + A34_*DG12 + A35_*DG13 + A36_*DG23;
        s->stnS_.rs12() += A14_*DS11 + A24_*DS22 + A34_*DS33 + A44_*DG12 + A45_*DG13 + A46_*DG23;
        s->stnS_.rs13() += A15_*DS11 + A25_*DS22 + A35_*DS33 + A45_*DG12 + A55_*DG13 + A56_*DG23;
        s->stnS_.rs23() += A16_*DS11 + A26_*DS22 + A36_*DS33 + A46_*DG12 + A56_*DG13 + A66_*DG23;
    }

    /* inverts matrix b(6,6) */
    /* returns False if singular */
    bool ModelAnisotropic::xmatinv(double b[6][6]) {
        bool flag = true;
        for (uint32 i = 0; i < 6; ++i) {
            double bii = b[i][i];
            if (std::abs(bii) <= limits<double>::epsilon())
                return false; // Singular matrix, inversion not possible
            double inv_bii = 1.0 / bii;
            // Normalize the current row and invert diagonal
            for (uint32 k = 0; k < 6; ++k) 
                b[i][k] *= inv_bii;
            b[i][i] = inv_bii;
            // Eliminate other rows
            for (uint32 j = 0; j < 6; ++j) {
                if (j == i) continue;
                double bji = b[j][i];
                for (uint32 k = 0; k < 6; ++k)
                    b[j][k] -= bji * b[i][k];
                b[j][i] = -bji * inv_bii;
            }
        }
        return flag;
    }

    /* multiplication of square matrices : B3 = B1*B2 */
    void ModelAnisotropic::xmatmul(double b1[6][6],double b2[6][6],double b3[6][6]) {
        for (uint32 i=0; i<6; ++i) {
            for (uint32 j=0; j<6; ++j) {
                double sum = 0.0;
                for (uint32 m=0; m<6; ++m) 
                    sum += b1[i][m] * b2[m][j];
                b3[i][j] = sum;
            }
        }
    }

    void ModelAnisotropic::buildMatrix(const double& e1, const double& e2, const double& e3,
                                       const double& v1, const double& v2, const double& v3, 
                                       const double& g1, const double& g2, const double& g3, 
                                       const double& dd, const double& dip, const double& rot) {
        double adBB[6][6] = {};
        double adTT[6][6] = {};
        double adTTa[6][6] = {};
        double adTTb[6][6] = {};
        double adTTc[6][6] = {};
        double adTTw[6][6] = {};

        // Initialize adBB with local compliance matrix values
        if (e1 > 0) adBB[0][0] = 1.0 / e1;
        if (e2 > 0) adBB[0][1] = -v3 / e2;
        if (e3 > 0) adBB[0][2] = -v2 / e3;
        adBB[1][0] = adBB[0][1];
        if (e2 > 0) adBB[1][1] = 1.0 / e2;
        if (e3 > 0) adBB[1][2] = -v1 / e3;
        adBB[2][0] = adBB[0][2];
        adBB[2][1] = adBB[1][2];
        if (e3 > 0) adBB[2][2] = 1.0 / e3;
        if (g3 > 0) adBB[3][3] = 1.0 / g3;
        if (g2 > 0) adBB[4][4] = 1.0 / g2;
        if (g1 > 0) adBB[5][5] = 1.0 / g1;

        // Invert the local compliance matrix
        if (xmatinv(adBB)) {
            // Strain transformation matrix 
            // (A) rotate about Znew to get Xnew along horiz. proj. of dip dir
            double ang = (90.0 - dd) * degrad;
            double C = cos(ang);   
            double S = sin(ang);
            double C2 = C * C;
            double S2 = S * S;
            double SC = S * C;
            // 6x6 rotational matrix Q
            adTTb[0][0] = C2; 
            adTTb[0][1] = S2;
            adTTb[0][3] = SC;
            adTTb[1][0] = S2;
            adTTb[1][1] = C2;
            adTTb[1][3] = -SC;
            adTTb[2][2] = 1.0;
            adTTb[3][0] = -2.0 * SC;
            adTTb[3][1] = 2.0 * SC;
            adTTb[3][3] = C2 - S2;
            adTTb[4][4] = C;
            adTTb[4][5] = S;
            adTTb[5][4] = -S;
            adTTb[5][5] = C;
            // (B) rotate about Ynew to get Xnew along dip dir
            ang = - dip * degrad;
            C = cos(ang);
            S = sin(ang);
            C2 = C * C;
            S2 = S * S;
            SC = S * C;
            adTTc[0][0] = C2;
            adTTc[0][2] = S2;
            adTTc[0][4] = SC;
            adTTc[2][0] = S2;
            adTTc[2][2] = C2;
            adTTc[2][4] = -SC;
            adTTc[1][1] = 1.0;
            adTTc[4][0] = -2.0 * SC;
            adTTc[4][2] = 2.0 * SC;
            adTTc[4][4] = C2 - S2;
            adTTc[3][3] = C;
            adTTc[3][5] = S;
            adTTc[5][3] = -S;
            adTTc[5][5] = C;

            // composition of transformation matrices 
            if (isD0(rot)) 
                xmatmul(adTTc, adTTb, adTT); //adTT = adTTc * adTTb
            else {
                // (C) rotate around Znew
                ang = rot * degrad;
                double C = cos(ang);
                double S = sin(ang);
                double C2 = C * C;
                double S2 = S * S;
                double SC = S * C;
                adTTa[0][0] = C2;
                adTTa[0][1] = S2;
                adTTa[0][3] = SC;
                adTTa[1][0] = S2;
                adTTa[1][1] = C2;
                adTTa[1][3] = -SC;
                adTTa[2][2] = 1.0;
                adTTa[3][0] = -2.0 * SC;
                adTTa[3][1] = 2.0 * SC;
                adTTa[3][3] = C2 - S2;
                adTTa[4][4] = C;
                adTTa[4][5] = S;
                adTTa[5][4] = -S;
                adTTa[5][5] = C;
                xmatmul(adTTa, adTTc, adTTw);
                xmatmul(adTTw, adTTb, adTT);
            }

            // global elasticity matrix 
            xmatmul(adBB,adTT,adTTw); //adTTw=K

            for (uint32 i=0;i<6;i++) {
                for (uint32 j=0;j<6;j++) 
                    adTTa[i][j] = adTT[j][i];
            }
            xmatmul(adTTa, adTTw, adBB);
            //
            A11_ = adBB[0][0];
            A12_ = adBB[0][1];
            A13_ = adBB[0][2];
            A14_ = adBB[0][3];
            A15_ = adBB[0][4];
            A16_ = adBB[0][5];
            A22_ = adBB[1][1];
            A23_ = adBB[1][2];
            A24_ = adBB[1][3];
            A25_ = adBB[1][4];
            A26_ = adBB[1][5];
            A33_ = adBB[2][2];
            A34_ = adBB[2][3];
            A35_ = adBB[2][4];
            A36_ = adBB[2][5];
            A44_ = adBB[3][3];
            A45_ = adBB[3][4];
            A46_ = adBB[3][5];
            A55_ = adBB[4][4];
            A56_ = adBB[4][5];
            A66_ = adBB[5][5];
        }
    }

}//namespace models
