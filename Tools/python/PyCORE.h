#include "../../Tools/JetCorrector.h"
#include "../../Tools/jetcorr/JetCorrectionUncertainty.h"
#include "../../Tools/jetcorr/SimpleJetCorrectionUncertainty.h"

// Want to be able to pass numpy arrays of jet information and get back
// a numpy array of JECs, b-tag SFs, ...
struct PyCOREBridge {
    std::vector<JetCorrectorParameters> JEC_vParam;
    FactorizedJetCorrector* JEC_corrector;
    void AddL1JECFile(const std::string& l1) { JEC_vParam.push_back(JetCorrectorParameters(l1)); }
    void AddL2JECFile(const std::string& l2) { JEC_vParam.push_back(JetCorrectorParameters(l2)); }
    void AddL3JECFile(const std::string& l3) { JEC_vParam.push_back(JetCorrectorParameters(l3)); }
    void MakeJetCorrector() { JEC_corrector = new FactorizedJetCorrector(JEC_vParam); }
    float GetJEC(float pt, float eta, float rho, float area) {
        JEC_corrector->setJetPt(pt);
        JEC_corrector->setJetEta(eta);
        JEC_corrector->setRho(rho);
        JEC_corrector->setJetA(area);
        return JEC_corrector->getCorrection();
    }
    // default numpy dtype is double
    void GetJECs(int N, Double_t* pts,Double_t* etas,Double_t* rhos,Double_t* areas,Double_t* outs) {
        for (int i = 0; i < N; i++) outs[i] = GetJEC(pts[i],etas[i],rhos[i],areas[i]);
    }
    // another overload in case dtype is np.float32
    void GetJECs(int N, Float_t* pts,Float_t* etas,Float_t* rhos,Float_t* areas,Float_t* outs) {
        for (int i = 0; i < N; i++) outs[i] = GetJEC(pts[i],etas[i],rhos[i],areas[i]);
    }
};
