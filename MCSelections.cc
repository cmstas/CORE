#include "MCSelections.h"
#include "Math/VectorUtil.h"

using namespace tas;

int dumpDocLines(){

  //nitialize particle database
  TDatabasePDG *pdg = new TDatabasePDG();

  //Print Header
  cout << "                " << "   pt    " << "  phi  " << "      eta   " << "    mass  " << "status " << "Mother  " << endl;     
  std::cout << "---------------------------------------------------------------------" << std::endl;

  //Loop over gen particles
  for (unsigned int j = 0; j < genps_id().size(); j++) {

    //mass
    float m2 = genps_p4().at(j).M2();
    float m  = m2 >= 0 ? sqrt(m2) : 0.0;

    //Print information about the jth gen particle
    const char* particle = (abs(genps_id().at(j)) == 4124) ? "Lambda_c*" : pdg->GetParticle(genps_id().at(j))->GetName();
    const char* mother_particle = (abs(genps_id_mother().at(j)) == 4124) ? "Lambda_c*" : pdg->GetParticle(genps_id_mother().at(j))->GetName();

    cout << setw(4)  << left  <<                    j                        << " "
         << setw(10) << left  <<                    particle                 << " "
         << setw(7)  << right << setprecision(4) << genps_p4().at(j).pt()    << "  "
         << setw(7)  << right << setprecision(4) << genps_p4().at(j).phi()   << "  "
         << setw(10) << right << setprecision(4) << genps_p4().at(j).eta()   << "  "
         << setw(7)  << right << setprecision(4) << m                        << "  "
         << setw(4)  << right <<                    genps_status().at(j)     << " "
         << setw(10) << left  <<                    mother_particle          << " " 
         << endl;

    //epton daughters of the jth gen particle
    if(genps_lepdaughter_id()[j].size() > 0) {
      cout << endl;
      cout << "  Daughters:" << endl;
      for(unsigned int i = 0; i < genps_lepdaughter_id()[j].size(); i++) {

        // mass
        float m2_daught = genps_lepdaughter_p4().at(j).at(i).M2();
        float m_daught  = m2_daught >= 0 ? sqrt(m2_daught) : 0.0;

        //Print information about the daughters
        const char* daughter_particle = abs(genps_lepdaughter_id().at(j).at(i)) == 4124 ? "Lambda_c*" : pdg->GetParticle(genps_lepdaughter_id().at(j).at(i))->GetName();
        cout << genps_lepdaughter_id().at(j).at(i) << endl;
        cout << setw(2)  << left  << "    "            << i                                                               << " "
             << setw(10) << left  << daughter_particle                                                                    << " "
             << setw(7)  << right << setprecision(4)   << genps_lepdaughter_p4().at(j).at(i).pt()                         << "  "
             << setw(7)  << right << setprecision(4)   << genps_lepdaughter_p4().at(j).at(i).phi()                        << "  "
             << setw(10) << right << setprecision(4)   << genps_lepdaughter_p4().at(j).at(i).eta()                        << "  "
             << setw(7)  << right << setprecision(4)   << m_daught                                                        << "  "  << endl;
      
      }//daughter-loop
      cout << endl;
    }
  
  }//gen-loop
  delete pdg;
  return 0;
}

bool isFromSUSY(int id, int idx){
  int mc_id         = abs(id) == 11          ? els_mc_id().at(idx)         : mus_mc_id().at(idx);
  int mc_motherid   = abs(id) == 11          ? els_mc_motherid().at(idx)   : mus_mc_motherid().at(idx);
  int mc3_motheridx = abs(id) == 11          ? els_mc3_motheridx().at(idx) : mus_mc3_motheridx().at(idx);
  int mc_grandmaid  = mc3_motheridx == -9999 ? -9999                       : genps_id_mother().at(mc3_motheridx);
  //For Leptons
  if ((abs(mc_id)==11 || abs(mc_id)==13)){
    if (abs(mc_motherid)==1000024) return true;
    if (abs(mc_motherid)==15 && abs(mc_grandmaid)==1000024) return true;
    return false;
  }
  return false;
}

bool isFromW(int id, int idx){
  int mc_id         = abs(id) == 11          ? els_mc_id().at(idx)         : mus_mc_id().at(idx);
  int mc_motherid   = abs(id) == 11          ? els_mc_motherid().at(idx)   : mus_mc_motherid().at(idx);
  int mc3_motheridx = abs(id) == 11          ? els_mc3_motheridx().at(idx) : mus_mc3_motheridx().at(idx);
  int mc3_id        = abs(id)==11            ? els_mc3_id().at(idx)        : mus_mc3_id().at(idx);
  int mc3_motherid  = abs(id)==11            ? els_mc3_motherid().at(idx)  : mus_mc3_motherid().at(idx);
  int mc_grandmaid  = mc3_motheridx == -9999 ? -9999                       : genps_id_mother().at(mc3_motheridx);
 
  //For Leptons
  if ((abs(mc_id)==11 || abs(mc_id)==13)){
    if (abs(mc_motherid)==24) return true;
    if (abs(mc_motherid)==15 && abs(mc_grandmaid)==24) return true;
    return false;
  }

  //For Photons
  if (abs(mc_id)==22){
    if (abs(mc_motherid)!=abs(id)) return false;
    if (abs(mc3_id)!=abs(id)) return false;
    if (abs(mc3_motherid)!=24) return false;
    return true;
  }

  //everything else
  return false;
}

bool isFromZ(int id, int idx){
  int mc_id         = abs(id) == 11          ? els_mc_id().at(idx)         : mus_mc_id().at(idx);
  int mc_motherid   = abs(id) == 11          ? els_mc_motherid().at(idx)   : mus_mc_motherid().at(idx);
  int mc3_motheridx = abs(id) == 11          ? els_mc3_motheridx().at(idx) : mus_mc3_motheridx().at(idx);
  int mc3_id        = abs(id) == 11          ? els_mc3_id().at(idx)        : mus_mc3_id().at(idx);
  int mc3_motherid  = abs(id) == 11          ? els_mc3_motherid().at(idx)  : mus_mc3_motherid().at(idx);
  int mc_grandmaid  = mc3_motheridx == -9999 ? -9999                       : genps_id_mother().at(mc3_motheridx);

  //For Leptons
  if ((abs(mc_id)==11 || abs(mc_id)==13)){
    if (abs(mc_motherid)==23) return true;
    if (abs(mc_motherid)==15 && abs(mc_grandmaid)==23) return true;
    return false;
  }

  //For Photons
  if (abs(mc_id)==22){
    if (abs(mc_motherid)!=abs(id)) return false;
    if (abs(mc3_id)!=abs(id)) return false;
    if (abs(mc3_motherid)!=23) return false;
    return true;
  }

  //For everything else
  return false;
}

bool isFromB(int id, int idx){
  int mc_id       = abs(id) == 11 ? els_mc_id().at(idx)       : mus_mc_id().at(idx);
  int mc_motherid = abs(id) == 11 ? els_mc_motherid().at(idx) : mus_mc_motherid().at(idx);
  return (abs(mc_id)==11 || abs(mc_id)==13) && idIsBeauty(mc_motherid);
}

bool isFromC(int id, int idx){
  int mc_id       = abs(id) == 11 ? els_mc_id().at(idx)       : mus_mc_id().at(idx);
  int mc_motherid = abs(id) == 11 ? els_mc_motherid().at(idx) : mus_mc_motherid().at(idx);
  return (abs(mc_id)==11 || abs(mc_id)==13) && idIsCharm(mc_motherid);
}

bool isFromLight(int id, int idx){
  int mc_id       = abs(id) == 11 ? els_mc_id().at(idx)       : mus_mc_id().at(idx);
  int mc_motherid = abs(id) == 11 ? els_mc_motherid().at(idx) : mus_mc_motherid().at(idx);
  if (abs(mc_id) != 11 && abs(mc_id) != 13) return false;
  if (abs(mc_motherid)>200 && abs(mc_motherid)<400) return true; 
  if (abs(mc_motherid)>0 && abs(mc_motherid)<4) return true;
  return false;
}

bool isFromLightFake(int id, int idx){
  int mc_id       = abs(id) == 11 ? els_mc_id().at(idx)       : mus_mc_id().at(idx);
  int mc_motherid = abs(id) == 11 ? els_mc_motherid().at(idx) : mus_mc_motherid().at(idx);

  //True if light hadron
  if (abs(mc_id)>200 && abs(mc_id)<400) return true;

  //True if non-quark whose mother is light hadron
  if (abs(mc_id) == 11) return false;
  if (abs(mc_id) == 13) return false;
  if (abs(mc_motherid)>200 && abs(mc_motherid)<400) return true;
  if (abs(mc_motherid)>0 && abs(mc_motherid)<4) return true;
 
  return false;
}

bool idIsCharm(int id) {
  id = abs(id);
  if (id == 4      ) return true;
  if (id == 411    ) return true;
  if (id == 421    ) return true;
  if (id == 10411  ) return true;
  if (id == 10421  ) return true;
  if (id == 413    ) return true;
  if (id == 423    ) return true;
  if (id == 10413  ) return true;
  if (id == 10423  ) return true;
  if (id == 20413  ) return true;
  if (id == 20423  ) return true;
  if (id == 415    ) return true;
  if (id == 425    ) return true;
  if (id == 431    ) return true;
  if (id == 10431  ) return true;
  if (id == 433    ) return true;
  if (id == 10433  ) return true;
  if (id == 20433  ) return true;
  if (id == 435    ) return true;
  if (id == 441    ) return true;
  if (id == 10441  ) return true;
  if (id == 100441 ) return true;
  if (id == 443    ) return true;
  if (id == 10443  ) return true;
  if (id == 20443  ) return true;
  if (id == 100443 ) return true;
  if (id == 30443  ) return true;
  if (id == 9000443) return true;
  if (id == 9010443) return true;
  if (id == 9020443) return true;
  if (id == 445    ) return true;
  if (id == 9000445) return true;
  if (id == 4122   ) return true;
  if (id == 4222   ) return true;
  if (id == 4212   ) return true;
  if (id == 4112   ) return true;
  if (id == 4224   ) return true;
  if (id == 4214   ) return true;
  if (id == 4114   ) return true;
  if (id == 4232   ) return true;
  if (id == 4132   ) return true;
  if (id == 4322   ) return true;
  if (id == 4312   ) return true;
  if (id == 4324   ) return true;
  if (id == 4314   ) return true;
  if (id == 4332   ) return true;
  if (id == 4334   ) return true;
  if (id == 4412   ) return true;
  if (id == 4422   ) return true;
  if (id == 4414   ) return true;
  if (id == 4424   ) return true;
  if (id == 4432   ) return true;
  if (id == 4434   ) return true;
  if (id == 4444   ) return true;
  return false;
}

bool idIsBeauty(int id) {
  id = abs(id);
  if (id == 5     ) return true;
  if (id == 511   ) return true;
  if (id == 521   ) return true;
  if (id == 10511 ) return true;
  if (id == 10521 ) return true;
  if (id == 513   ) return true;
  if (id == 523   ) return true;
  if (id == 10513 ) return true;
  if (id == 10523 ) return true;
  if (id == 20513 ) return true;
  if (id == 20523 ) return true;
  if (id == 515   ) return true;
  if (id == 525   ) return true;
  if (id == 531   ) return true;
  if (id == 10531 ) return true;
  if (id == 533   ) return true;
  if (id == 10533 ) return true;
  if (id == 20533 ) return true;
  if (id == 535   ) return true;
  if (id == 541   ) return true;
  if (id == 10541 ) return true;
  if (id == 543   ) return true;
  if (id == 10543 ) return true;
  if (id == 20543 ) return true;
  if (id == 545   ) return true;
  if (id == 551   ) return true;
  if (id == 10551 ) return true;
  if (id == 100551) return true;
  if (id == 110551) return true;
  if (id == 200551) return true;
  if (id == 210551) return true;
  if (id == 553   ) return true;
  if (id == 10553 ) return true;
  if (id == 20553 ) return true;
  if (id == 30553 ) return true;
  if (id == 100553) return true;
  if (id == 110553) return true;
  if (id == 120553) return true;
  if (id == 130553) return true;
  if (id == 200553) return true;
  if (id == 210553) return true;
  if (id == 220553) return true;
  if (id == 300553) return true;
  if (id == 900055) return true;
  if (id == 901055) return true;
  if (id == 555   ) return true;
  if (id == 10555 ) return true;
  if (id == 20555 ) return true;
  if (id == 100555) return true;
  if (id == 110555) return true;
  if (id == 120555) return true;
  if (id == 200555) return true;
  if (id == 557   ) return true;
  if (id == 100557) return true;
  if (id == 5122  ) return true; 
  if (id == 5112  ) return true;
  if (id == 5212  ) return true;
  if (id == 5222  ) return true;
  if (id == 5114  ) return true;
  if (id == 5214  ) return true;
  if (id == 5224  ) return true;
  if (id == 5132  ) return true;
  if (id == 5232  ) return true;
  if (id == 5312  ) return true;
  if (id == 5322  ) return true;
  if (id == 5314  ) return true;
  if (id == 5324  ) return true;
  if (id == 5332  ) return true;
  if (id == 5334  ) return true;
  if (id == 5142  ) return true;
  if (id == 5242  ) return true;
  if (id == 5412  ) return true;
  if (id == 5422  ) return true;
  if (id == 5414  ) return true;
  if (id == 5424  ) return true;
  if (id == 5342  ) return true;
  if (id == 5432  ) return true;
  if (id == 5434  ) return true;
  if (id == 5442  ) return true;
  if (id == 5444  ) return true;
  if (id == 5512  ) return true;
  if (id == 5522  ) return true;
  if (id == 5514  ) return true;
  if (id == 5524  ) return true;
  if (id == 5532  ) return true;
  if (id == 5534  ) return true;
  if (id == 5542  ) return true;
  if (id == 5544  ) return true;
  if (id == 5554  ) return true;
  return false;
}

//________________________________________________________________
int getSourceId(int genpIdx){

  // look for earliest mother that is: top, W, Z, or H
  //   terminate search if top or H is found
  //   continue if W, Z is found, to see if there's a top or H
  //   can also add SUSY particles to logic, not yet done

  int currentIdx = genpIdx;
  int motherIdx = cms3.genps_idx_simplemother().at(currentIdx);
  int motherId = abs(cms3.genps_id_simplemother().at(currentIdx));
  int sourceId = -1;
  if (useSourceId(motherId)) sourceId = motherId;

  while ( (currentIdx >= 0) && (motherIdx >= 0) && !terminateSourceId(sourceId) && !terminateMotherId(motherId) ) {
    currentIdx = motherIdx;
    motherIdx = cms3.genps_idx_simplemother().at(currentIdx);
    motherId = abs(cms3.genps_id_simplemother().at(currentIdx));
    if (useSourceId(motherId)) sourceId = motherId;
  }

  return sourceId;

}

//________________________________________________________________
// list of sourceIds to terminate search on
//   can add susy particles here
bool terminateSourceId(int sourceId){
  int id = abs(sourceId);
  if (id == 6 || id == 25) return true;
  return false;
}

//________________________________________________________________
// list of sourceIds to allow
//   top, Z, W, H
//   can add susy particles here
bool useSourceId(int motherId){
  int id = abs(motherId);
  if (id == 6 || (id >= 23 && id <= 25)) return true;
  return false;
}

//________________________________________________________________
// list of mother PdgIds to terminate search on
//   quarks, gluons, hadrons
bool terminateMotherId(int motherId){
  int id = abs(motherId);
  if (id <= 5 || id == 21 || (id > 100 && id < 1000000)) return true;
  return false;
}

//________________________________________________________________
// top pt reweighting from the TOP PAG
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
//  using the 8 TeV "All Combined" weights
float topPtWeight(float pt_top, float pt_tbar) {
  if (pt_top < 0. || pt_tbar < 0.) return 1.;
  float weight_top = exp(0.156 - 0.00137 * pt_top);
  float weight_tbar = exp(0.156 - 0.00137 * pt_tbar);
  return sqrt( weight_top * weight_tbar );
}

//________________________________________________________________
// top pt reweighting from the TOP PAG
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
//  using the 8 TeV "All Combined" weights
//  adding the recommendation to cutoff the reweighting at pt of 400
float topPtWeight_cutoff(float pt_top, float pt_tbar) {
  if (pt_top < 0. || pt_tbar < 0.) return 1.;
  float pt_top_cutoff = std::min(pt_top, float(400.));
  float pt_tbar_cutoff = std::min(pt_tbar, float(400.));
  float weight_top = exp(0.156 - 0.00137 * pt_top_cutoff);
  float weight_tbar = exp(0.156 - 0.00137 * pt_tbar_cutoff);
  return sqrt( weight_top * weight_tbar );
}

//________________________________________________________________
// code to get number of reco jets matched to ISR jets
// inputs are jets with pt > 30, |eta| < 2.4, PF Loose ID, leptons removed
// https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
int get_nisrMatch ( const std::vector<LorentzVector>& clean_jets ) {

  int nisr(0);
  for (size_t ijet(0); ijet<clean_jets.size(); ijet++){
    if (clean_jets.at(ijet).pt() < 30.0) continue;
    if (fabs(clean_jets.at(ijet).eta()) > 2.4) continue;
    bool matched=false;
    for (size_t imc = 0; imc < cms3.genps_id().size(); imc++) {
      if (matched) break;
      if (cms3.genps_status().at(imc)!=23 || abs(cms3.genps_id().at(imc))>5) continue;
      int momid = abs(cms3.genps_id_simplemother().at(imc));
      if(!(momid==6 || momid==23 || momid==24 || momid==25 || momid>1e6)) continue; 
      // check against daughter in case of hard initial splitting
      // - since we don't save daughters, only mothers, we do a double loop here...
      for (size_t jmc = 0; jmc < cms3.genps_id().size(); jmc++) {
	if (cms3.genps_idx_simplemother().at(jmc) != int(imc) ) continue;
	float dR = ROOT::Math::VectorUtil::DeltaR(clean_jets.at(ijet), cms3.genps_p4().at(jmc));
        if(dR<0.3){
	  matched = true;
	  break;
        }
      }
    } // Loop over MC particles
    if(!matched) {
      nisr++;
    }
  } // Loop over jets

  return nisr;
}

//________________________________________________________________
// returns year dependent ISR weight based on nisrMatch
float get_isrWeight( int nisrMatch, int year ) {
  if (year >= 2016) {
    // Moriond 2017 values for ISR weight based on nisrMatch
    if      (nisrMatch == 0) return 1.000;
    else if (nisrMatch == 1) return 0.920;
    else if (nisrMatch == 2) return 0.821;
    else if (nisrMatch == 3) return 0.715;
    else if (nisrMatch == 4) return 0.662;
    else if (nisrMatch == 5) return 0.561;
    else if (nisrMatch >= 6) return 0.511;
  }
  return 1.;
}

//________________________________________________________________
// returns year dependent ISR weight uncertainty based on nisrMatch
float get_isrUnc( int nisrMatch, int year ) {
  if (year >= 2016) {
    // Moriond 2017 values
    if      (nisrMatch == 0) return 0.000;
    else if (nisrMatch == 1) return 0.040;
    else if (nisrMatch == 2) return 0.090;
    else if (nisrMatch == 3) return 0.143;
    else if (nisrMatch == 4) return 0.170;
    else if (nisrMatch == 5) return 0.221;
    else if (nisrMatch >= 6) return 0.258;
  }
  return 0.;
}

//________________________________________________________________
// returns 2016 ICHEP ISR weight based on nisrMatch
// s13 of https://indico.cern.ch/event/557678/contributions/2247944/attachments/1311994/1963568/16-07-19_ana_manuelf_isr.pdf
float get_isrWeight( int nisrMatch ) {

  if      (nisrMatch == 0) return 1.;
  else if (nisrMatch == 1) return 0.882;
  else if (nisrMatch == 2) return 0.792;
  else if (nisrMatch == 3) return 0.702;
  else if (nisrMatch == 4) return 0.648;
  else if (nisrMatch == 5) return 0.601;
  else if (nisrMatch >= 6) return 0.515;

  // shouldn't get here
  std::cout << "ERROR: get_isrWeight: bad input for nisrMatch: " << nisrMatch << std::endl;
  return 1.;
}

//________________________________________________________________
// returns 2016 ICHEP ISR weight uncertainty based on nisrMatch
// s13 of https://indico.cern.ch/event/557678/contributions/2247944/attachments/1311994/1963568/16-07-19_ana_manuelf_isr.pdf
float get_isrUnc( int nisrMatch ) {

  if      (nisrMatch == 0) return 0.;
  else if (nisrMatch == 1) return 0.059;
  else if (nisrMatch == 2) return 0.104;
  else if (nisrMatch == 3) return 0.149;
  else if (nisrMatch == 4) return 0.176;
  else if (nisrMatch == 5) return 0.199;
  else if (nisrMatch >= 6) return 0.242;

  // shouldn't get here
  std::cout << "ERROR: get_isrUnc: bad input for nisrMatch: " << nisrMatch << std::endl;
  return 0.;
}

//________________________________________________________________
// returns 2017 Moriond ISR weight based on ISR pt, for ewkino models
// s6 of https://www.dropbox.com/s/rv79odbk1plltm1/17-02-22_ana_isr_ewk.pdf?dl=1 (link to be updated once in a more permanent place)
float get_isrWeight_ewk( float isr_pt ) {

  if      (isr_pt < 50. ) return 1.000;
  else if (isr_pt < 100.) return 1.052;
  else if (isr_pt < 150.) return 1.179;
  else if (isr_pt < 200.) return 1.150;
  else if (isr_pt < 300.) return 1.057;
  else if (isr_pt < 400.) return 1.000;
  else if (isr_pt < 600.) return 0.912;
  else                    return 0.783;

  // shouldn't get here
  std::cout << "ERROR: get_isrWeight_ewk: bad input for isr_pt: " << isr_pt << std::endl;
  return 1.;
}

