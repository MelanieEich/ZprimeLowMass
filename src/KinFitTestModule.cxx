#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/CommonModules.h>
#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>

#include <UHH2/ZprimeLowMass/include/KinFitTestSelections.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestSetup.h>
//#include <UHH2/ZprimeLowMass/include/KinFitTestFit.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestHists.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestHistsMuonUnc.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestHistsPulls.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestHistsFitted.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestRECOHists.h>
#include <UHH2/ZprimeLowMass/include/KinFitTestFITHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/AdditionalSelections.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/LuminosityHists.h>

#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicUtils.h>
#include <UHH2/ZprimeSemiLeptonic/include/ZprimeSemiLeptonicSelections.h>

//KinFit includes
#include <UHH2/KinFitter/interface/TKinFitter.h>
#include <UHH2/KinFitter/interface/TFitParticleEMomDev.h>
#include <UHH2/KinFitter/interface/TFitParticleEtEtaPhi.h>
#include <UHH2/KinFitter/interface/TFitParticleEScaledMomDev.h>
#include <UHH2/KinFitter/interface/TFitConstraintEp.h>
//#include <UHH2/KinFitter/interface/TFitConstraintMGaus.h>
#include <UHH2/KinFitter/interface/TFitConstraintM.h>




using namespace std;
using namespace uhh2;

namespace uhh2examples {


class KinFitTestModule: public AnalysisModule {
public:
  explicit KinFitTestModule(Context & ctx);
  bool process(Event & event) override;
  void setupJets();
  void setupLeptons();
  void setupConstraints();


//  enum Param{ kEMom, kEtEtaPhi};
  /// supported constraints
  enum Constraint { kWHadMass = 1, kWLepMass, kTopHadMass, kTopLepMass, kSumPt};
  //enum Constraint { kWHadMass = 1, kWLepMass, kEqualTopMasses, kSumPt};//,   kSumPt,  kEqualTopMasses,  kEqualWMasses, ,kMET,
  //  enum ObjectType{kUdscJet, kBJet, kMuon, kElectron, kMet};
  void setBestFit(Event & event, std::vector<TLorentzVector> solution_vec, std::vector<TLorentzVector> input_vector, KinFitTestSetup::ObjectType leptonType,  float chi2_lowest,float fit_prob, float fit_status, bool fit_found, bool is_input, std::vector<int> jet_combi);
  void fit(Event & event, int bl, int bh, int h1, int h2, KinFitTestSetup::ObjectType leptonType, KinFitTestSetup::Param jetParam_, KinFitTestSetup::Param lepParam_, KinFitTestSetup::Param metParam_);
  bool dm_check(Event & event, int bl, int bh, int h1, int h2);
  //  TMatrixD setupMatrix(const TLorentzVector& object, const ObjectType objType, const Param param);

  void setPull(int jet_variable_size, const TTbarGen * ttgen, int all_best, bool is_MC_);

  bool fit_found, met_cut, isMC_, high_met, high_blep;
  bool matchable, b_jet_matching, match_ttbar_decay, correct_match;
  std::vector<TLorentzVector> solution;
  std::vector<TLorentzVector> solution2;
  std::vector<TLorentzVector> input_vectors;
  std::vector<TLorentzVector> input_vectors_fit_found;
  std::vector<TLorentzVector> input_vectors_all;
  std::vector<float> sol_met_pz;
  std::vector<int> jet_combi_best;
  std::vector<int> jet_combi_all;
  float chi2_lowest, fit_prob, dR_tlep, dR_thad, dR_wlep, dR_jetq1, dR_jetq2, dR_jetbh, dR_jetbl, dR_lepton, dR_neutrino, dR_combined, fit_status, fit_status_each;
  double dm_w_lep, dm_w_had, dm_w_lep_all, dm_w_had_all, dm_top_lep, dm_top_had, dm_top_lep_all, dm_top_had_all;

  //Fit particles
  TAbsFitParticle* jethad1_;
  TAbsFitParticle* jethad2_;
  TAbsFitParticle* jethadb_;
  TAbsFitParticle* jetlepb_;
  TAbsFitParticle* lepton_;
  TAbsFitParticle* met_;


private:
  //histograms
  std::unique_ptr<Hists> h_nocuts, h_trigger, h_lepton_size, h_lepton_jet_sel, h_before_metcut, h_after_metcut, h_aftersel, h_fitbest, h_fitbest_chi2_cut;
  std::unique_ptr<MuonHists> h_mu_nocuts, h_mu_trigger, h_mu_lepton_size, h_mu_lepton_jet_sel,  h_mu_before_metcut, h_mu_after_sel;
  std::unique_ptr<EventHists> h_event_nocuts, h_event_trigger, h_event_lepton_size, h_event_lepton_jet_sel, h_event_before_metcut, h_event_after_sel;
  std::unique_ptr<JetHists> h_jet_nocuts, h_jet_trigger, h_jet_lepton_size, h_jet_lepton_jet_sel, h_jet_before_metcut, h_jet_after_sel;
  //  std::unique_ptr<BTagMCEfficiencyHists> h_btag_hists;
  std::unique_ptr<KinFitTestFITHists> h_before_fit, h_fitted, h_fitted_allcomb, h_fitted_allcomb_fit_converged, h_fitted_0btag,h_fitted_1btag,h_fitted_2btag, h_fitted_matchable, h_fitted_matchable_0btag, h_fitted_matchable_1btag, h_fitted_matchable_2btag, h_fitted_not_matchable, h_fitted_chi2_matchable, h_fitted_chi2_not_matchable, h_fitted_b_jet_matching, h_fitted_b_jet_matching_0btag, h_fitted_b_jet_matching_1btag, h_fitted_b_jet_matching_2btag, h_fitted_match_ttbar_decay, h_fitted_match_ttbar_decay_b_jet_matching, h_fitted_match_ttbar_decay_b_jet_matching_0btag, h_fitted_match_ttbar_decay_b_jet_matching_1btag, h_fitted_match_ttbar_decay_b_jet_matching_2btag, h_fitted_correct_match, h_fitted_correct_match_chi2_cut, h_fitted_correct_match_0btag, h_fitted_correct_match_1btag, h_fitted_correct_match_2btag, h_fitted_dR_cut, h_fitted_chi2_cut;
  std::unique_ptr<KinFitTestRECOHists> h_reco,h_reco_matchable,h_reco_not_matchable,h_reco_matchable_chi2_30,h_reco_matchable_chi2_6, h_reco_not_matchable_chi2_6, h_reco_ttbar_m, h_reco_ttbar_m_chi2_30, h_reco_ttbar_m_chi2_6, h_reco_chi2_30, h_reco_chi2_0btag, h_reco_chi2_1btag, h_reco_chi2_2btag, h_reco_chi2_20, h_reco_chi2_10, h_reco_chi2_6;
  std::unique_ptr<KinFitTestHistsMuonUnc> h_muon_unc;
  std::unique_ptr<KinFitTestHistsPulls> h_pulls;

  //lumi hist
  std::unique_ptr<Hists> lumi_h;
  std::unique_ptr<Hists> lumi_after_selection_h;

   //KinFit object
  std::unique_ptr<TKinFitter> myKinFitter;

  /// jet parametrization
  KinFitTestSetup::Param  jetParam_;
  /// lepton parametrization
  KinFitTestSetup::Param lepParam_;
  /// met parametrization
  KinFitTestSetup::Param metParam_;
  //Lepton type
  KinFitTestSetup::ObjectType leptonType;

  //Objects for fit
  TFitConstraintEp* sumPxConstr_;
  TFitConstraintEp* sumPyConstr_;
  std::vector<Constraint> constrList_;
  Bool_t constrainSumPt_;

  //Fitted particles
  TLorentzVector jetLepB__p4;
  TLorentzVector jetHadB__p4;
  TLorentzVector jetHad1__p4;
  TLorentzVector jetHad2__p4;
  TLorentzVector lepton__p4;
  TLorentzVector met__p4;

  //input particles
  TLorentzVector p4Had1;
  TLorentzVector p4Had2;
  TLorentzVector p4HadB;
  TLorentzVector p4LepB;
  TLorentzVector p4Lep;
  TLorentzVector p4Neutrino;


  const TMatrixD* covMat;


protected:
  enum lepton { muon, elec, both};
  lepton channel_;

  bool do_fit, debug, pulls, unc_test;
  int matchable_;
  double mW_, decw_w, mTop_, decw_top, maxDeltaS, maxF;
  double chi2_curr, F_curr, chi2_prev, F_prev, fitprob_curr, fitprob_prev, fitprob_best;
  int nrjetcomb_before_dm_cut, nrjetcomb_after_dm_cut, maxNrIter, verbosity, nr_suitable_events;

  //cleaners
  std::unique_ptr<MuonCleaner>     muoSR_cleaner;
  std::unique_ptr<ElectronCleaner> eleSR_cleaner;
  std::unique_ptr<JetCleaner>      jet_IDcleaner;
  std::unique_ptr<JetCorrector>                    jet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer>     jetER_smearer;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> jetlepton_cleaner;
  std::unique_ptr<JetCleaner>                      jet_cleaner2;

  std::unique_ptr<uhh2::AnalysisModule> lumiweight;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;

  //selections
  std::unique_ptr<uhh2::Selection> lumi_sel;
  std::unique_ptr<uhh2::AndSelection> metfilters_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> jet4_sel;
  std::unique_ptr<uhh2::Selection> jet3_sel;
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> jet1_sel;
  std::unique_ptr<uhh2::Selection> ttagevt_sel;
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::Selection> trigger_sel;
  std::unique_ptr<uhh2::Selection> trigger2_sel;

  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> muonID_SF;
  std::unique_ptr<uhh2::AnalysisModule> muonHLT_SF;
  std::unique_ptr<uhh2::AnalysisModule> muonTRK_SF;
  std::unique_ptr<uhh2::AnalysisModule> btagSF;
 
  JetId btag_ID_;
  CSVBTag::wp b_working_point;
  std::vector<int> nrbtagedjet;

  //constraint list
  //  std::map<Constraint, TFitConstraintMGaus*> massConstr_;
  std::map<Constraint, TFitConstraintM*> massConstr_;

  Event::Handle<float> h_wgtMC__GEN;

  //fitted objects
  Event::Handle<int> parametrisation;
  Event::Handle<float> rec_fit_chi2;
  Event::Handle<float> fit_probab;
  Event::Handle<float> rec_fit_status;
  Event::Handle<bool> rec_fit_fit_found;
  Event::Handle<bool> rec_fit_is_input;
  Event::Handle<TLorentzVector> rec_fit_whad_v4;
  Event::Handle<TLorentzVector> rec_fit_wlep_v4;
  Event::Handle<TLorentzVector> rec_fit_thad_v4;
  Event::Handle<TLorentzVector> rec_fit_tlep_v4;
  Event::Handle<TLorentzVector> rec_fit_lepton_v4;
  Event::Handle<TLorentzVector> rec_fit_neutrino_v4;
  Event::Handle<TLorentzVector> rec_fit_bhad_v4;
  Event::Handle<TLorentzVector> rec_fit_blep_v4;
  Event::Handle<TLorentzVector> rec_fit_jetq1_v4;
  Event::Handle<TLorentzVector> rec_fit_jetq2_v4;
  Event::Handle<float> rec_fit_dR_combined; 

  Event::Handle<TLorentzVector> input_jetq1_v4; 
  Event::Handle<TLorentzVector> input_jetq2_v4; 
  Event::Handle<TLorentzVector> input_jetbh_v4; 
  Event::Handle<TLorentzVector> input_jetbl_v4; 
  Event::Handle<TLorentzVector> input_lepton_v4; 
  Event::Handle<TLorentzVector> input_neutrino_v4; 

  Event::Handle<int> jet_permutation_before_dm_cut;
  Event::Handle<int> jet_permutation_after_dm_cut;

  Event::Handle<int> jet_combi_1;
  Event::Handle<int> jet_combi_2;
  Event::Handle<int> jet_combi_3;
  Event::Handle<int> jet_combi_4;

  Event::Handle<TMatrixD> pull_had1;
  Event::Handle<TMatrixD> pull_had2;
  Event::Handle<TMatrixD> pull_hadb;
  Event::Handle<TMatrixD> pull_lepb;
  Event::Handle<TMatrixD> pull_lepton;
  Event::Handle<TMatrixD> pull_neu;

  Event::Handle<float> pull_gtr_had1_et;
  Event::Handle<float> pull_gtr_had1_eta;
  Event::Handle<float> pull_gtr_had1_phi;
  Event::Handle<float> pull_gtr_had2_et;
  Event::Handle<float> pull_gtr_had2_eta;
  Event::Handle<float> pull_gtr_had2_phi;
  Event::Handle<float> pull_gtr_hadb_et;
  Event::Handle<float> pull_gtr_hadb_eta;
  Event::Handle<float> pull_gtr_hadb_phi;
  Event::Handle<float> pull_gtr_lepb_et;
  Event::Handle<float> pull_gtr_lepb_eta;
  Event::Handle<float> pull_gtr_lepb_phi;
  Event::Handle<float> pull_gtr_lepton_et;
  Event::Handle<float> pull_gtr_lepton_eta;
  Event::Handle<float> pull_gtr_lepton_phi;
  Event::Handle<float> pull_gtr_neu_et;
  Event::Handle<float> pull_gtr_neu_eta;
  Event::Handle<float> pull_gtr_neu_phi;

  Event::Handle<TMatrixD> unc_had1;
  Event::Handle<TMatrixD> unc_had2;
  Event::Handle<TMatrixD> unc_hadb;
  Event::Handle<TMatrixD> unc_lepb;
  Event::Handle<TMatrixD> unc_lepton;
  Event::Handle<TMatrixD> unc_neu;

  bool use_ttagging_;
  TopJetId ttag_ID_;
  float ttag_minDR_jet_;

  TMatrixD pull_had1_best, pull_had2_best, pull_hadb_best, pull_lepb_best, pull_lepton_best, pull_neu_best;
  TMatrixD pull_had1_all, pull_had2_all, pull_hadb_all, pull_lepb_all, pull_lepton_all, pull_neu_all;
  TMatrixD pull_gtr_had1_best, pull_gtr_had2_best, pull_gtr_hadb_best, pull_gtr_lepb_best, pull_gtr_lepton_best, pull_gtr_neu_best;

  // // ttbar reconstruction
  std::unique_ptr<uhh2::Selection> chi2_sel_30;
  std::unique_ptr<uhh2::Selection> chi2_sel_20;
  std::unique_ptr<uhh2::Selection> chi2_sel_10;
  std::unique_ptr<uhh2::Selection> chi2_sel_6;
  std::unique_ptr<uhh2::AnalysisModule>  reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule>  ttbar_reco__ttag0;
  std::unique_ptr<uhh2::AnalysisModule>  ttbar_reco__ttag1;
  std::unique_ptr<Chi2Discriminator>     ttbar_chi2__ttag0;
  std::unique_ptr<Chi2DiscriminatorTTAG> ttbar_chi2__ttag1;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  uhh2::Event::Handle<TTbarGen> h_ttbar_gen;
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis> > h_ttbar_hyps;
};


  KinFitTestModule::KinFitTestModule(Context & ctx)
{
  debug = false;   //if debug is true, lots of prints are done

  const std::string& keyword = ctx.get("keyword");
  if (keyword== "pt")           jetParam_ = metParam_ = lepParam_ = KinFitTestSetup::kEMom;
  else if (keyword== "etetaphi")     jetParam_ = metParam_ = lepParam_ = KinFitTestSetup::kEtEtaPhi;
  else {
    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- invalid argument for 'KEYW' key in xml file (must be 'pt' or 'etetaphi'): \""+keyword+"\"";
    throw std::runtime_error(log);
  } 
  const std::string& do_fit_ = ctx.get("do_fit");
  if(do_fit_=="true") do_fit = true;
  if(do_fit_=="false") do_fit = false;

  const std::string& debug_ = ctx.get("debugging");
  if(debug_=="true") debug = true;
  if(debug_=="false") debug = false;

  const bool isMC = (ctx.get("dataset_type") == "MC");
  isMC_ = isMC;

  const std::string& channel = ctx.get("channel", "");
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else if(channel == "both") channel_ = both;
  else {
    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- invalid argument for 'channel' key in xml file (must be 'muon', 'elec' or 'both'): \""+channel+"\"";
    throw std::runtime_error(log);
  }

  const std::string& unc_test_ = ctx.get("unc_test");
  if(unc_test_=="true") unc_test = true;
  if(unc_test_=="false") unc_test = false;

  const std::string& pulls_ = ctx.get("pulls");
  if(pulls_=="true") pulls = true;
  if(pulls_=="false") pulls = false;

  matchable_ = 0;
  const std::string& matchable_tt = ctx.get("matchable");
  if     (matchable_tt == "yes")  matchable_ = 1;
  else if(matchable_tt == "no")   matchable_ = 2;
  else if(matchable_tt == "both") matchable_ = 0;
  else {
    std::string log("TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- ");
    log += "TTbarLJAnalysisLiteModule::TTbarLJAnalysisLiteModule -- invalid argument for 'matchable' key in xml file (must be 'yes', 'no' or 'both'): \""+matchable_tt+"\"";
    throw std::runtime_error(log);
  }

  const std::string ttbar_gen_label ("ttbargen");
  const std::string ttbar_hyps_label ("TTbarReconstruction");
  const std::string ttbar_chi2_label ("Chi2");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  h_ttbar_gen = ctx.get_handle<TTbarGen> (ttbar_gen_label);

  /* t-tagging */
  use_ttagging_ = true;
  const std::string& ttag_wp = ctx.get("ttag_wp");
  ttag_ID_ = TopTagID_SoftDrop(ttag_wp);
  ttag_minDR_jet_ = 1.2;
  ttagevt_sel.reset(new TopTagEventSelection(ttag_ID_, ttag_minDR_jet_));

  //  Reconstuction
  ttbar_reco__ttag0.reset(new HighMassTTbarReconstruction(ctx, NeutrinoReconstruction, ttbar_hyps_label));
  ttbar_chi2__ttag0.reset(new Chi2Discriminator(ctx, ttbar_hyps_label));
  ttbar_chi2__ttag0->set_Mtlep_mean (175.);
  ttbar_chi2__ttag0->set_Mtlep_sigma( 19.);
  ttbar_chi2__ttag0->set_Mthad_mean (177.);
  ttbar_chi2__ttag0->set_Mthad_sigma( 16.);
  reco_primlep.reset(new PrimaryLepton(ctx));

  h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis> >(ttbar_hyps_label);
  chi2_sel_30.reset(new HypothesisDiscriminatorCut(ctx, 0., 30., ttbar_hyps_label, ttbar_chi2_label));
  chi2_sel_20.reset(new HypothesisDiscriminatorCut(ctx, 0., 20., ttbar_hyps_label, ttbar_chi2_label));
  chi2_sel_10.reset(new HypothesisDiscriminatorCut(ctx, 0., 10., ttbar_hyps_label, ttbar_chi2_label));
  chi2_sel_6.reset(new HypothesisDiscriminatorCut(ctx, 0., 6., ttbar_hyps_label, ttbar_chi2_label));

  if(ctx.get("dataset_version") == "TTbar_Mtt0000to0700") genmttbar_sel.reset(new MttbarGenSelection( 0., 700.));
  else                                                    genmttbar_sel.reset(new uhh2::AndSelection(ctx));


  const std::string& btag_wp = ctx.get("btag_wp");

  if     (btag_wp == "CSVL") btag_ID_ = CSVBTag(CSVBTag::WP_LOOSE);
  else if(btag_wp == "CSVM") btag_ID_ = CSVBTag(CSVBTag::WP_MEDIUM);
  else if(btag_wp == "CSVT") btag_ID_ = CSVBTag(CSVBTag::WP_TIGHT);

  if     (btag_wp == "CSVL") b_working_point = CSVBTag::WP_LOOSE;
  else if(btag_wp == "CSVM") b_working_point = CSVBTag::WP_MEDIUM;
  else if(btag_wp == "CSVT") b_working_point = CSVBTag::WP_TIGHT;

  ElectronId eleID;
  float ele_pt(-1.), muon_pt(-1.), jet1_pt(-1.), jet2_pt(-1.),jet3_pt(-1.), jet4_pt(-1.), MET(-1.);

  //KinFit
  myKinFitter.reset(new TKinFitter());

  const std::string& trigger = ctx.get("trigger", "NULL");
  const std::string& trigger2 = ctx.get("trigger2", "NULL");

  //values
  ele_pt = 24.;
  if(trigger=="HLT_IsoMu22_v*"){
    muon_pt = 24.;
  }
  else if(trigger=="HLT_IsoMu24_v*"){
    muon_pt = 26.;
  }
  else if(trigger=="HLT_IsoMu27_v*"){
    muon_pt = 29.;
  }
  //TODO: Change electron ID to 2016 one
  eleID = ElectronID_Spring15_25ns_tight;
  jet1_pt = 25.;
  jet2_pt =  25.;
  jet3_pt =  25.;
  jet4_pt =  25.;
  MET     =  60.;


  //COMMON MODULES
  if(isMC){
    pileup_SF.reset(new MCPileupReweight(ctx)); 
    lumiweight.reset(new MCLumiWeight(ctx));
  }
  if(!isMC) lumi_sel.reset(new LumiSelection(ctx));
  
  //  muon-ID
  const std::string& muonID_SFac    = ctx.get("muonID_SF_file");
  const std::string& muonID_directory    = ctx.get("muonID_SF_directory");
  // //  muon-HLT
  const std::string& muonHLT_SFac   = ctx.get("muonHLT_SF_file");
  const std::string& muonHLT_directory   = ctx.get("muonHLT_SF_directory");
  //  muon-Trk
  const std::string& muonTRK_SFac = ctx.get("muonTRK_SF_file");

  //muon ID scale factors
  muonID_SF.reset(new MCMuonScaleFactor(ctx, muonID_SFac, muonID_directory, 1.12, "ID")); 
  muonHLT_SF.reset(new MCMuonScaleFactor(ctx, muonHLT_SFac, muonHLT_directory, 0.0, "HLT",true));
  muonTRK_SF.reset(new MCMuonTrkScaleFactor(ctx, muonTRK_SFac, 0.0, "TRK"));
  //b-tagging scale factors
  btagSF.reset(new MCBTagScaleFactor(ctx, b_working_point,"jets","central","comb","incl","MCBtagEfficiencies"));//CSV

  //FIXME: Include more MET filters, when Ntuples produced
  /* MET filters */
  PrimaryVertexId pvid=StandardPrimaryVertexId();
  metfilters.emplace_back(new PrimaryVertexCleaner(pvid));
  metfilters_sel.reset(new uhh2::AndSelection(ctx, "metfilters"));
  metfilters_sel->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_sel->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_sel->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  // metfilters_sel->add<TriggerSelection>("1-good-vtx", "Flag_goodVertices");
  metfilters_sel->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");

  //// OBJ CLEANING
  const     MuonId muoSR(AndId<Muon>(MuonIso(), (AndId<Muon>    (PtEtaCut  (muon_pt,    2.1), MuonIDTight()))));
  const ElectronId eleSR(AndId<Electron>(PtEtaSCCut(ele_pt,     2.1), eleID));

  muoSR_cleaner.reset(new     MuonCleaner(muoSR));
  eleSR_cleaner.reset(new ElectronCleaner(eleSR));

    //

  const JetId jetID(JetPFID(JetPFID::WP_LOOSE));

  std::vector<std::string> JEC_AK4, JEC_AK8;
  if(isMC){
    JEC_AK4 = JERFiles::Spring16_25ns_L123_AK4PFchs_MC;
    JEC_AK8 = JERFiles::Spring16_25ns_L123_AK8PFchs_MC;
  }
  else {
    JEC_AK4 = JERFiles::Spring16_25ns_L123_AK4PFchs_DATA;
    JEC_AK8 = JERFiles::Spring16_25ns_L123_AK8PFchs_DATA;
  }

  jet_IDcleaner.reset(new JetCleaner(ctx, jetID));
  jet_corrector.reset(new JetCorrector(ctx, JEC_AK4));
  if(isMC) jetER_smearer.reset(new GenericJetResolutionSmearer(ctx));
  jetlepton_cleaner.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_AK4));
  jet_cleaner2.reset(new JetCleaner(ctx, 25., 2.4));

  //EVENT SELECTION
  jet4_sel.reset(new NJetSelection(4, -1, JetId(PtEtaCut(jet4_pt, 2.4))));
  jet3_sel.reset(new NJetSelection(3, -1, JetId(PtEtaCut(jet3_pt, 2.4))));
  jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut(jet2_pt, 2.4))));
  jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut(jet1_pt, 2.4))));

  met_sel.reset(new METCut  (MET   , uhh2::infinity));

  if(trigger != "NULL") trigger_sel.reset(new TriggerSelection(trigger));
  else                  trigger_sel.reset(new uhh2::AndSelection(ctx));
  if(trigger2 != "NULL" && !isMC) trigger2_sel.reset(new TriggerSelection(trigger2));
  else                  trigger2_sel.reset(new uhh2::AndSelection(ctx));

  sol_met_pz.clear();
  solution.clear();
  solution2.clear();
  input_vectors.clear();
  input_vectors_fit_found.clear();
  input_vectors_all.clear();

  pull_gtr_had1_best.ResizeTo(3,1);
  pull_gtr_had2_best.ResizeTo(3,1);
  pull_gtr_hadb_best.ResizeTo(3,1);
  pull_gtr_lepb_best.ResizeTo(3,1);
  pull_gtr_lepton_best.ResizeTo(3,1);
  pull_gtr_neu_best.ResizeTo(3,1); 


  chi2_prev = 100.;
  F_prev = 100.;
  nrjetcomb_before_dm_cut = 0;
  nrjetcomb_after_dm_cut  = 0;
  nr_suitable_events = 0;


  //Constraints
  //Default: const int maxNrIter=200, const double maxDeltaS=5e-5, const double maxF=1e-4, const double mW=80.4, const double mTop=173.
  mW_ = 80.4;
  decw_w = 2.;
  mTop_ = 173.;
  decw_top = 1.7;
  maxNrIter = 300;
  maxDeltaS = 5e-5;
  maxF=1e-4;
  verbosity=0; // 0=no print of matrices, 1=print matrices, 2=print lorentz vectors, 3= for high met print lorentz vectors

  if(debug){
    std::cout<<"Electron pt cut "<<ele_pt<<std::endl;
    std::cout<<"Muon pt cut "<<muon_pt<<std::endl;
    std::cout<<"Jets pt cuts "<<jet1_pt<<" "<<jet2_pt<<" "<<jet3_pt<<" "<<jet4_pt<<std::endl;
    std::cout<<"MET cut "<<MET<<std::endl;
    std::cout<<"Do the fit "<<do_fit<<std::endl;
    std::cout<<"Mass, width W "<<mW_<<", "<<decw_w<<std::endl;
    std::cout<<"Mass, width top"<<mTop_<<", "<<decw_top<<std::endl;
    std::cout<<"Number max iterations "<<maxNrIter<<std::endl;
    std::cout<<"Delta S max "<<maxDeltaS<<std::endl;
    std::cout<<"F max "<<maxF<<std::endl;
  }

  // //Values: WMass, TopMass   kWHadMass = 1, kWLepMass, kTopHadMass, kTopLepMass, kMET, kEqualTopMasses, kSumPt
  constrList_.push_back(kWHadMass);
  constrList_.push_back(kWLepMass);
  constrList_.push_back(kTopHadMass);
  constrList_.push_back(kTopLepMass);
  //  constrList_.push_back(kMET);
  //constrList_.push_back(kEqualTopMasses);
  // //  constrList_.push_back(kEqualWMasses);
  constrList_.push_back(kSumPt); 




  lumi_h.reset(new LuminosityHists(ctx,"lumi"));
  lumi_after_selection_h.reset(new LuminosityHists(ctx,"lumi_after_selection"));
  h_wgtMC__GEN = ctx.declare_event_output<float>("wgtMC__GEN");  parametrisation = ctx.declare_event_output<int>("parametrisation");
  rec_fit_chi2 = ctx.declare_event_output<float>("rec_fit_chi2");
  fit_probab = ctx.declare_event_output<float>("fit_probab");
  rec_fit_status = ctx.declare_event_output<float>("rec_fit_status");
  rec_fit_fit_found = ctx.declare_event_output<bool>("rec_fit_fit_found");
  rec_fit_is_input = ctx.declare_event_output<bool>("rec_fit_is_input");
  rec_fit_whad_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_whad_v4");
  rec_fit_wlep_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_wlep_v4");
  rec_fit_thad_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_thad_v4");
  rec_fit_tlep_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_tlep_v4");
  rec_fit_lepton_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_lepton_v4");
  rec_fit_neutrino_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_neutrino_v4");
  rec_fit_bhad_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_bhad_v4");
  rec_fit_blep_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_blep_v4");
  rec_fit_jetq1_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_jetq2_v4");
  rec_fit_jetq2_v4 = ctx.declare_event_output<TLorentzVector>("rec_fit_jetq1_v4");
  rec_fit_dR_combined = ctx.declare_event_output<float>("rec_fit_dR_combined"); 

  input_jetq1_v4 = ctx.declare_event_output<TLorentzVector>("input_jetq1_v4");
  input_jetq2_v4 = ctx.declare_event_output<TLorentzVector>("input_jetq2_v4");
  input_jetbh_v4 = ctx.declare_event_output<TLorentzVector>("input_jetbh_v4");
  input_jetbl_v4 = ctx.declare_event_output<TLorentzVector>("input_jetbl_v4");
  input_lepton_v4 = ctx.declare_event_output<TLorentzVector>("input_lepton_v4");
  input_neutrino_v4 = ctx.declare_event_output<TLorentzVector>("input_neutrino_v4");

  jet_permutation_before_dm_cut = ctx.declare_event_output<int>("jet_permutation_before_dm_cut");
  jet_permutation_after_dm_cut = ctx.declare_event_output<int>("jet_permutation_after_dm_cut");

  jet_combi_1 = ctx.declare_event_output<int>("jet_combi_1");
  jet_combi_2 = ctx.declare_event_output<int>("jet_combi_2");
  jet_combi_3 = ctx.declare_event_output<int>("jet_combi_3");
  jet_combi_4 = ctx.declare_event_output<int>("jet_combi_4");

  pull_had1 = ctx.declare_event_output<TMatrixD>("pull_had1");
  pull_had2 = ctx.declare_event_output<TMatrixD>("pull_had2");
  pull_hadb = ctx.declare_event_output<TMatrixD>("pull_hadb");
  pull_lepb = ctx.declare_event_output<TMatrixD>("pull_lepb");
  pull_lepton = ctx.declare_event_output<TMatrixD>("pull_lepton");
  pull_neu = ctx.declare_event_output<TMatrixD>("pull_neu");

  pull_gtr_had1_et = ctx.declare_event_output<float>("pull_gtr_had1_et");
  pull_gtr_had1_eta = ctx.declare_event_output<float>("pull_gtr_had1_eta");
  pull_gtr_had1_phi = ctx.declare_event_output<float>("pull_gtr_had1_phi");
  pull_gtr_had2_et = ctx.declare_event_output<float>("pull_gtr_had2_et");
  pull_gtr_had2_eta = ctx.declare_event_output<float>("pull_gtr_had2_eta");
  pull_gtr_had2_phi = ctx.declare_event_output<float>("pull_gtr_had2_phi");
  pull_gtr_hadb_et = ctx.declare_event_output<float>("pull_gtr_hadb_et");
  pull_gtr_hadb_eta = ctx.declare_event_output<float>("pull_gtr_hadb_eta");
  pull_gtr_hadb_phi = ctx.declare_event_output<float>("pull_gtr_hadb_phi");
  pull_gtr_lepb_et = ctx.declare_event_output<float>("pull_gtr_lepb_et");
  pull_gtr_lepb_eta = ctx.declare_event_output<float>("pull_gtr_lepb_eta");
  pull_gtr_lepb_phi = ctx.declare_event_output<float>("pull_gtr_lepb_phi");
  pull_gtr_lepton_et = ctx.declare_event_output<float>("pull_gtr_lepton_et");
  pull_gtr_lepton_eta = ctx.declare_event_output<float>("pull_gtr_lepton_eta");
  pull_gtr_lepton_phi = ctx.declare_event_output<float>("pull_gtr_lepton_phi");
  pull_gtr_neu_et = ctx.declare_event_output<float>("pull_gtr_neu_et");
  pull_gtr_neu_eta = ctx.declare_event_output<float>("pull_gtr_neu_eta");
  pull_gtr_neu_phi = ctx.declare_event_output<float>("pull_gtr_neu_phi");

  unc_had1 = ctx.declare_event_output<TMatrixD>("unc_had1");
  unc_had2 = ctx.declare_event_output<TMatrixD>("unc_had2");
  unc_hadb = ctx.declare_event_output<TMatrixD>("unc_hadb");
  unc_lepb = ctx.declare_event_output<TMatrixD>("unc_lepb");
  unc_lepton = ctx.declare_event_output<TMatrixD>("unc_lepton");
  unc_neu = ctx.declare_event_output<TMatrixD>("unc_neu");


  h_nocuts.reset(new KinFitTestHists(ctx, "NoCuts"));
  h_trigger.reset(new KinFitTestHists(ctx, "Trigger"));
  h_lepton_size.reset(new KinFitTestHists(ctx, "LeptonSize"));
  h_lepton_jet_sel.reset(new KinFitTestHists(ctx, "LeptonJetSel"));
  h_before_metcut.reset(new KinFitTestHists(ctx, "BeforeMETcut"));
  h_after_metcut.reset(new KinFitTestHists(ctx, "AfterMETcut"));
  h_aftersel.reset(new KinFitTestHists(ctx, "AfterSel"));
  h_mu_nocuts.reset(new MuonHists(ctx, "MuNoCuts"));
  h_mu_trigger.reset(new MuonHists(ctx, "MuTrigger"));
  h_mu_lepton_size.reset(new MuonHists(ctx, "MuLeptonSize"));
  h_mu_lepton_jet_sel.reset(new MuonHists(ctx, "MuLeptonJetSel"));
  h_mu_before_metcut.reset(new MuonHists(ctx, "MuBeforeMETcut"));
  h_mu_after_sel.reset(new MuonHists(ctx, "MuAfterSelection"));
  h_event_nocuts.reset(new EventHists(ctx, "EventNoCuts"));
  h_event_trigger.reset(new EventHists(ctx, "EventTrigger"));
  h_event_lepton_size.reset(new EventHists(ctx, "EventLeptonSize"));
  h_event_lepton_jet_sel.reset(new EventHists(ctx, "EventLeptonJetSel"));
  h_event_before_metcut.reset(new EventHists(ctx, "EventBeforeMETcut"));
  h_event_after_sel.reset(new EventHists(ctx, "EventAfterSel"));  
  h_jet_nocuts.reset(new JetHists(ctx, "JetNoCuts"));
  h_jet_trigger.reset(new JetHists(ctx, "Trigger"));
  h_jet_lepton_size.reset(new JetHists(ctx, "LeptonSize"));
  h_jet_lepton_jet_sel.reset(new JetHists(ctx, "LeptonJetSel"));
  h_jet_before_metcut.reset(new JetHists(ctx, "JetBeforeMETcut"));
  h_jet_after_sel.reset(new JetHists(ctx, "JetAfterSel"));  
  //  h_btag_hists.reset(new BTagMCEfficiencyHists(ctx, "BTagEff", b_working_point, "jets"));  
  if(do_fit){
    h_fitbest.reset(new KinFitTestHistsFitted(ctx, "FitBest"));
    h_fitbest_chi2_cut.reset(new KinFitTestHistsFitted(ctx, "FitBestChi2cut"));
  }
  h_before_fit.reset(new KinFitTestFITHists(ctx, "Before_fit_ttbar", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));

    h_reco.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_ttbar_m.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_matchable_ttbar", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_ttbar_m_chi2_30.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_matchable_ttbar_chi2_30", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_ttbar_m_chi2_6.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_matchable_ttbar_chi2_6", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_matchable.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_matchable", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_not_matchable.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_not_matchable", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_matchable_chi2_30.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_matchable_chi2_30", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_matchable_chi2_6.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_matchable_chi2_6", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_not_matchable_chi2_6.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_not_matchable_chi2_6", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_30.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_30", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_0btag.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_0btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_1btag.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_1btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_2btag.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_2btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_20.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_20", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_10.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_10", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_reco_chi2_6.reset(new KinFitTestRECOHists(ctx, "Reco_ttbar_chi2_6", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));

  if(do_fit){
    h_fitted.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_allcomb.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_allcombinations", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_allcomb_fit_converged.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_allcombinations_fit_converges", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_chi2_cut.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_chi2_cut", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_0btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_0b", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_1btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_1b", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_2btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_2b", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_dR_cut.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_dR_cut", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_matchable.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_matchable", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_matchable_0btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_matchable_0btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_matchable_1btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_matchable_1btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_matchable_2btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_matchable_2btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_not_matchable.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_not_matchable", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_chi2_matchable.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_chi2_matchable", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_chi2_not_matchable.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_chi2_not_matchable", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_b_jet_matching.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_b_jet_matching", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_b_jet_matching_0btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_b_jet_matching_0btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_b_jet_matching_1btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_b_jet_matching_1btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_b_jet_matching_2btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_b_jet_matching_2btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_match_ttbar_decay.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_match_ttbar_decay", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_match_ttbar_decay_b_jet_matching.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_match_ttbar_decay_b_jet_matching", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_match_ttbar_decay_b_jet_matching_0btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_match_ttbar_decay_b_jet_matching_0btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_match_ttbar_decay_b_jet_matching_1btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_match_ttbar_decay_b_jet_matching_1btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_match_ttbar_decay_b_jet_matching_2btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_match_ttbar_decay_b_jet_matching_2btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_correct_match.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_correct_match", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_correct_match_chi2_cut.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_correct_match_chi2_cut", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_correct_match_0btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_correct_match_0btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_correct_match_1btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_correct_match_1btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    h_fitted_correct_match_2btag.reset(new KinFitTestFITHists(ctx, "Fitted_ttbar_correct_match_2btag", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));

    if(unc_test)    h_muon_unc.reset(new KinFitTestHistsMuonUnc(ctx, "Muon_Uncertainty", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));
    if(pulls)    h_pulls.reset(new KinFitTestHistsPulls(ctx, "Pull_hists", ttbar_gen_label, ttbar_hyps_label, ttbar_chi2_label));

  }

}



bool KinFitTestModule::process(Event & event) {
  //  if(!(event.event==14697293))  return false;
  if(verbosity==1 || verbosity==2 || verbosity==3 || debug==true)    cout << " %%%%% KinFitTestModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  float w_GEN(1.);
  if(!event.isRealData){
    // /* GEN M-ttbar selection */
    ttgenprod->process(event);
    w_GEN = event.weight;
    //Selects only ttbar events form 0 to 700 GeV
    // if(!genmttbar_sel->passes(event)) return false;
  }
  event.set(h_wgtMC__GEN , w_GEN);

  const TTbarGen* ttgen_test(0);
  if(isMC_){
    const auto& ttbargen_test = event.get(h_ttbar_gen);
    ttgen_test = &ttbargen_test;
   }


  /* MET filters */
  if(!metfilters_sel->passes(event)) return false;
  for(auto & m : metfilters){
    m->process(event);
  }
  if(debug)    std::cout<<"met filters passed"<<std::endl;
  /* CMS-certified luminosity sections */
  if(event.isRealData && !lumi_sel->passes(event)) return false;
  /* Lumi weight and pilup SF */
  if(!event.isRealData){ 
    pileup_SF->process(event);
    lumiweight->process(event);
  }
  btagSF->process(event);


  h_nocuts->fill(event);
  h_mu_nocuts->fill(event);
  h_event_nocuts->fill(event);
  h_jet_nocuts->fill(event);

 
  //KinFit
  //TODO: 0. Add all objects as in bool TTbarLJAnalysisLiteModule::process
  
  //// LEPTON selection
  muoSR_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);
  
  eleSR_cleaner->process(event);
  sort_by_pt<Electron>(*event.electrons);

  jet_IDcleaner->process(event);
  jet_corrector->process(event);
  if(jetER_smearer.get()) jetER_smearer->process(event);
  jetlepton_cleaner->process(event);

  jet_cleaner2->process(event);
  sort_by_pt<Jet>(*event.jets);

  h_lepton_jet_sel->fill(event);
  h_mu_lepton_jet_sel->fill(event);
  h_event_lepton_jet_sel->fill(event);
  h_jet_lepton_jet_sel->fill(event);

  //HLT selection
  if(channel_ == muon || channel_ == elec){
    const bool pass_trigger = trigger_sel->passes(event);
     bool pass_trigger2 = true;
    pass_trigger2 = trigger2_sel->passes(event);
    if(!pass_trigger || !pass_trigger2) return false;
  }
   
  muonID_SF->process(event);
  muonHLT_SF->process(event);
  muonTRK_SF->process(event);

  h_trigger->fill(event);
  h_mu_trigger->fill(event);
  h_event_trigger->fill(event);
  h_jet_trigger->fill(event);

  lumi_h->fill(event);

 
  if(channel_ == elec){
    const bool pass_lep1_ = (event.electrons->size() == 1);
    if(!pass_lep1_ || event.muons->size()>0) return false;
    leptonType = KinFitTestSetup::kElectron;
  }
  else if(channel_ == muon){
    const bool pass_lep1_ = (event.muons->size() == 1);
    if(!pass_lep1_ || event.electrons->size()>0) return false;
    leptonType = KinFitTestSetup::kMuon;
  }
  else if(channel_ == both){
    const bool pass_lep1_ = ((event.muons->size() == 1) || (event.electrons->size() == 1));
    if((event.muons->size() == 1) && (event.electrons->size() == 1)) return false;
    if(!pass_lep1_) return false;
    if(event.muons->size() == 1) leptonType = KinFitTestSetup::kMuon;
    if(event.electrons->size() == 1) leptonType = KinFitTestSetup::kElectron;
  }
  
  h_lepton_size->fill(event);
  h_mu_lepton_size->fill(event);
  h_event_lepton_size->fill(event);
  h_jet_lepton_size->fill(event);
 
  //// JET selection
  /*  4th AK4 jet selection */
  const bool pass_jet4 = jet4_sel->passes(event);
  if(!pass_jet4) return false;
  /*  3rd AK4 jet selection */
  const bool pass_jet3 = jet3_sel->passes(event);
  if(!pass_jet3) return false;
  /*  2nd AK4 jet selection */
  const bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;
  /*  1st AK4 jet selection */
  const bool pass_jet1 = jet1_sel->passes(event);
  if(!pass_jet1) return false;


  h_before_metcut->fill(event);
  h_mu_before_metcut->fill(event);
  h_event_before_metcut->fill(event);
  h_jet_before_metcut->fill(event);


  // /*  MET selection */
  const bool pass_met = met_sel->passes(event);
  if(!pass_met) return false;
  
  if(debug){
    std::cout<<"Numbers of electrons "<<event.electrons->size()<<std::endl;
    std::cout<<"Numbers of muons "<<event.muons->size()<<std::endl;
    std::cout<<"Lepton type "<<leptonType<<std::endl;
    std::cout<<"Numbers of jets "<<event.jets->size()<<std::endl;
  }
  h_after_metcut->fill(event);
  h_mu_after_sel->fill(event);
  h_event_after_sel->fill(event);
  h_jet_after_sel->fill(event);
  //  h_btag_hists->fill(event);
  h_aftersel->fill(event);
  lumi_after_selection_h->fill(event);

  //Fill kinematic reconstruction

  /* TOPTAG-EVENT boolean */
  reco_primlep->process(event);
  ttbar_reco__ttag0->process(event);
  ttbar_chi2__ttag0->process(event);
  
  int jetbtagN(0);
  std::vector<ReconstructionHypothesis>& ttbar_hyps = event.get(h_ttbar_hyps);
  const ReconstructionHypothesis* rec_ttbar = get_best_hypothesis(ttbar_hyps, "Chi2");
  const bool pass_chi2_30 = chi2_sel_30->passes(event);
  const bool pass_chi2_20 = chi2_sel_20->passes(event);
  const bool pass_chi2_10 = chi2_sel_10->passes(event);
  const bool pass_chi2_6 = chi2_sel_6->passes(event);
  for(const auto& j : *event.jets){
    if(btag_ID_(j, event)) ++jetbtagN;
  }
  //  const float rec_chi2 = rec_ttbar->discriminator("Chi2");
  h_reco->fill(event);
  if(pass_chi2_30){
    h_reco_chi2_30->fill(event);
    if(jetbtagN==0)     h_reco_chi2_0btag->fill(event);
    if(jetbtagN==1)     h_reco_chi2_1btag->fill(event);
    if(jetbtagN>=2)     h_reco_chi2_2btag->fill(event);
  }
  if(pass_chi2_20)    h_reco_chi2_20->fill(event);
  if(pass_chi2_10)    h_reco_chi2_10->fill(event);
  if(pass_chi2_6)     h_reco_chi2_6->fill(event);

  if(!do_fit) return false;

  bool matchable_ttbar(false), match_ttbar_decay_reco(false);
  const TTbarGen* ttgen(0);
  if(isMC_){
    const auto& ttbargen = event.get(h_ttbar_gen);
    ttgen = &ttbargen;
    if(ttgen->IsSemiLeptonicDecay()){
      const LorentzVector gen_jetq1_v4 = ttgen->Q1().v4();
      const LorentzVector gen_jetq2_v4 = ttgen->Q2().v4();
      const LorentzVector gen_bhad_v4 = ttgen->BHad().v4();
      const LorentzVector gen_blep_v4 = ttgen->BLep().v4();
      const LorentzVector gen_lepton_v4 = ttgen->ChargedLepton().v4();
      //check if event is "matchable": pt and eta of objects inside the selection criteria
      bool pt_jets(false), eta_jets(false), pt_lepton(false), eta_lepton(false);
      if(gen_jetq1_v4.Pt()>25. && gen_jetq2_v4.Pt()>25. && gen_bhad_v4.Pt()>25. && gen_blep_v4.Pt()>25.) pt_jets = true;
      if(fabs(gen_jetq1_v4.Eta())<2.4 && fabs(gen_jetq2_v4.Eta())<2.4 && fabs(gen_bhad_v4.Eta())<2.4 && fabs(gen_blep_v4.Eta())<2.4 ) eta_jets = true;
      
      switch(leptonType){
      case KinFitTestSetup::kElectron :       if(gen_lepton_v4.Pt()>24.) pt_lepton = true;	break;
      case KinFitTestSetup::kMuon :           if(gen_lepton_v4.Pt()>24.) pt_lepton = true; break;
      case KinFitTestSetup::kMet : break;
      case KinFitTestSetup::kBJet : break;
      case KinFitTestSetup::kUdscJet : break;
      }
      if(fabs(gen_lepton_v4.Eta())<2.1) eta_lepton = true;
      
      if(pt_jets && eta_jets && pt_lepton && eta_lepton) matchable_ttbar = true;

      if(matchable_==1 && matchable_ttbar==false) return false;
      else if(matchable_==2 && matchable_ttbar==true) return false;

      std::vector<Jet> toplep_jets_reco;
      for(unsigned int jet_lep=0; jet_lep<rec_ttbar->toplep_jets().size();jet_lep++){
	toplep_jets_reco.push_back(rec_ttbar->toplep_jets().at(jet_lep));
      }
      std::vector<Jet> tophad_jets_reco;
      for(unsigned int jet_had=0; jet_had<rec_ttbar->tophad_jets().size();jet_had++){
	tophad_jets_reco.push_back(rec_ttbar->tophad_jets().at(jet_had));
      }

      LorentzVector lepton_reco = rec_ttbar->lepton().v4();
      float dR_lep = h_fitted_dR_cut->deltaR_fit(lepton_reco, gen_lepton_v4);

      std::vector<LorentzVector> gen_jets;
      gen_jets.clear();
      gen_jets.push_back(gen_jetq1_v4);
      gen_jets.push_back(gen_jetq2_v4);
      gen_jets.push_back(gen_bhad_v4);
      //     gen_jets.push_back(gen_blep_v4);
      
      int nr_matched_reco_jets(0);
      
      for(unsigned int gj=0; gj<gen_jets.size(); gj++){
        int matches(0);  
        for(unsigned int fj=0; fj<toplep_jets_reco.size(); fj++){
          float dR = h_fitted_dR_cut->deltaR_fit_new(toplep_jets_reco[fj].eta(), gen_blep_v4.eta(),toplep_jets_reco[fj].phi(), gen_blep_v4.phi());
          if(dR<0.3){
            matches++;
          }}
        for(unsigned int fjh=0; fjh<tophad_jets_reco.size(); fjh++){
          float dR = h_fitted_dR_cut->deltaR_fit_new(tophad_jets_reco[fjh].eta(), gen_jets[gj].eta(),tophad_jets_reco[fjh].phi(), gen_jets[gj].phi());

          if(dR<0.3){
            matches++;
         }
        }
        if(matches>=1) nr_matched_reco_jets++;
      }
      
      if(nr_matched_reco_jets>=4 && dR_lep<0.3) match_ttbar_decay_reco = true;
    }
  }

  if(matchable_ttbar){
    h_reco_matchable->fill(event);
    if(pass_chi2_30)    h_reco_matchable_chi2_30->fill(event);
    if(pass_chi2_6)    h_reco_matchable_chi2_6->fill(event);
    if(match_ttbar_decay_reco){
      h_reco_ttbar_m->fill(event);
      if(pass_chi2_30) h_reco_ttbar_m_chi2_30->fill(event);
      if(pass_chi2_6) h_reco_ttbar_m_chi2_6->fill(event);
    }
  }else{
    h_reco_not_matchable->fill(event);
    if(pass_chi2_6)    h_reco_not_matchable_chi2_6->fill(event);
  }


  //TODO: 1. Fill kinfit as in void TtSemiLepKinFitter::setupFitter() 
  myKinFitter->reset();
  
  setupJets();
  setupLeptons();
  setupConstraints();
 
  myKinFitter->addMeasParticle(jethad1_);
  myKinFitter->addMeasParticle(jethad2_);
  myKinFitter->addMeasParticle(jethadb_);
  myKinFitter->addMeasParticle(jetlepb_);
  myKinFitter->addMeasParticle(lepton_);
  myKinFitter->addMeasParticle(met_);
 

  // add constraints
  for(unsigned int j=0; j<constrList_.size(); j++){
    if(constrList_[j]!=kSumPt){
      myKinFitter->addConstraint(massConstr_[constrList_[j]]);
      if(debug) std::cout<<"Constraint list "<<constrList_[j]<<std::endl;
    }
  }

  if(constrainSumPt_) {
    myKinFitter->addConstraint(sumPxConstr_);
    myKinFitter->addConstraint(sumPyConstr_);
  }

  myKinFitter->setMaxDeltaS(maxDeltaS);
  myKinFitter->setMaxF(maxF);
  myKinFitter->setMaxNbIter(maxNrIter);
  myKinFitter->setVerbosity(verbosity);


  
  //2. Call fit as in int TtSemiLepKinFitter::fit(const TLorentzVector& p4HadP...) 
  //NB: TLorentzVector from the event for the particles as in bool TTbarLJAnalysisLiteModule::process(uhh2::Event& event)
  
  
  //b-tag
  //  unsigned int bjet;
  nrbtagedjet.clear();
  for (unsigned int bjet=0;bjet<event.jets->size();bjet++){
    if(btag_ID_(event.jets->at(bjet), event)){
      nrbtagedjet.push_back(bjet);
    }
  }

  if(debug) std::cout<<"Number of b tagged jets "<<nrbtagedjet.size()<<std::endl;
  
  chi2_lowest = 10000.;
  fit_prob = 0.;
  fit_found = false;
  fitprob_curr = -1.;
  fitprob_prev = -1.;
  fitprob_best = 0.;
  chi2_curr = -1.;
  F_curr = -1.;
  dR_tlep = infinity;
  dR_thad = infinity;
  dR_wlep = infinity;
  dR_jetq1 = infinity;
  dR_jetq2 = infinity;
  dR_jetbh = infinity;
  dR_jetbl = infinity;
  dR_lepton = infinity;
  dR_neutrino = infinity;
  dR_combined = infinity;
  fit_status = -100.;
  fit_status_each = -100.;
  jet_combi_all.clear();
  jet_combi_best.clear();
  input_vectors.clear();
  dm_w_lep = -100.;    
  dm_w_had = -100.;
  dm_top_lep = -100.;    
  dm_top_had = -100.;
  nrjetcomb_before_dm_cut = 0;
  nrjetcomb_after_dm_cut = 0;
  pull_gtr_had1_best.Zero();
  pull_gtr_had2_best.Zero();
  pull_gtr_hadb_best.Zero();
  pull_gtr_lepb_best.Zero();
  pull_gtr_lepton_best.Zero();
  pull_gtr_neu_best.Zero(); 

  if(debug){
    std::cout<<"number of jets "<<event.jets->size()<<std::endl;
    for(unsigned int i=0; i<event.jets->size(); i++){
      std::cout<<"pt of jet "<<i<<" "<<event.jets->at(i).v4().Pt()<<std::endl;
    }
  }
  //  more than two btaged jets
  if(event.jets->size()>=4 && nrbtagedjet.size()>=2){
    // protect against reading beyond array boundaries
    int nrCombJets = event.jets->size();
    nr_suitable_events++;
    for (const auto& bl : nrbtagedjet) {
  	for (const auto& bh : nrbtagedjet) {
  	  for (int h1=0; h1<nrCombJets; h1++) {
  	    if (!(bl==bh) && !(bh==h1 || bl==h1)) {
  	      for (int h2=h1+1; h2<nrCombJets; h2++) {
  		if (!(h2==h1 || h2==bl || h2==bh)) {
  		  nrjetcomb_before_dm_cut++;
		  bool dm_obj = dm_check(event, bl, bh, h1, h2);
		  if(dm_obj)    
		    fit(event, bl, bh, h1, h2, leptonType, jetParam_, lepParam_, metParam_);
		  else{
		    if(debug)   std::cout<<"no fit"<<std::endl;
		  }
		  nrjetcomb_after_dm_cut++;
  		}
  	      }
  	    }
  	  }
  	}
    }
  }


  //one btaged jet
  else if(event.jets->size()>=4 && nrbtagedjet.size()==1){
    // protect against reading beyond array boundaries
    int nrCombJets = event.jets->size();
    nr_suitable_events++;
    //tagged b jet is on leptonic side
    int bl = nrbtagedjet[0];
    for (int bh=0;bh<nrCombJets; bh++) {
  	for (int h1=0; h1<nrCombJets; h1++) {
  	  if (!(bl==bh) && !(bh==h1 || bl==h1)) {
  	    for (int h2=h1+1; h2<nrCombJets; h2++) {
  	      if (!(h2==h1 || h2==bl || h2==bh)) {
		nrjetcomb_before_dm_cut++;
		bool dm_obj = dm_check(event, bl, bh, h1, h2);
		if(dm_obj)    
		  fit(event, bl, bh, h1, h2, leptonType, jetParam_, lepParam_, metParam_);
		else{
		  if(debug)   std::cout<<"no fit"<<std::endl;
		}
		nrjetcomb_after_dm_cut++;
  	      }
  	    }
  	  }
  	}
    }
  
    //tagged b jet is on hadronic side
  for (int bl=0; bl<nrCombJets; bl++) {
    int bh = nrbtagedjet[0];
    for (int h1=0; h1<nrCombJets; h1++) {
  	if (!(bl==bh) && !(bh==h1 || bl==h1)) {
  	  for (int h2=h1+1; h2<nrCombJets; h2++) {
  	    if (!(h2==h1 || h2==bl || h2==bh)) {
	      nrjetcomb_before_dm_cut++;
	      bool dm_obj = dm_check(event, bl, bh, h1, h2);
	      if(dm_obj)    
		fit(event, bl, bh, h1, h2, leptonType, jetParam_, lepParam_, metParam_);
	      else{
		if(debug)   std::cout<<"no fit"<<std::endl;
	      }
	      nrjetcomb_after_dm_cut++;
  	    }
  	  }
  	}
    }
  }
  }


  // const TTbarGen* ttgen_input(0);
  // if(isMC_){
  //   const auto& ttbargen_input = event.get(h_ttbar_gen);
  //   ttgen_input = &ttbargen_input;

 
  //no btaged jets
  else if(event.jets->size()>=4 && nrbtagedjet.size()==0){
    int nrCombJets = event.jets->size();
    // int nrCombJets = 4;
    nr_suitable_events++;

    for (int bl=0; bl<nrCombJets; bl++) {
      for (int bh=0;bh<nrCombJets; bh++) {
    	for (int h1=0; h1<nrCombJets; h1++) {
    	  if (!(bl==bh) && !(bh==h1 || bl==h1)) {
    	    for (int h2=h1+1; h2<nrCombJets; h2++) {
    	      if (!(h2==h1 || h2==bl || h2==bh)) {
		nrjetcomb_before_dm_cut++;
		bool dm_obj = dm_check(event, bl, bh, h1, h2);
		if(dm_obj)    
		  fit(event, bl, bh, h1, h2, leptonType, jetParam_, lepParam_, metParam_);
		else{
		  if(debug)   std::cout<<"no fit"<<std::endl;
		}
		nrjetcomb_after_dm_cut++;
	      }
	    }
	  }
	}
      }
    }
  }


  if(!fit_found) return false;
  else{
    if(debug){
      std::cout<<"**** Best fit result *****"<<std::endl;
      std::cout<<"Lowest chi2 value "<<chi2_lowest<<std::endl;
      std::cout<<"Highest fit probability "<<fit_prob<<std::endl;
      std::cout<<"Status of fit "<<fit_status<<std::endl;
      std::cout<<"SOLUTION:"<<std::endl;
      for(unsigned int i=0; i<solution.size(); i++){
	std::cout<<"Vector "<<i<<std::endl;
	solution[i].Print();
	std::cout<<"Et "<<solution[i].Et()<<std::endl;
	std::cout<<"Pt "<<solution[i].Pt()<<std::endl;
	std::cout<<"Mass "<<solution[i].M()<<std::endl;
      }
      std::cout<<"INPUT:"<<std::endl;
      for(unsigned int i=0; i<input_vectors_fit_found.size(); i++){
	std::cout<<"Vector "<<i<<std::endl;
	input_vectors_fit_found[i].Print();
	std::cout<<"Et "<<input_vectors_fit_found[i].Et()<<std::endl;
	std::cout<<"Pt "<<input_vectors_fit_found[i].Pt()<<std::endl;
	std::cout<<"Mass "<<input_vectors_fit_found[i].M()<<std::endl;
      }
    }

    setBestFit(event, solution, input_vectors_fit_found, leptonType, chi2_lowest, fit_prob, fit_status, fit_found, false, jet_combi_best);

    event.set(pull_had1, pull_had1_best);
    event.set(pull_had2, pull_had2_best);
    event.set(pull_hadb, pull_hadb_best);
    event.set(pull_lepb, pull_lepb_best);
    event.set(pull_lepton, pull_lepton_best);
    event.set(pull_neu, pull_neu_best);
    if(isMC_){
      if(ttgen_test->IsSemiLeptonicDecay()){
	event.set(pull_gtr_had1_et, pull_gtr_had1_best(0,0));
	event.set(pull_gtr_had1_eta, pull_gtr_had1_best(1,0));
	event.set(pull_gtr_had1_phi, pull_gtr_had1_best(2,0));
	event.set(pull_gtr_had2_et, pull_gtr_had2_best(0,0));
	event.set(pull_gtr_had2_eta, pull_gtr_had2_best(1,0));
	event.set(pull_gtr_had2_phi, pull_gtr_had2_best(2,0));
	event.set(pull_gtr_hadb_et, pull_gtr_hadb_best(0,0));
	event.set(pull_gtr_hadb_eta, pull_gtr_hadb_best(1,0));
	event.set(pull_gtr_hadb_phi, pull_gtr_hadb_best(2,0));
	event.set(pull_gtr_lepb_et, pull_gtr_lepb_best(0,0));
	event.set(pull_gtr_lepb_eta, pull_gtr_lepb_best(1,0));
	event.set(pull_gtr_lepb_phi, pull_gtr_lepb_best(2,0));
	event.set(pull_gtr_lepton_et, pull_gtr_lepton_best(0,0));
	event.set(pull_gtr_lepton_eta, pull_gtr_lepton_best(1,0));
	event.set(pull_gtr_lepton_phi, pull_gtr_lepton_best(2,0));
	event.set(pull_gtr_neu_et, pull_gtr_neu_best(0,0));
	event.set(pull_gtr_neu_eta, pull_gtr_neu_best(1,0));
	event.set(pull_gtr_neu_phi, pull_gtr_neu_best(2,0));
      }
    }
    else{
       event.set(pull_gtr_had1_et, -100.);
       event.set(pull_gtr_had1_eta, -100.);
       event.set(pull_gtr_had1_phi, -100.);
       event.set(pull_gtr_had2_et, -100.);
       event.set(pull_gtr_had2_eta, -100.);
       event.set(pull_gtr_had2_phi, -100.);
       event.set(pull_gtr_hadb_et, -100.);
       event.set(pull_gtr_hadb_eta, -100.);
       event.set(pull_gtr_hadb_phi, -100.);
       event.set(pull_gtr_lepb_et, -100.);
       event.set(pull_gtr_lepb_eta, -100.);
       event.set(pull_gtr_lepb_phi, -100.);
       event.set(pull_gtr_lepton_et, -100.);
       event.set(pull_gtr_lepton_eta, -100.);
       event.set(pull_gtr_lepton_phi, -100.);
       event.set(pull_gtr_neu_et, -100.);
       event.set(pull_gtr_neu_eta, -100.);
       event.set(pull_gtr_neu_phi, -100.);

    }

    event.set(jet_permutation_before_dm_cut, nrjetcomb_before_dm_cut);
    event.set(jet_permutation_after_dm_cut, nrjetcomb_after_dm_cut);

    h_fitted->fill(event); 
    if(unc_test)    h_muon_unc->fill(event);
    if(pulls && correct_match)    h_muon_unc->fill(event);
    h_fitbest->fill(event);

    if(debug){
      std::cout<<"!!!!!!!!! plot"<<std::endl;
      std::cout<<"wlep  "<<dm_w_lep<<std::endl;
      std::cout<<"whad  "<<dm_w_had<<std::endl;
      std::cout<<"top lep  "<<dm_top_lep<<std::endl;
      std::cout<<"top had  "<<dm_top_had<<std::endl;
    }

    if(nrbtagedjet.size()==0) h_fitted_0btag->fill(event);
    if(nrbtagedjet.size()==1) h_fitted_1btag->fill(event);
    if(nrbtagedjet.size()>=2) h_fitted_2btag->fill(event);
    
    if(matchable){
      h_fitted_matchable->fill(event);
      if(nrbtagedjet.size()==0)	h_fitted_matchable_0btag->fill(event);
      if(nrbtagedjet.size()==1)	h_fitted_matchable_1btag->fill(event);
      if(nrbtagedjet.size()>=2)	h_fitted_matchable_2btag->fill(event);
    }
    else h_fitted_not_matchable->fill(event);
    if(matchable && b_jet_matching){
      h_fitted_b_jet_matching->fill(event);
      if(nrbtagedjet.size()==0)	h_fitted_b_jet_matching_0btag->fill(event);
      if(nrbtagedjet.size()==1)	h_fitted_b_jet_matching_1btag->fill(event);
      if(nrbtagedjet.size()>=2)	h_fitted_b_jet_matching_2btag->fill(event);
    }

    if(fit_prob>0.2){
      h_fitted_chi2_cut->fill(event);
      h_fitbest_chi2_cut->fill(event);
      if(matchable)      h_fitted_chi2_matchable->fill(event);
      else               h_fitted_chi2_not_matchable->fill(event);
      if(matchable && match_ttbar_decay){
	h_fitted_match_ttbar_decay->fill(event);
	if(b_jet_matching){
	  h_fitted_match_ttbar_decay_b_jet_matching->fill(event);
	  if(nrbtagedjet.size()==0)	h_fitted_match_ttbar_decay_b_jet_matching_0btag->fill(event);
	  if(nrbtagedjet.size()==1)	h_fitted_match_ttbar_decay_b_jet_matching_1btag->fill(event);
	  if(nrbtagedjet.size()>=2)	h_fitted_match_ttbar_decay_b_jet_matching_2btag->fill(event);
	}
    }
      if(matchable && correct_match){
	h_fitted_correct_match->fill(event);
	if(nrbtagedjet.size()==0)	h_fitted_correct_match_0btag->fill(event);
	if(nrbtagedjet.size()==1)	h_fitted_correct_match_1btag->fill(event);
	if(nrbtagedjet.size()>=2)	h_fitted_correct_match_2btag->fill(event);
	if(fit_prob>0.2) h_fitted_correct_match_chi2_cut->fill(event);
      }
    }
  }


  if(debug){
    pull_gtr_had1_best.Print();
    pull_gtr_had2_best.Print();
    pull_gtr_hadb_best.Print();
    pull_gtr_lepb_best.Print();
    pull_gtr_lepton_best.Print();
    pull_gtr_neu_best.Print();
  }

  //3. Compare fit results with GEN info
  
  delete jethad1_; 
  delete jethad2_;
  delete jethadb_;
  delete jetlepb_;
  delete lepton_;
  delete met_;
  if(constrainSumPt_){
    delete sumPxConstr_;
    delete sumPyConstr_;
  }
  for(unsigned int j=0; j<constrList_.size(); j++){
    delete massConstr_[constrList_[j]];
  }

  pull_had1_best.Zero();
  pull_had2_best.Zero();
  pull_hadb_best.Zero();
  pull_lepb_best.Zero();
  pull_lepton_best.Zero();
  pull_neu_best.Zero();
  pull_had1_all.Zero();
  pull_had2_all.Zero();
  pull_hadb_all.Zero();
  pull_lepb_all.Zero();
  pull_lepton_all.Zero();
  pull_neu_all.Zero();
  pull_gtr_had1_best.Zero();
  pull_gtr_had2_best.Zero();
  pull_gtr_hadb_best.Zero();
  pull_gtr_lepb_best.Zero();
  pull_gtr_lepton_best.Zero();
  pull_gtr_neu_best.Zero();
  // pull_gtr_had1_all.Zero();
  // pull_gtr_had2_all.Zero();
  // pull_gtr_hadb_all.Zero();
  // pull_gtr_lepb_all.Zero();
  // pull_gtr_lepton_all.Zero();
  // pull_gtr_neu_all.Zero();
  
  massConstr_.clear(); 
  //  if(!fit_found) return false;

  return true;
}

  void KinFitTestModule::fit(Event & event, int bl, int bh, int h1, int h2, KinFitTestSetup::ObjectType leptonType, KinFitTestSetup::Param jetParam_, KinFitTestSetup::Param lepParam_, KinFitTestSetup::Param metParam_){

    dm_w_lep_all = -100.;    
    dm_w_had_all = -100.;
    dm_top_lep_all = -100.;    
    dm_top_had_all = -100.;
  
    p4Had1 = TLorentzVector(event.jets->at(h1).v4().Px(),event.jets->at(h1).v4().Py(),event.jets->at(h1).v4().Pz(),event.jets->at(h1).v4().E());
    p4Had2 = TLorentzVector(event.jets->at(h2).v4().Px(),event.jets->at(h2).v4().Py(),event.jets->at(h2).v4().Pz(),event.jets->at(h2).v4().E());
    p4HadB = TLorentzVector(event.jets->at(bh).v4().Px(),event.jets->at(bh).v4().Py(),event.jets->at(bh).v4().Pz(),event.jets->at(bh).v4().E());
    p4LepB = TLorentzVector(event.jets->at(bl).v4().Px(),event.jets->at(bl).v4().Py(),event.jets->at(bl).v4().Pz(),event.jets->at(bl).v4().E());
   
    switch(leptonType){
    case KinFitTestSetup::kElectron :       p4Lep  = TLorentzVector(event.electrons->at(0).v4().Px(),event.electrons->at(0).v4().Py(),event.electrons->at(0).v4().Pz(),event.electrons->at(0).v4().E()); break;
    case KinFitTestSetup::kMuon :           p4Lep  = TLorentzVector(event.muons->at(0).v4().Px(),event.muons->at(0).v4().Py(),event.muons->at(0).v4().Pz(),event.muons->at(0).v4().E()); break;
    case KinFitTestSetup::kMet : break;
    case KinFitTestSetup::kBJet : break;
    case KinFitTestSetup::kUdscJet : break;
    }

    p4Neutrino  = TLorentzVector(event.met->v4().Px(),event.met->v4().Py(),0,sqrt(pow(event.met->v4().M(),2) + pow(event.met->v4().Pt(), 2)));  
 
    // set covariance matrices of the objects to be fitted
    KinFitTestSetup setup;
    TMatrixD covHad1 = setup.setupMatrix(p4Had1, KinFitTestSetup::kUdscJet, jetParam_);
    TMatrixD covHad2 = setup.setupMatrix(p4Had2, KinFitTestSetup::kUdscJet, jetParam_);
    TMatrixD covHadB = setup.setupMatrix(p4HadB, KinFitTestSetup::kBJet, jetParam_);
    TMatrixD covLepB = setup.setupMatrix(p4LepB, KinFitTestSetup::kBJet, jetParam_);
    TMatrixD covLep  = setup.setupMatrix(p4Lep,  leptonType, lepParam_);
    TMatrixD covMET  = setup.setupMatrix(p4Neutrino, KinFitTestSetup::kMet, metParam_);

    event.set(unc_had1, covHad1);
    event.set(unc_had2, covHad2);
    event.set(unc_hadb, covHadB);
    event.set(unc_lepb, covLepB);
    event.set(unc_lepton, covLep);
    event.set(unc_neu, covMET);

    // set the kinematics of the objects to be fitted
    jethad1_->setIni4Vec( &p4Had1 );
    jethad2_->setIni4Vec( &p4Had2 );
    jethadb_->setIni4Vec( &p4HadB );
    jetlepb_->setIni4Vec( &p4LepB );
    lepton_ ->setIni4Vec( &p4Lep  );
    met_    ->setIni4Vec( &p4Neutrino  );


    jethad1_->setCovMatrix( &covHad1 );
    jethad2_->setCovMatrix( &covHad2 );
    jethadb_->setCovMatrix( &covHadB );
    jetlepb_->setCovMatrix( &covLepB );
    lepton_ ->setCovMatrix( &covLep  );
    met_    ->setCovMatrix( &covMET  );
    
    if(constrainSumPt_){
      // setup Px and Py constraint for curent event configuration so that sum Pt will be conserved
      sumPxConstr_->setConstraint( p4Had1.Px() + p4Had2.Px() + p4HadB.Px() + p4LepB.Px() + p4Lep.Px() + p4Neutrino.Px() );
      sumPyConstr_->setConstraint( p4Had1.Py() + p4Had2.Py() + p4HadB.Py() + p4LepB.Py() + p4Lep.Py() + p4Neutrino.Py() );
    }

    //do the fit
    myKinFitter->fit();

    //get chi^2 of fit
    chi2_prev = chi2_curr;
    F_prev = F_curr;
    chi2_curr = 1000.;
    F_curr = 0.;

    chi2_curr = myKinFitter->getS();
    F_curr = myKinFitter->getF();

    fitprob_curr = 0.;
    fitprob_curr =  TMath::Prob(myKinFitter->getS(), myKinFitter->getNDF());
    fit_status_each = myKinFitter->getStatus();

    jet_combi_all.push_back(h1);
    jet_combi_all.push_back(h2);
    jet_combi_all.push_back(bh);
    jet_combi_all.push_back(bl);

    //TEST: Set values for all combinations
    solution2.clear();
    input_vectors_all.clear();

    TLorentzVector sol_had1_2 = TLorentzVector(jethad1_->getCurr4Vec()->Px(), jethad1_->getCurr4Vec()->Py(), jethad1_->getCurr4Vec()->Pz(),jethad1_->getCurr4Vec()->E());
    TLorentzVector sol_had2_2 = TLorentzVector(jethad2_->getCurr4Vec()->Px(), jethad2_->getCurr4Vec()->Py(),jethad2_->getCurr4Vec()->Pz(),jethad2_->getCurr4Vec()->E());
    TLorentzVector sol_hadb_2 = TLorentzVector(jethadb_->getCurr4Vec()->Px(), jethadb_->getCurr4Vec()->Py(), jethadb_->getCurr4Vec()->Pz(), jethadb_->getCurr4Vec()->E());
    TLorentzVector sol_lepb_2 = TLorentzVector(jetlepb_->getCurr4Vec()->Px(), jetlepb_->getCurr4Vec()->Py(), jetlepb_->getCurr4Vec()->Pz(), jetlepb_->getCurr4Vec()->E());
    TLorentzVector sol_lepton_2 = TLorentzVector(lepton_->getCurr4Vec()->Px(), lepton_->getCurr4Vec()->Py(), lepton_->getCurr4Vec()->Pz(), lepton_->getCurr4Vec()->E());
    TLorentzVector sol_met_2 =  TLorentzVector(met_->getCurr4Vec()->Px(), met_->getCurr4Vec()->Py(), met_->getCurr4Vec()->Pz(), met_->getCurr4Vec()->E());
 

    solution2.push_back(sol_had1_2);
    solution2.push_back(sol_had2_2);
    solution2.push_back(sol_hadb_2);
    solution2.push_back(sol_lepb_2);
    solution2.push_back(sol_lepton_2);
    solution2.push_back(sol_met_2);

    input_vectors_all.push_back(p4Had1);
    input_vectors_all.push_back(p4Had2);
    input_vectors_all.push_back(p4HadB);
    input_vectors_all.push_back(p4LepB);
    input_vectors_all.push_back(p4Lep);
    input_vectors_all.push_back(p4Neutrino);

    dm_w_lep_all = (p4Lep+p4Neutrino).M() - (sol_lepton_2+sol_met_2).M();
    dm_w_had_all = (p4Had1+p4Had2).M() - (sol_had1_2+sol_had2_2).M();
    dm_top_lep_all = (p4Lep+p4Neutrino+p4LepB).M() - (sol_lepton_2+sol_met_2+sol_lepb_2).M();
    dm_top_had_all = (p4Had1+p4Had2+p4HadB).M() - (sol_had1_2+sol_had2_2+sol_hadb_2).M();
    if(debug){
      std::cout<<"wlep all "<<dm_w_lep_all<<std::endl;
      std::cout<<"whad all "<<dm_w_had_all<<std::endl;
      std::cout<<"top lep all "<<dm_top_lep_all<<std::endl;
      std::cout<<"top had all "<<dm_top_had_all<<std::endl;
    }
    //Status of Fitter: 0 converged, 1 not converged (reaches e.g. max iteration), -1 no fit performed, -10 aborted, 10 running
    if(fit_status_each==0 && fitprob_curr>fit_prob){
  	fit_prob = fitprob_curr;
  	chi2_lowest = chi2_curr;
  	fit_status = fit_status_each;

  	jet_combi_best.clear();
  	jet_combi_best.push_back(h1);
  	jet_combi_best.push_back(h2);
  	jet_combi_best.push_back(bh);
  	jet_combi_best.push_back(bl);

  	solution.clear();
  	input_vectors_fit_found.clear();

        TLorentzVector sol_had1 = TLorentzVector(jethad1_->getCurr4Vec()->Px(), jethad1_->getCurr4Vec()->Py(), jethad1_->getCurr4Vec()->Pz(),jethad1_->getCurr4Vec()->E());
  	TLorentzVector sol_had2 = TLorentzVector(jethad2_->getCurr4Vec()->Px(), jethad2_->getCurr4Vec()->Py(),jethad2_->getCurr4Vec()->Pz(),jethad2_->getCurr4Vec()->E());
  	TLorentzVector sol_hadb = TLorentzVector(jethadb_->getCurr4Vec()->Px(), jethadb_->getCurr4Vec()->Py(), jethadb_->getCurr4Vec()->Pz(), jethadb_->getCurr4Vec()->E());
  	TLorentzVector sol_lepb = TLorentzVector(jetlepb_->getCurr4Vec()->Px(), jetlepb_->getCurr4Vec()->Py(), jetlepb_->getCurr4Vec()->Pz(), jetlepb_->getCurr4Vec()->E());
  	TLorentzVector sol_lepton = TLorentzVector(lepton_->getCurr4Vec()->Px(), lepton_->getCurr4Vec()->Py(), lepton_->getCurr4Vec()->Pz(), lepton_->getCurr4Vec()->E());
  	TLorentzVector sol_met =  TLorentzVector(met_->getCurr4Vec()->Px(), met_->getCurr4Vec()->Py(), met_->getCurr4Vec()->Pz(), met_->getCurr4Vec()->E());

  	solution.push_back(sol_had1);
  	solution.push_back(sol_had2);
  	solution.push_back(sol_hadb);
  	solution.push_back(sol_lepb);
  	solution.push_back(sol_lepton);
  	solution.push_back(sol_met);

  	input_vectors_fit_found.push_back(p4Had1);
  	input_vectors_fit_found.push_back(p4Had2);
  	input_vectors_fit_found.push_back(p4HadB);
  	input_vectors_fit_found.push_back(p4LepB);
  	input_vectors_fit_found.push_back(p4Lep);
  	input_vectors_fit_found.push_back(p4Neutrino);

  	int jet_variables_size(0);
  	switch(jetParam_){
  	case KinFitTestSetup::kEMom :         
  	  jet_variables_size = 4;
  	  event.set(parametrisation, 0);
  	  break;
  	case KinFitTestSetup::kEtEtaPhi :     
  	  jet_variables_size = 3;
  	  event.set(parametrisation, 1);
  	  break;
  	}
  	const TTbarGen* ttgen_test(0);
  	if(isMC_){
  	  const auto& ttbargen_test = event.get(h_ttbar_gen);
  	  ttgen_test = &ttbargen_test;
  	}
  	setPull(jet_variables_size, ttgen_test, 0, isMC_);

  	if(debug){
  	  std::cout<<"fit found"<<std::endl;
  	  std::cout<<"SOLUTION:"<<std::endl;
  	  for(unsigned int j=0; j<solution.size();j++){
  	    solution[j].Print();
  	    std::cout<<"ET "<<solution[j].Et()<<std::endl;
  	    std::cout<<"pT "<<solution[j].Pt()<<std::endl;
  	    std::cout<<"Mass "<<solution[j].M()<<std::endl;
  	  }
  	  std::cout<<"chi2 "<<chi2_curr<<std::endl;
  	  std::cout<<"Fit probability "<<fitprob_curr<<std::endl;
  	  std::cout<<"INPUT:"<<std::endl;
  	  for(unsigned int j=0; j<input_vectors_fit_found.size();j++){
  	    input_vectors_fit_found[j].Print();
  	    std::cout<<"Et "<<input_vectors_fit_found[j].Et()<<std::endl;
  	    std::cout<<"Pt "<<input_vectors_fit_found[j].Pt()<<std::endl;
  	    std::cout<<"Mass "<<input_vectors_fit_found[j].M()<<std::endl;
  	  }
  	}

  	fit_found = true;

  	dm_w_lep = (p4Lep+p4Neutrino).M() - (sol_lepton+sol_met).M();
  	dm_w_had = (p4Had1+p4Had2).M() - (sol_had1+sol_had2).M();
  	dm_top_lep = (p4Lep+p4Neutrino+p4LepB).M() - (sol_lepton+sol_met+sol_lepb).M();
  	dm_top_had = (p4Had1+p4Had2+p4HadB).M() - (sol_had1+sol_had2+sol_hadb).M();
    if(debug){
      std::cout<<"%%%%%%%%%%% fit found"<<std::endl;
      std::cout<<"wlep  "<<dm_w_lep<<std::endl;
      std::cout<<"whad  "<<dm_w_had<<std::endl;
      std::cout<<"top lep  "<<dm_top_lep<<std::endl;
      std::cout<<"top had  "<<dm_top_had<<std::endl;
    }
    }
    

    int jet_variables_size_all(0);
    switch(jetParam_){
    case KinFitTestSetup::kEMom :         
      jet_variables_size_all = 4; 	  
      event.set(parametrisation, 0);
     break;
    case KinFitTestSetup::kEtEtaPhi :     
      jet_variables_size_all = 3; 
      event.set(parametrisation, 1);
      break;
    }

    const TTbarGen* ttgen_test1(0);
    if(isMC_){
      const auto& ttbargen_test1 = event.get(h_ttbar_gen);
      ttgen_test1 = &ttbargen_test1;
    }

    setPull(jet_variables_size_all, ttgen_test1, 1, isMC_);
 
    event.set(pull_had1, pull_had1_all);
    event.set(pull_had2, pull_had2_all);
    event.set(pull_hadb, pull_hadb_all);
    event.set(pull_lepb, pull_lepb_all);
    event.set(pull_lepton, pull_lepton_all);
    event.set(pull_neu, pull_neu_all); 
    event.set(pull_gtr_had1_et, -100.);
    event.set(pull_gtr_had1_eta, -100.);
    event.set(pull_gtr_had1_phi, -100.);
    event.set(pull_gtr_had2_et, -100.);
    event.set(pull_gtr_had2_eta, -100.);
    event.set(pull_gtr_had2_phi, -100.);
    event.set(pull_gtr_hadb_et, -100.);
    event.set(pull_gtr_hadb_eta, -100.);
    event.set(pull_gtr_hadb_phi, -100.);
    event.set(pull_gtr_lepb_et, -100.);
    event.set(pull_gtr_lepb_eta, -100.);
    event.set(pull_gtr_lepb_phi, -100.);
    event.set(pull_gtr_lepton_et, -100.);
    event.set(pull_gtr_lepton_eta, -100.);
    event.set(pull_gtr_lepton_phi, -100.);
    event.set(pull_gtr_neu_et, -100.);
    event.set(pull_gtr_neu_eta, -100.);
    event.set(pull_gtr_neu_phi, -100.);

    if((verbosity==2 || verbosity==3 ) ){
      std::cout<<"##############################"<<std::endl;
      std::cout<<"input: "<<std::endl;
      std::cout<<"lepton type "<<leptonType<<std::endl;
      std::cout<<"p4LepB = "<<std::endl;
      p4LepB.Print();
      std::cout<<"p4LepB pt = "<<p4LepB.Pt()<<std::endl;
      std::cout<<"p4LepB mass = "<<p4LepB.M()<<std::endl;
      std::cout<<"p4HadB = "<<std::endl;
      p4HadB.Print();
      std::cout<<"p4HadB pt = "<<p4HadB.Pt()<<std::endl;
      std::cout<<"p4HadB mass = "<<p4HadB.M()<<std::endl;
      std::cout<<"p4Had1 = "<<std::endl;
      p4Had1.Print();
      std::cout<<"p4Had1 pt = "<<p4Had1.Pt()<<std::endl;
      std::cout<<"p4Had1 mass = "<<p4Had1.M()<<std::endl;
      std::cout<<"p4Had2 = "<<std::endl;
      p4Had2.Print();
      std::cout<<"p4Had2 pt = "<<p4Had2.Pt()<<std::endl;
      std::cout<<"p4Had2 mass = "<<p4Had2.M()<<std::endl;
      std::cout<<"p4Lep  = "<<std::endl;
      p4Lep.Print();
      std::cout<<"p4Lep pt = "<<p4Lep.Pt()<<std::endl;
      std::cout<<"p4Lep mass = "<<p4Lep.M()<<std::endl;
      std::cout<<"p4Neutrino  = "<<std::endl;
      p4Neutrino.Print();
      std::cout<<"p4Neutrino pt = "<<p4Neutrino.Pt()<<std::endl;
      std::cout<<"p4Neutrino mass = "<<p4Neutrino.M()<<std::endl;
      std::cout<<"WHad mass = "<<(p4Had1+p4Had2).M()<<std::endl;
      std::cout<<"WLep mass = "<<(p4Lep+p4Neutrino).M()<<std::endl;
      std::cout<<"TLep mass = "<<(p4Lep+p4Neutrino+p4LepB).M()<<std::endl;
      std::cout<<"THad mass = "<<(p4Had1+p4Had2+p4HadB).M()<<std::endl;
      std::cout<<"##############################"<<std::endl;
    }    

    event.set(jet_permutation_before_dm_cut, 0);
    event.set(jet_permutation_after_dm_cut, 0);


    setBestFit(event, solution2, input_vectors_all, leptonType, chi2_curr, fitprob_curr, fit_status_each, fit_found, false, jet_combi_all);
    h_fitted_allcomb->fill(event);
    if(fit_status_each==0) h_fitted_allcomb_fit_converged->fill(event);

  }



void KinFitTestModule::setBestFit(Event & event, std::vector<TLorentzVector> solution_vec, std::vector<TLorentzVector> input_vector, KinFitTestSetup::ObjectType leptonType, float chi2_lowest, float fit_prob, float fit_status, bool fit_found, bool is_input, std::vector<int> jet_combi){

    jetHad1__p4 = solution_vec[0];
    jetHad2__p4 = solution_vec[1];
    jetHadB__p4 = solution_vec[2];
    jetLepB__p4 = solution_vec[3];
    lepton__p4 =  solution_vec[4];
    met__p4 =     solution_vec[5];

    TLorentzVector input_p4Had1 = input_vector[0];
    TLorentzVector input_p4Had2 = input_vector[1];
    TLorentzVector input_p4HadB = input_vector[2];
    TLorentzVector input_p4LepB = input_vector[3];
    TLorentzVector input_p4Lep  = input_vector[4];
    TLorentzVector input_p4MET  = input_vector[5];

    if(debug){
      std::cout<<"Set best fit"<<std::endl;
      std::cout<<"SOLUTION:"<<std::endl;
      std::cout<<"WHad mass = "<<(jetHad1__p4+jetHad2__p4).M()<<std::endl;
      std::cout<<"WLep mass = "<<(lepton__p4+met__p4).M()<<std::endl;
      std::cout<<"TLep mass = "<<(lepton__p4+met__p4+jetLepB__p4).M()<<std::endl;
      std::cout<<"THad mass = "<<(jetHad1__p4+jetHad2__p4+jetHadB__p4).M()<<std::endl;  
      std::cout<<"chi2 "<<chi2_lowest<<std::endl;
      std::cout<<"INPUT:"<<std::endl;
      std::cout<<"WHad mass = "<<(input_p4Had1+input_p4Had2).M()<<std::endl;
      std::cout<<"WLep mass = "<<(input_p4Lep+input_p4MET).M()<<std::endl;
      std::cout<<"TLep mass = "<<(input_p4Lep+input_p4MET+input_p4LepB).M()<<std::endl;
      std::cout<<"THad mass = "<<(input_p4Had1+input_p4Had2+input_p4HadB).M()<<std::endl;    
    }

    float chi2(chi2_lowest);
    float fit_prob_(fit_prob);
    TLorentzVector whad_v4 = jetHad1__p4 + jetHad2__p4;
    TLorentzVector wlep_v4 = lepton__p4 + met__p4;
    TLorentzVector toplep_v4 = wlep_v4 + jetLepB__p4;
    TLorentzVector tophad_v4 = jetHadB__p4 + jetHad1__p4 + jetHad2__p4 ;
    float neu_pt(met__p4.Pt());
    float status(fit_status);
    bool found(fit_found);

    met_cut = false;
    matchable = false;
    b_jet_matching = false;
    match_ttbar_decay = false;
    correct_match = false;
    high_met = false;

    if(neu_pt>50.) met_cut = true;

    dR_tlep = infinity;
    dR_thad = infinity;
    dR_wlep = infinity;
    dR_jetq1 = infinity;
    dR_jetq2 = infinity;
    dR_jetbh = infinity;
    dR_jetbl = infinity;
    dR_lepton = infinity;
    dR_neutrino = infinity;
    dR_combined = infinity;

    high_blep = false;

    const TTbarGen* ttgen(0);
    if(isMC_){
      const auto& ttbargen = event.get(h_ttbar_gen);
      ttgen = &ttbargen;

    if(ttgen->IsSemiLeptonicDecay()){
      if(ttgen->Neutrino().v4().Pt()>100.){
       	high_met=true;
      }

      high_blep = fabs(input_p4LepB.Pt() - ttgen->BLep().v4().Pt())>50.;

      if((verbosity==3 && high_blep) || debug){
	std::cout<<"q1 pT gen "<<ttgen->Q1().v4().Pt()<<std::endl;
	std::cout<<"q2 pT gen "<<ttgen->Q2().v4().Pt()<<std::endl;
	std::cout<<"bhad pT gen "<<ttgen->BHad().v4().Pt()<<std::endl;
	std::cout<<"blep pT gen "<<ttgen->BLep().v4().Pt()<<std::endl;
	std::cout<<"Lepton pT gen "<<ttgen->ChargedLepton().v4().Pt()<<std::endl;
	std::cout<<"MET pT gen "<<ttgen->Neutrino().v4().Pt()<<std::endl;

      }

      const LorentzVector gen_Tlep_v4 = ttgen->TopLep().v4();
      TLorentzVector rec_fit_Tlep(lepton__p4 + met__p4 + jetLepB__p4);
      LorentzVector rec_fit_Tlep_v4;
      rec_fit_Tlep_v4.SetXYZT(rec_fit_Tlep.X(), rec_fit_Tlep.Y(), rec_fit_Tlep.Z(), rec_fit_Tlep.T());
      dR_tlep = h_fitted_dR_cut->deltaR_fit(rec_fit_Tlep_v4, gen_Tlep_v4);
      
      const LorentzVector gen_Thad_v4 = ttgen->TopHad().v4();
      TLorentzVector rec_fit_Thad(jetHad1__p4 + jetHad2__p4 + jetHadB__p4);
      LorentzVector rec_fit_Thad_v4;
      rec_fit_Thad_v4.SetXYZT(rec_fit_Thad.X(), rec_fit_Thad.Y(), rec_fit_Thad.Z(), rec_fit_Thad.T());
      dR_thad = h_fitted_dR_cut->deltaR_fit(rec_fit_Thad_v4, gen_Thad_v4);
      
      const LorentzVector gen_Wlep_v4 = ttgen->WLep().v4();
      TLorentzVector rec_fit_Wlep(lepton__p4 + met__p4);
      LorentzVector rec_fit_Wlep_v4;
      rec_fit_Wlep_v4.SetXYZT(rec_fit_Wlep.X(), rec_fit_Wlep.Y(), rec_fit_Wlep.Z(), rec_fit_Wlep.T());
      dR_wlep = h_fitted_dR_cut->deltaR_fit(rec_fit_Wlep_v4, gen_Wlep_v4);

      //matched and unmatched ttbar events
      const LorentzVector gen_jetq1_v4 = ttgen->Q1().v4();
      const LorentzVector gen_jetq2_v4 = ttgen->Q2().v4();
      LorentzVector rec_fit_jetq1_lv;
      rec_fit_jetq1_lv.SetXYZT(jetHad1__p4.X(), jetHad1__p4.Y(), jetHad1__p4.Z(), jetHad1__p4.T());
      LorentzVector rec_fit_jetq2_lv;
      rec_fit_jetq2_lv.SetXYZT(jetHad2__p4.X(), jetHad2__p4.Y(), jetHad2__p4.Z(), jetHad2__p4.T());
      float dR_jetq1_1 = h_fitted_dR_cut->deltaR_fit(rec_fit_jetq1_lv, gen_jetq1_v4);
      float dR_jetq1_2 = h_fitted_dR_cut->deltaR_fit(rec_fit_jetq1_lv, gen_jetq2_v4);
      float dR_jetq2_1 = h_fitted_dR_cut->deltaR_fit(rec_fit_jetq2_lv, gen_jetq1_v4);
      float dR_jetq2_2 = h_fitted_dR_cut->deltaR_fit(rec_fit_jetq2_lv, gen_jetq2_v4);
      if(fabs(dR_jetq1_1)<fabs(dR_jetq1_2)){
	dR_jetq1 = dR_jetq1_1;
	dR_jetq2 = dR_jetq2_2;
      }
      else{
	dR_jetq1 = dR_jetq1_2;
	dR_jetq2 = dR_jetq2_1;
      }
 
      const LorentzVector gen_bhad_v4 = ttgen->BHad().v4();
      LorentzVector rec_fit_jetbh_lv;
      rec_fit_jetbh_lv.SetXYZT(jetHadB__p4.X(), jetHadB__p4.Y(), jetHadB__p4.Z(), jetHadB__p4.T());
      dR_jetbh= h_fitted_dR_cut->deltaR_fit(rec_fit_jetbh_lv, gen_bhad_v4);

      const LorentzVector gen_blep_v4 = ttgen->BLep().v4();
      LorentzVector rec_fit_jetbl_lv;
      rec_fit_jetbl_lv.SetXYZT(jetLepB__p4.X(), jetLepB__p4.Y(), jetLepB__p4.Z(), jetLepB__p4.T());
      dR_jetbl= h_fitted_dR_cut->deltaR_fit(rec_fit_jetbl_lv, gen_blep_v4);

      const LorentzVector gen_lepton_v4 = ttgen->ChargedLepton().v4();
      LorentzVector rec_fit_lepton_lv;
      rec_fit_lepton_lv.SetXYZT(lepton__p4.X(), lepton__p4.Y(), lepton__p4.Z(), lepton__p4.T());
      dR_lepton= h_fitted_dR_cut->deltaR_fit(rec_fit_lepton_lv, gen_lepton_v4);

      const LorentzVector gen_neutrino_v4 = ttgen->Neutrino().v4();
      LorentzVector rec_fit_neutrino_lv;
      rec_fit_neutrino_lv.SetXYZT(met__p4.X(), met__p4.Y(), met__p4.Z(), met__p4.T());
      dR_neutrino= h_fitted_dR_cut->deltaR_fit(rec_fit_neutrino_lv, gen_neutrino_v4);


      //check if event is "matchable": pt and eta of objects inside the selection criteria
      bool pt_jets(false), eta_jets(false), pt_lepton(false), eta_lepton(false);
      if(gen_jetq1_v4.Pt()>25. && gen_jetq2_v4.Pt()>25. && gen_bhad_v4.Pt()>25. && gen_blep_v4.Pt()>25.) pt_jets = true;
      if(fabs(gen_jetq1_v4.Eta())<2.4 && fabs(gen_jetq2_v4.Eta())<2.4 && fabs(gen_bhad_v4.Eta())<2.4 && fabs(gen_blep_v4.Eta())<2.4 ) eta_jets = true;

      switch(leptonType){
      case KinFitTestSetup::kElectron :       if(gen_lepton_v4.Pt()>24.) pt_lepton = true;	break;
      case KinFitTestSetup::kMuon :           if(gen_lepton_v4.Pt()>24.) pt_lepton = true; break;
      case KinFitTestSetup::kMet : break;
      case KinFitTestSetup::kBJet : break;
      case KinFitTestSetup::kUdscJet : break;
      }
      if(fabs(gen_lepton_v4.Eta())<2.1) eta_lepton = true;

      if(pt_jets && eta_jets && pt_lepton && eta_lepton) matchable = true;


      //check if b-jets of fit matched to b-quarks
      float dR_b_lep_had = h_fitted_dR_cut->deltaR_fit(rec_fit_jetbl_lv, gen_bhad_v4);
      float dR_b_had_lep = h_fitted_dR_cut->deltaR_fit(rec_fit_jetbh_lv, gen_blep_v4);
      if((dR_jetbh<0.3 || dR_b_had_lep<0.3) && (dR_jetbl<0.3 || dR_b_lep_had<0.3)) b_jet_matching = true;

      //check if objects of fit matched to (any) ttbar decay products
      std::vector<LorentzVector> gen_jets;
      gen_jets.clear();
      gen_jets.push_back(gen_jetq1_v4);
      gen_jets.push_back(gen_jetq2_v4);
      gen_jets.push_back(gen_bhad_v4);
      //     gen_jets.push_back(gen_blep_v4);

      std::vector<LorentzVector> fit_jets;
      fit_jets.clear();
      fit_jets.push_back(rec_fit_jetq1_lv);
      fit_jets.push_back(rec_fit_jetq2_lv);
      fit_jets.push_back(rec_fit_jetbh_lv);
      //     fit_jets.push_back(rec_fit_jetbl_lv);

      //Setup dpt_Dpt = (pt_fit - pt_gen)/pt_gen for the jets
      float dpt_Dpt_had1_1 = (jetHad1__p4.Pt() - gen_jetq1_v4.Pt())/gen_jetq1_v4.Pt();
      float dpt_Dpt_had1_2 = (jetHad1__p4.Pt() - gen_jetq2_v4.Pt())/gen_jetq2_v4.Pt();
      float dpt_Dpt_had2_1 = (jetHad2__p4.Pt() - gen_jetq1_v4.Pt())/gen_jetq1_v4.Pt();
      float dpt_Dpt_had2_2 = (jetHad2__p4.Pt() - gen_jetq2_v4.Pt())/gen_jetq2_v4.Pt();
      float dpt_Dpt_hadb, dpt_Dpt_lepb;
      dpt_Dpt_hadb  = (jetHadB__p4.Pt() - gen_bhad_v4.Pt())/gen_bhad_v4.Pt();
      dpt_Dpt_lepb  = (jetLepB__p4.Pt() - gen_blep_v4.Pt())/gen_blep_v4.Pt();

      int nr_matched_jets(0);

      for(unsigned int gj=0; gj<gen_jets.size(); gj++){
	int matches(0);
	for(unsigned int fj=0; fj<fit_jets.size(); fj++){
	  float dR = h_fitted_dR_cut->deltaR_fit(fit_jets[fj], gen_jets[gj]);
	  if(dR<0.3){
	    matches++;
	  }
	}
	if(matches>=1) nr_matched_jets++;
      }
      float dR_blep_fit =  h_fitted_dR_cut->deltaR_fit(rec_fit_jetbl_lv, gen_blep_v4);

      float dR_lep_fit =  h_fitted_dR_cut->deltaR_fit(rec_fit_lepton_lv, gen_lepton_v4);
      if(nr_matched_jets>=4 && dR_blep_fit<0.3 && dR_lep_fit<0.3) match_ttbar_decay = true;

      //check if correct match: jets of fit are exact the quarks
      if(((dR_jetq1_1<0.3 && dR_jetq2_2<0.3) || (dR_jetq2_1<0.3 && dR_jetq1_2<0.3)) && dR_jetbh<0.3 && dR_jetbl<0.3) {
	dR_combined = dR_jetq1 + dR_jetq2 + dR_jetbh + dR_jetbl + dR_lepton + dR_neutrino;
	//	if(((fabs(dpt_Dpt_had1_1)<0.15 && fabs(dpt_Dpt_had2_2)<0.15) || (fabs(dpt_Dpt_had2_1)<0.15 && fabs(dpt_Dpt_had1_2)<0.15)) && fabs(dpt_Dpt_hadb)<0.15 && fabs(dpt_Dpt_lepb)<0.15){
	  correct_match = true;
	  if(debug){
	    std::cout<<"dpt_Dpt_had1_1 "<<dpt_Dpt_had1_1<<std::endl;
	    std::cout<<"dpt_Dpt_had1_2 "<<dpt_Dpt_had1_2<<std::endl;
	    std::cout<<"dpt_Dpt_had2_1 "<<dpt_Dpt_had2_1<<std::endl;
	    std::cout<<"dpt_Dpt_had2_2 "<<dpt_Dpt_had2_2<<std::endl;
	    std::cout<<"dpt_Dpt_hadb "<<dpt_Dpt_hadb<<std::endl;
	    std::cout<<"dpt_Dpt_lepb "<<dpt_Dpt_lepb<<std::endl;
	  }

	if(debug)	std::cout<<"correctly matched event "<<dR_combined<<std::endl;
      }

    }
    }

    if(debug){
     std::cout<<"combinations: "<<std::endl;
     std::cout<<"status of fit "<<fit_status<<std::endl;
     std::cout<<"jetLepB__p4 = "<<std::endl;
     jetLepB__p4.Print();
     std::cout<<"jetLepB pt = "<<jetLepB__p4.Pt()<<std::endl;
     std::cout<<"jetHadB__p4 = "<<std::endl;
     jetHadB__p4.Print();
     std::cout<<"jetHadB pt = "<<jetHadB__p4.Pt()<<std::endl;
     std::cout<<"jetHad1__p4 = "<<std::endl;
     jetHad1__p4.Print();
     std::cout<<"jetHad1 pt = "<<jetHad1__p4.Pt()<<std::endl;
     std::cout<<"jetHad2__p4 = "<<std::endl;
     jetHad2__p4.Print();
     std::cout<<"jetHad2 pt = "<<jetHad2__p4.Pt()<<std::endl;
     std::cout<<"lepton__p4 = "<<std::endl;
     lepton__p4.Print();
     std::cout<<"lepton pt = "<<lepton__p4.Pt()<<std::endl;
     std::cout<<"met__p4 = "<<std::endl;
     met__p4.Print();
     std::cout<<"met pt = "<<met__p4.Pt()<<std::endl;
     std::cout<<"wlep sol mass = "<<(lepton__p4 + met__p4).M()<<std::endl;
     std::cout<<"whad sol mass = "<<(jetHad1__p4 + jetHad2__p4).M()<<std::endl;
     std::cout<<"tlep sol mass = "<<(lepton__p4 + met__p4 + jetLepB__p4).M()<<std::endl;
     std::cout<<"thad sol mass = "<<(jetHad1__p4 + jetHad2__p4 + jetHadB__p4).M()<<std::endl;
     std::cout<<"fitprob_curr = "<<fit_prob_<<std::endl;
     std::cout<<"chi_curr = "<<chi2<<std::endl;
     std::cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;

  } 

    event.set(rec_fit_chi2, chi2);
    event.set(fit_probab, fit_prob_);
    event.set(rec_fit_whad_v4, whad_v4);
    event.set(rec_fit_wlep_v4, wlep_v4);
    event.set(rec_fit_blep_v4,jetLepB__p4);
    event.set(rec_fit_bhad_v4,jetHadB__p4);
    event.set(rec_fit_thad_v4, tophad_v4);
    event.set(rec_fit_tlep_v4, toplep_v4);
    event.set(rec_fit_lepton_v4,  lepton__p4);
    event.set(rec_fit_neutrino_v4,  met__p4);
    event.set(rec_fit_jetq1_v4,jetHad1__p4);
    event.set(rec_fit_jetq2_v4,jetHad2__p4);
    event.set(rec_fit_status, status);
    event.set(rec_fit_fit_found, found);
    event.set(rec_fit_is_input, is_input);

    event.set(input_lepton_v4, input_p4Lep);
    event.set(input_neutrino_v4, input_p4MET);
    event.set(input_jetq1_v4, input_p4Had1);
    event.set(input_jetq2_v4, input_p4Had2);
    event.set(input_jetbh_v4, input_p4HadB);
    event.set(input_jetbl_v4, input_p4LepB);
    event.set(rec_fit_dR_combined, dR_combined);
    event.set(jet_combi_1, jet_combi[0]);
    event.set(jet_combi_2, jet_combi[1]);
    event.set(jet_combi_3, jet_combi[2]);
    event.set(jet_combi_4, jet_combi[3]);
  }



  void KinFitTestModule::setupJets()
{
  TMatrixD empty3x3(3,3); 
  TMatrixD empty4x4(4,4);
  switch(jetParam_){ // setup jets according to parameterization
  case KinFitTestSetup::kEMom :
    jethad1_= new TFitParticleEMomDev   ("Jet1_pt", "Jet1_pt", 0, &empty4x4);
    jethad2_= new TFitParticleEMomDev   ("Jet2_pt", "Jet2_pt", 0, &empty4x4);
    jethadb_= new TFitParticleEMomDev   ("Jet3_pt", "Jet3_pt", 0, &empty4x4);
    jetlepb_= new TFitParticleEMomDev   ("Jet4_pt", "Jet4_pt", 0, &empty4x4);
     break;
  case KinFitTestSetup::kEtEtaPhi :
    jethad1_= new TFitParticleEtEtaPhi  ("Jet1_eta", "Jet1_eta", 0, &empty3x3);
    jethad2_= new TFitParticleEtEtaPhi  ("Jet2_eta", "Jet2_eta", 0, &empty3x3);
    jethadb_= new TFitParticleEtEtaPhi  ("Jet3_eta", "Jet3_eta", 0, &empty3x3);
    jetlepb_= new TFitParticleEtEtaPhi  ("Jet4_eta", "Jet4_eta", 0, &empty3x3);
    break;
  }
}
void KinFitTestModule::setupLeptons()
{
  TMatrixD empty3x3(3,3); 
  switch(lepParam_){ // setup lepton according to parameterization
  case KinFitTestSetup::kEMom       : 
    lepton_  = new TFitParticleEScaledMomDev("Lepton",   "Lepton",   0, &empty3x3);
    break;
  case KinFitTestSetup::kEtEtaPhi   : 
    lepton_  = new TFitParticleEtEtaPhi     ("Lepton",   "Lepton",   0, &empty3x3); 
    break;
  }
  switch(metParam_){ // setup MET according to parameterization
  case KinFitTestSetup::kEMom       : 
    met_= new TFitParticleEScaledMomDev("MET", "MET", 0, &empty3x3);
    break;
  case KinFitTestSetup::kEtEtaPhi   : 
    met_= new TFitParticleEtEtaPhi     ("MET", "MET", 0, &empty3x3);
    break;
  }
}




void KinFitTestModule::setupConstraints() 
{
  //Decay width top, approximately
  // float decw_top = 0.175*pow((mTop_/mW_),3); //Result with mTop=173. and mW_=80.4 is 1.74 GeV
  // float decw_w = 2.085;
  

  massConstr_[kWHadMass] = new TFitConstraintM("WMassHad",      "WMassHad",      0, 0, mW_);
  massConstr_[kWLepMass      ] = new TFitConstraintM("WMassLep",      "WMassLep",      0, 0, mW_);
  massConstr_[kTopHadMass    ] = new TFitConstraintM("TopMassHad",    "TopMassHad",    0, 0, mTop_);
  massConstr_[kTopLepMass    ] = new TFitConstraintM("TopMassLep",    "TopMassLep",    0, 0, mTop_);
  //  massConstr_[kMET  ] = new TFitConstraintM("MET",  "MET",  0, 0,    0.);
  //massConstr_[kEqualTopMasses] = new TFitConstraintM("EqualTopMasses","EqualTopMasses",0, 0,    0.);
  // //  massConstr_[kEqualWMasses] = new TFitConstraintM("EqualWMasses","EqualWMasses",0, 0,    0.);
  sumPxConstr_                 = new TFitConstraintEp("SumPx",        "SumPx", 0, TFitConstraintEp::pX, 0.);
  sumPyConstr_                 = new TFitConstraintEp("SumPy",        "SumPy", 0, TFitConstraintEp::pY, 0.);


  massConstr_[kWHadMass      ]->addParticles1(jethad1_,   jethad2_    );
  massConstr_[kWLepMass      ]->addParticles1(lepton_, met_);
  massConstr_[kTopHadMass    ]->addParticles1(jethad1_, jethad2_, jethadb_);
  massConstr_[kTopLepMass    ]->addParticles1(lepton_, met_, jetlepb_);
  //  massConstr_[kMET           ]->addParticle1 (met_);
  //massConstr_[kEqualTopMasses]->addParticles1(jethad1_, jethad2_, jethadb_);
  //massConstr_[kEqualTopMasses]->addParticles2(lepton_, met_, jetlepb_);
  // // massConstr_[kEqualWMasses]->addParticles1(jethad1_, jethad2_);
  // // massConstr_[kEqualWMasses]->addParticles2(lepton_, met_);
  sumPxConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);
  sumPyConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);


  // float decw_top = 1.7;
  // float decw_w = 2.;
  // float decw_top = 60.;
  // float decw_w = 30.;

  // massConstr_[kWHadMass] = new TFitConstraintMGaus("WMassHad",      "WMassHad",      0, 0, mW_  , decw_w);
  // massConstr_[kWLepMass      ] = new TFitConstraintMGaus("WMassLep",      "WMassLep",      0, 0, mW_ , decw_w);
  // massConstr_[kTopHadMass    ] = new TFitConstraintMGaus("TopMassHad",    "TopMassHad",    0, 0, mTop_, decw_top);
  // massConstr_[kTopLepMass    ] = new TFitConstraintMGaus("TopMassLep",    "TopMassLep",    0, 0, mTop_, decw_top );
  // // massConstr_[kEqualTopMasses] = new TFitConstraintMGaus("EqualTopMasses","EqualTopMasses",0, 0,    0.1, 1.);
  // //massConstr_[kEqualWMasses] = new TFitConstraintMGaus("EqualWMasses","EqualWMasses",0, 0,    .1, 1.);
  // sumPxConstr_                 = new TFitConstraintEp("SumPx",        "SumPx", 0, TFitConstraintEp::pX, 0.);
  // sumPyConstr_                 = new TFitConstraintEp("SumPy",        "SumPy", 0, TFitConstraintEp::pY, 0.);

  // massConstr_[kWHadMass      ]->addParticles1(jethad1_,   jethad2_    );
  // massConstr_[kWLepMass      ]->addParticles1(lepton_, met_);
  // massConstr_[kTopHadMass    ]->addParticles1(jethad1_, jethad2_, jethadb_);
  // massConstr_[kTopLepMass    ]->addParticles1(lepton_, met_, jetlepb_);
  // // massConstr_[kEqualTopMasses]->addParticles1(jethad1_, jethad2_, jethadb_);
  // // massConstr_[kEqualTopMasses]->addParticles2(lepton_, met_, jetlepb_);
  // // massConstr_[kEqualWMasses]->addParticles1(jethad1_, jethad2_);
  // // massConstr_[kEqualWMasses]->addParticles2(lepton_, met_);
  // sumPxConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);
  // sumPyConstr_->addParticles(lepton_, met_, jethad1_, jethad2_, jethadb_, jetlepb_);



  // if(std::find(constrList_.begin(), constrList_.end(), kSumPt)!=constrList_.end())    constrainSumPt_ = true;
  //else  constrainSumPt_ = false;
  constrainSumPt_ = false;
  //constrainSumPt_ = true;

  }


void KinFitTestModule::setPull(int jet_variables_size, const TTbarGen * ttgen, int all_best, bool is_MC_){
  TMatrixD pull_had1_, pull_had2_, pull_hadb_, pull_lepb_, pull_lepton_, pull_neu_;
  TMatrixD pull_gtr_had1_, pull_gtr_had2_, pull_gtr_hadb_, pull_gtr_lepb_, pull_gtr_lepton_, pull_gtr_neu_;

  const TMatrixD* pull_ = jethad1_->getPull();
  pull_had1_.ResizeTo(jet_variables_size,1);
  for(int i=0; i<jet_variables_size; i++){
    pull_had1_(i,0) = (*pull_)(i,0);
  }

  const TMatrixD* pull_2 = jethad2_->getPull();
  pull_had2_.ResizeTo(jet_variables_size,1);
  for(int i=0; i<jet_variables_size; i++){
    pull_had2_(i,0) = (*pull_2)(i,0);
  }

  const TMatrixD* pull_3 = jethadb_->getPull();
  pull_hadb_.ResizeTo(jet_variables_size,1);
  for(int i=0; i<jet_variables_size; i++){
    pull_hadb_(i,0) = (*pull_3)(i,0);
  }

  const TMatrixD* pull_4 = jetlepb_->getPull();
  pull_lepb_.ResizeTo(jet_variables_size,1);
  for(int i=0; i<jet_variables_size; i++){
    pull_lepb_(i,0) = (*pull_4)(i,0);
  }
  const TMatrixD* pull_5 = lepton_->getPull();
  pull_lepton_.ResizeTo(3,1);
  for(int i=0; i<3; i++){
    pull_lepton_(i,0) = (*pull_5)(i,0);
  }

  const TMatrixD* pull_6 = met_->getPull();
  pull_neu_.ResizeTo(3,1);
  for(int i=0; i<3; i++){
    pull_neu_(i,0) = (*pull_6)(i,0);
  }
  
  if(all_best==0){
    pull_had1_best.ResizeTo(jet_variables_size, 1);
    pull_had2_best.ResizeTo(jet_variables_size,1);
    pull_hadb_best.ResizeTo(jet_variables_size,1);
    pull_lepb_best.ResizeTo(jet_variables_size,1);
    pull_lepton_best.ResizeTo(3,1);
    pull_neu_best.ResizeTo(3,1);

    pull_had1_best = pull_had1_;
    pull_had2_best = pull_had2_;
    pull_hadb_best = pull_hadb_;
    pull_lepb_best = pull_lepb_;
    pull_lepton_best = pull_lepton_;
    pull_neu_best = pull_neu_;
  }
  else if (all_best==1){
    pull_had1_all.ResizeTo(jet_variables_size, 1);
    pull_had2_all.ResizeTo(jet_variables_size,1);
    pull_hadb_all.ResizeTo(jet_variables_size,1);
    pull_lepb_all.ResizeTo(jet_variables_size,1);
    pull_lepton_all.ResizeTo(3,1);
    pull_neu_all.ResizeTo(3,1);

    pull_had1_all = pull_had1_;
    pull_had2_all = pull_had2_;
    pull_hadb_all = pull_hadb_;
    pull_lepb_all = pull_lepb_;
    pull_lepton_all = pull_lepton_;
    pull_neu_all = pull_neu_;
  }

  pull_gtr_had1_.ResizeTo(jet_variables_size,1);
  pull_gtr_had2_.ResizeTo(jet_variables_size,1);
  pull_gtr_hadb_.ResizeTo(jet_variables_size,1);
  pull_gtr_lepb_.ResizeTo(jet_variables_size,1);
  pull_gtr_lepton_.ResizeTo(3,1);
  pull_gtr_neu_.ResizeTo(3,1);
  
  if(is_MC_){
    if(ttgen->IsSemiLeptonicDecay()){
      switch(jetParam_){
      case KinFitTestSetup::kEMom :
	break;
      case KinFitTestSetup::kEtEtaPhi :{
      const TMatrixD* had1_sigmaf = jethad1_->getCovMatrixFit();
      float had1_et1 = (jethad1_->getCurr4Vec()->Et() - ttgen->Q1().v4().Et())/(*had1_sigmaf)(0,0);
      float had1_et2 = (jethad1_->getCurr4Vec()->Et() - ttgen->Q2().v4().Et())/(*had1_sigmaf)(0,0);
      float had1_eta1 = (jethad1_->getCurr4Vec()->Eta() - ttgen->Q1().v4().Eta())/(*had1_sigmaf)(1,1);
      float had1_eta2 = (jethad1_->getCurr4Vec()->Eta() - ttgen->Q2().v4().Eta())/(*had1_sigmaf)(1,1);
      float had1_phi1 = (jethad1_->getCurr4Vec()->Phi() - ttgen->Q1().v4().Phi())/(*had1_sigmaf)(2,2);
      float had1_phi2 = (jethad1_->getCurr4Vec()->Phi() - ttgen->Q2().v4().Phi())/(*had1_sigmaf)(2,2);
      const TMatrixD* had2_sigmaf = jethad2_->getCovMatrixFit();
      float had2_et1 = (jethad2_->getCurr4Vec()->Et() - ttgen->Q1().v4().Et())/(*had2_sigmaf)(0,0);
      float had2_et2 = (jethad2_->getCurr4Vec()->Et() - ttgen->Q2().v4().Et())/(*had2_sigmaf)(0,0);
      float had2_eta1 = (jethad2_->getCurr4Vec()->Eta() - ttgen->Q1().v4().Eta())/(*had2_sigmaf)(1,1);
      float had2_eta2 = (jethad2_->getCurr4Vec()->Eta() - ttgen->Q2().v4().Eta())/(*had2_sigmaf)(1,1);
      float had2_phi1 = (jethad2_->getCurr4Vec()->Phi() - ttgen->Q1().v4().Phi())/(*had2_sigmaf)(2,2);
      float had2_phi2 = (jethad2_->getCurr4Vec()->Phi() - ttgen->Q2().v4().Phi())/(*had2_sigmaf)(2,2);
      if(fabs(had1_et1)<fabs(had2_et1)){
	pull_gtr_had1_(0,0) = had1_et1;
	pull_gtr_had2_(0,0) = had2_et2;
      }else{
	pull_gtr_had1_(0,0) = had1_et2;
	pull_gtr_had2_(0,0) = had2_et1;
      }
      if(fabs(had1_eta1)<fabs(had2_eta1)){
	pull_gtr_had1_(1,0) = had1_eta1;
	pull_gtr_had2_(1,0) = had2_eta2;
      } else{
	pull_gtr_had1_(1,0) = had1_eta2;
	pull_gtr_had2_(1,0) = had2_eta1;
      }
      if(fabs(had1_phi1)<fabs(had2_phi1)){
	pull_gtr_had1_(2,0) = had1_phi1;
	pull_gtr_had2_(2,0) = had2_phi2;
      }else{
	pull_gtr_had1_(2,0) = had1_phi2;
	pull_gtr_had2_(2,0) = had2_phi1;
      }
      const TMatrixD* hadb_sigmaf = jethadb_->getCovMatrixFit();
      float hadb_et = (jethadb_->getCurr4Vec()->Et() - ttgen->BHad().v4().Et())/(*hadb_sigmaf)(0,0);
      float hadb_eta = (jethadb_->getCurr4Vec()->Eta() - ttgen->BHad().v4().Eta())/(*hadb_sigmaf)(1,1);
      float hadb_phi = (jethadb_->getCurr4Vec()->Phi() - ttgen->BHad().v4().Phi())/(*hadb_sigmaf)(2,2);
      pull_gtr_hadb_(0,0) = hadb_et;
      pull_gtr_hadb_(1,0) = hadb_eta;
      pull_gtr_hadb_(2,0) = hadb_phi;

      const TMatrixD* lepb_sigmaf = jetlepb_->getCovMatrixFit();
      float lepb_et = (jetlepb_->getCurr4Vec()->Et() - ttgen->BLep().v4().Et())/(*lepb_sigmaf)(0,0);
      float lepb_eta = (jetlepb_->getCurr4Vec()->Eta() - ttgen->BLep().v4().Eta())/(*lepb_sigmaf)(1,1);
      float lepb_phi = (jetlepb_->getCurr4Vec()->Phi() - ttgen->BLep().v4().Phi())/(*lepb_sigmaf)(2,2);
      pull_gtr_lepb_(0,0) = lepb_et;
      pull_gtr_lepb_(1,0) = lepb_eta;
      pull_gtr_lepb_(2,0) = lepb_phi;

      const TMatrixD* lepton_sigmaf = lepton_->getCovMatrixFit();
      if(debug){
	std::cout<<"lepton et fit "<<lepton_->getCurr4Vec()->Et()<<std::endl;
	std::cout<<"lepton et gen "<<ttgen->ChargedLepton().v4().Et()<<std::endl;
	std::cout<<"lepton et matrix "<<(*lepton_sigmaf)(0,0)<<std::endl;
	std::cout<<"lepton eta fit "<<lepton_->getCurr4Vec()->Eta()<<std::endl;
	std::cout<<"lepton eta gen "<<ttgen->ChargedLepton().v4().Eta()<<std::endl;
	std::cout<<"lepton eta matrix "<<(*lepton_sigmaf)(1,1)<<std::endl;
	std::cout<<"lepton phi fit "<<lepton_->getCurr4Vec()->Phi()<<std::endl;
	std::cout<<"lepton phi gen "<<ttgen->ChargedLepton().v4().Phi()<<std::endl;
	std::cout<<"lepton phi matrix "<<(*lepton_sigmaf)(2,2)<<std::endl;
	lepton_sigmaf->Print();
      }
      float lepton_et = (lepton_->getCurr4Vec()->Et() - ttgen->ChargedLepton().v4().Et())/(*lepton_sigmaf)(0,0);
      float lepton_eta = (lepton_->getCurr4Vec()->Eta() - ttgen->ChargedLepton().v4().Eta())/(*lepton_sigmaf)(1,1);
      float lepton_phi = (lepton_->getCurr4Vec()->Phi() - ttgen->ChargedLepton().v4().Phi())/(*lepton_sigmaf)(2,2);
      pull_gtr_lepton_(0,0) = lepton_et;
      pull_gtr_lepton_(1,0) = lepton_eta;
      pull_gtr_lepton_(2,0) = lepton_phi;
      if(debug){
	std::cout<<"lepton pull gtr et "<<lepton_et<<std::endl;
	std::cout<<"lepton pull gtr eta "<<lepton_eta<<std::endl;
	std::cout<<"lepton pull gtr phi "<<lepton_phi<<std::endl;
      }
      const TMatrixD* met_sigmaf = met_->getCovMatrixFit();
      float met_et = (met_->getCurr4Vec()->Et() - ttgen->Neutrino().v4().Et())/(*met_sigmaf)(0,0);
      float met_eta = (met_->getCurr4Vec()->Eta() - ttgen->Neutrino().v4().Eta())/(*met_sigmaf)(1,1);
      float met_phi = (met_->getCurr4Vec()->Phi() - ttgen->Neutrino().v4().Phi())/(*met_sigmaf)(2,2);
      pull_gtr_neu_(0,0) = met_et;
      pull_gtr_neu_(1,0) = met_eta;
      pull_gtr_neu_(2,0) = met_phi;

      if(all_best==0){
	pull_gtr_had1_best.ResizeTo(jet_variables_size,1);
	pull_gtr_had2_best.ResizeTo(jet_variables_size,1);
	pull_gtr_hadb_best.ResizeTo(jet_variables_size,1);
	pull_gtr_lepb_best.ResizeTo(jet_variables_size,1);
	pull_gtr_lepton_best.ResizeTo(3,1);
	pull_gtr_neu_best.ResizeTo(3,1);
	
	pull_gtr_had1_best = pull_gtr_had1_;
	pull_gtr_had2_best = pull_gtr_had2_;
	pull_gtr_hadb_best = pull_gtr_hadb_;
	pull_gtr_lepb_best = pull_gtr_lepb_;
	pull_gtr_lepton_best = pull_gtr_lepton_;
	pull_gtr_neu_best = pull_gtr_neu_;
      }
      break;}
      }}}else{
      if(all_best==0){
	pull_gtr_had1_best.Zero();
	pull_gtr_had2_best.Zero();
	pull_gtr_hadb_best.Zero();
	pull_gtr_lepb_best.Zero();
	pull_gtr_lepton_best.Zero();
	pull_gtr_neu_best.Zero();
      }
    }
}

  bool KinFitTestModule::dm_check(Event & event, int bl, int bh, int h1, int h2){
    double dm_w_had_check(-1.), dm_w_lep_check(-1.), dm_top_had_check(-1.), dm_top_lep_check(-1.);
    bool check_dm;
    dm_w_had_check = (event.jets->at(h1).v4()+event.jets->at(h2).v4()).M() - mW_;
    dm_top_had_check = (event.jets->at(h1).v4()+event.jets->at(h2).v4()+event.jets->at(bh).v4()).M() - mTop_;
   
    switch(leptonType){
    case KinFitTestSetup::kElectron :      
      dm_w_lep_check = (event.electrons->at(0).v4()+event.met->v4()).M() - mW_;
      dm_top_lep_check = (event.electrons->at(0).v4()+event.met->v4()+event.jets->at(bl).v4()).M() - mTop_;
      break;
    case KinFitTestSetup::kMuon :          
      dm_w_lep_check = (event.muons->at(0).v4()+event.met->v4()).M() - mW_;
      dm_top_lep_check = (event.muons->at(0).v4()+event.met->v4()+event.jets->at(bl).v4()).M() - mTop_;
      break;
    case KinFitTestSetup::kMet : break;
    case KinFitTestSetup::kBJet : break;
    case KinFitTestSetup::kUdscJet : break;
    }

    //    check_dm = (fabs(dm_w_had_check)<35. && fabs(dm_w_lep_check)<70. && fabs(dm_top_had_check)<50. && fabs(dm_top_lep_check)<100.);
    check_dm = (fabs(dm_w_had_check)<35. && fabs(dm_top_had_check)<50. );

    if(debug){
      std::cout<<"////////// check"<<std::endl;
      std::cout<<"dm wlep  "<<dm_w_lep_check<<std::endl;
      std::cout<<"dm whad  "<<dm_w_had_check<<std::endl;
      std::cout<<"dm top lep  "<<dm_top_lep_check<<std::endl;
      std::cout<<"dm top had  "<<dm_top_had_check<<std::endl;
      std::cout<<"dm check "<<check_dm<<std::endl;
    }
    return check_dm;
  }


UHH2_REGISTER_ANALYSIS_MODULE(KinFitTestModule)
}
