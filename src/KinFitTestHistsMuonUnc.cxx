#include "UHH2/KinFitTest/include/KinFitTestHistsMuonUnc.h"
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

#include <iostream>


using namespace std;
using namespace uhh2;
using namespace uhh2examples;


KinFitTestHistsMuonUnc::KinFitTestHistsMuonUnc(uhh2::Context& ctx, const std::string& dir, const std::string& ttgen, const std::string& hyps, const std::string& disc):
  HistsBASE(ctx, dir) {

  h_ttbar_gen  = ctx.get_handle<TTbarGen>(ttgen);
  h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis> >(hyps);
  disc_name_ = disc;
  input_jetq1_v4 = ctx.get_handle<TLorentzVector>("input_jetq1_v4");
  input_jetq2_v4 = ctx.get_handle<TLorentzVector>("input_jetq2_v4");
  input_jetbh_v4 = ctx.get_handle<TLorentzVector>("input_jetbh_v4");
  input_jetbl_v4 = ctx.get_handle<TLorentzVector>("input_jetbl_v4");
  input_lepton_v4 = ctx.get_handle<TLorentzVector>("input_lepton_v4");
  input_neutrino_v4 = ctx.get_handle<TLorentzVector>("input_neutrino_v4");
  unc_had1 = ctx.get_handle<TMatrixD>("unc_had1");
  unc_had2 = ctx.get_handle<TMatrixD>("unc_had2");
  unc_hadb = ctx.get_handle<TMatrixD>("unc_hadb");
  unc_lepb = ctx.get_handle<TMatrixD>("unc_lepb");
  unc_lepton = ctx.get_handle<TMatrixD>("unc_lepton");
  unc_neu = ctx.get_handle<TMatrixD>("unc_neu");
  init();
}


void KinFitTestHistsMuonUnc::init(){

  book_TH1F("eta_muon_eta_bin_1_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_1_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_1_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_1_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_1_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_2_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_2_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_2_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_2_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_2_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_3_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_3_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_3_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_3_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_3_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_4_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_4_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_4_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_4_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_4_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_5_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_5_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_5_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_5_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_5_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_6_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_6_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_6_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_6_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_6_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_7_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_7_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_7_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_7_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_7_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_8_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_8_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_8_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_8_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_8_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_9_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_9_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_9_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_9_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_9_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_10_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_10_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_10_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_10_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_10_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_11_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_11_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_11_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_11_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_11_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_12_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_12_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_12_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_12_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_12_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_13_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_13_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_13_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_13_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_13_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_14_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_14_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_14_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_14_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_14_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_15_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_15_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_15_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_15_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_15_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_16_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_16_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_16_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_16_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_16_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_17_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_17_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_17_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_17_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_17_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_18_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_18_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_18_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_18_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_18_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_19_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_19_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_19_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_19_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_19_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_20_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_20_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_20_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_20_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_20_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_21_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_21_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_21_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_21_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_21_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_22_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_22_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_22_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_22_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_22_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_23_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_23_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_23_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_23_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_23_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_24_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_24_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_24_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_24_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("eta_muon_eta_bin_24_pt_bin_5"    , 1000,-0.005,0.005);

  book_TH1F("phi_muon_eta_bin_1_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_1_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_1_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_1_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_1_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_2_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_2_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_2_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_2_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_2_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_3_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_3_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_3_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_3_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_3_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_4_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_4_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_4_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_4_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_4_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_5_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_5_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_5_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_5_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_5_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_6_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_6_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_6_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_6_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_6_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_7_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_7_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_7_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_7_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_7_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_8_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_8_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_8_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_8_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_8_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_9_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_9_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_9_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_9_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_9_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_10_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_10_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_10_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_10_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_10_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_11_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_11_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_11_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_11_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_11_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_12_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_12_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_12_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_12_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_12_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_13_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_13_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_13_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_13_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_13_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_14_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_14_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_14_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_14_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_14_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_15_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_15_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_15_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_15_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_15_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_16_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_16_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_16_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_16_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_16_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_17_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_17_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_17_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_17_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_17_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_18_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_18_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_18_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_18_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_18_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_19_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_19_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_19_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_19_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_19_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_20_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_20_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_20_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_20_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_20_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_21_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_21_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_21_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_21_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_21_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_22_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_22_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_22_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_22_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_22_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_23_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_23_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_23_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_23_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_23_pt_bin_5"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_24_pt_bin_1"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_24_pt_bin_2"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_24_pt_bin_3"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_24_pt_bin_4"    , 1000,-0.005,0.005);
  book_TH1F("phi_muon_eta_bin_24_pt_bin_5"    , 1000,-0.005,0.005);

  book_TH1F("pt_neu_pt_bin_1"    , 1000,-50,50);
  book_TH1F("pt_neu_pt_bin_2"    , 1000,-50,50);
  book_TH1F("pt_neu_pt_bin_3"    , 1000,-50,50);
  book_TH1F("pt_neu_pt_bin_4"    , 1000,-50,50);
  book_TH1F("pt_neu_pt_bin_5"    , 1000,-50,50);
  book_TH1F("eta_neu_pt_bin_1"    , 1000,-0.5,0.5);
  book_TH1F("eta_neu_pt_bin_2"    , 1000,-0.5,0.5);
  book_TH1F("eta_neu_pt_bin_3"    , 1000,-0.5,0.5);
  book_TH1F("eta_neu_pt_bin_4"    , 1000,-0.5,0.5);
  book_TH1F("eta_neu_pt_bin_5"    , 1000,-0.5,0.5);
  book_TH1F("phi_neu_pt_bin_1"    , 1000,-0.5,0.5);
  book_TH1F("phi_neu_pt_bin_2"    , 1000,-0.5,0.5);
  book_TH1F("phi_neu_pt_bin_3"    , 1000,-0.5,0.5);
  book_TH1F("phi_neu_pt_bin_4"    , 1000,-0.5,0.5);
  book_TH1F("phi_neu_pt_bin_5"    , 1000,-0.5,0.5);

}


void KinFitTestHistsMuonUnc::fill(const uhh2::Event& event){

  const float weight(event.weight);
  const TTbarGen* ttgen(0);


  if(!event.isRealData){
    const auto& ttbargen = event.get(h_ttbar_gen);
    ttgen = &ttbargen;
   }

 const ReconstructionHypothesis* hyp(0);
  float hyp_val(0.);
  if(event.is_valid(h_ttbar_hyps)){

    const auto& ttbar_hyps = event.get(h_ttbar_hyps);

    hyp = get_best_hypothesis(ttbar_hyps, disc_name_);
    if(!hyp) throw std::runtime_error("EffyTTbarHists::fill -- null pointer to ReconstructionHypothesis object (from \"get_best_hypothesis\")");

    hyp_val = hyp->discriminator(disc_name_);
  }



  fill(event, hyp, hyp_val, ttgen, weight);



  return;
}

void  KinFitTestHistsMuonUnc::fill(const uhh2::Event& event, const ReconstructionHypothesis* hyp, const float hyp_val, const TTbarGen* ttgen, const double weight){

   bool ttljets(false);

 if(ttgen){
   ttljets = (ttgen->DecayChannel() == TTbarGen::e_muhad || ttgen->DecayChannel() == TTbarGen::e_ehad);
  }

 TLorentzVector input_jetq1 = event.get(input_jetq1_v4);
 TLorentzVector input_jetq2 = event.get(input_jetq2_v4);
 TLorentzVector input_jetbh = event.get(input_jetbh_v4);
 TLorentzVector input_jetbl = event.get(input_jetbl_v4);
 TLorentzVector input_lepton = event.get(input_lepton_v4);
 TLorentzVector input_neutrino = event.get(input_neutrino_v4);
 TLorentzVector input_wlep = (input_lepton + input_neutrino);
 TLorentzVector input_whad = (input_jetq1 + input_jetq2);
 TLorentzVector input_tlep = (input_lepton + input_neutrino + input_jetbl);
 TLorentzVector input_thad = (input_jetq1 + input_jetq2 + input_jetbh);
 TLorentzVector input_ttbar = (input_tlep + input_thad);


 if(ttljets){
   float dR_gen_meas = deltaR_fit_new(ttgen->ChargedLepton().eta(), input_lepton.Eta(), ttgen->ChargedLepton().phi(), input_lepton.Phi());
   if(dR_gen_meas<0.3){
     if(0.000<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.100){
       if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_1_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_1_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	 H1("eta_muon_eta_bin_1_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_1_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
       else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	 H1("eta_muon_eta_bin_1_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_1_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	 H1("eta_muon_eta_bin_1_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_1_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	 H1("eta_muon_eta_bin_1_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_1_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
     }
     else if(0.100<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.200){
       if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	 H1("eta_muon_eta_bin_2_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_2_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	 H1("eta_muon_eta_bin_2_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_2_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	 H1("eta_muon_eta_bin_2_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_2_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	 H1("eta_muon_eta_bin_2_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_2_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
       else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	 H1("eta_muon_eta_bin_2_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	 H1("phi_muon_eta_bin_2_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
       }
     }
	else if(0.200<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.300){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_3_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_3_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_3_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_3_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_3_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_3_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_3_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_3_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_3_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_3_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.300<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.400){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_4_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_4_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_4_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_4_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_4_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_4_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_4_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_4_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_4_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_4_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.400<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.500){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_5_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_5_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_5_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_5_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_5_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_5_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_5_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_5_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_5_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_5_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.500<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.600){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_6_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_6_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_6_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_6_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_6_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_6_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_6_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_6_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_6_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_6_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.600<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.700){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_7_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_7_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_7_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_7_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_7_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_7_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_7_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_7_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_7_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_7_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.700<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.800){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_8_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_8_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_8_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_8_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_8_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_8_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_8_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_8_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_8_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_8_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.800<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<0.900){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_9_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_9_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_9_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_9_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_9_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_9_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_9_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_9_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_9_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_9_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(0.900<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.000){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_10_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_10_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_10_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_10_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_10_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_10_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_10_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_10_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_10_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_10_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.000<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.100){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_11_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_11_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_11_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_11_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_11_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_11_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_11_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_11_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_11_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_11_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.100<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.200){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_12_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_12_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_12_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_12_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_12_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_12_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_12_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_12_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_12_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_12_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.200<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.300){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_13_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_13_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_13_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_13_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_13_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_13_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_13_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_13_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_13_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_13_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.300<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.400){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_14_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_14_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_14_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_14_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_14_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_14_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_14_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_14_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_14_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_14_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.400<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.500){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_15_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_15_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_15_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_15_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_15_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_15_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_15_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_15_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_15_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_15_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.500<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.600){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_16_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_16_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_16_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_16_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_16_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_16_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_16_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_16_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_16_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_16_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.600<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.700){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_17_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_17_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_17_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_17_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_17_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_17_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_17_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_17_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_17_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_17_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.700<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.800){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_18_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_18_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_18_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_18_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_18_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_18_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_18_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_18_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_18_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_18_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.800<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<1.900){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_19_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_19_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_19_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_19_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_19_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_19_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_19_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_19_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_19_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_19_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(1.900<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<2.000){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_20_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_20_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_20_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_20_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_20_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_20_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_20_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_20_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_20_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_20_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(2.000<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<2.100){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_21_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_21_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_21_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_21_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_21_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_21_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_21_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_21_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_21_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_21_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(2.100<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<2.200){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_22_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_22_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_22_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_22_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_22_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_22_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_22_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_22_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_22_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_22_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(2.200<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<2.300){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_23_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_23_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_23_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_23_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_23_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_23_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_23_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_23_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_23_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_23_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}
	else if(2.300<=fabs(input_lepton.Eta()) && fabs(input_lepton.Eta())<2.400){
	  if(30<=input_lepton.Pt() && input_lepton.Pt()<40){
	    H1("eta_muon_eta_bin_24_pt_bin_1")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_24_pt_bin_1")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(50<=input_lepton.Pt() && input_lepton.Pt()<60){
	    H1("eta_muon_eta_bin_24_pt_bin_2")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_24_pt_bin_2")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(70<=input_lepton.Pt() && input_lepton.Pt()<80){
	    H1("eta_muon_eta_bin_24_pt_bin_3")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_24_pt_bin_3")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(90<=input_lepton.Pt() && input_lepton.Pt()<100){
	    H1("eta_muon_eta_bin_24_pt_bin_4")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_24_pt_bin_4")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	  else if(110<=input_lepton.Pt() && input_lepton.Pt()<120){
	    H1("eta_muon_eta_bin_24_pt_bin_5")->Fill((input_lepton.Eta() - ttgen->ChargedLepton().eta()));
	    H1("phi_muon_eta_bin_24_pt_bin_5")->Fill((input_lepton.Phi() - ttgen->ChargedLepton().phi()));
	  }
	}    
      }

      float dR_gen_meas_neu = deltaR_fit_new(ttgen->Neutrino().eta(), input_neutrino.Eta(), ttgen->Neutrino().phi(), input_neutrino.Phi());
      if(dR_gen_meas_neu<0.3){
	if(30<=input_neutrino.Pt() && input_neutrino.Pt()<40){   
	  H1("pt_neu_pt_bin_1")->Fill((input_neutrino.Pt() - ttgen->Neutrino().pt() ));
	  H1("eta_neu_pt_bin_1")->Fill((input_neutrino.Eta() - ttgen->Neutrino().eta()));
	  H1("phi_neu_pt_bin_1")->Fill((input_neutrino.Phi() - ttgen->Neutrino().phi()));
	}
	if(50<=input_neutrino.Pt() && input_neutrino.Pt()<60){   	
	  H1("pt_neu_pt_bin_2")->Fill((input_neutrino.Pt() - ttgen->Neutrino().pt()));
	  H1("eta_neu_pt_bin_2")->Fill((input_neutrino.Eta() - ttgen->Neutrino().eta()));
	  H1("phi_neu_pt_bin_2")->Fill((input_neutrino.Phi() - ttgen->Neutrino().phi()));
	}
	if(70<=input_neutrino.Pt() && input_neutrino.Pt()<80){   	
	  H1("pt_neu_pt_bin_3")->Fill((input_neutrino.Pt() - ttgen->Neutrino().pt()));
	  H1("eta_neu_pt_bin_3")->Fill((input_neutrino.Eta() - ttgen->Neutrino().eta()));
	  H1("phi_neu_pt_bin_3")->Fill((input_neutrino.Phi() - ttgen->Neutrino().phi()));
	}
	if(90<=input_neutrino.Pt() && input_neutrino.Pt()<100){ 	        
	  H1("pt_neu_pt_bin_4")->Fill((input_neutrino.Pt() - ttgen->Neutrino().pt()));
	  H1("eta_neu_pt_bin_4")->Fill((input_neutrino.Eta() - ttgen->Neutrino().eta()));
	  H1("phi_neu_pt_bin_4")->Fill((input_neutrino.Phi() - ttgen->Neutrino().phi()));
	}
	if(110<=input_neutrino.Pt() && input_neutrino.Pt()<120){ 	
	  H1("pt_neu_pt_bin_5")->Fill((input_neutrino.Pt() - ttgen->Neutrino().pt()));
	  H1("eta_neu_pt_bin_5")->Fill((input_neutrino.Eta() - ttgen->Neutrino().eta()));
	  H1("phi_neu_pt_bin_5")->Fill((input_neutrino.Phi() - ttgen->Neutrino().phi()));
	}
      }
    }
}


float KinFitTestHistsMuonUnc::deltaR_fit_new(const float eta1, const float eta2, const float phi1, const float phi2){
    float deltaeta = eta1 - eta2;
    float dphi = delta_phi(phi1, phi2);
    return sqrt(deltaeta * deltaeta + dphi * dphi);
}

float KinFitTestHistsMuonUnc::delta_phi(const float phi1, const float phi2){

  float a_dphi = fabs(phi1-phi2);
  if(a_dphi > M_PI) a_dphi = 2*M_PI - a_dphi;

  return a_dphi;
}



KinFitTestHistsMuonUnc::~KinFitTestHistsMuonUnc(){}

