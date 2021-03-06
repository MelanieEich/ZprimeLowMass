#include "TMatrixD.h"

#include <UHH2/core/include/Utils.h>

#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/KinFitTest/include/HistsBASE.h>

#include <TLorentzVector.h>

namespace uhh2examples {


class KinFitTestHistsPulls: public HistsBASE {


 public:
  KinFitTestHistsPulls(uhh2::Context&, const std::string&, const std::string& ttgen="ttbargen", const std::string& hyps="TTbarReconstruction", const std::string& disc="Chi2");
  virtual ~KinFitTestHistsPulls();

  virtual void fill(const uhh2::Event&) override;
  virtual void fill(const uhh2::Event&, const ReconstructionHypothesis* hyp, const float, const TTbarGen* ttgen=0, const double wgt=1.);
  float deltaR_fit_new(const float eta1, const float eta2, const float phi1, const float phi2);

  protected:
  uhh2::Event::Handle<TTbarGen> h_ttbar_gen;
  uhh2::Event::Handle<std::vector<ReconstructionHypothesis>> h_ttbar_hyps;
  std::string disc_name_;

  virtual void init() override;

  uhh2::Event::Handle<int> parametrisation_;

  uhh2::Event::Handle<TMatrixD> pull_had1;
  uhh2::Event::Handle<TMatrixD> pull_had2;
  uhh2::Event::Handle<TMatrixD> pull_hadb;
  uhh2::Event::Handle<TMatrixD> pull_lepb;
  uhh2::Event::Handle<TMatrixD> pull_lepton;
  uhh2::Event::Handle<TMatrixD> pull_neu;

  uhh2::Event::Handle<float> pull_gtr_had1_et;
  uhh2::Event::Handle<float> pull_gtr_had1_eta;
  uhh2::Event::Handle<float> pull_gtr_had1_phi;
  uhh2::Event::Handle<float> pull_gtr_had2_et;
  uhh2::Event::Handle<float> pull_gtr_had2_eta;
  uhh2::Event::Handle<float> pull_gtr_had2_phi;
  uhh2::Event::Handle<float> pull_gtr_hadb_et;
  uhh2::Event::Handle<float> pull_gtr_hadb_eta;
  uhh2::Event::Handle<float> pull_gtr_hadb_phi;
  uhh2::Event::Handle<float> pull_gtr_lepb_et;
  uhh2::Event::Handle<float> pull_gtr_lepb_eta;
  uhh2::Event::Handle<float> pull_gtr_lepb_phi;
  uhh2::Event::Handle<float> pull_gtr_lepton_et;
  uhh2::Event::Handle<float> pull_gtr_lepton_eta;
  uhh2::Event::Handle<float> pull_gtr_lepton_phi;
  uhh2::Event::Handle<float> pull_gtr_neu_et;
  uhh2::Event::Handle<float> pull_gtr_neu_eta;
  uhh2::Event::Handle<float> pull_gtr_neu_phi;

  float delta_phi(const float, const float); 

};
}
