//Resolution of the Jets to set up covariance matrix
class HelperJet {
  
 public:
  enum Flavor {kUds, kB};
  
  HelperJet(){};
  ~HelperJet(){};
  
  inline double pt (double pt, double eta,  Flavor flav);
  inline double eta(double pt, double eta,  Flavor flav);
  inline double phi(double pt, double eta,  Flavor flav);
  inline double a(double pt, double eta,  Flavor flav);
  inline double b(double pt, double eta,  Flavor flav);
  inline double c(double pt, double eta,  Flavor flav);
  inline double d(double pt, double eta,  Flavor flav);
  inline double et(double pt, double eta,  Flavor flav);
  inline double theta(double pt, double eta,  Flavor flav);
};

inline double HelperJet::pt(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);
  if(fabs(eta)<1.4) res+= 6.68;
  else res+=-3.14*fabs(eta)+11.89;
  if(flav==kB){
    res=0.33*sqrt(pt);
    if(fabs(eta)<1.4) res+= 6.57;
    else res+=-1.09*fabs(eta)+8.50;    
  }
  return res;
}

inline double HelperJet::a(double pt, double eta,  Flavor flav)
{

  double res = 0.29*sqrt(pt);
  if(flav==kB){
    if(fabs(eta)<0.17) res = 0.142188+0.288157*exp(-(0.0275628*pt));
    else if(fabs(eta)<0.35) res = 0.141215+0.277502*exp(-(0.0282327*pt));
    else if(fabs(eta)<0.5) res = 0.140493+0.26136*exp(-(0.0254784*pt));
    else if(fabs(eta)<0.7) res = 0.142608+0.371256*exp(-(0.0299132*pt));
    else if(fabs(eta)<0.9) res = 0.138817+0.237167*exp(-(0.0268413*pt));
    else if(fabs(eta)<1.15) res = 0.137322+0.293834*exp(-(0.0264771*pt));
    else if(fabs(eta)<1.4) res = 0.147184+0.40834*exp(-(0.0366769*pt));
    else if(fabs(eta)<1.7) res = 0.141529+0.195563*exp(-(0.023933*pt));
    else if(fabs(eta)<2.1) res = 0.129132+0.261012*exp(-(0.0323658*pt));
    else if(fabs(eta)<2.5) res = 0.120065+0.239351*exp(-(0.0356176*pt));
  }else{
    if(fabs(eta)<0.17) res = 0.142145+0.274541*exp(-(0.0286487*pt));
    else if(fabs(eta)<0.35) res = 0.146112+0.284777*exp(-(0.031039*pt));
    else if(fabs(eta)<0.5) res = 0.146703+0.343982*exp(-(0.0333917*pt));
    else if(fabs(eta)<0.7) res = 0.149067+0.241197*exp(-(0.0303415*pt));
    else if(fabs(eta)<0.9) res = 0.141233+0.304052*exp(-(0.0314644*pt));
    else if(fabs(eta)<1.15) res = 0.144531+0.266871*exp(-(0.0283385*pt));
    else if(fabs(eta)<1.4) res = 0.156593+0.41251*exp(-(0.0411896*pt));
    else if(fabs(eta)<1.7) res = 0.13452+0.269*exp(-(0.0230171*pt));
    else if(fabs(eta)<2.1) res = 0.0969416+0.211636*exp(-(0.0166117*pt));
    else if(fabs(eta)<2.5) res = 0.0727398+0.205132*exp(-(0.01661*pt));
  }
  return res;
}

inline double HelperJet::b(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);
  
  if(flav==kB){
    if(fabs(eta)<0.17) res = 2.02521+0.174301*exp(-(-0.0136646*pt));
    else if(fabs(eta)<0.35) res = 0.690872+1.32174*exp(-(-0.00459989*pt));
    else if(fabs(eta)<0.5) res = -956.522+958.486*exp(-(-7.57551e-06*pt));
    else if(fabs(eta)<0.7) res = 2.01017+0.0670861*exp(-(-0.0194879*pt));
    else if(fabs(eta)<0.9) res = 1.92701+0.124384*exp(-(-0.0144246*pt));
    else if(fabs(eta)<1.15) res = 1.94854+0.034138*exp(-(-0.0227995*pt));
    else if(fabs(eta)<1.4) res = 2.07599+0.0208406*exp(-(-0.0214439*pt));
    else if(fabs(eta)<1.7) res = 1.99002+0.0077067*exp(-(-0.0278379*pt));
    else if(fabs(eta)<2.1) res = 1.93181+0.00196985*exp(-(-0.0390155*pt));
    else if(fabs(eta)<2.5) res = 1.75611+0.0287099*exp(-(-0.0202162*pt));
  }else{
    if(fabs(eta)<0.17) res = -644.18+645.342*exp(-(-1.7462e-05*pt));
    else if(fabs(eta)<0.35) res = -1.34986+2.48937*exp(-(-0.00396448*pt));
    else if(fabs(eta)<0.5) res = -4.8078+5.98195*exp(-(-0.0015*pt));
    else if(fabs(eta)<0.7) res = -9.55949+10.7419*exp(-(-0.000891802*pt));
    else if(fabs(eta)<0.9) res = 0.374264+0.884391*exp(-(-0.00642501*pt));
    else if(fabs(eta)<1.15) res = -9.51561+10.7429*exp(-(-0.000639173*pt));
    else if(fabs(eta)<1.4) res = -0.158886+1.43841*exp(-(-0.00347282*pt));
    else if(fabs(eta)<1.7) res = 0.28876+0.992659*exp(-(-0.00454673*pt));
    else if(fabs(eta)<2.1) res = 0.292201+0.922609*exp(-(-0.00536289*pt));
    else if(fabs(eta)<2.5) res = -657.479+658.533*exp(-(-1.40481e-05*pt));
  }
  
  return res;
}

inline double HelperJet::c(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);
  
  if(flav==kB){
    if(fabs(eta)<0.17) res = -2048.11+2051.27*exp(-(1.81337e-06*pt));
    else if(fabs(eta)<0.35) res = 2.86314+0.598484*exp(-(0.0322582*pt));
    else if(fabs(eta)<0.5) res = 2.54785+0.77807*exp(-(0.0100979*pt));
    else if(fabs(eta)<0.7) res = 2.78865+6.82558*exp(-(0.0677953*pt));
    else if(fabs(eta)<0.9) res = 2.63261+2.66562*exp(-(0.044099*pt));
    else if(fabs(eta)<1.15) res = 2.51691+2.75456*exp(-(0.0425201*pt));
    else if(fabs(eta)<1.4) res = 2.60458+0.623776*exp(-(0.0148415*pt));
    else if(fabs(eta)<1.7) res = 2.39005+1.60119*exp(-(0.0154536*pt));
    else if(fabs(eta)<2.1) res = 2.49783+0.408896*exp(-(0.04005*pt));
    else if(fabs(eta)<2.5) res = 2.02143+0.000354016*exp(-(-0.0476218*pt));
  }else{
    if(fabs(eta)<0.17) res = 1.72736+0.212526*exp(-(-0.00672611*pt));
    else if(fabs(eta)<0.35) res = -3480.32+3482.17*exp(-(-7.08939e-07*pt));
    else if(fabs(eta)<0.5) res = 1.89173+0.0381423*exp(-(-0.0151248*pt));
    else if(fabs(eta)<0.7) res = 2+0.00575323*exp(-(-0.0233459*pt));
    else if(fabs(eta)<0.9) res = 1.92618+0.00541047*exp(-(-0.0269683*pt));
    else if(fabs(eta)<1.15) res = 1.95535+0.0194328*exp(-(-0.0175058*pt));
    else if(fabs(eta)<1.4) res = 1.99855+0.00403398*exp(-(-0.027822*pt));
    else if(fabs(eta)<1.7) res = 2.02356+1.33745e-15*exp(-(-0.201028*pt));
    else if(fabs(eta)<2.1) res = 1.70305+0.0428122*exp(-(-0.0206051*pt));
    else if(fabs(eta)<2.5) res = 0.964128+0.370407*exp(-(-0.0119056*pt));
  }
  
  return res;
}

inline double HelperJet::d(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);

  if(flav==kB){
    if(fabs(eta)<0.17) res = 0.144545+0.374307*exp(-(0.0314642*pt));
    else if(fabs(eta)<0.35) res = 0.141366+0.241457*exp(-(0.0263273*pt));
    else if(fabs(eta)<0.5) res = 0.142509+0.324478*exp(-(0.0288208*pt));
    else if(fabs(eta)<0.7) res = 0.139809+0.292392*exp(-(0.0257797*pt));
    else if(fabs(eta)<0.9) res = 0.137762+0.258213*exp(-(0.0265974*pt));
    else if(fabs(eta)<1.15) res = 0.135831+0.338346*exp(-(0.0277574*pt));
    else if(fabs(eta)<1.4) res = 0.148253+0.53999*exp(-(0.0415757*pt));
    else if(fabs(eta)<1.7) res = 0.145785+0.297049*exp(-(0.0319928*pt));
    else if(fabs(eta)<2.1) res = 0.126895+0.291297*exp(-(0.0328976*pt));
    else if(fabs(eta)<2.5) res = 0.119259+0.228403*exp(-(0.0341018*pt));
  }else{
    if(fabs(eta)<0.17) res = 0.139738+0.310387*exp(-(0.0296992*pt));
    else if(fabs(eta)<0.35) res = 0.142229+0.274556*exp(-(0.0286575*pt));
    else if(fabs(eta)<0.5) res = 0.14344+0.313425*exp(-(0.0294567*pt));
    else if(fabs(eta)<0.7) res = 0.147284+0.252005*exp(-(0.0303533*pt));
    else if(fabs(eta)<0.9) res = 0.13541+0.261376*exp(-(0.0260026*pt));
    else if(fabs(eta)<1.15) res = 0.144023+0.27422*exp(-(0.0282282*pt));
    else if(fabs(eta)<1.4) res = 0.158694+0.352571*exp(-(0.0385597*pt));
    else if(fabs(eta)<1.7) res = 0.135936+0.24773*exp(-(0.0226282*pt));
    else if(fabs(eta)<2.1) res = 0.0957325+0.213108*exp(-(0.016424*pt));
    else if(fabs(eta)<2.5) res = 0.0622281+0.185045*exp(-(0.0129321*pt));
  } 
  

  return res;
}
	
inline double HelperJet::theta(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);
  if(flav==kB){ 
    if(fabs(eta)<0.17) res = 0.0348657+0.0966067*exp(-(0.0388496*pt));
    else if(fabs(eta)<0.35) res = 0.0327521+0.0901303*exp(-(0.0354149*pt));
    else if(fabs(eta)<0.5) res = 0.030283+0.0970732*exp(-(0.036806*pt));
    else if(fabs(eta)<0.7) res = 0.0266898+0.0994089*exp(-(0.0379773*pt));
    else if(fabs(eta)<0.9) res = 0.022614+0.091974*exp(-(0.0384986*pt));
    else if(fabs(eta)<1.15) res = 0.0173219+0.0923928*exp(-(0.0395886*pt));
    else if(fabs(eta)<1.4) res = 0.0129851+0.077473*exp(-(0.0363236*pt));
    else if(fabs(eta)<1.7) res = 0.00967747+0.0667463*exp(-(0.0383101*pt));
    else if(fabs(eta)<2.1) res = 0.00744593+0.050027*exp(-(0.0432323*pt));
    else if(fabs(eta)<2.5) res = 0.00589309+0.0346904*exp(-(0.0484705*pt));
  }else{
    if(fabs(eta)<0.17) res = 0.0384999+0.105382*exp(-(0.0440247*pt));
    else if(fabs(eta)<0.35) res = 0.0384758+0.086093*exp(-(0.0436464*pt));
    else if(fabs(eta)<0.5) res = 0.0325808+0.0860212*exp(-(0.0371074*pt));
    else if(fabs(eta)<0.7) res = 0.0304022+0.094116*exp(-(0.0402575*pt));
    else if(fabs(eta)<0.9) res = 0.0246649+0.0883167*exp(-(0.0386519*pt));
    else if(fabs(eta)<1.15) res = 0.0197668+0.080217*exp(-(0.0368397*pt));
    else if(fabs(eta)<1.4) res = 0.0148531+0.0692569*exp(-(0.0354383*pt));
    else if(fabs(eta)<1.7) res = 0.0113007+0.0601861*exp(-(0.0373138*pt));
    else if(fabs(eta)<2.1) res = 0.00860103+0.0522148*exp(-(0.0450405*pt));
    else if(fabs(eta)<2.5) res = 0.00692992+0.031671*exp(-(0.0490789*pt));
  }
  return res;
}

inline double HelperJet::phi(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);
  if(flav==kB){ 
    if(fabs(eta)<0.17) res = 0.0217379+0.157228*exp(-(0.0237049*pt));
    else if(fabs(eta)<0.35) res = 0.0246735+0.165694*exp(-(0.0269728*pt));
    else if(fabs(eta)<0.5) res = 0.0226496+0.162875*exp(-(0.0250225*pt));
    else if(fabs(eta)<0.7) res = 0.0232004+0.167887*exp(-(0.0259427*pt));
    else if(fabs(eta)<0.9) res = 0.0220083+0.183889*exp(-(0.0273017*pt));
    else if(fabs(eta)<1.15) res = 0.0216601+0.18482*exp(-(0.027148*pt));
    else if(fabs(eta)<1.4) res = 0.0225191+0.19558*exp(-(0.0271747*pt));
    else if(fabs(eta)<1.7) res = 0.0208365+0.19241*exp(-(0.0262627*pt));
    else if(fabs(eta)<2.1) res = 0.0233749+0.180684*exp(-(0.032052*pt));
    else if(fabs(eta)<2.5) res = 0.0231891+0.165484*exp(-(0.0389916*pt));
  }else{
    if(fabs(eta)<0.17) res = 0.025945+0.16529*exp(-(0.027385*pt));
    else if(fabs(eta)<0.35) res = 0.0257294+0.155366*exp(-(0.0268639*pt));
    else if(fabs(eta)<0.5) res = 0.0231924+0.166508*exp(-(0.0264153*pt));
    else if(fabs(eta)<0.7) res = 0.0252736+0.176977*exp(-(0.0281594*pt));
    else if(fabs(eta)<0.9) res = 0.0269299+0.174505*exp(-(0.0295757*pt));
    else if(fabs(eta)<1.15) res = 0.0250191+0.178321*exp(-(0.0284006*pt));
    else if(fabs(eta)<1.4) res = 0.0242345+0.174682*exp(-(0.027202*pt));
    else if(fabs(eta)<1.7) res = 0.0242234+0.199735*exp(-(0.0291262*pt));
    else if(fabs(eta)<2.1) res = 0.0300017+0.175617*exp(-(0.0357093*pt));
    else if(fabs(eta)<2.5) res = 0.0282978+0.141994*exp(-(0.0383663*pt));
  }
  return res;
}

inline double HelperJet::et(double pt, double eta,  Flavor flav)
{
  double res = 0.29*sqrt(pt);
  if(flav==kB){ 
    if(fabs(eta)<0.17) res = 10.1445+0.0445196*pt;
    else if(fabs(eta)<0.35) res = 8.944+0.0576391*pt;
    else if(fabs(eta)<0.5) res = 8.8462+0.0630703*pt;
    else if(fabs(eta)<0.7) res = 8.90812+0.057605*pt;
    else if(fabs(eta)<0.9) res = 8.42307+0.0623024*pt;
    else if(fabs(eta)<1.15) res = 8.64735+0.0619689*pt;
    else if(fabs(eta)<1.4) res = 8.90574+0.0595526*pt;
    else if(fabs(eta)<1.7) res = 8.74415+0.0612412*pt;
    else if(fabs(eta)<2.1) res = 7.80525+0.0511843*pt;
    else if(fabs(eta)<2.5) res = 6.07877+0.0569574*pt;
  }else{
    if(fabs(eta)<0.17) res = 8.04345+0.0705188*pt;
    else if(fabs(eta)<0.35) res = 7.77565+0.0722247*pt;
    else if(fabs(eta)<0.5) res = 8.53558+0.0647725*pt;
    else if(fabs(eta)<0.7) res = 7.99374+0.0681428*pt;
    else if(fabs(eta)<0.9) res = 7.51695+0.0755128*pt;
    else if(fabs(eta)<1.15) res = 7.70804+0.0748981*pt;
    else if(fabs(eta)<1.4) res = 7.2996+0.0852397*pt;
    else if(fabs(eta)<1.7) res = 7.95923+0.0715722*pt;
    else if(fabs(eta)<2.1) res = 7.68012+0.0543317*pt;
    else if(fabs(eta)<2.5) res = 6.59141+0.0452384*pt;
  }
  return res;
}

inline double HelperJet::eta(double pt, double eta,  Flavor flav)
{
  double res=-1.53e-4*pt+0.05;
  if(flav==kB){
    if(fabs(eta)<0.17) res = 0.0354411+0.130211*exp(-(0.0423321*pt));
    else if(fabs(eta)<0.35) res = 0.0341742+0.115784*exp(-(0.0375787*pt));
    else if(fabs(eta)<0.5) res = 0.0336537+0.143109*exp(-(0.0416702*pt));
    else if(fabs(eta)<0.7) res = 0.031894+0.140572*exp(-(0.0405564*pt));
    else if(fabs(eta)<0.9) res = 0.0301507+0.126005*exp(-(0.0384231*pt));
    else if(fabs(eta)<1.15) res = 0.0268018+0.140956*exp(-(0.03872*pt));
    else if(fabs(eta)<1.4) res = 0.0242941+0.131887*exp(-(0.0337869*pt));
    else if(fabs(eta)<1.7) res = 0.0223985+0.142949*exp(-(0.0351412*pt));
    else if(fabs(eta)<2.1) res = 0.0229694+0.145993*exp(-(0.0391759*pt));
    else if(fabs(eta)<2.5) res = 0.0243171+0.145419*exp(-(0.0430824*pt));
  }else{
    if(fabs(eta)<0.17) res = 0.038575+0.129026*exp(-(0.0427157*pt));
    else if(fabs(eta)<0.35) res = 0.0400221+0.127222*exp(-(0.0458212*pt));
    else if(fabs(eta)<0.5) res = 0.0358472+0.119962*exp(-(0.0395586*pt));
    else if(fabs(eta)<0.7) res = 0.0359612+0.136094*exp(-(0.04212*pt));
    else if(fabs(eta)<0.9) res = 0.0327581+0.127339*exp(-(0.0393077*pt));
    else if(fabs(eta)<1.15) res = 0.0305095+0.116238*exp(-(0.0347537*pt));
    else if(fabs(eta)<1.4) res = 0.027529+0.117893*exp(-(0.032858*pt));
    else if(fabs(eta)<1.7) res = 0.0253065+0.112751*exp(-(0.0307028*pt));
    else if(fabs(eta)<2.1) res = 0.0266957+0.131998*exp(-(0.0380667*pt));
    else if(fabs(eta)<2.5) res = 0.0298413+0.13651*exp(-(0.0458475*pt));
  }
  return res;
}

//Resolution for muons
  class HelperMuon {

  public:
    HelperMuon(){};
    ~HelperMuon(){};

    inline double pt_new(double pt);
    inline double pt (double pt, double eta);
    inline double eta(double pt, double eta);
    inline double phi(double pt, double eta);
    inline double a(double pt, double eta);
    inline double b(double pt, double eta);
    inline double c(double pt, double eta);
    inline double d(double pt, double eta);
    inline double et(double pt, double eta);
    inline double theta(double pt, double eta);
  };

inline double HelperMuon::pt_new(double pt){
  double relres(0.);
  if (pt<=100) relres = 0.01;
  else if (pt>100 && pt<=200) relres = 0.02;
  else if (pt>200 && pt<=400) relres = 0.035;
  else if (pt>400 && pt<=1000) relres = 0.07;
  double res = pt*relres;
  return res;
}

inline double HelperMuon::pt(double pt, double eta)
{
  double res = 1.5e-4*(pt*pt)+0.534*fabs(eta)+1.9e-2;
  return res;
}

inline double HelperMuon::a(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = -0.00163044+0.00921744*exp(-(-0.00517804*pt));
  else if(fabs(eta)<0.35) res = -38.9811+38.9892*exp(-(-1.58728e-06*pt));
  else if(fabs(eta)<0.5) res = -17.254+17.2634*exp(-(-2.86961e-06*pt));
  else if(fabs(eta)<0.7) res = 0.00651163+0.0038473*exp(-(-0.00716166*pt));
  else if(fabs(eta)<0.9) res = -14.9098+14.9207*exp(-(-3.61147e-06*pt));
  else if(fabs(eta)<1.15) res = -0.0130723+0.028881*exp(-(-0.00143687*pt));
  else if(fabs(eta)<1.3) res = 0.0102039+0.00629269*exp(-(-0.00659229*pt));
  else if(fabs(eta)<1.6) res = -26.1501+26.1657*exp(-(-2.75489e-06*pt));
  else if(fabs(eta)<1.9) res = -26.1006+26.1168*exp(-(-2.65457e-06*pt));
  else if(fabs(eta)<2.5) res = -110.342+110.361*exp(-(-1.3011e-06*pt));
  return res;
}

inline double HelperMuon::b(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = -25.2075+25.2104*exp(-(-1.10483e-05*pt));
  else if(fabs(eta)<0.35) res = -117.663+117.666*exp(-(-1.71416e-06*pt));
  else if(fabs(eta)<0.5) res = -21.8718+21.8742*exp(-(-8.13178e-06*pt));
  else if(fabs(eta)<0.7) res = -35.0557+35.0573*exp(-(-5.75421e-06*pt));
  else if(fabs(eta)<0.9) res = -6.37721+6.38003*exp(-(-2.51376e-05*pt));
  else if(fabs(eta)<1.15) res = -73.9844+73.9873*exp(-(-1.75066e-06*pt));
  else if(fabs(eta)<1.3) res = -32.7368+32.7402*exp(-(-3.28819e-06*pt));
  else if(fabs(eta)<1.6) res = -46.9103+46.9132*exp(-(-2.64771e-06*pt));
  else if(fabs(eta)<1.9) res = -63.3183+63.3218*exp(-(-1.56186e-06*pt));
  else if(fabs(eta)<2.5) res = -13.913+13.9174*exp(-(-6.62559e-06*pt));
  return res;
}

inline double HelperMuon::c(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = -0.00784191+0.0107731*exp(-(-0.00262573*pt));
  else if(fabs(eta)<0.35) res = -4.16489+4.16781*exp(-(-8.18221e-06*pt));
  else if(fabs(eta)<0.5) res = -0.00758502+0.0109898*exp(-(-0.0023199*pt));
  else if(fabs(eta)<0.7) res = -0.0190548+0.0222179*exp(-(-0.00140706*pt));
  else if(fabs(eta)<0.9) res = -7.74332+7.74699*exp(-(-3.61138e-06*pt));
  else if(fabs(eta)<1.15) res = 0.00239644+0.00222721*exp(-(-0.00729116*pt));
  else if(fabs(eta)<1.3) res = -50.7222+50.7266*exp(-(-5.72416e-07*pt));
  else if(fabs(eta)<1.6) res = -9.91368+9.91795*exp(-(-3.8653e-06*pt));
  else if(fabs(eta)<1.9) res = 0.000103356+0.00464674*exp(-(-0.00750142*pt));
  else if(fabs(eta)<2.5) res = -0.00425536+0.00971869*exp(-(-0.00546016*pt));
  return res;
}

inline double HelperMuon::d(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = -0.00146069+0.00904233*exp(-(-0.00524103*pt));
  else if(fabs(eta)<0.35) res = -34.455+34.4631*exp(-(-1.78448e-06*pt));
  else if(fabs(eta)<0.5) res = -15.3417+15.3511*exp(-(-3.23796e-06*pt));
  else if(fabs(eta)<0.7) res = 0.00639493+0.00393752*exp(-(-0.00707984*pt));
  else if(fabs(eta)<0.9) res = -12.855+12.866*exp(-(-4.16994e-06*pt));
  else if(fabs(eta)<1.15) res = -0.00526993+0.0211877*exp(-(-0.00186356*pt));
  else if(fabs(eta)<1.3) res = 0.0105021+0.00600376*exp(-(-0.00677709*pt));
  else if(fabs(eta)<1.6) res = -23.5742+23.5898*exp(-(-3.0526e-06*pt));
  else if(fabs(eta)<1.9) res = -27.317+27.3331*exp(-(-2.56587e-06*pt));
  else if(fabs(eta)<2.5) res = -112.151+112.17*exp(-(-1.25452e-06*pt));
  return res;
}

inline double HelperMuon::theta(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.000327129+0.000101907*exp(-(0.0172489*pt));
  else if(fabs(eta)<0.35) res = 0.00025588+9.99049e-05*exp(-(0.0189315*pt));
  else if(fabs(eta)<0.5) res = 0.000215701+0.000109968*exp(-(0.0406069*pt));
  else if(fabs(eta)<0.7) res = 0.000196151+5.66998e-05*exp(-(0.0263501*pt));
  else if(fabs(eta)<0.9) res = 0.000106579+0.000113065*exp(-(0.00377145*pt));
  else if(fabs(eta)<1.15) res = 0.000120697+0.0001408*exp(-(0.0520758*pt));
  else if(fabs(eta)<1.3) res = 9.57227e-05+0.000205436*exp(-(0.0633277*pt));
  else if(fabs(eta)<1.6) res = 8.56706e-05+0.000117908*exp(-(0.0582652*pt));
  else if(fabs(eta)<1.9) res = 6.62861e-05+0.000110841*exp(-(0.0640963*pt));
  else if(fabs(eta)<2.5) res = 6.3783e-05+0.000122656*exp(-(0.0974097*pt));
  return res;
}

inline double HelperMuon::phi(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 7.21523e-05+0.000293781*exp(-(0.0518546*pt));
  else if(fabs(eta)<0.35) res = 7.15456e-05+0.000290324*exp(-(0.0496431*pt));
  else if(fabs(eta)<0.5) res = 7.25417e-05+0.000322288*exp(-(0.0497559*pt));
  else if(fabs(eta)<0.7) res = 7.24273e-05+0.000301504*exp(-(0.0480936*pt));
  else if(fabs(eta)<0.9) res = 7.23791e-05+0.0003355*exp(-(0.0454216*pt));
  else if(fabs(eta)<1.15) res = 8.13896e-05+0.000432844*exp(-(0.0480919*pt));
  else if(fabs(eta)<1.3) res = 7.93329e-05+0.000333341*exp(-(0.0367028*pt));
  else if(fabs(eta)<1.6) res = 9.34279e-05+0.000372581*exp(-(0.0429296*pt));
  else if(fabs(eta)<1.9) res = 0.000112312+0.000479423*exp(-(0.0513205*pt));
  else if(fabs(eta)<2.5) res = 0.000144398+0.000432592*exp(-(0.0400788*pt));
  return res;
}

inline double HelperMuon::et(double pt, double eta)
{
  double res = 1000.; 
  if(fabs(eta)<0.17) res = -0.0552605+0.0115814*pt;
  else if(fabs(eta)<0.35) res = -0.05039+0.0122729*pt;
  else if(fabs(eta)<0.5) res = -0.0435167+0.0128949*pt;
  else if(fabs(eta)<0.7) res = -0.038473+0.0129088*pt;
  else if(fabs(eta)<0.9) res = -0.0333693+0.0140788*pt;
  else if(fabs(eta)<1.15) res = -0.0102406+0.0180217*pt;
  else if(fabs(eta)<1.3) res = -0.0244845+0.0190696*pt;
  else if(fabs(eta)<1.6) res = -0.055785+0.0205605*pt;
  else if(fabs(eta)<1.9) res = -0.0457006+0.0204167*pt;
  else if(fabs(eta)<2.5) res = -0.0399952+0.027388*pt;
  return res;
}

inline double HelperMuon::eta(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.000322451+0.000107167*exp(-(0.0156347*pt));
  else if(fabs(eta)<0.35) res = 0.000271279+0.000106197*exp(-(0.0235732*pt));
  else if(fabs(eta)<0.5) res = 0.000235783+0.000146536*exp(-(0.0491182*pt));
  else if(fabs(eta)<0.7) res = 0.000241657+0.00012571*exp(-(0.0681321*pt));
  else if(fabs(eta)<0.9) res = 0.000218419+0.000380006*exp(-(0.0749053*pt));
  else if(fabs(eta)<1.15) res = 0.000184903+0.000217257*exp(-(0.0520574*pt));
  else if(fabs(eta)<1.3) res = 0.000169888+0.000490542*exp(-(0.0747289*pt));
  else if(fabs(eta)<1.6) res = 0.000179318+0.000290048*exp(-(0.0613212*pt));
  else if(fabs(eta)<1.9) res = 0.000144481+0.000260476*exp(-(0.0392984*pt));
  else if(fabs(eta)<2.5) res = 0.000158888+0.000442249*exp(-(0.0537515*pt));
  return res;
}

//Resolution for electrons
  class HelperElectron {

  public:
    HelperElectron(){};
    ~HelperElectron(){};

    inline double pt_new (double pt, double eta);
    inline double pt (double pt, double eta);
    inline double eta(double pt, double eta);
    inline double phi(double pt, double eta);
    inline double a(double pt, double eta);
    inline double b(double pt, double eta);
    inline double c(double pt, double eta);
    inline double d(double pt, double eta);
    inline double et(double pt, double eta);
    inline double theta(double pt, double eta);
  };


inline double HelperElectron::pt_new(double pt, double eta)
{
  double res= 1000.;
  if( fabs(eta)<=1.1){
    res = 0.02*pt;
  }
  else{
    if(pt<10.) res = 0.105*pt;
    if(pt<15.) res = 0.085*pt;
    if(pt<20.) res = 0.07*pt;
    if(pt<25.) res = 0.06*pt;
    if(pt<30.) res = 0.055*pt;
    if(pt<35.) res = 0.05*pt;
    if(pt<40.) res = 0.045*pt;
    if(pt<50.) res = 0.04*pt;
    else res = 0.035*pt;
  }
  return res;
}


inline double HelperElectron::pt(double pt, double eta)
{
  double res=0.2*sqrt(pt);
  if( fabs(eta)<=0.8 ) 
    res+=-0.28*fabs(eta)+0.54;
  else if( 0.8<fabs(eta) && fabs(eta)<=1.4 ) 
    res+= 1.52*fabs(eta)-1.07;
  else  
    res+=-0.158*eta*eta +0.97;
  return res;
}

inline double HelperElectron::a(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.0114228+0.1135*exp(-(0.111787*pt));
  else if(fabs(eta)<0.35) res = 0.010297+0.0163219*exp(-(0.0248655*pt));
  else if(fabs(eta)<0.5) res = 0.0081634+0.0230718*exp(-(0.0250963*pt));
  else if(fabs(eta)<0.7) res = 0.0122627+0.0228103*exp(-(0.0420525*pt));
  else if(fabs(eta)<0.9) res = 0.0124118+0.0384051*exp(-(0.0545988*pt));
  else if(fabs(eta)<1.15) res = 0.0133398+0.0307574*exp(-(0.0316605*pt));
  else if(fabs(eta)<1.3) res = 0.0144956+0.0355239*exp(-(0.0273916*pt));
  else if(fabs(eta)<1.6) res = -13.9017+13.9315*exp(-(7.01823e-06*pt));
  else if(fabs(eta)<1.9) res = 0.0106309+0.0230149*exp(-(0.00792017*pt));
  else if(fabs(eta)<2.5) res = 0.0187822+0.119922*exp(-(0.140598*pt));
  return res;
}

inline double HelperElectron::b(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = -36.5147+36.516*exp(-(-9.36847e-06*pt));
  else if(fabs(eta)<0.35) res = -0.132294+0.134289*exp(-(-0.0021212*pt));
  else if(fabs(eta)<0.5) res = -96.9689+96.9704*exp(-(-2.83158e-06*pt));
  else if(fabs(eta)<0.7) res = -122.755+122.757*exp(-(-1.90859e-06*pt));
  else if(fabs(eta)<0.9) res = -226.455+226.457*exp(-(-8.91395e-07*pt));
  else if(fabs(eta)<1.15) res = -249.279+249.281*exp(-(-7.57645e-07*pt));
  else if(fabs(eta)<1.3) res = -16.7465+16.7481*exp(-(-1.17848e-05*pt));
  else if(fabs(eta)<1.6) res = -128.535+128.537*exp(-(-1.52236e-06*pt));
  else if(fabs(eta)<1.9) res = -66.1731+66.1762*exp(-(-2.49121e-06*pt));
  else if(fabs(eta)<2.5) res = -0.0199509+0.0237796*exp(-(-0.00591733*pt));
  return res;
}

inline double HelperElectron::c(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = -68.6101+68.6161*exp(-(-7.88064e-07*pt));
  else if(fabs(eta)<0.35) res = -130.219+130.224*exp(-(-8.39696e-07*pt));
  else if(fabs(eta)<0.5) res = -262.825+262.83*exp(-(-4.26008e-07*pt));
  else if(fabs(eta)<0.7) res = -47.898+47.9041*exp(-(-2.07568e-06*pt));
  else if(fabs(eta)<0.9) res = -178.795+178.799*exp(-(-5.78263e-07*pt));
  else if(fabs(eta)<1.15) res = -249.26+249.267*exp(-(-3.16408e-07*pt));
  else if(fabs(eta)<1.3) res = -173.603+173.606*exp(-(-1.58982e-06*pt));
  else if(fabs(eta)<1.6) res = -98.726+98.7326*exp(-(-2.67151e-06*pt));
  else if(fabs(eta)<1.9) res = -82.6028+82.6122*exp(-(-2.79483e-06*pt));
  else if(fabs(eta)<2.5) res = -119.94+119.95*exp(-(-1.69882e-06*pt));
  return res;
}

inline double HelperElectron::d(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.0114197+0.081753*exp(-(0.0969625*pt));
  else if(fabs(eta)<0.35) res = 0.0097562+0.0137737*exp(-(0.0187112*pt));
  else if(fabs(eta)<0.5) res = 0.00844899+0.021595*exp(-(0.0241679*pt));
  else if(fabs(eta)<0.7) res = 0.0124758+0.0243678*exp(-(0.0459914*pt));
  else if(fabs(eta)<0.9) res = 0.0117518+0.0357218*exp(-(0.0463595*pt));
  else if(fabs(eta)<1.15) res = 0.0133048+0.0292511*exp(-(0.0305812*pt));
  else if(fabs(eta)<1.3) res = 0.0144949+0.0337369*exp(-(0.0269872*pt));
  else if(fabs(eta)<1.6) res = 0.0137634+0.0343427*exp(-(0.0208025*pt));
  else if(fabs(eta)<1.9) res = 0.00572644+0.0270719*exp(-(0.00536847*pt));
  else if(fabs(eta)<2.5) res = 0.0189177+0.169591*exp(-(0.152597*pt));
  return res;
}

inline double HelperElectron::theta(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.000282805+0.000157786*exp(-(0.0343273*pt));
  else if(fabs(eta)<0.35) res = 0.000184362+4.34076e-05*exp(-(-0.0131909*pt));
  else if(fabs(eta)<0.5) res = 0.000249332+5.83114e-05*exp(-(0.0508729*pt));
  else if(fabs(eta)<0.7) res = -6.56357e-05+0.000325051*exp(-(0.00177319*pt));
  else if(fabs(eta)<0.9) res = 0.000182277+0.000125324*exp(-(0.0581923*pt));
  else if(fabs(eta)<1.15) res = 0.000140771+0.000407914*exp(-(0.0971668*pt));
  else if(fabs(eta)<1.3) res = 0.000125551+0.001266*exp(-(0.180176*pt));
  else if(fabs(eta)<1.6) res = 0.000107631+101920*exp(-(1.17024*pt));
  else if(fabs(eta)<1.9) res = 8.33927e-05+158936*exp(-(1.20127*pt));
  else if(fabs(eta)<2.5) res = 6.55271e-05+0.12459*exp(-(0.437044*pt));
  return res;
}

inline double HelperElectron::phi(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.000175676+0.000471783*exp(-(0.0383161*pt));
  else if(fabs(eta)<0.35) res = 0.000202185+0.00048635*exp(-(0.0373331*pt));
  else if(fabs(eta)<0.5) res = 0.000150868+0.000444216*exp(-(0.0268835*pt));
  else if(fabs(eta)<0.7) res = 0.000243624+0.00182347*exp(-(0.0850746*pt));
  else if(fabs(eta)<0.9) res = 0.000254463+0.000431233*exp(-(0.0446507*pt));
  else if(fabs(eta)<1.15) res = 0.000309592+0.000918965*exp(-(0.0555677*pt));
  else if(fabs(eta)<1.3) res = 0.000502204+0.000277996*exp(-(0.076721*pt));
  else if(fabs(eta)<1.6) res = 0.000361181+0.000655126*exp(-(0.0238519*pt));
  else if(fabs(eta)<1.9) res = 0.000321587+0.00155721*exp(-(0.0337709*pt));
  else if(fabs(eta)<2.5) res = 0.000819101+0.00205336*exp(-(0.0992806*pt));
  return res;
}

inline double HelperElectron::et(double pt, double eta)
{
  double res = 1000.; 
  if(fabs(eta)<0.17) res = 0.326238+0.00760789*pt;
  else if(fabs(eta)<0.35) res = 0.40493+0.00659958*pt;
  else if(fabs(eta)<0.5) res = 0.369785+0.00690331*pt;
  else if(fabs(eta)<0.7) res = 0.437539+0.00703808*pt;
  else if(fabs(eta)<0.9) res = 0.456138+0.0078252*pt;
  else if(fabs(eta)<1.15) res = 0.518685+0.00907836*pt;
  else if(fabs(eta)<1.3) res = 0.733672+0.00953255*pt;
  else if(fabs(eta)<1.6) res = 1.02678+0.0116056*pt;
  else if(fabs(eta)<1.9) res = 0.948368+0.00977619*pt;
  else if(fabs(eta)<2.5) res = 0.418302+0.0127816*pt;
  return res;
}

inline double HelperElectron::eta(double pt, double eta)
{
  double res = 1000.;
  if(fabs(eta)<0.17) res = 0.000266154+0.000104322*exp(-(0.0140464*pt));
  else if(fabs(eta)<0.35) res = -0.251539+0.251791*exp(-(-7.37147e-07*pt));
  else if(fabs(eta)<0.5) res = 0.000290074+1.54664e-12*exp(-(-0.115541*pt));
  else if(fabs(eta)<0.7) res = 8.37182e-05+0.000233453*exp(-(0.00602386*pt));
  else if(fabs(eta)<0.9) res = 0.000229422+0.000114253*exp(-(0.0188935*pt));
  else if(fabs(eta)<1.15) res = 0.000191525+0.000404238*exp(-(0.0554545*pt));
  else if(fabs(eta)<1.3) res = 0.000195461+1.43699e-07*exp(-(-0.0315088*pt));
  else if(fabs(eta)<1.6) res = 0.000223422+2.05169e+07*exp(-(1.41408*pt));
  else if(fabs(eta)<1.9) res = -0.867114+0.867336*exp(-(4.08639e-07*pt));
  else if(fabs(eta)<2.5) res = -0.866567+0.866746*exp(-(-8.62871e-07*pt));
  return res;
}

//Resolution for MET
  class HelperMET {
    
  public:
    HelperMET(){};
    ~HelperMET(){};
    
    inline double met(double met);
    inline double a(double pt);
    inline double b(double pt);
    inline double c(double pt);
    inline double d(double pt);
    inline double theta();
    inline double phi(double pt);
    inline double et(double pt);
    inline double eta();
  };


inline double HelperMET::met(double met)
{
  return 1.14*exp(-2.16e-3*met)+0.258;
}

inline double HelperMET::a(double pt)
{
  double res = 0.241096+0.790046*exp(-(0.0248773*pt));
  return res;
}

inline double HelperMET::b(double pt)
{
  double res = -141945+141974*exp(-(-1.20077e-06*pt));
  return res;
}

inline double HelperMET::c(double pt)
{
  double res = 21.5615+1.13958*exp(-(-0.00921408*pt));
  return res;
}

inline double HelperMET::d(double pt)
{
	double res = 0.376192+15.2485*exp(-(0.116907*pt));
  return res;
}

inline double HelperMET::theta()
{
	double res = 1000000.;
  return res;
}

inline double HelperMET::phi(double pt)
{
	double res = 0.201336+1.53501*exp(-(0.0216707*pt));
  return res;
}

inline double HelperMET::et(double pt)
{
	double res = 11.7801+0.145218*pt;
  return res;
}

inline double HelperMET::eta()
{
	double res = 1000000.;
  return res;
}

