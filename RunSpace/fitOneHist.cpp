
Double_t myfunction(Double_t *xx, Double_t *p)
{
   Double_t x =xx[0];
   //Double_t f = p[0]*exp(-0.5*pow((x-p[1])/p[2],2)) + p[3] + p[4]*x + p[5]*pow(x,2);
   Double_t f = p[0]*exp(-0.5*pow((x-p[1])/p[2],2)) + p[3]*exp(-0.5*pow((x-p[4])/p[5],2));
   return f;
}




void fitOneHist(){
gStyle->SetOptFit(1111);

std::string mdcHitType = "MdcMcHits";
std::string detType = "baseline";
std::string magType = "simUniField";
std::string angle= "0";
std::string momentum="125";

std::string fileName = "perf/upgrade/beamPosition0/uniformField/PixelR35mm/MdcNoise/ana_momentum_scan/CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.recUniField/pt/mu-/CosThetaDeg_-0.8_0.8_momentumMev_400_400/tracksummary_ckf.root";

//std::string fileName = "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/" + detType + "/ana_momentum_angle_scan/CKF.smeared." + mdcHitType + ".PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss." + magType + ".recUniField/pt/mu-/absCosThetaDeg_" + angle + "_" + angle + "_momentumMev_" + momentum + "_" + momentum + "/tracksummary_ckf.root";

std::cout<<"Reading file "<< fileName << std::endl;

TFile* file = TFile::Open(Form("%s",fileName.c_str()),"READ");

TTree * tree = (TTree*)file->Get("tracksummary");

//TH1F* hist = new TH1F("hist", "", 200, -0.02, 0.02);
//tree->Draw("-1./eQOP_fit*sin(eTHETA_fit)-t_pT>>hist", "nMeasurements>=10&&nHoles<=20&&nOutliers<=1000");

TH1F* hist = new TH1F("hist", "", 200, -3, 3);
tree->Draw("res_eLOC0_fit>>hist", "nMeasurements>=6&&nHoles<=20&&nOutliers<=1000");

hist->Draw();

//TF1 *fa = new TF1("fa","gaus(0)", -0.004, 0.004);
TF1 *fa = new TF1("fa","gaus(0)", hist->GetBinLowEdge(1),hist->GetBinLowEdge(1) +  hist->GetBinWidth(1)*hist->GetNbinsX());
TFitResultPtr rt = hist->Fit("fa","S");
//fa->Draw("same");

//TF1* fb = new TF1("fb", myfunction, hist->GetBinLowEdge(1),hist->GetBinLowEdge(1) +  hist->GetBinWidth(1)*hist->GetNbinsX(), 6);
TF1* fb = new TF1("fb", myfunction, -rt->Parameter(2)*5, rt->Parameter(2)*5, 6);
fb->SetParameter(0,rt->Parameter(0)*0.9);
fb->SetParameter(1,rt->Parameter(1));
fb->SetParameter(2,rt->Parameter(2)*0.9);
fb->SetParameter(3,rt->Parameter(0)*0.1);
fb->SetParameter(4,rt->Parameter(1));
fb->SetParameter(5,rt->Parameter(2)*1.5);


TFitResultPtr r = hist->Fit("fb","S");
if(not r.Get()){
  std::cout<<"Fit failed"  << std::endl;
}
fb->Draw("same");


}
