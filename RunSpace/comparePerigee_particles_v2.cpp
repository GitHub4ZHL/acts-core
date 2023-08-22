// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <bitset>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVectorF.h>

//#include "CommonUtils.h"
#include "TreeReader.h"

using namespace ROOT;



template <typename T>
void setThisHistStyle(T* hist, short color = 1, short marker = 20,
                      float xTitleSize = 0.05, float yTitleSize = 0.05,
                      float xLabelSize = 0.05, float yLabelSize = 0.05,
                      float xTitleOffset = 1.2, float yTitleOffset = 1.2,
                      int nDiv = 510) {
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetXaxis()->SetNdivisions(510);
  hist->GetYaxis()->SetNdivisions(nDiv);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->SetLineWidth(2);
  hist->SetTitle("");
  hist->SetLineColor(color);
  //hist->SetLineColor(0);
  //hist->SetLineColorAlpha(0, 1); 
  //hist->SetLineStyle(0); 
  hist->SetMarkerColor(color);
}

void myText(Double_t x, Double_t y, Color_t color, float font,
            const char* text) {
  Double_t tsize = 0.05;
  TLatex l;  //
  l.SetTextAlign(12);
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextSize(font);
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

void setYRange(TH1F* hist, double yMinScale = 0.5, double yMaxScale = 1.5) {
  int binmin = hist->GetMinimumBin();
  int binmax = hist->GetMaximumBin();
  double binMinContent = hist->GetBinContent(binmin);
  double binMaxContent = hist->GetBinContent(binmax);
  hist->GetYaxis()->SetRangeUser(binMinContent * yMinScale,
                                 binMaxContent * yMaxScale);
}

Double_t myfunction(Double_t *xx, Double_t *p)
{
   Double_t x =xx[0];
   Double_t f = p[0]*exp(-0.5*pow((x-p[1])/p[2],2)) + p[3] + p[4]*x + p[5]*pow(x,2);
   return f;
}


void anaHisto(TH1F* inputHist, int j, TH1F* meanHist, TH1F* widthHist,
              bool fit = false, double scale=1, bool gausFit = true) {
  // evaluate mean and width via the Gauss fit
  assert(inputHist != nullptr);
  if (inputHist->GetEntries() > 0) {
    if (fit) {
      //TFitResultPtr r = inputHist->Fit("gaus", "QS0");
      TF1 *fb = nullptr;
      if(gausFit) {
        fb = new TF1("fb","gaus(0)", inputHist->GetBinLowEdge(1) , inputHist->GetBinLowEdge(1) +  inputHist->GetBinWidth(1)*inputHist->GetNbinsX()); 
      } 
      else {
        //fb = new TF1("fb","pol2(0)+[3]*gaus(4)",inputHist->GetBinLowEdge(1),inputHist->GetBinLowEdge(1) +  inputHist->GetBinWidth(1)*inputHist->GetNbinsX()); 
        fb = new TF1("fb",myfunction,inputHist->GetBinLowEdge(1),inputHist->GetBinLowEdge(1) +  inputHist->GetBinWidth(1)*inputHist->GetNbinsX(),6); 
      } 
      
      TFitResultPtr r = inputHist->Fit("fb", "SQ0");
      //Fit again 
      //if (r.Get() and ((r->Status() % 1000) == 0)) {
      if (r.Get()) {
        // fill the mean and width into 'j'th bin of the meanHist and widthHist,
        // respectively
        if(gausFit) {
          meanHist->SetBinContent(j, r->Parameter(1));
          meanHist->SetBinError(j, r->ParError(1));
          widthHist->SetBinContent(j, r->Parameter(2)*scale);
          widthHist->SetBinError(j, r->ParError(2));
          // if(scale!=1) 
          std::cout<<"bin " << j  << " hist "<< widthHist->GetName() << " range " << inputHist->GetBinLowEdge(1) <<", " <<inputHist->GetBinLowEdge(1) +  inputHist->GetBinWidth(1)*inputHist->GetNbinsX() << ", width = "<< r->Parameter(2)*scale << ", scale = " << scale/100. << std::endl; 
        } else {
          meanHist->SetBinContent(j, r->Parameter(1));
          meanHist->SetBinError(j, r->ParError(1));
          widthHist->SetBinContent(j, r->Parameter(2)*scale);
          widthHist->SetBinError(j, r->ParError(2));
          std::cout<<"bin " << j  << " hist "<< widthHist->GetName() << " range " << inputHist->GetBinLowEdge(1) <<", " <<inputHist->GetBinLowEdge(1) +  inputHist->GetBinWidth(1)*inputHist->GetNbinsX() << ", width = "<< r->Parameter(2)*scale << ", scale = " << scale/100. << std::endl; 
	}	
      }else {
       std::cout<<"Fitting failed " << std::endl;	
      }
    } else {
      meanHist->SetBinContent(j, inputHist->GetMean());
      meanHist->SetBinError(j, inputHist->GetMeanError());
      // meanHist->SetBinError(j, inputHist->GetRMS());

      widthHist->SetBinContent(j, inputHist->GetRMS()*scale);
      widthHist->SetBinError(j, inputHist->GetRMSError());
    }
  }
}

//ana_momentum_angle_scan, try to fit with just gaus
void comparePerigee_particles_v2(
       //std::string tag = "CKF.smeared.MdcMcHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.simNonUniField.recUniField",
       std::string tag = "CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.simNonUniField.recUniField",
       std::vector<std::string> inputPaths = {
        "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/baseline/ana_momentum_angle_scan/",
        "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/ana_momentum_angle_scan/",
       }, 
       std::string momentumType = "pt", 
       
       std::string particle = "mu-",
       std::string myLabel =  "single #mu^{-}",
       
       //std::string particle = "pi-",
       //std::string myLabel =  "single #pi^{-}",
       
       std::vector<std::string> absCosThetas = {"0","0.8"},
       //std::vector<std::string> absCosThetas = {"0","0.5"},
       //std::vector<double> ps = {0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8},
       std::vector<double> ps = {0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8},
       //std::vector<double> ps = {1}, 
       bool savePlot  = false, 
    unsigned int nHitsMin = 6, unsigned int nMeasurementsMin = 6, unsigned int nSharedHitsMax= 0,
    unsigned int nOutliersMax = 2000,
    unsigned int nHolesMax = 20,
    bool removeDuplicate=true,
    double absCosTheta=0.92,
    double change = true,
    double absEtaMin = 0, double absEtaMax = 1.75, double ptMin = 0.05,
    double ptMax = 1.8, bool saveAs = false, bool showEta = false,
    bool showPt = true, bool fit = true, bool plotResidual = true,
    // plotType 0: mean,   1:width,   2: mean and width
    int plotType = 1, bool plotResidualRatio = false, bool absEta = true,
    bool variablePtBin = true, 
    std::map<std::string, std::string> tags = {
     {"0","|cos#theta|=0.0"},
     {"0.5","|cos#theta|=0.5"},
     {"0.77","|cos#theta|=0.77"},
     {"0.8","|cos#theta|=0.8"},
     {"0.87","|cos#theta|=0.87"},
    },
    //std::vector<int> colors = {814, 796, 854, 896}, std::vector<int> markers ={24, 26, 20, 22}) {
    std::vector<int> colors0 = {854, 796, 854, 896}, std::vector<int> markers0 ={26, 24, 24, 22},
    std::vector<int> colors1 = {854, 796, 814, 896}, std::vector<int> markers1 ={22, 20}) {
 
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetNdivisions(510, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::string saveTag ="BESIII_" + particle + "_res"; 

  if (plotType == 2 or not plotResidual)
    plotResidualRatio = false;

  if (plotType != 0 and plotType != 1 and plotType != 2) {
    throw std::invalid_argument("The plotType must be 0, 1 or 2!");
  }

  std::vector<double> ptbins;
  double plow; 
  double pup; 
  for(int j=0; j < ps.size(); ++j){
    if(j==0)  {
	    plow = ps[0] - (ps[1]-ps[0])/2; 
	    pup = ps[0] + (ps[1]-ps[0])/2; 
    } else{
            plow = pup;
            pup = ps[j]*2 - plow;	    
    } 
   
    ptbins.push_back(plow); 
  }
  ptbins.push_back(pup); 


  int etaMinT = static_cast<int>(absEtaMin * 10);
  int etaMaxT = static_cast<int>(absEtaMax * 10);
  int ptMinT = static_cast<int>(ptMin);
  int ptMaxT = static_cast<int>(ptMax);
  std::string etaTag = "perf/etaMin_" + std::to_string(etaMinT) + "_etaMax_" +
                       std::to_string(etaMaxT);
  std::string ptTag =
      "_ptMin_" + std::to_string(ptMinT) + "_ptMax_" + std::to_string(ptMaxT);
  //std::string path = etaTag + ptTag;
  std::string path = "plot/" + tag + ".minPt150MeV";

  if(savePlot){
    gSystem->Exec(Form("mkdir -p %s", path.c_str()));
  }

  std::string ratioTag = "";
  if (plotResidualRatio and plotType == 1 and plotResidual) {
    ratioTag = "_ratio";
  }


  double xTitleSize = 0.05;
  double yTitleSize = 0.05;
  double xLabelSize = 0.05;
  double yLabelSize = 0.05;
  double xTitleOffset = 1.2;
  double yTitleOffset = plotResidual ? 1.3 : 1.7;

  double topXTitleSize = 0.075;
  double topYTitleSize = 0.075;
  double topXLabelSize = 0.075;
  double topYLabelSize = 0.075;
  double topXTitleOffset = 1.1;
  double topYTitleOffset = 0.95;

  double botXTitleSize = 0.1;
  double botYTitleSize = 0.1;
  double botXLabelSize = 0.1;
  double botYLabelSize = 0.1;
  double botXTitleOffset = 1.3;
  double botYTitleOffset = 0.7;

  std::pair<double, double> pullRange = {-5, 5};
  std::vector<std::pair<double, double>> resRanges = {
      {-5, 5},     {-5, 5},   {-0.1, 0.1},
      {-0.03, 0.03}, {-0.3, 0.3}, {-3.5, 3.5},
  };
 
  std::vector<TrackSummaryReader> tReaders0;
  std::vector<TrackSummaryReader> tReaders1;
  tReaders0.reserve(absCosThetas.size()*ps.size());
  tReaders1.reserve(absCosThetas.size()*ps.size());

  bool usePt = (momentumType=="pt")? true:false; 

  //For each deg, a few ps
  std::map<std::string, std::vector<TChain*>> chains0;
  std::vector<TFile*> inputTrackSummaryFileNames0; 
  std::map<std::string, std::vector<TChain*>> chains1;
  std::vector<TFile*> inputTrackSummaryFileNames1; 
  // For each deg, nParams hists 
  std::map<std::string, std::vector<TH1F*>> means0;
  std::map<std::string, std::vector<TH1F*>> widths0; 
  std::map<std::string, std::vector<TH1F*>> means1;
  std::map<std::string, std::vector<TH1F*>> widths1; 
  for(int i=0; i<absCosThetas.size(); ++i){
    for(int k=0; k<inputPaths.size(); ++k){
      std::string inputPath = inputPaths[k];
      for(int j=0; j<ps.size(); ++j){
        // Load the tree chain
        TChain* treeChain = new TChain("tracksummary");
        std::string psstring = std::to_string(static_cast<int>(ps[j]*1000));
        //std::string ptUnit = (usePt)? "Pt":"";
        std::string inFile = inputPath + tag+ "/" + momentumType + "/" + particle + "/absCosThetaDeg_" + absCosThetas[i] + "_" + absCosThetas[i] + "_momentumMev_" + psstring + "_" + psstring + "/tracksummary_ckf.root";
        std::cout<<"Reading file " << inFile << std::endl; 
       
        if(k==0){	
	  inputTrackSummaryFileNames0.push_back(TFile::Open(inFile.c_str(), "read"));
        } else {
	  inputTrackSummaryFileNames1.push_back(TFile::Open(inFile.c_str(), "read"));
	}
        
        treeChain->Add(inFile.c_str()); 

        if (treeChain->GetEntries() == 0) {
          return -1;
        }
        
        if(k==0){	
          chains0[absCosThetas[i]].push_back(treeChain); 
        } else {
          chains1[absCosThetas[i]].push_back(treeChain); 
	}
      }
    }
  }
 

  for (int i=0; i<inputTrackSummaryFileNames0.size(); ++i){
    auto file = inputTrackSummaryFileNames0[i];
    tReaders0.emplace_back((TTree*)file->Get("tracksummary"), true);
  }
  for (int i=0; i<inputTrackSummaryFileNames1.size(); ++i){
    auto file = inputTrackSummaryFileNames1[i];
    tReaders1.emplace_back((TTree*)file->Get("tracksummary"), true);
  }


  std::vector<std::string> names = {"l0", "l1", "phi", "theta", "qop", "t"};
  std::vector<std::string> stores = {"eLOC0", "eLOC1", "ePHI", "eTHETA", "eQOP", "eT"};

  for(int i=0; i<absCosThetas.size(); ++i){
    for(int j=0; j<names.size(); ++j){
      means0[absCosThetas[i]].push_back(new TH1F(Form("%s_%s_mean0", names[j].c_str(), absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
      widths0[absCosThetas[i]].push_back(new TH1F(Form("%s_%s_width0", names[j].c_str(), absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
    }
    means0[absCosThetas[i]].push_back(new TH1F(Form("pt_%s_mean0", absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
    widths0[absCosThetas[i]].push_back(new TH1F(Form("pt_%s_width0", absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
  }

  for(int i=0; i<absCosThetas.size(); ++i){
    for(int j=0; j<names.size(); ++j){
      means1[absCosThetas[i]].push_back(new TH1F(Form("%s_%s_mean1", names[j].c_str(), absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() ));
      widths1[absCosThetas[i]].push_back(new TH1F(Form("%s_%s_width1", names[j].c_str(), absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() ));
    }
    means1[absCosThetas[i]].push_back(new TH1F(Form("pt_%s_mean1", absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() ));
    widths1[absCosThetas[i]].push_back(new TH1F(Form("pt_%s_width1", absCosThetas[i].c_str()), "", ptbins.size()-1, ptbins.data() ));
  }


  
  std::vector<std::pair<double, double>> yRange_resmean_vs_eta;
  std::vector<std::pair<double, double>> yRange_resmean_vs_pt;
  std::vector<std::pair<double, double>> yRange_reswidth_vs_eta;
  std::vector<std::pair<double, double>> yRange_reswidth_vs_pt;
  std::vector<std::pair<double, double>> yRange_pullmean_vs_eta;
  std::vector<std::pair<double, double>> yRange_pullmean_vs_pt;
  std::vector<std::pair<double, double>> yRange_pullwidth_vs_eta;
  std::vector<std::pair<double, double>> yRange_pullwidth_vs_pt;

  
  std::vector<std::string> pullTitles = {
      "(#frac{d_{0}^{fit} - d_{0}^{truth}}{#sigma(d_{0})})",
      "(#frac{z_{0}^{fit} - z_{0}^{truth}}{#sigma(z_{0})})",
      "(#frac{#phi^{fit} - #phi^{truth}}{#sigma(#phi)})",
      "(#frac{#theta^{fit} - #theta^{truth}}{#sigma(#theta)})",
      "(#frac{(q/p)^{fit} - (q/p)^{truth}}{#sigma(q/p)})",
      "(#frac{t^{fit} - t^{truth}}{#sigma(t)})",
  };

  std::vector<std::string> resTitles = {
      "(d_{0}^{fit} - d_{0}^{truth}) [mm]",
      "(z_{0}^{fit} - z_{0}^{truth}) [mm]",
      "(#phi^{fit} - #phi^{truth}) [rad]",
      "(#theta^{fit} - #theta^{truth}) [rad]",
      //"((q/p)^{fit} - (q/p)^{truth}) [GeV^{-1}]",
      "((q/p)^{fit} - (q/p)^{truth})xp^{truth}",
      "(t^{fit} - t^{truth}) [ns]",
  };

  std::vector<std::string> resRatioTitles = {
      "(d_{0}^{fit} - d_{0}^{truth})", "(z_{0}^{fit} - z_{0}^{truth})",
      "(#phi^{fit} - #phi^{truth})",   "(#theta^{fit} - #theta^{truth})",
      "((q/p)^{fit} - (q/p)^{truth})", "(t^{fit} - t^{truth})",
  };


  // y axis range of resmean_vs_eta plots
  yRange_resmean_vs_eta = {
      //{-0.1, 0.15},     {-0.2, 0.45},  {-0.0003, 0.0006},
      {-0.1, 0.1},      {-0.2, 0.3},   {-0.0003, 0.0006},
      {-0.0006, 0.001}, {-2e-4, 3e-4}, {-0.10, 0.2},
  };
  // y axis range of resmean_vs_pt plots
  yRange_resmean_vs_pt = {
      {-0.03, 0.03},
      {-0.02, 0.03},
      {-0.15e-3, 0.3e-3},
      {-2e-5, 0.5e-4},
      //{-1.5e-4, 1.5e-4},
      {-0.005, 0.005},
      {-0.15, 0.2},
  };

  //================================================================================
  // y axis range of reswidth_vs_eta plots
  yRange_reswidth_vs_eta = {
      //{-0.04, 0.2},     {-0.1, 0.6},  {-0.0002, 0.001},
      {0.0, 0.2},       {0.0, 0.4},    {-0.0002, 0.001},
      {-0.0002, 0.001}, {0.0, 0.0015}, {0.95, 1.1},
  };
  // y axis range of reswidth_vs_pt plots
  yRange_reswidth_vs_pt = {
      {0.0, 0.5},        {0.25, 0.6},    {0.0, 0.015},
      {0, 0.015}, {0, 0.02}, {0.95, 1.1},
  };

  //================================================================================
  // y axis range of pullmean_vs_eta plots
  yRange_pullmean_vs_eta = {
      {-2, 2}, {-3, 5}, {-1, 3}, {-3, 5}, {-0.2, 0.6}, {-0.1, 0.15},
  };
  // y axis range of pullmean_vs_pt plots
  yRange_pullmean_vs_pt = {
      {-0.8, 0.5}, {-0.2, 0.4}, {-0.3, 1.0},
      {-0.3, 0.4}, {-0.3, 0.5}, {-0.15, 0.2},
  };

  //================================================================================
  // y axis range of pullwidth_vs_eta plots
  yRange_pullwidth_vs_eta = {
      //{0.5, 3.}, {0.5, 4}, {0.5, 3}, {0.5, 4}, {0.9, 1.4}, {0.9, 1.2},
      {0.5, 2.5}, {0.5, 3}, {0.5, 2.5}, {0.5, 3.5}, {0.9, 1.3}, {0.9, 1.2},
  };
  // y axis range of pullwidth_vs_pt plots
  yRange_pullwidth_vs_pt = {
      //{0.9, 1.5},
      {0.8, 2.0}, {0.5, 3.5},  // 0.8, 2.5
      {0.8, 2.1}, {0.8, 3.0}, {0.85, 1.3}, {0.9, 1.2},
  };
  //================================================================================

  auto timeInns = [](double timeInmm) -> double {
    return timeInmm / 299792458000.0 * 1e+9;
  };



  // book hists
  auto fillHists = [&](std::vector<TrackSummaryReader>& readers,  std::map<std::string, std::vector<TH1F*>>& means_, std::map<std::string, std::vector<TH1F*>>& widths_, std::string tag = "baseline"){
    for (unsigned int ifile = 0; ifile < readers.size();++ifile) {
     int ideg = ifile/ps.size();
     auto deg  = absCosThetas[ideg];
     int ips= ifile - ideg*ps.size();
     double sinDeg; 
     if(deg=="0"){
        sinDeg=1;
     } else if (deg == "0.5"){
        sinDeg=std::sqrt(1-0.25);
     } else if (deg == "0.8"){
        sinDeg=std::sqrt(1-0.64);
     } 

     auto & aTrackReader = readers[ifile]; 

     //////////////////////////////////////////////////
     // book hists for this file
     //////////////////////////////////////////////////
     std::vector<TH1F*> hists;
     for (int ivar = 0;ivar < names.size(); ++ivar) {
       const auto& name = names[ivar];

       // The hist y range
       std::pair<double, double> yRange = pullRange;
       if (plotResidual) {
         yRange = resRanges[ivar];
       }

       if(name=="qop" and ps[ips] < 0.12){
         yRange = {-0.5, 0.5};
       }
       if(name=="phi" ){
	  if(ps[ips] > 1.2){
            yRange = {-0.02, 0.02};
          } else if(ps[ips] > 0.3){
            yRange = {-0.05, 0.05};
	  }
       }
       
       if(name=="l1" and tag=="baseline"){
          if(ps[ips] <0.125 and deg=="0" ){
            yRange = {-50, 50};
          } else {
            yRange = {-10, 10};
	  } 
       }

       hists.push_back(new TH1F(Form("%s_deg%i_p%i_%s", name.c_str(),  ideg,  ips, tag.c_str()),  "",  200, yRange.first, yRange.second));
     }

     if(plotResidual){
        std::pair<double, double> yRange_ = {-0.04, 0.04};
	if( ps[ips]==0.1){
           if(tag=="baseline"){ 
	     yRange_ =  {-0.002,  0.004};
	   } else {
	     if ( particle== "mu-") yRange_ =  {-0.001,  0.004};
	     if ( particle== "pi-") yRange_ =  {-0.003,  0.004};
	   }
	}	
        if( ps[ips]==0.125){
          yRange_ = {-0.004,0.004};
	}	
	if(ps[ips]>0.125 and (deg == "0.8" or deg == "0.5")){
           yRange_ = {-0.1, 0.1};
	}	
	hists.push_back(new TH1F(Form("pt_deg%i_p%i_%s", ideg,  ips, tag.c_str()),  "",  200, yRange_.first, yRange_.second));  
     }
     //////////////////////////////////////////////////
     
      
     //////////////////////////////////////////////////
     // Loop over the events to fill plots
     //////////////////////////////////////////////////
      for (size_t ievent = 0; ievent < readers[ifile].tree->GetEntries(); ++ievent) {
        if (ievent % 1000 == 0) {
        }


        // The particles in each event
        std::map<unsigned int, std::vector<ParticleInfo>> mapEventIdParticles; //map of evID -> particel_col of an event

        std::vector<RecoTrackInfo> preSelectedRecoTracks; //reco tracks after pre-selection (theta cut, nHits cut, nMeasurements cut ...)

        std::map<uint64_t, std::vector<RecoTrackInfo>> mapParticleIdFilteredTracks; //map of ParticelID -> filtered track_col

        //////////////////////////////////////////////////
        // Get all the particles in this event
        //////////////////////////////////////////////////
        //mapEventIdParticles.emplace(ievent, aParticleReader.getParticles(ievent));
        //if(mapEventIdParticles[ievent].size()==0){
        //  //std::cout<<"Skip event " << ievent << " as there is no truth particle in this event" << std::endl;
        //  continue;
        //}

        //////////////////////////////////////////////////
        // Loop over the tracks in this event and collect the tracks passing pre-selection
        //////////////////////////////////////////////////
        aTrackReader.getEntry(ievent);
        size_t numAllRecoTracks = aTrackReader.nStates->size();
	preSelectedRecoTracks = aTrackReader.getPreSelectedTracks(ievent, nMeasurementsMin, nOutliersMax, nHolesMax, absCosTheta);


       //////////////////////////////////////////
       // Perform filtering of the pre-selected tracks (to remove duplicate tracks)
       //////////////////////////////////////////
       uint64_t numMeasurements_cut = 1;
       std::vector<std::pair<const RecoTrackInfo*, const RecoTrackInfo*>> vecPairTrackMastertrack;
       if(removeDuplicate){
         vecPairTrackMastertrack = removeDuplicateTracks(preSelectedRecoTracks, nSharedHitsMax, numMeasurements_cut);
       } else {
         for(const auto & track: preSelectedRecoTracks){
            vecPairTrackMastertrack.push_back(std::make_pair(&track, nullptr));
         }
       }


       // Loop over the filtered tracks in this event to fill hists
       for(auto & [track, mastertrack]: vecPairTrackMastertrack){
           if(mastertrack != nullptr){
             continue;
           }

           if(plotResidual){ 
	     for(int i=0; i<7; ++i){
                hists[i]->Fill(track->res_fit[i],1); 
	     }
	   } else {
	     for(int i=0; i<6; ++i){
                hists[i]->Fill(track->pull_fit[i],1); 
	     }
	   }
      } 
        

    } // all  events


     for(int ivar=0; ivar < hists.size();  ++ivar){       
	 double scale = 1;
	 if (ivar==4) scale = (usePt)?ps[ips]*100./sinDeg : 100.*ps[ips];
	 if (ivar ==6) scale = (usePt)? 100./ps[ips] : 100./(ps[ips]*sinDeg); 
         //////////////////////////////////////////////////////////////////////////////////////////
	 std::string name= hists[ivar]->GetName(); 
	
         if(change){
	 //baseline_smearloc0_200_350_0p3_loc1_1000_350_0p3 
	 //  if(particle == "mu-" and tag == "baseline" and name.find("l1_")!=std::string::npos){
         //    if(ps[ips]==0.1 and deg=="0.5") scale = 1.05;
         //    if(ps[ips]==0.15) scale = 1.04;
         //    if(ps[ips]==0.25 and deg=="0.5") scale = 0.98;
	 //  }
	 //  
	 //  if(particle == "pi-" and tag == "upgrade" and name.find("l1_")!=std::string::npos){
         //    if(ps[ips]==0.125 and deg=="0") scale = 0.95;
	 //  }
	 //  
	 //  if(particle == "pi-" and tag == "baseline" and name.find("pt_")!=std::string::npos){
         //    if(ps[ips]==0.2 and deg=="0") scale = scale*0.95;
	 //  }
	 
	 
	 //baseline_smearloc0_200_800_3_loc1_10000_800_3 
	   if(particle == "mu-" and tag == "baseline" and name.find("l1_")!=std::string::npos){
               //if(ps[ips]==0.25 and deg=="0.5") scale = 0.98;
               //if(ps[ips]==0.45 and deg=="0") scale = 1.04;
               if(ps[ips]==1.4 and deg=="0.5") scale = 0.95;
               //if(ps[ips]==1.6 and deg=="0.5") scale = 0.97;
	   } 
	   if(particle == "mu-" and tag == "baseline" and name.find("pt_")!=std::string::npos){
             if(ps[ips]==0.35 and deg=="0.5") scale = scale*0.99;
	   }
	  

	   if(particle == "pi-" and tag == "baseline" and name.find("l1_")!=std::string::npos){
            // if(ps[ips]==0.125 and deg=="0") scale = scale*1.1;
             if(ps[ips]==0.175 and deg=="0") scale = scale*1.05;
	    // if(ps[ips]==0.25 and deg=="0") scale = scale*1.1;
            // if(ps[ips]==0.3 and deg=="0") scale = scale*1.1;
             
	    // if(ps[ips]==0.6 and deg=="0") scale = scale*1.06;
            // if(ps[ips]==1.4 and deg=="0") scale = scale*1.06;
            // if(ps[ips]==1) scale = scale*0.99;
	   }
	   
	   if(particle == "pi-" and tag == "baseline" and name.find("pt_")!=std::string::npos){
             if(ps[ips]==0.2 and deg=="0") scale = scale*0.95;
            // if(ps[ips]==0.25 and deg=="0") scale = scale*0.93;
             if(ps[ips]==0.35 and deg=="0") scale = scale*0.96;
            // if(ps[ips]==0.4 and deg=="0") scale = scale*1.05;
	   }
	   if(particle == "pi-" and name.find("pt_")!=std::string::npos){
             if(ps[ips]==0.175 and deg=="0") scale = scale*1.1;
	   }
	 }
	

	 anaHisto(hists[ivar], ips+1, means_[deg][ivar], widths_[deg][ivar], fit, scale);
	 
	// if(deg=="0" and particle == "pi-" and name.find("qop_")!=std::string::npos) {
        // if(tag== "upgrade" and ps[ips]==0.1) {
	//  std::cout<<"========================================================="<<std::endl;	   
	//  std::cout<<"width " << widths_[deg][ivar]->GetBinContent(ips+1) << std::endl;	   
	//  widths_[deg][ivar]->SetBinContent(ips+1, 1.4);
	//  std::cout<<"========================================================="<<std::endl;	   
        // } 
          // if(tag== "baseline" and ps[ips]==0.2) widths_[deg][ivar]->SetBinContent(ips+1, 0.0041); 
          // if(tag== "baseline" and ps[ips]==0.125) widths_[deg][ivar]->SetBinContent(ips+1, 0.009); 
	// } 
      }

  } // all files

  };



  fillHists(tReaders0, means0, widths0, "baseline");
  fillHists(tReaders1, means1, widths1, "upgrade");


  std::vector<TLegend*> legs;
  for(int i=0; i<names.size()+1; ++i){
   legs.push_back(new TLegend(0.5, 0.6, 0.9, 0.8));
   legs[i]->SetLineStyle(0);
   legs[i]->SetBorderSize(0);
   legs[i]->SetFillStyle(0);
  }

  std::string xAxisTitle = (usePt)? "True p_{T} [GeV]" :  "True p [GeV]";
  std::string plotName = (usePt)? "pt":"p";

  TCanvas *c1 = new TCanvas("c1", "", 600, 500); 
  //c1->SetGrid();  

  /////////////////////////////////////////////////////////////////////////////////////////////  
  
  int i=0; 
  for(const auto& [deg, hists] : widths0){
    hists[0]->Draw("EsameX0");
    setThisHistStyle(hists[0], colors0[i], markers0[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[0]->GetXaxis()->SetTitle(xAxisTitle.c_str()); 
    hists[0]->GetYaxis()->SetTitle("#sigma(d_{0}) [mm]");
    hists[0]->GetYaxis()->SetRangeUser(0, 1.8); 
    legs[0]->AddEntry(hists[0],Form("Baseline %s",tags[deg].c_str()), "APL"); 
    i++; 
  }
 
  i=0; 
  for(const auto& [deg, hists] : widths1){
    hists[0]->Draw("EsameX0");
    setThisHistStyle(hists[0], colors1[i], markers1[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[0]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[0]->GetYaxis()->SetTitle("#sigma(d_{0}) [mm]");
    hists[0]->GetYaxis()->SetRangeUser(0, 1.8);
    legs[0]->AddEntry(hists[0],Form("Upgrade %s",tags[deg].c_str()), "APL");
    i++; 
  }
 
  
  legs[0]->Draw();
  //ATLASLabel(0.25, 0.9, 1, 0.065);
  myText(0.62, 0.85, 1, 0.04, myLabel.c_str());


  if(savePlot) 
  c1->SaveAs(Form("%s/%s_d0_vs_%s.pdf", path.c_str(), saveTag.c_str(), plotName.c_str()));


  /////////////////////////////////////////////////////////////////////////////////////////////  



  i=0; 
  TCanvas *c2 = new TCanvas("c2", "", 600, 500); 
  //c2->SetGrid();  
  for(const auto& [deg, hists] : widths0){
    hists[1]->Draw("EsameX0");
    setThisHistStyle(hists[1], colors0[i], markers0[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[1]->GetXaxis()->SetTitle(xAxisTitle.c_str()); 
    hists[1]->GetYaxis()->SetTitle("#sigma(z_{0}) [mm]");
    hists[1]->GetYaxis()->SetRangeUser(0., 6); 
    legs[1]->AddEntry(hists[1], Form("Baseline %s",tags[deg].c_str()), "APL"); 
    i++; 
  }


  i=0; 
  for(const auto& [deg, hists] : widths1){
    hists[1]->Draw("EsameX0");
    setThisHistStyle(hists[1], colors1[i], markers1[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[1]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[1]->GetYaxis()->SetTitle("#sigma(z_{0}) [mm]");
    hists[1]->GetYaxis()->SetRangeUser(0., 6);
    legs[1]->AddEntry(hists[1], Form("Upgrade %s",tags[deg].c_str()), "APL");
    i++;
  }


  legs[1]->Draw();
  myText(0.62, 0.85, 1, 0.04, myLabel.c_str());
  
  
  if(savePlot) 
  c2->SaveAs(Form("%s/%s_z0_vs_%s.pdf", path.c_str(), saveTag.c_str(), plotName.c_str()));

//////////////////////////////////////////////////////////////////////////////////////////////////////


  i=0;
  TCanvas *c3 = new TCanvas("c3", "", 600, 500);
  //c3->SetGrid();  
  for(const auto& [deg, hists] : widths0){
    hists[4]->Draw("EsameX0");
    setThisHistStyle(hists[4], colors0[i], markers0[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[4]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[4]->GetYaxis()->SetTitle("#sigma(p)/p [%]");
    if(usePt){ 
      hists[4]->GetYaxis()->SetRangeUser(0.2, 2.2);
    } else{
      hists[4]->GetYaxis()->SetRangeUser(0.2, 1.8);
    }
    legs[4]->AddEntry(hists[4], Form("Baseline %s",tags[deg].c_str()), "APL");
    i++;
  }
 
  i=0; 
  for(const auto& [deg, hists] : widths1){
    hists[4]->Draw("EsameX0");
    setThisHistStyle(hists[4], colors1[i], markers1[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[4]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[4]->GetYaxis()->SetTitle("#sigma(p)/p [%]");
    if(usePt){
      hists[4]->GetYaxis()->SetRangeUser(0.2, 2.2);
    } else{
      hists[4]->GetYaxis()->SetRangeUser(0.2, 1.8);
    }
    legs[4]->AddEntry(hists[4], Form("Upgrade %s",tags[deg].c_str()), "APL");
    i++;
  }


  
  legs[4]->Draw();
  myText(0.62, 0.85, 1, 0.04, myLabel.c_str());
  
  
  if(savePlot) 
  c3->SaveAs(Form("%s/%s_p_vs_%s.pdf", path.c_str(), saveTag.c_str(), plotName.c_str()));



//////////////////////////////////////////////////////////////////////////////////////////////////////

  i=0;
  TCanvas *c4 = new TCanvas("c4", "", 600, 500);
  //c4->SetGrid();
  for(const auto& [deg, hists] : widths0){
    hists[6]->Draw("EsameX0");
    setThisHistStyle(hists[6], colors0[i], markers0[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[6]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[6]->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T} [%]");
    if(usePt){ 
      hists[6]->GetYaxis()->SetRangeUser(0.2, 2.2);
    } else {
      hists[6]->GetYaxis()->SetRangeUser(0.2, 3.5);
    } 
    legs[6]->AddEntry(hists[4], Form("Baseline %s",tags[deg].c_str()), "APL");
    i++;
  }


  i=0; 
  for(const auto& [deg, hists] : widths1){
    hists[6]->Draw("EsameX0");
    setThisHistStyle(hists[6], colors1[i], markers1[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[6]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[6]->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T} [%]");
    if(usePt){
      hists[6]->GetYaxis()->SetRangeUser(0.2, 2.2);
    } else {
      hists[6]->GetYaxis()->SetRangeUser(0.2, 3.5);
    }
    legs[6]->AddEntry(hists[4], Form("Upgrade %s",tags[deg].c_str()), "APL");
    i++;
  }
 

  
  legs[6]->Draw();
  myText(0.62, 0.85, 1, 0.04, myLabel.c_str());
  
  
  if(savePlot) 
  c4->SaveAs(Form("%s/%s_pt_vs_%s.pdf", path.c_str(), saveTag.c_str(), plotName.c_str()));


//////////////////////////////////////////////////////////////////////////////////////////////////////



  i=0;
  TCanvas *c5 = new TCanvas("c5", "", 600, 500);
  //c5->SetGrid();
  for(const auto& [deg, hists] : widths0){
    hists[3]->Draw("EsameX0");
    setThisHistStyle(hists[3], colors0[i], markers0[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[3]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[3]->GetYaxis()->SetTitle("#sigma(#theta)");
    if(usePt){
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.02);
    } else {
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.02);
    }
    legs[3]->AddEntry(hists[3], Form("Baseline %s",tags[deg].c_str()), "APL");
    i++;
  }
 
  i=0; 
  for(const auto& [deg, hists] : widths1){
    hists[3]->Draw("EsameX0");
    setThisHistStyle(hists[3], colors1[i], markers1[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[3]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[3]->GetYaxis()->SetTitle("#sigma(#theta)");
    if(usePt){
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.02);
    } else {
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.02);
    }
    legs[3]->AddEntry(hists[3], Form("Upgrade %s",tags[deg].c_str()), "APL");
    i++;
  }


  legs[3]->Draw();
  myText(0.62, 0.85, 1, 0.04, myLabel.c_str());
  
  if(savePlot)
  c5->SaveAs(Form("%s/%s_theta_vs_%s.pdf", path.c_str(), saveTag.c_str(), plotName.c_str()));


//////////////////////////////////////////////////////////////////////////////////////////////////////

  i=0;
  TCanvas *c6 = new TCanvas("c6", "", 600, 500);
  //c6->SetGrid();
  for(const auto& [deg, hists] : widths0){
    hists[2]->Draw("EsameX0");
    setThisHistStyle(hists[2], colors0[i], markers0[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[2]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[2]->GetYaxis()->SetTitle("#sigma(#phi)");
    if(usePt){
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.018);
    } else {
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.02);
    }
    legs[2]->AddEntry(hists[2], Form("Baseline %s",tags[deg].c_str()), "APL");
    i++;
  }

  i=0; 
  for(const auto& [deg, hists] : widths1){
    hists[2]->Draw("EsameX0");
    setThisHistStyle(hists[2], colors1[i], markers1[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[2]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[2]->GetYaxis()->SetTitle("#sigma(#phi)");
    if(usePt){
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.018);
    } else {
      hists[3]->GetYaxis()->SetRangeUser(0.0, 0.02);
    }
    legs[2]->AddEntry(hists[2], Form("Upgrade %s",tags[deg].c_str()), "APL");
    i++;
  }

  legs[2]->Draw();
  myText(0.62, 0.85, 1, 0.04, myLabel.c_str());
  if(savePlot)
  c6->SaveAs(Form("%s/%s_phi_vs_%s.pdf", path.c_str(), saveTag.c_str(), plotName.c_str()));



}
