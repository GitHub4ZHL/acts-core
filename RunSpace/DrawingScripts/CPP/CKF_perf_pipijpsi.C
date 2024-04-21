// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>

#include "CommonUtils.h"
#include "TreeReader.h"
#include "SetStyleATLAS.hpp"

/// This script/function reads all the reconstructed tracks from e.g.
/// 'tracksummary_ckf.root' and the (possibly selected) truth particles from
/// e.g. 'track_finder_particles.root' (which contains the info of 'nHits'), and
/// defines the efficiency, fake rate and duplicaiton rate. It aims to make
/// custom definition and tuning of the reconstruction performance easier.
/// Multiple files for the reconstructed tracks are allowed.
/// 
/// NB: It's very likely that fiducal cuts are already imposed on the truth
/// particles. Please check the selection criteria in the truth fitting example
/// which writes out the 'track_finder_particles.root'. For instance, if the
/// truth particles are already required to have pT > 1 GeV, it does not make
/// sense to have ptMin = 0.5 GeV here.
///
//

//std::string tag = "CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps100000.betheLoss.recUniField";
//std::string tag = "CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps100000.betheLoss.recUniField.extToZero";
//std::string tag = "CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps100000.betheLoss.recUniField.extToVertex";
//std::string tag = "CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps100000.betheLoss.recUniField.extToVertex.TestSeed";
std::string tag = "CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps100000.betheLoss.recUniField.extToVertex.TestSeed5";
//std::string tag = "test";

void CKF_perf_pipijpsi(
  const std::vector<std::string>& inputSimParticleFileNames =
  {
    // tree: particles, event per entry 
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/baseline/beamPosition0/uniformField/pipijpsi/" +tag  + "/particles.root",
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/beamPosition0/uniformField/PixelR35mm/MdcNoise/pipijpsi/" +tag  + "/particles.root",
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/beamPosition0/uniformField/PixelR45mm/MdcNoise/pipijpsi/" +tag  + "/particles.root",
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/beamPosition0/uniformField/PixelR55mm/MdcNoise/pipijpsi/" +tag  + "/particles.root",
  },
  const std::vector<std::string>& inputTrackSummaryFileNames =
  {
    // tree: track_reconstructed, event per entry
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/baseline/beamPosition0/uniformField/pipijpsi/" +tag  + "/tracksummary_ckf.root",
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/beamPosition0/uniformField/PixelR35mm/MdcNoise/pipijpsi/" +tag  + "/tracksummary_ckf.root",
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/beamPosition0/uniformField/PixelR45mm/MdcNoise/pipijpsi/" +tag  + "/tracksummary_ckf.root",
    "/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/upgrade/beamPosition0/uniformField/PixelR55mm/MdcNoise/pipijpsi/" +tag  + "/tracksummary_ckf.root",
  },
  bool useBESIIITracking = true,
  std::string myLabel = "#psi(3686)#rightarrow #pi^{+}#pi^{-}J/#psi(#rightarrow #mu^{+}#mu^{-}), #mu", 
  const std::vector<std::string>& trackSummaryFileLegends_plus =
  {
    //"MDC, #pi^{+}",
    //"MDC + Pixel (r=35 mm), #pi^{+}",
    //"MDC + Pixel (r=45 mm), #pi^{+}",
    //"MDC + Pixel (r=55 mm), #pi^{+}",
    "MDC",
    "MDC + Pixel (r=35 mm)",
    "MDC + Pixel (r=45 mm)",
    "MDC + Pixel (r=55 mm)",
  },
  const std::vector<std::string>& trackSummaryFileLegends_minus =
  {
    //"MDC, #pi^{-}",
    //"MDC + Pixel (r=35 mm), #pi^{-}",
    //"MDC + Pixel (r=45 mm), #pi^{-}",
    //"MDC + Pixel (r=55 mm), #pi^{-}",
    "MDC",
    "MDC + Pixel (r=35 mm)",
    "MDC + Pixel (r=45 mm)",
    "MDC + Pixel (r=55 mm)",
  },
  std::vector<int> colors_plus={
    1,
    832, 
    //854,
    874,
    796,
    //896,
  },
  std::vector<int> colors_minus={
    1, 
    832, 
    //854,
    874,
    796,
  //  896,
  },
  std::vector<int> markers_plus ={24, 20, 21, 22},
  std::vector<int> markers_minus ={24, 20, 21, 22},
  const std::string& simParticleTreeName = "particles",
  const std::string& trackSummaryTreeName = "tracksummary",
  bool  savePlot= false,
  bool  removeDuplicate= false,
  unsigned int nHitsMin = 5, unsigned int nMeasurementsMin = 1, unsigned int nSharedHitsMax= 1, unsigned int nOutliersMax = 1000, unsigned int nHolesMax = 1000,
  double ptMin = 0.15, double absCosTheta=0.93, double truthMatchProbMin = 0.5, int absPdgId=13) {
  gStyle->SetOptFit(0011);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.05);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.05, "xy");
  gStyle->SetLabelSize(0.05, "xy");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleOffset(1.5, "y");
  gStyle->SetNdivisions(510, "x");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  double XTitleSize = 0.05;
  double YTitleSize = 0.05;
  double XLabelSize = 0.05;
  double YLabelSize = 0.05;
  double XTitleOffset = 1.;
  double YTitleOffset = 0.8;

  /*
  double topXTitleSize = 0.06;
  double topYTitleSize = 0.06;
  double topXLabelSize = 0.06;
  double topYLabelSize = 0.06;
  double topXTitleOffset = 1.;
  double topYTitleOffset = 0.9;

  double botXTitleSize = 0.15;
  double botYTitleSize = 0.13;
  double botXLabelSize = 0.15;
  double botYLabelSize = 0.15;
  double botXTitleOffset = 1.3;
  double botYTitleOffset = 0.35;
*/

  double topXTitleSize = 0.065;
  double topYTitleSize = 0.065;
  double topXLabelSize = 0.065;
  double topYLabelSize = 0.065;
  double topXTitleOffset = 0.9;
  double topYTitleOffset = 0.8;

  double botXTitleSize = 0.14;
  double botYTitleSize = 0.12;
  double botXLabelSize = 0.13;
  double botYLabelSize = 0.13;
  double botXTitleOffset = 1.4;
  double botYTitleOffset = 0.4;

  //float myLabelSize=0.05; //0.7:0.3
  float myLabelSize=0.055; //0.7:0.3


  //std::vector<double> ptRanges = {40, 0.05, 0.45};
  std::vector<double> ptRanges;
  std::string particle; 
  if(absPdgId==211){
    //ptRanges  = {20, 0.05, 0.45};
    //ptRanges  = {8, 0.05, 0.45};  // 0.05 GeV
    //ptRanges  = {6, 0.1, 0.4};  // 0.05 GeV
    ptRanges  = {5, 0.15, 0.4};  // 0.05 GeV
    particle ="pi"; 
  } else {
    //ptRanges  = {14, 0.4, 1.8}; //0.1 GeV/
    ptRanges  = {13, 0.5, 1.8}; //0.1 GeV/
    particle ="mu"; 
  }


  // Check the inputs are valid
  if (inputSimParticleFileNames.size()!=inputTrackSummaryFileNames.size() or inputTrackSummaryFileNames.size() != trackSummaryFileLegends_plus.size()) {
    throw std::invalid_argument(
      "Please specify the legends you want to show for all the track files");
  }

  // The number of track files to read
  unsigned int nTrackFiles = inputTrackSummaryFileNames.size();
  std::vector<TFile*> particleFiles;
  std::vector<TFile*> trackFiles;
  particleFiles.reserve(nTrackFiles);
  trackFiles.reserve(nTrackFiles);
  for (const auto& fileName : inputSimParticleFileNames) {
    particleFiles.push_back(TFile::Open(fileName.c_str(), "read")); }
  for (const auto& fileName : inputTrackSummaryFileNames) {
    trackFiles.push_back(TFile::Open(fileName.c_str(), "read"));
  }

  // Define variables for tree reading (turn on the events sorting since we have more than one root files to read)
  std::vector<RootParticleReader> pReaders;
  std::vector<TrackSummaryReader> tReaders;
  pReaders.reserve(nTrackFiles);
  tReaders.reserve(nTrackFiles);


  for (const auto& particleFile : particleFiles) {
    pReaders.emplace_back((TTree*)particleFile->Get(simParticleTreeName.c_str()), true);
  }
  for (const auto& trackFile : trackFiles) {
    tReaders.emplace_back((TTree*)trackFile->Get(trackSummaryTreeName.c_str()), true);
  }

  
  std::vector<size_t> nEvents;
  nEvents.reserve(nTrackFiles);
  for (const auto& tReader : tReaders) {
    size_t entries = tReader.tree->GetEntries();
    //size_t entries = 1000; //10000
    
    nEvents.push_back(entries);
  }

  // Define the efficiency plots
  //////////////////////////////////////////////////////////////////////////
  // plus + minus 
  std::vector<TEfficiency*> trackEff_vs_theta;
  std::vector<TProfile*> trackPurity_vs_theta;
  std::vector<TProfile*> nParticleHits_vs_theta;
  std::vector<TEfficiency*> fakeRate_vs_theta;
  std::vector<TEfficiency*> duplicateRate_vs_theta;

  std::vector<TEfficiency*> trackEff_vs_pt;
  std::vector<TProfile*> trackPurity_vs_pt;
  std::vector<TProfile*> nParticleHits_vs_pt;
  std::vector<TEfficiency*> fakeRate_vs_pt;
  std::vector<TEfficiency*> duplicateRate_vs_pt;

  std::vector<TEfficiency*> trackEff_2d;
  std::vector<TProfile2D*> trackPurity_2d; 
  
  std::vector<TH1F*> ratio_trackEff_vs_theta;
  std::vector<TH1F*> ratio_trackEff_vs_pt;
  
  
  //////////////////////////////////////////////////////////////////////////
  // plus
  std::vector<TEfficiency*> trackEff_vs_theta_plus;
  std::vector<TProfile*> trackPurity_vs_theta_plus;
  std::vector<TProfile*> nParticleHits_vs_theta_plus;
  std::vector<TEfficiency*> fakeRate_vs_theta_plus;
  std::vector<TEfficiency*> duplicateRate_vs_theta_plus;

  std::vector<TEfficiency*> trackEff_vs_pt_plus;
  std::vector<TProfile*> trackPurity_vs_pt_plus;
  std::vector<TProfile*> nParticleHits_vs_pt_plus;
  std::vector<TEfficiency*> fakeRate_vs_pt_plus;
  std::vector<TEfficiency*> duplicateRate_vs_pt_plus;
  
  std::vector<TEfficiency*> trackEff_2d_plus;
  std::vector<TProfile2D*> trackPurity_2d_plus;
 
  std::vector<TH1F*> ratio_trackEff_vs_theta_plus;
  std::vector<TH1F*> ratio_trackEff_vs_pt_plus;
  
  //////////////////////////////////////////////////////////////////////////
  // minus 
  std::vector<TEfficiency*> trackEff_vs_theta_minus;
  std::vector<TProfile*> trackPurity_vs_theta_minus;
  std::vector<TProfile*> nParticleHits_vs_theta_minus;
  std::vector<TEfficiency*> fakeRate_vs_theta_minus;
  std::vector<TEfficiency*> duplicateRate_vs_theta_minus;
  
  std::vector<TEfficiency*> trackEff_vs_pt_minus;
  std::vector<TProfile*> trackPurity_vs_pt_minus;
  std::vector<TProfile*> nParticleHits_vs_pt_minus;
  std::vector<TEfficiency*> fakeRate_vs_pt_minus;
  std::vector<TEfficiency*> duplicateRate_vs_pt_minus;
  
  std::vector<TEfficiency*> trackEff_2d_minus;
  std::vector<TProfile2D*> trackPurity_2d_minus;

  std::vector<TH1F*> ratio_trackEff_vs_theta_minus;
  std::vector<TH1F*> ratio_trackEff_vs_pt_minus;

  std::cout<<"Start booking " << std::endl;

  for (int i = 0; i < nTrackFiles; ++i) {
    trackEff_vs_theta.push_back(new TEfficiency(
                                       Form("trackeff_vs_theta_%i", i), ";Truth cos#theta;Efficiency", 20, -1, 1));
    trackPurity_vs_theta.push_back(new TProfile(
                                       Form("hiteff_vs_theta_%i", i), ";Truth cos#theta;Track purity", 20, -1, 1, 0, 1.4));
    nParticleHits_vs_theta.push_back(new TProfile(
                                       Form("nParticleHits_vs_theta_%i", i), ";Truth cos#theta;nParticleHits", 20, -1, 1, 0, 300));
    fakeRate_vs_theta.push_back(new TEfficiency(
                                       Form("fakerate_vs_theta_%i", i), ";cos#theta;Fake rate", 20, -1, 1));
    duplicateRate_vs_theta.push_back(new TEfficiency(
                                            Form("duplicaterate_vs_theta_%i", i), ";cos#theta;Duplicate rate", 20, -1, 1));
    trackEff_vs_pt.push_back(new TEfficiency(
                                    Form("trackeff_vs_pt_%i", i), ";Truth p_{T} [GeV/c];Efficiency", ptRanges[0], ptRanges[1], ptRanges[2]));
    trackPurity_vs_pt.push_back(new TProfile(
                                    Form("hiteff_vs_pt_%i", i), ";Truth p_{T} [GeV/c];Track purity", ptRanges[0], ptRanges[1], ptRanges[2], 0, 1.4));
    nParticleHits_vs_pt.push_back(new TProfile(
                                    Form("nParticleHits_vs_pt_%i", i), ";Truth p_{T} [GeV/c];nParticleHits", ptRanges[0], ptRanges[1], ptRanges[2], 0, 300));
    fakeRate_vs_pt.push_back(new TEfficiency(
                                    Form("fakerate_vs_pt_%i", i), ";p_{T} [GeV/c];Fake rate", ptRanges[0], ptRanges[1], ptRanges[2]));
    duplicateRate_vs_pt.push_back(new TEfficiency(
                                         Form("duplicaterate_vs_pt_%i", i), ";p_{T} [GeV/c];Duplicate rate", ptRanges[0], ptRanges[1], ptRanges[2]));
    trackEff_2d.push_back(new TEfficiency(
                                    Form("trackeff_2d_%i", i), ";Truth p_{T} [GeV/c];Truth cos#theta;Efficiency", ptRanges[0], ptRanges[1], ptRanges[2], 20, -1, 1));
    trackPurity_2d.push_back(new TProfile2D(
                                    Form("hiteff_2d_%i", i), ";Truth p_{T} [GeV/c];Truth cos#theta;Track purity", ptRanges[0], ptRanges[1], ptRanges[2], 20, -1, 1, 0, 1.));


    trackEff_vs_theta_plus.push_back(new TEfficiency(
                                       Form("trackeff_vs_theta_plus_%i", i), ";Truth cos#theta;Efficiency", 20, -1, 1));
    trackPurity_vs_theta_plus.push_back(new TProfile(
                                       Form("hiteff_vs_theta_plus_%i", i), ";Truth cos#theta;Track purity", 20, -1, 1, 0, 1.4));
    nParticleHits_vs_theta_plus.push_back(new TProfile(
                                       Form("nParticleHits_vs_theta_plus_%i", i), ";Truth cos#theta;nParticleHits", 20, -1, 1, 0, 300));
    fakeRate_vs_theta_plus.push_back(new TEfficiency(
                                       Form("fakerate_vs_theta_plus_%i", i), ";cos#theta;Fake rate", 20, -1, 1));
    duplicateRate_vs_theta_plus.push_back(new TEfficiency(
                                            Form("duplicaterate_vs_theta_plus_%i", i), ";cos#theta;Duplicate rate", 20, -1, 1));
    trackEff_vs_pt_plus.push_back(new TEfficiency(
                                    Form("trackeff_vs_pt_plus_%i", i), ";Truth p_{T} [GeV/c];Efficiency", ptRanges[0], ptRanges[1], ptRanges[2]));
    trackPurity_vs_pt_plus.push_back(new TProfile(
                                    Form("hiteff_vs_pt_plus_%i", i), ";Truth p_{T} [GeV/c];Track purity", ptRanges[0], ptRanges[1], ptRanges[2], 0, 1.4));
    nParticleHits_vs_pt_plus.push_back(new TProfile(
                                    Form("nParticleHits_vs_pt_plus_%i", i), ";Truth p_{T} [GeV/c];nParticleHits", ptRanges[0], ptRanges[1], ptRanges[2], 0, 300));
    fakeRate_vs_pt_plus.push_back(new TEfficiency(
                                    Form("fakerate_vs_pt_plus_%i", i), ";p_{T} [GeV/c];Fake rate", ptRanges[0], ptRanges[1], ptRanges[2]));
    duplicateRate_vs_pt_plus.push_back(new TEfficiency(
                                         Form("duplicaterate_vs_pt_plus_%i", i), ";p_{T} [GeV/c];Duplicate rate", ptRanges[0], ptRanges[1], ptRanges[2]));
    
    trackEff_2d_plus.push_back(new TEfficiency(
                                    Form("trackeff_2d_plus_%i", i), ";Truth p_{T} [GeV/c];Truth cos#theta;Efficiency", ptRanges[0], ptRanges[1], ptRanges[2], 20, -1, 1));
    trackPurity_2d_plus.push_back(new TProfile2D(
                                    Form("hiteff_2d_plus_%i", i), ";Truth p_{T} [GeV/c];Truth cos#theta;Track purity", ptRanges[0], ptRanges[1], ptRanges[2], 20, -1, 1, 0, 1.));
   


    trackEff_vs_theta_minus.push_back(new TEfficiency(
                                        Form("trackeff_vs_theta_minus_%i", i), ";Truth cos#theta;Efficiency", 20, -1, 1)); // cosin
    trackPurity_vs_theta_minus.push_back(new TProfile(
                                        Form("hiteff_vs_theta_minus_%i", i), ";Truth cos#theta;Track purity", 20, -1, 1, 0, 1.4)); // cosin
    nParticleHits_vs_theta_minus.push_back(new TProfile(
                                        Form("nParticleHits_vs_theta_minus_%i", i), ";Truth cos#theta;nParticleHits", 20, -1, 1, 0, 300)); // cosin
    fakeRate_vs_theta_minus.push_back(new TEfficiency(
                                        Form("fakerate_vs_theta_minus_%i", i), ";cos#theta;Fake rate", 20, -1, 1));
    duplicateRate_vs_theta_minus.push_back(new TEfficiency(
                                             Form("duplicaterate_vs_theta_minus_%i", i), ";cos#theta;Duplicate rate", 20, -1, 1));
    trackEff_vs_pt_minus.push_back(new TEfficiency(
                                     Form("trackeff_vs_pt_minus_%i", i), ";Truth p_{T} [GeV/c];Efficiency", ptRanges[0], ptRanges[1], ptRanges[2]));
    trackPurity_vs_pt_minus.push_back(new TProfile(
                                     Form("hiteff_vs_pt_minus_%i", i), ";Truth p_{T} [GeV/c];Track purity", ptRanges[0], ptRanges[1], ptRanges[2], 0, 1.4));
    nParticleHits_vs_pt_minus.push_back(new TProfile(
                                     Form("nParticleHits_vs_pt_minus_%i", i), ";Truth p_{T} [GeV/c];nParticleHits", ptRanges[0], ptRanges[1], ptRanges[2], 0, 300));
    fakeRate_vs_pt_minus.push_back(new TEfficiency(
                                     Form("fakerate_vs_pt_minus_%i", i), ";p_{T} [GeV/c];Fake rate", ptRanges[0], ptRanges[1], ptRanges[2]));
    duplicateRate_vs_pt_minus.push_back(new TEfficiency(
                                          Form("duplicaterate_vs_pt_minus_%i", i), ";p_{T} [GeV/c];Duplicate rate", ptRanges[0], ptRanges[1], ptRanges[2]));
 
    trackEff_2d_minus.push_back(new TEfficiency(
                                    Form("trackeff_2d_minus_%i", i), ";Truth p_{T} [GeV/c];Truth cos#theta;Efficiency", ptRanges[0], ptRanges[1], ptRanges[2], 20, -1, 1));
    trackPurity_2d_minus.push_back(new TProfile2D(
                                    Form("hiteff_2d_minus_%i", i), ";Truth p_{T} [GeV/c];Truth cos#theta;Track purity", ptRanges[0], ptRanges[1], ptRanges[2], 20, -1, 1, 0, 1.));
  
  }

  
  for (int i = 1; i < nTrackFiles; ++i) {
     ratio_trackEff_vs_theta.push_back(new TH1F(
                                        Form("ratio_trackeff_vs_theta_%i", i), "", 20, -1, 1)); // cosin
     ratio_trackEff_vs_pt.push_back(new TH1F(
                                        Form("ratio_trackeff_vs_pt_%i", i), "", ptRanges[0], ptRanges[1], ptRanges[2])); // cosin
     ratio_trackEff_vs_theta[i-1] -> GetXaxis() -> SetTitle("Truth cos#theta"); 
     ratio_trackEff_vs_pt[i-1] -> GetXaxis() -> SetTitle("p_{T} [GeV/c]"); 
    

     ratio_trackEff_vs_theta_minus.push_back(new TH1F(
                                        Form("ratio_trackeff_vs_theta_minus_%i", i), "", 20, -1, 1)); // cosin
     ratio_trackEff_vs_pt_minus.push_back(new TH1F(
                                        Form("ratio_trackeff_vs_pt_minus_%i", i), "", ptRanges[0], ptRanges[1], ptRanges[2])); // cosin
     ratio_trackEff_vs_theta_minus[i-1] -> GetXaxis() -> SetTitle("Truth cos#theta"); 
     ratio_trackEff_vs_pt_minus[i-1] -> GetXaxis() -> SetTitle("p_{T} [GeV/c]"); 
    

     ratio_trackEff_vs_theta_plus.push_back(new TH1F(
                                        Form("ratio_trackeff_vs_theta_plus_%i", i), "", 20, -1, 1)); // cosin
     ratio_trackEff_vs_pt_plus.push_back(new TH1F(
                                        Form("ratio_trackeff_vs_pt_plus_%i", i), "", ptRanges[0], ptRanges[1], ptRanges[2])); // cosin
     ratio_trackEff_vs_theta_plus[i-1] -> GetXaxis() -> SetTitle("Truth cos#theta"); 
     ratio_trackEff_vs_pt_plus[i-1] -> GetXaxis() -> SetTitle("p_{T} [GeV/c]"); 
  }


  std::cout<<"Finish booking "<< std::endl;

  // Set styles
  for (int i = 0; i < nTrackFiles; ++i) {
    auto color_plus = colors_plus[i];
    auto marker_plus  = markers_plus[i];
    auto color_minus = colors_minus[i];
    auto marker_minus  = markers_minus[i];
    auto color = colors_minus[i];
    auto marker  = markers_minus[i];
   
    setEffStyle(trackEff_vs_theta[i], color, marker);
    setEffStyle(trackPurity_vs_theta[i], color, marker);
    setEffStyle(nParticleHits_vs_theta[i], color, marker);
    setEffStyle(fakeRate_vs_theta[i], color, marker);
    setEffStyle(duplicateRate_vs_theta[i], color, marker);
    setEffStyle(trackEff_vs_pt[i], color, marker);
    setEffStyle(trackPurity_vs_pt[i], color, marker);
    setEffStyle(nParticleHits_vs_pt[i], color, marker);
    setEffStyle(fakeRate_vs_pt[i], color, marker);
    setEffStyle(duplicateRate_vs_pt[i], color, marker);
    
    
    setEffStyle(trackEff_vs_theta_plus[i], color_plus, marker_plus);
    setEffStyle(trackPurity_vs_theta_plus[i], color_plus, marker_plus);
    setEffStyle(nParticleHits_vs_theta_plus[i], color_plus, marker_plus);
    setEffStyle(fakeRate_vs_theta_plus[i], color_plus, marker_plus);
    setEffStyle(duplicateRate_vs_theta_plus[i], color_plus, marker_plus);
    setEffStyle(trackEff_vs_pt_plus[i], color_plus, marker_plus);
    setEffStyle(trackPurity_vs_pt_plus[i], color_plus, marker_plus);
    setEffStyle(nParticleHits_vs_pt_plus[i], color_plus, marker_plus);
    setEffStyle(fakeRate_vs_pt_plus[i], color_plus, marker_plus);
    setEffStyle(duplicateRate_vs_pt_plus[i], color_plus, marker_plus);
 
    setEffStyle(trackEff_vs_theta_minus[i], color_minus, marker_minus);
    setEffStyle(trackPurity_vs_theta_minus[i], color_minus, marker_minus);
    setEffStyle(nParticleHits_vs_theta_minus[i], color_minus, marker_minus);
    setEffStyle(fakeRate_vs_theta_minus[i], color_minus, marker_minus);
    setEffStyle(duplicateRate_vs_theta_minus[i], color_minus, marker_minus);
    setEffStyle(trackEff_vs_pt_minus[i], color_minus, marker_minus);
    setEffStyle(trackPurity_vs_pt_minus[i], color_minus, marker_minus);
    setEffStyle(nParticleHits_vs_pt_minus[i], color_minus, marker_minus);
    setEffStyle(fakeRate_vs_pt_minus[i], color_minus, marker_minus);
    setEffStyle(duplicateRate_vs_pt_minus[i], color_minus, marker_minus);
  
    if(i>0){
      setThisHistStyle(ratio_trackEff_vs_theta[i-1], color, marker, botXTitleSize, botYTitleSize,
                         botXLabelSize, botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
      setThisHistStyle(ratio_trackEff_vs_pt[i-1], color, marker, botXTitleSize, botYTitleSize,
                         botXLabelSize, botYLabelSize, botXTitleOffset, botYTitleOffset, 510);
      
      setThisHistStyle(ratio_trackEff_vs_theta_minus[i-1], color_minus, marker_minus, botXTitleSize, botYTitleSize,
                         botXLabelSize, botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
      setThisHistStyle(ratio_trackEff_vs_pt_minus[i-1], color_minus, marker_minus, botXTitleSize, botYTitleSize,
                         botXLabelSize, botYLabelSize, botXTitleOffset, botYTitleOffset, 510);

      setThisHistStyle(ratio_trackEff_vs_theta_plus[i-1], color_plus, marker_plus, botXTitleSize, botYTitleSize,
                         botXLabelSize, botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
      setThisHistStyle(ratio_trackEff_vs_pt_plus[i-1], color_plus, marker_plus, botXTitleSize, botYTitleSize,
                         botXLabelSize, botYLabelSize, botXTitleOffset, botYTitleOffset, 510);
    } 
  }


  // Loop over the track files
  for (unsigned int ifile = 0; ifile < nTrackFiles; ++ifile) {
    //std::cout << "Processing track file: " << inputTrackSummaryFileNames[ifile]
    //          << std::endl;

    auto & aTrackReader = tReaders[ifile];
    auto & aParticleReader = pReaders[ifile];

    // Loop over the events to fill plots
    //for (size_t ievent = 0; ievent < 1000; ++ievent) {
    for (size_t ievent = 0; ievent < nEvents[ifile]; ++ievent) {
      if (ievent>0 and ievent % 10000 == 0) {
        std::cout << "Processed events: " << ievent << std::endl;
      }

      // The particles in each event
      std::map<unsigned int, std::vector<ParticleInfo>> mapEventIdParticles; //map of evID -> particel_col of an event

      std::vector<RecoTrackInfo> preSelectedRecoTracks; //reco tracks after pre-selection (theta cut, nHits cut, nMeasurements cut ...)
      
      std::map<uint64_t, std::vector<RecoTrackInfo>> mapParticleIdFilteredTracks; //map of ParticelID -> filtered track_col

      //////////////////////////////////////////////////
      // Get all the particles in this event
      //////////////////////////////////////////////////
      mapEventIdParticles.emplace(ievent, aParticleReader.getParticles(ievent));
      if(mapEventIdParticles[ievent].size()==0){
        //std::cout<<"Skip event " << ievent << " as there is no truth particle in this event" << std::endl; 
        continue; 
      } 


      //////////////////////////////////////////////////
      // Loop over the tracks in this event and collect the tracks passing pre-selection
      //////////////////////////////////////////////////
      aTrackReader.getEntry(ievent);
      size_t numAllRecoTracks = aTrackReader.nStates->size();
      preSelectedRecoTracks = aTrackReader.getPreSelectedTracks(ievent, nMeasurementsMin, nOutliersMax, nHolesMax, absCosTheta, false, 10, 100);   
      //preSelectedRecoTracks = aTrackReader.getPreSelectedTracks(ievent, nMeasurementsMin, nOutliersMax, nHolesMax, absCosTheta, false);   


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

      uint64_t numMastertracks = 0;
      for(auto & [track, mastertrack]: vecPairTrackMastertrack){
        if(mastertrack == nullptr){
          numMastertracks ++;
          mapParticleIdFilteredTracks[track->majorityParticleId].push_back(
            *track);
	}
      }

      if(numMastertracks < vecPairTrackMastertrack.size()){
        std::cout<< "all reco tracks num: "<< numAllRecoTracks<<std::endl;
        std::cout<< "pre-selected tracks num: "<< vecPairTrackMastertrack.size()<<std::endl;
        std::cout<< "filtered tracks num : "<< numMastertracks <<std::endl;
      }


      //////////////////////////////////////////
      //Loop over all filtered reco tracks for filling fake rate plots 
      // The fake rate is defined as the ratio of selected truth-matched tracks
      // over all filtered tracks
      //////////////////////////////////////////
      for(auto & [track, mastertrack]: vecPairTrackMastertrack){
          if(mastertrack != nullptr){
            continue; 
	  } 

          auto nMajorityHits = track->nMajorityHits;	  
          auto nMeasurements = track->nMeasurements;	  
          auto nSharedHits = track->nSharedHits;	  
          auto qop = track->qop;	  
          auto pt  = track->pt;	  
          auto eta = track->eta;	  
          double theta = std::atan(std::exp(-eta))*2;
            
	  trackPurity_vs_theta[ifile]->Fill(std::cos(theta), nMajorityHits*1.0/nMeasurements, 1.);	    
          trackPurity_vs_pt[ifile]->Fill(pt, nMajorityHits*1.0/nMeasurements, 1.);	    
          trackPurity_2d[ifile]->Fill(pt, std::cos(theta), nMajorityHits*1.0/nMeasurements); 
	  
	  if(qop>0){ 
	    trackPurity_vs_theta_plus[ifile]->Fill(std::cos(theta), nMajorityHits*1.0/nMeasurements, 1.);	    
            trackPurity_vs_pt_plus[ifile]->Fill(pt, nMajorityHits*1.0/nMeasurements, 1.);	    
            trackPurity_2d_plus[ifile]->Fill(pt, std::cos(theta), nMajorityHits*1.0/nMeasurements); 
          } else {
	    trackPurity_vs_theta_minus[ifile]->Fill(std::cos(theta), nMajorityHits*1.0/nMeasurements, 1.);	    
            trackPurity_vs_pt_minus[ifile]->Fill(pt, nMajorityHits*1.0/nMeasurements, 1.);	    
            trackPurity_2d_minus[ifile]->Fill(pt, std::cos(theta), nMajorityHits*1.0/nMeasurements); 
	  }

          // // Fill the fake rate plots
          if (nMajorityHits * 1. / nMeasurements > truthMatchProbMin) {

             fakeRate_vs_theta[ifile]->Fill(false, std::cos(theta));
             fakeRate_vs_pt[ifile]->Fill(false, pt);
            
	     if(qop>0){
               fakeRate_vs_theta_plus[ifile]->Fill(false, std::cos(theta));
               fakeRate_vs_pt_plus[ifile]->Fill(false, pt);
             } else {
               fakeRate_vs_theta_minus[ifile]->Fill(false, std::cos(theta));
               fakeRate_vs_pt_minus[ifile]->Fill(false, pt);
             }

          } else {
             //std::cout<<"Find fake track in event "<< ievent << ", nMeasurements " << nMeasurements <<" nSharedHits " <<  nSharedHits << std::endl; 

	     fakeRate_vs_theta[ifile]->Fill(true, std::cos(theta));
             fakeRate_vs_pt[ifile]->Fill(true, pt);
             if(qop>0){ 
	       fakeRate_vs_theta_plus[ifile]->Fill(true, std::cos(theta));
               fakeRate_vs_pt_plus[ifile]->Fill(true, pt);
             } else {
               fakeRate_vs_theta_minus[ifile]->Fill(true, std::cos(theta));
               fakeRate_vs_pt_minus[ifile]->Fill(true, pt);
             }
          }
      }

      //std::cout<<"Finish fake rate plots" << std::endl;

      /////////////////////////////////////////
      // Loop over all selected and truth-matched tracks
      // The duplicate rate is defined as the ratio of duplicate tracks among
      // all the selected truth-matched tracks (only one track is 'real'; others
      // are 'duplicated')
      //////////////////////////////////////////
      for (auto& [pid, tracks] : mapParticleIdFilteredTracks) {
        
	// Sort all tracks matched to this particle according to majority prob
        // and track quality
        std::sort(tracks.begin(), tracks.end(),
                  [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                    //if (lhs.nMajorityHits > rhs.nMajorityHits) {
                    if (lhs.nMajorityHits/lhs.nMeasurements > rhs.nMajorityHits/rhs.nMeasurements) {
                      return true;
                    }
                    if (lhs.nMeasurements > rhs.nMeasurements) {
                      return true;
                    }
                    return false;
                  });
        // Fill the duplication rate plots
        for (size_t k = 0; k < tracks.size(); ++k) {
          auto eta = tracks[k].eta;
          auto pt = tracks[k].pt;
          auto qop = tracks[k].qop;
          double theta = std::atan(std::exp(-eta))*2;
          if (k == 0) {
            duplicateRate_vs_theta[ifile]->Fill(false, std::cos(theta));
            duplicateRate_vs_pt[ifile]->Fill(false, pt);
	    if(qop>0){
              duplicateRate_vs_theta_plus[ifile]->Fill(false, std::cos(theta));
              duplicateRate_vs_pt_plus[ifile]->Fill(false, pt);
            } else {
              duplicateRate_vs_theta_minus[ifile]->Fill(false, std::cos(theta));
              duplicateRate_vs_pt_minus[ifile]->Fill(false, pt);
	    } 
	  } else {
            //std::cout<<"Find duplicate track for event " << ievent << ", pt " << pt << std::endl; 
	    duplicateRate_vs_theta[ifile]->Fill(true, std::cos(theta));
            duplicateRate_vs_pt[ifile]->Fill(true, pt);
	    if(qop>0){
	      duplicateRate_vs_theta_plus[ifile]->Fill(true, std::cos(theta));
              duplicateRate_vs_pt_plus[ifile]->Fill(true, pt);
            } else {
	      duplicateRate_vs_theta_minus[ifile]->Fill(true, std::cos(theta));
              duplicateRate_vs_pt_minus[ifile]->Fill(true, pt);
	    } 
          }
        }
      }  // end of all selected truth-matched tracks


      //std::cout<<"Finish duplicat rate plots" << std::endl;

      //////////////////////////////////////////
      // Loop over all truth particles in this event
      // The effiency is define as the ratio of selected particles that have
      // been matched with reco
      //////////////////////////////////////////
      for (const auto& particle : mapEventIdParticles[ievent]) {
        auto nHits = particle.nHits;
        auto eta = particle.eta;
        double theta = std::atan(std::exp(-eta))*2;
        auto pt = particle.pt;
        float q = particle.q;
	if (absPdgId!=999 and abs(particle.particlePdg) != absPdgId){
          continue;	
	}	
        if (abs(std::cos(theta)) > absCosTheta){
          continue;	
	}	
	if (nHits < nHitsMin or pt < ptMin ) {
          continue;
        }
        uint64_t pid = particle.particleId;
	

        // Fill the efficiency plots
        bool found = false; 
	auto ip = mapParticleIdFilteredTracks.find(pid);
        //if (ip != mapParticleIdFilteredTracks.end()) found = true; 
        if (ip != mapParticleIdFilteredTracks.end()) {
           found = true;
        	
           auto tracks = (*ip).second;		
           //std::cout<<"tracks size "<< tracks.size() << std::endl;  
	   double nOutliers = tracks[0].nOutliers; 
	   auto nMajorityHits = tracks[0].nMajorityHits; 
	   //std::cout<<"nMeasurements = " << nMeasurements << std::endl; 

            if(abs(particle.particlePdg)==211 and std::abs(std::cos(theta)) >0.85){
              if(nMajorityHits<6) found = false;
	    }

           //if(nMajorityHits*1./nHits<0.5) found = false;
        } 


	trackEff_vs_theta[ifile]->Fill(found, std::cos(theta));
        trackEff_vs_pt[ifile]->Fill(found, pt);
        trackEff_2d[ifile]->Fill(found, pt, std::cos(theta)); 

	nParticleHits_vs_theta[ifile]->Fill(std::cos(theta), nHits, 1.);	    
        nParticleHits_vs_pt[ifile]->Fill(pt, nHits, 1.);	    
	
	if(q>0){ 
	  trackEff_vs_theta_plus[ifile]->Fill(found, std::cos(theta));
          trackEff_vs_pt_plus[ifile]->Fill(found, pt);
          trackEff_2d_plus[ifile]->Fill(found, pt, std::cos(theta)); 

	  nParticleHits_vs_theta_plus[ifile]->Fill(std::cos(theta), nHits, 1.);	    
          nParticleHits_vs_pt_plus[ifile]->Fill(pt, nHits, 1.);	    
	
	} else {
	  trackEff_vs_theta_minus[ifile]->Fill(found, std::cos(theta));
          trackEff_vs_pt_minus[ifile]->Fill(found, pt);
          trackEff_2d_minus[ifile]->Fill(found, pt, std::cos(theta)); 
          
	  nParticleHits_vs_theta_minus[ifile]->Fill(std::cos(theta), nHits, 1.);	    
          nParticleHits_vs_pt_minus[ifile]->Fill(pt, nHits, 1.);	    
	
	}
//        } else {
//          //std::cout<<"preSelectedRecoTracks.size() " << preSelectedRecoTracks.size() << std::endl; 
//	  //std::cout<<" mapParticleIdFilteredTracks.size() " <<  mapParticleIdFilteredTracks.size() << std::endl; 
//	  //std::cout<<"Particle in event " << ievent <<" with pt " << pt <<" theta = " << theta*180/M_PI << ", particle_id = " << pid << " not reconstructed " << std::endl;
//          
//          trackEff_vs_theta[ifile]->Fill(false, std::cos(theta));
//          trackEff_vs_pt[ifile]->Fill(false, pt);
//          trackEff_2d[ifile]->Fill(false, pt, std::cos(theta)); 
//            
//	  nParticleHits_vs_theta[ifile]->Fill(std::cos(theta), nHits, 1.);	    
//          nParticleHits_vs_pt[ifile]->Fill(pt, nHits, 1.);	    
//          
//	  if(q>0){ 
//            trackEff_vs_theta_plus[ifile]->Fill(false, std::cos(theta));
//            trackEff_vs_pt_plus[ifile]->Fill(false, pt);
//            trackEff_2d_plus[ifile]->Fill(false, pt, std::cos(theta)); 
//            
//	    nParticleHits_vs_theta_plus[ifile]->Fill(std::cos(theta), nHits, 1.);	    
//            nParticleHits_vs_pt_plus[ifile]->Fill(pt, nHits, 1.);	    
//          } else {
//            trackEff_vs_theta_minus[ifile]->Fill(false, std::cos(theta));
//            trackEff_vs_pt_minus[ifile]->Fill(false, pt);
//            trackEff_2d_minus[ifile]->Fill(false, pt, std::cos(theta)); 
//	    
//	    nParticleHits_vs_theta_minus[ifile]->Fill(std::cos(theta), nHits, 1.);	    
//            nParticleHits_vs_pt_minus[ifile]->Fill(pt, nHits, 1.);	    
//	  }
      }  // end of all particles

    }  // end of all events
  }    // end of all track files


  if(useBESIIITracking){
    //TFile* file = TFile::Open("perf_besIII_baseline_realBeamPosition_nonUniformField_NodrdzCuts_nMeasurementsMin0.5.root","READ");
    TFile* file = TFile::Open("perf_besIII_baseline_realBeamPosition_nonUniformField.root","READ");
    //TFile* file = TFile::Open("perf_besIII_baseline_beamPosition0_uniformField.root","READ");
    trackEff_vs_theta[0] = (TEfficiency*)file->Get(Form("track_eff_vs_costhta_%s", particle.c_str())); 
    trackEff_vs_pt[0] = (TEfficiency*)file->Get(Form("track_eff_vs_pt_%s", particle.c_str())); 

    setEffStyle(trackEff_vs_theta[0], colors_minus[0], markers_minus[0]);
    setEffStyle(trackEff_vs_pt[0], colors_minus[0], markers_minus[0]);

    trackEff_vs_pt[0]->SetTitle(";Truth p_{T} [GeV/c];Efficiency");
    trackEff_vs_theta[0]->SetTitle(";Truth cos#theta [GeV/c];Efficiency");

    
    //This does not work well 
    //replaceEff(trackEff_vs_theta[0], "BESIII_pipijpsi_"+particle+"_eff_vs_costheta.txt");
    //replaceEff(trackEff_vs_pt[0], "BESIII_pipijpsi_"+particle+"_eff_vs_pt.txt");
  }



  std::cout<<"Start fill ratio hists" << std::endl;
  for(int i=1; i < nTrackFiles; ++i){
    getRatio(trackEff_vs_theta[0], trackEff_vs_theta[i], ratio_trackEff_vs_theta[i-1]);
    getRatio(trackEff_vs_pt[0], trackEff_vs_pt[i], ratio_trackEff_vs_pt[i-1]);
    
    getRatio(trackEff_vs_theta_minus[0], trackEff_vs_theta_minus[i], ratio_trackEff_vs_theta_minus[i-1]);
    getRatio(trackEff_vs_pt_minus[0], trackEff_vs_pt_minus[i], ratio_trackEff_vs_pt_minus[i-1]);
    
    getRatio(trackEff_vs_theta_plus[0], trackEff_vs_theta_plus[i], ratio_trackEff_vs_theta_plus[i-1]);
    getRatio(trackEff_vs_pt_plus[0], trackEff_vs_pt_plus[i], ratio_trackEff_vs_pt_plus[i-1]);
  } 

  std::cout << "All good. Now plotting..." << std::endl;

  // The legends
  std::vector<TLegend*> legs;
  for (int i = 0; i < 12; ++i) {
    //TLegend* legend = new TLegend(0.6, 0.65, 0.9, 0.8);
    TLegend* legend = new TLegend(0.5, 0.58, 0.9, 0.8);
    //legend->SetName("#psi(2S)#rightarrow#pi^{+}#pi^{-}J/#psi, J/#psi#rightarrow#mu^{+}#mu^{-}"); 
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legs.push_back(legend);
  }
  // Add entry for the legends
  for (int i = 0; i < nTrackFiles; ++i) {
/*

    legs[0]->AddEntry(trackEff_vs_theta_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[1]->AddEntry(fakeRate_vs_theta_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[2]->AddEntry(duplicateRate_vs_theta_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[3]->AddEntry(trackEff_vs_pt_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[4]->AddEntry(fakeRate_vs_pt_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[5]->AddEntry(duplicateRate_vs_pt_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[6]->AddEntry(trackPurity_vs_pt_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[7]->AddEntry(trackPurity_vs_theta_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[8]->AddEntry(nParticleHits_vs_pt_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
    legs[9]->AddEntry(nParticleHits_vs_theta_plus[i],
                      Form("%s", trackSummaryFileLegends_plus[i].c_str()), "lp");
 
    legs[0]->AddEntry(trackEff_vs_theta_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[1]->AddEntry(fakeRate_vs_theta_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[2]->AddEntry(duplicateRate_vs_theta_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[3]->AddEntry(trackEff_vs_pt_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[4]->AddEntry(fakeRate_vs_pt_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[5]->AddEntry(duplicateRate_vs_pt_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[6]->AddEntry(trackPurity_vs_pt_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[7]->AddEntry(trackPurity_vs_theta_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[8]->AddEntry(nParticleHits_vs_pt_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[9]->AddEntry(nParticleHits_vs_theta_minus[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
 */ 
 

    legs[0]->AddEntry(trackEff_vs_theta[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[1]->AddEntry(fakeRate_vs_theta[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[2]->AddEntry(duplicateRate_vs_theta[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[3]->AddEntry(trackEff_vs_pt[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[4]->AddEntry(fakeRate_vs_pt[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[5]->AddEntry(duplicateRate_vs_pt[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[6]->AddEntry(trackPurity_vs_pt[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[7]->AddEntry(trackPurity_vs_theta[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[8]->AddEntry(nParticleHits_vs_pt[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");
    legs[9]->AddEntry(nParticleHits_vs_theta[i],
                      Form("%s", trackSummaryFileLegends_minus[i].c_str()), "lp");

  }


  // Define the canvas etc.
  std::vector<TCanvas*> cans(2);
  std::vector<TPad*> tpads(2);
  std::vector<TPad*> bpads(2);
  std::vector<TPad*> pads(2);
  std::vector<TLine*> lines(2);

  for(int i=0; i<2; ++i){
    cans[i] = new TCanvas(Form("%i",i), "", 700, 600);
	  
    tpads[i] = new TPad(Form("tpad_%i",i), "", 0, 0.36, 1, 1.0);
    bpads[i] = new TPad(Form("bpad_%i",i), "", 0, 0.05, 1, 0.35);
    pads[i] = new TPad(Form("pad_%i",i), "", 0, 0.05, 1, 1.0);

    tpads[i]->SetRightMargin(0.05);
    tpads[i]->SetBottomMargin(0.005);
    bpads[i]->SetTopMargin(0.005);
    bpads[i]->SetRightMargin(0.05);
    bpads[i]->SetBottomMargin(0.4);
    pads[i]->SetBottomMargin(0.15);

    if(i==0)
      lines[i] = new TLine(-1, 1., 1, 1.);
    if(i==1)
      lines[i] = new TLine(ptRanges[1], 1., ptRanges[2], 1.);
    lines[i]->SetLineWidth(2);
    lines[i]->SetLineStyle(kDashed);
  } 

  // eff_vs_costheta
  cans[0]->cd();
  tpads[0]->Draw();
  bpads[0]->Draw();
  tpads[0]->cd();
  for (int i = 0; i < nTrackFiles; ++i) {
    std::string mode = (i == 0) ? "" : "same";
    trackEff_vs_theta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[0]->Draw();
    }
    //myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
    myText(0.5, 0.85, 1, myLabelSize, myLabel.c_str());
    if(particle=="pi"){ 
      //setEffRange(trackEff_vs_theta[i], -1, 1, 0.801,  1.20, topXTitleSize, topYTitleSize,
      setEffRange(trackEff_vs_theta[i], -1, 1, 0.701,  1.20, topXTitleSize, topYTitleSize,
                         topXLabelSize, topYLabelSize, topXTitleOffset, topYTitleOffset);
    } else {
      //setEffRange(trackEff_vs_theta[i], -1, 1, 0.921,  1.06, topXTitleSize, topYTitleSize,
      setEffRange(trackEff_vs_theta[i], -1, 1, 0.91,  1.09, topXTitleSize, topYTitleSize,
                         topXLabelSize, topYLabelSize, topXTitleOffset, topYTitleOffset);
    }
  }
  bpads[0]->cd();
  for(int i=1; i<nTrackFiles; ++i){
    ratio_trackEff_vs_theta[i-1]->Draw("E1same");
    ratio_trackEff_vs_theta[i-1]->GetYaxis()->SetRangeUser(0.99,1.099);
    ratio_trackEff_vs_theta[i-1]->GetYaxis()->SetTitle("(MDC+Pixel)/MDC");
  } 
  lines[0]->Draw();
  //cans[0]->Update();

  //eff_vs_pt 
  cans[1]->cd();
  tpads[1]->Draw();
  bpads[1]->Draw();
  tpads[1]->cd();
  for (int i = 0; i < nTrackFiles; ++i) {
    std::string mode = (i == 0) ? "" : "same";
    trackEff_vs_pt[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[1]->Draw();
    }
    //myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
    myText(0.5, 0.85, 1, myLabelSize, myLabel.c_str());
    if(particle=="pi"){
      setEffRange(trackEff_vs_pt[i], ptRanges[1], ptRanges[2], 0.925,  1.06, topXTitleSize, topYTitleSize,
                         topXLabelSize, topYLabelSize, topXTitleOffset, topYTitleOffset);
    } else {
      setEffRange(trackEff_vs_pt[i], ptRanges[1], ptRanges[2], 0.839,  1.15, topXTitleSize, topYTitleSize,
                         topXLabelSize, topYLabelSize, topXTitleOffset, topYTitleOffset);
    }
  }
  bpads[1]->cd();
  for(int i=1; i<nTrackFiles; ++i){
    ratio_trackEff_vs_pt[i-1]->Draw("E1same");
    if(particle == "pi"){
      //ratio_trackEff_vs_pt[i-1]->GetYaxis()->SetRangeUser(0.99,1.07); //pt>100 MeV
      ratio_trackEff_vs_pt[i-1]->GetYaxis()->SetRangeUser(0.9,1.45); //pt>150 MeV
    } else {
      ratio_trackEff_vs_pt[i-1]->GetYaxis()->SetRangeUser(0.95,1.19);
    } 
    ratio_trackEff_vs_pt[i-1]->GetYaxis()->SetTitle("(MDC+Pixel)/MDC");
  }
  lines[1]->Draw();


// Now draw the plots (no ratio plots)
//  std::vector<TCanvas*> cs;
//  for(int i=0; i < 14; ++i){
//    cs.push_back(new TCanvas(Form("c_%i", i), "", 600, 500)); 
//    cs[i]->SetGrid(); 
//  }
//
//  float scaleRangeMax = 1.2;
//  for (int i = 0; i < nTrackFiles; ++i) {
//    std::string mode = (i == 0) ? "" : "same";
//    cs[0]->cd();
//    trackEff_vs_theta_minus[i]->Draw(mode.c_str());
//    trackEff_vs_theta_plus[i]->Draw("same");
//    if (i == nTrackFiles - 1) {
//      legs[0]->Draw();
//    }
//    myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//    //adaptEffRange(trackEff_vs_theta_minus[i], 1, scaleRangeMax);
//    if(particle=="pi"){ 
//      setEffRange(trackEff_vs_theta_minus[i], 0.92,  1.04);
//    } else {
//      setEffRange(trackEff_vs_theta_minus[i], 0.94,  1.06);
//    }
//
//    cs[1]->cd();
//    fakeRate_vs_theta_minus[i]->Draw(mode.c_str());
//    fakeRate_vs_theta_plus[i]->Draw("same");
//    if (i == nTrackFiles - 1) {
//      legs[1]->Draw();
//    }
//    myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//    //adaptEffRange(fakeRate_vs_theta_minus[i], 1, scaleRangeMax);
//    if(particle=="pi"){ 
//      setEffRange(fakeRate_vs_theta_minus[i], 0, 0.04);
//    } else {
//      setEffRange(fakeRate_vs_theta_minus[i], 0, 0.02);
//    }

//    cs[2]->cd();
//    duplicateRate_vs_theta_minus[i]->Draw(mode.c_str());
//    duplicateRate_vs_theta_plus[i]->Draw("same");
//    //duplicateRate_vs_theta_plus[i]->GetYAxis()->SetRangeUser(0.,0.0035); 
//    //duplicateRate_vs_theta_minus[i]->GetYAxis()->SetRangeUser(0.,0.0035); 
//    if (i == nTrackFiles - 1) {
//      legs[2]->Draw();
//    }
//    myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//    //adaptEffRange(duplicateRate_vs_theta_minus[i], 1, scaleRangeMax);
//    if(particle=="pi"){ 
//      setEffRange(duplicateRate_vs_theta_minus[i], 0, 0.05);
//    } else {
//      setEffRange(duplicateRate_vs_theta_minus[i], 0, 0.02);
//    }
//
//    cs[3]->cd();
//    trackEff_vs_pt_minus[i]->Draw(mode.c_str());
//    trackEff_vs_pt_plus[i]->Draw("same");
//    //trackEff_vs_pt_plus[i]->GetYAxis()->SetRangeUser(0.92,1.09); 
//    //trackEff_vs_pt_minus[i]->GetYAxis()->SetRangeUser(0.92,1.09); 
//    if (i == nTrackFiles - 1) {
//      legs[3]->Draw();
//    }
//    myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//    //adaptEffRange(trackEff_vs_pt_minus[i], 1, scaleRangeMax);
//    if(particle=="pi"){ 
//      setEffRange(trackEff_vs_pt_minus[i], 0.92, 1.09);
//    } else {
//      setEffRange(trackEff_vs_pt_minus[i], 0.94, 1.06);
//    }
//
//    cs[4]->cd();
//    fakeRate_vs_pt_minus[i]->Draw(mode.c_str());
//    fakeRate_vs_pt_plus[i]->Draw("same");
//    //fakeRate_vs_pt_plus[i]->GetYAxis()->SetRangeUser(0.,0.002); 
//    //fakeRate_vs_pt_minus[i]->GetYAxis()->SetRangeUser(0.,0.002); 
//    if (i == nTrackFiles - 1) {
//      legs[4]->Draw();
//    }
//    myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//    //adaptEffRange(fakeRate_vs_pt_minus[i], 1, scaleRangeMax);
//    if(particle=="pi"){ 
//      setEffRange(fakeRate_vs_pt_minus[i], 0.005, 0.035);
//    } else {
//      setEffRange(fakeRate_vs_pt_minus[i], 0, 0.02);
//    }
//
//    cs[5]->cd();
//    duplicateRate_vs_pt_minus[i]->Draw(mode.c_str());
//    duplicateRate_vs_pt_plus[i]->Draw("same");
//    //duplicateRate_vs_pt_plus[i]->GetYAxis()->SetRangeUser(0.,0.01); 
//    //duplicateRate_vs_pt_minus[i]->GetYAxis()->SetRangeUser(0.,0.01); 
//    if (i == nTrackFiles - 1) {
//      legs[5]->Draw();
//    }
//    myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//    //adaptEffRange(duplicateRate_vs_pt_minus[i], 1, scaleRangeMax);
//    if(particle=="pi"){ 
//      setEffRange(duplicateRate_vs_pt_minus[i], 0, 0.04);
//    } else {
//      setEffRange(duplicateRate_vs_pt_minus[i], 0, 0.02);
//    } 


//     cs[6]->cd();
//     trackPurity_vs_pt_minus[i]->Draw("same"); 
//     trackPurity_vs_pt_plus[i]->Draw("same"); 
//     if (i == nTrackFiles - 1) {
//      legs[6]->Draw();
//     }
//     myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//     trackPurity_vs_pt_minus[i]->GetYaxis()->SetRangeUser(0.93, 1.0); 
//    
//
//
//     cs[7]->cd();
//     trackPurity_vs_theta_minus[i]->Draw("same"); 
//     trackPurity_vs_theta_plus[i]->Draw("same"); 
//     if (i == nTrackFiles - 1) {
//      legs[7]->Draw();
//     }
//     myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//     trackPurity_vs_theta_minus[i]->GetYaxis()->SetRangeUser(0.94, 1.0); 


//     cs[8]->cd();
//     nParticleHits_vs_pt_minus[i]->Draw("same");
//     nParticleHits_vs_pt_plus[i]->Draw("same");
//     if (i == nTrackFiles - 1) {
//      legs[8]->Draw();
//     }
//     myText(0.5, 0.85, 1, 0.04, myLabel.c_str());
//
//
//
//     cs[9]->cd();
//     nParticleHits_vs_theta_minus[i]->Draw("same");
//     nParticleHits_vs_theta_plus[i]->Draw("same");
//     if (i == nTrackFiles - 1) {
//      legs[9]->Draw();
//     }
//     myText(0.5, 0.85, 1, 0.04, myLabel.c_str());

//2d plots
     
//     cs[10]->cd();
//     trackEff_2d_plus[i]->Draw("colz");
//     myText(0.45, 0.85, 1, 0.04, myLabel.c_str());
//     myText(0.7, 0.8, 1, 0.04, trackSummaryFileLegends_plus[0].c_str());
//   
//     cs[11]->cd();
//     trackEff_2d_minus[i]->Draw("colz");
//     myText(0.45, 0.85, 1, 0.04, myLabel.c_str());
//     myText(0.7, 0.8, 1, 0.04, trackSummaryFileLegends_minus[0].c_str());
//
//
//     cs[12]->cd();
//     trackPurity_2d_plus[i]->Draw("colz");
//     myText(0.45, 0.85, 1, 0.04, myLabel.c_str());
//     myText(0.7, 0.8, 1, 0.04, trackSummaryFileLegends_plus[0].c_str());
//
//     cs[13]->cd();
//     trackPurity_2d_minus[i]->Draw("colz");
//     myText(0.45, 0.85, 1, 0.04, myLabel.c_str());
//     myText(0.7, 0.8, 1, 0.04, trackSummaryFileLegends_minus[0].c_str());
    
//  }





  std::string path_str=gSystem->pwd();
  TTimeStamp tts;
  std::string tts_str(100, '\n');
  std::snprintf(tts_str.data(), tts_str.size(), "D%.8u-T%.6u" , tts.GetDate(),tts.GetTime());
  tts_str=tts_str.c_str();
  path_str = path_str+"/plot/"+tts_str;
  std::string cmd_mkdir_str = "mkdir -p ";
  cmd_mkdir_str = cmd_mkdir_str+path_str;
  //gSystem->Exec(cmd_mkdir_str.c_str());

  std::string path = "plot/" + tag;
  if(savePlot){
    gSystem->Exec(Form("mkdir -p %s", path.c_str()));
  }

  //Save to text
  std::string fileName = path + "/BESIII_pipijpsi_" + particle + "_eff_vs_costheta.txt";
  std::cout<<"Generating file to " << fileName << std::endl; 
  std::ofstream file_eff_vs_theta;
  file_eff_vs_theta.open(fileName.c_str());
  file_eff_vs_theta.setf(std::ios_base::fixed, std::ios_base::floatfield);
  file_eff_vs_theta << std::setprecision(4);

  file_eff_vs_theta <<"#Efficiency_vs_costheta"<< std::endl;
  file_eff_vs_theta <<"bin,pt,MDC,MDC+Pixel(r=35mm),MDC+Pixel(r=45mm),MDC+Pixel(r=55mm)"<< std::endl;
  for(int i=0; i<20; ++i){
    double binCenter = -1 + 0.1*i + 0.05;
    double baseline = trackEff_vs_theta[0]->GetEfficiency(i+1);
    double r1 = trackEff_vs_theta[1]->GetEfficiency(i+1);
    double r2 = trackEff_vs_theta[2]->GetEfficiency(i+1);
    double r3 = trackEff_vs_theta[3]->GetEfficiency(i+1);
    file_eff_vs_theta <<i <<","<< binCenter << "," << baseline <<","<< r1 <<"," << r2 <<"," << r3 << std::endl;
  }
  file_eff_vs_theta.close();


  fileName = path + "/BESIII_pipijpsi_" + particle + "_eff_vs_pt.txt";
  std::cout<<"Generating file to " << fileName << std::endl; 
  std::ofstream file_eff_vs_pt;
  file_eff_vs_pt.open(fileName.c_str());
  file_eff_vs_pt.setf(std::ios_base::fixed, std::ios_base::floatfield);
  file_eff_vs_pt << std::setprecision(4);

  file_eff_vs_pt <<"#Efficiency_vs_pt"<< std::endl;
  file_eff_vs_pt <<"bin,pt,MDC,MDC+Pixel(r=35mm),MDC+Pixel(r=45mm),MDC+Pixel(r=55mm)"<< std::endl;
  for(int i=0; i<ptRanges[0]; ++i){
    double binCenter = ptRanges[1] + (ptRanges[2]-ptRanges[1])/ptRanges[0]*i + (ptRanges[2]-ptRanges[1])/ptRanges[0]/2;
    double baseline = trackEff_vs_pt[0]->GetEfficiency(i+1);
    double r1 = trackEff_vs_pt[1]->GetEfficiency(i+1);
    double r2 = trackEff_vs_pt[2]->GetEfficiency(i+1);
    double r3 = trackEff_vs_pt[3]->GetEfficiency(i+1);
    file_eff_vs_pt <<i <<","<< binCenter << "," << baseline <<","<< r1 <<"," << r2 <<"," << r3 << std::endl;
  }
  file_eff_vs_pt.close();



 if(savePlot){

   cans[0]->SaveAs(Form("%s/BESIII_pipijpsi_%s_eff_vs_costheta.pdf", path.c_str(), particle.c_str()));
   cans[1]->SaveAs(Form("%s/BESIII_pipijpsi_%s_eff_vs_pt.pdf", path.c_str(), particle.c_str()));

//     cs[0]->SaveAs(Form("%s/BESIII_pipijpsi_%s_eff_vs_costheta.pdf", path.c_str(), particle.c_str()));
//     cs[1]->SaveAs(Form("%s/BESIII_pipijpsi_%s_fakerate_vs_costheta.pdf", path.c_str(),particle.c_str()));
//     cs[2]->SaveAs(Form("%s/BESIII_pipijpsi_%s_duplirate_vs_costheta.pdf", path.c_str(),particle.c_str()));
//     cs[3]->SaveAs(Form("%s/BESIII_pipijpsi_%s_eff_vs_pt.pdf", path.c_str(),particle.c_str()));
//     cs[4]->SaveAs(Form("%s/BESIII_pipijpsi_%s_fakerate_vs_pt.pdf", path.c_str(),particle.c_str()));
//     cs[5]->SaveAs(Form("%s/BESIII_pipijpsi_%s_duplirate_vs_pt.pdf", path.c_str(),particle.c_str()));
//     cs[6]->SaveAs(Form("%s/BESIII_pipijpsi_%s_trackPurity_vs_pt.pdf", path.c_str(),particle.c_str()));
//     cs[7]->SaveAs(Form("%s/BESIII_pipijpsi_%s_trackPurity_vs_theta.pdf", path.c_str(),particle.c_str()));
//     cs[8]->SaveAs(Form("%s/BESIII_pipijpsi_%s_nParticleHits_vs_pt.pdf", path.c_str(),particle.c_str()));
//     cs[9]->SaveAs(Form("%s/BESIII_pipijpsi_%s_nParticleHits_vs_theta.pdf", path.c_str(),particle.c_str()));
//     cs[10]->SaveAs(Form("%s/BESIII_pipijpsi_%s_eff_2D_plus.pdf", path.c_str(),particle.c_str()));
//     cs[11]->SaveAs(Form("%s/BESIII_pipijpsi_%s_eff_2D_minus.pdf", path.c_str(),particle.c_str()));
//     cs[12]->SaveAs(Form("%s/BESIII_pipijpsi_%s_trackPurity_2D_plus.pdf", path.c_str(),particle.c_str()));
//     cs[13]->SaveAs(Form("%s/BESIII_pipijpsi_%s_trackPurity_2D_minus.pdf", path.c_str(),particle.c_str()));
  }
}
