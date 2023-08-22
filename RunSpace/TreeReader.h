// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TTree.h>

struct ParticleInfo {
  ULong64_t particleId = 0;
  int particlePdg = 0;
  double eta = 0;
  double p = 0;
  double pt = 0;
  double q = 0;
  UShort_t nHits = 0;

};



struct RecoTrackInfo {
  double eta;
  double pt;
  double qop; 
  unsigned int nMajorityHits;
  unsigned int nMeasurements;
  unsigned int nOutliers;
  unsigned int nSharedHits;

  
  uint64_t majorityParticleId;
  int majorityParticlePdg;
  double majorityParticlePt;
  double majorityParticleEta;
  std::vector<uint64_t> measurementIndex;
  std::vector<uint64_t> twinsMeasurementIndex;
  std::vector<double> res_fit = {0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> pull_fit = {0, 0, 0, 0, 0, 0};
};


/// Helper for reading tree
///
struct TreeReader {
  // The constructor
  TreeReader(TTree* tree_) : tree(tree_){};

  // Get entry
  void getEntry(unsigned int i) const {
    if (entryNumbers.size() > i) {
      tree->GetEntry(entryNumbers[i]);
    } else {
      tree->GetEntry(i);
    }
  };

  // The tree being read
  TTree* tree = nullptr;

 protected:
  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> entryNumbers = {};
};

/// Struct used for reading track states written out by the
/// RootTrajectoryStatesWriter
///
struct TrackStatesReader : public TreeReader {
  // Delete the default constructor
  TrackStatesReader() = delete;

  // The constructor
  TrackStatesReader(TTree* tree_, bool sortEvents) : TreeReader(tree_) {
    tree->SetBranchAddress("event_nr", &eventId);
    tree->SetBranchAddress("eLOC0_prt", &LOC0_prt);
    tree->SetBranchAddress("eLOC1_prt", &LOC1_prt);
    tree->SetBranchAddress("ePHI_prt", &PHI_prt);
    tree->SetBranchAddress("eTHETA_prt", &THETA_prt);
    tree->SetBranchAddress("eQOP_prt", &QOP_prt);
    tree->SetBranchAddress("eT_prt", &T_prt);
    tree->SetBranchAddress("eLOC0_flt", &LOC0_flt);
    tree->SetBranchAddress("eLOC1_flt", &LOC1_flt);
    tree->SetBranchAddress("ePHI_flt", &PHI_flt);
    tree->SetBranchAddress("eTHETA_flt", &THETA_flt);
    tree->SetBranchAddress("eQOP_flt", &QOP_flt);
    tree->SetBranchAddress("eT_flt", &T_flt);
    tree->SetBranchAddress("eLOC0_smt", &LOC0_smt);
    tree->SetBranchAddress("eLOC1_smt", &LOC1_smt);
    tree->SetBranchAddress("ePHI_smt", &PHI_smt);
    tree->SetBranchAddress("eTHETA_smt", &THETA_smt);
    tree->SetBranchAddress("eQOP_smt", &QOP_smt);
    tree->SetBranchAddress("eT_smt", &T_smt);

    tree->SetBranchAddress("res_eLOC0_prt", &res_LOC0_prt);
    tree->SetBranchAddress("res_eLOC1_prt", &res_LOC1_prt);
    tree->SetBranchAddress("res_ePHI_prt", &res_PHI_prt);
    tree->SetBranchAddress("res_eTHETA_prt", &res_THETA_prt);
    tree->SetBranchAddress("res_eQOP_prt", &res_QOP_prt);
    tree->SetBranchAddress("res_eT_prt", &res_T_prt);
    tree->SetBranchAddress("res_eLOC0_flt", &res_LOC0_flt);
    tree->SetBranchAddress("res_eLOC1_flt", &res_LOC1_flt);
    tree->SetBranchAddress("res_ePHI_flt", &res_PHI_flt);
    tree->SetBranchAddress("res_eTHETA_flt", &res_THETA_flt);
    tree->SetBranchAddress("res_eQOP_flt", &res_QOP_flt);
    tree->SetBranchAddress("res_eT_flt", &res_T_flt);
    tree->SetBranchAddress("res_eLOC0_smt", &res_LOC0_smt);
    tree->SetBranchAddress("res_eLOC1_smt", &res_LOC1_smt);
    tree->SetBranchAddress("res_ePHI_smt", &res_PHI_smt);
    tree->SetBranchAddress("res_eTHETA_smt", &res_THETA_smt);
    tree->SetBranchAddress("res_eQOP_smt", &res_QOP_smt);
    tree->SetBranchAddress("res_eT_smt", &res_T_smt);

    tree->SetBranchAddress("pull_eLOC0_prt", &pull_LOC0_prt);
    tree->SetBranchAddress("pull_eLOC1_prt", &pull_LOC1_prt);
    tree->SetBranchAddress("pull_ePHI_prt", &pull_PHI_prt);
    tree->SetBranchAddress("pull_eTHETA_prt", &pull_THETA_prt);
    tree->SetBranchAddress("pull_eQOP_prt", &pull_QOP_prt);
    tree->SetBranchAddress("pull_eT_prt", &pull_T_prt);
    tree->SetBranchAddress("pull_eLOC0_flt", &pull_LOC0_flt);
    tree->SetBranchAddress("pull_eLOC1_flt", &pull_LOC1_flt);
    tree->SetBranchAddress("pull_ePHI_flt", &pull_PHI_flt);
    tree->SetBranchAddress("pull_eTHETA_flt", &pull_THETA_flt);
    tree->SetBranchAddress("pull_eQOP_flt", &pull_QOP_flt);
    tree->SetBranchAddress("pull_eT_flt", &pull_T_flt);
    tree->SetBranchAddress("pull_eLOC0_smt", &pull_LOC0_smt);
    tree->SetBranchAddress("pull_eLOC1_smt", &pull_LOC1_smt);
    tree->SetBranchAddress("pull_ePHI_smt", &pull_PHI_smt);
    tree->SetBranchAddress("pull_eTHETA_smt", &pull_THETA_smt);
    tree->SetBranchAddress("pull_eQOP_smt", &pull_QOP_smt);
    tree->SetBranchAddress("pull_eT_smt", &pull_T_smt);

    tree->SetBranchAddress("g_x_prt", &g_x_prt);
    tree->SetBranchAddress("g_y_prt", &g_y_prt);
    tree->SetBranchAddress("g_z_prt", &g_z_prt);
    tree->SetBranchAddress("g_x_flt", &g_x_flt);
    tree->SetBranchAddress("g_y_flt", &g_y_flt);
    tree->SetBranchAddress("g_z_flt", &g_z_flt);
    tree->SetBranchAddress("g_x_smt", &g_x_smt);
    tree->SetBranchAddress("g_y_smt", &g_y_smt);
    tree->SetBranchAddress("g_z_smt", &g_z_smt);

    tree->SetBranchAddress("nStates", &nStates);
    tree->SetBranchAddress("nMeasurements", &nMeasurements);
    tree->SetBranchAddress("volume_id", &volume_id);
    tree->SetBranchAddress("layer_id", &layer_id);
    tree->SetBranchAddress("module_id", &module_id);
    tree->SetBranchAddress("predicted", &predicted);
    tree->SetBranchAddress("filtered", &filtered);
    tree->SetBranchAddress("smoothed", &smoothed);

    // It's not necessary if you just need to read one file, but please do it to
    // synchronize events if multiple root files are read
    if (sortEvents) {
      entryNumbers.resize(tree->GetEntries());
      tree->Draw("event_nr", "", "goff");
      // Sort to get the entry numbers of the ordered events
      TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(),
                  false);
    }
  }

  // The variables
  uint32_t eventId;
  std::vector<float>* LOC0_prt =
      new std::vector<float>;  ///< predicted parameter local x
  std::vector<float>* LOC1_prt =
      new std::vector<float>;  ///< predicted parameter local y
  std::vector<float>* PHI_prt =
      new std::vector<float>;  ///< predicted parameter phi
  std::vector<float>* THETA_prt =
      new std::vector<float>;  ///< predicted parameter theta
  std::vector<float>* QOP_prt =
      new std::vector<float>;  ///< predicted parameter q/p
  std::vector<float>* T_prt =
      new std::vector<float>;  ///< predicted parameter t
  std::vector<float>* LOC0_flt =
      new std::vector<float>;  ///< filtered parameter local x
  std::vector<float>* LOC1_flt =
      new std::vector<float>;  ///< filtered parameter local y
  std::vector<float>* PHI_flt =
      new std::vector<float>;  ///< filtered parameter phi
  std::vector<float>* THETA_flt =
      new std::vector<float>;  ///< filtered parameter theta
  std::vector<float>* QOP_flt =
      new std::vector<float>;  ///< filtered parameter q/p
  std::vector<float>* T_flt = new std::vector<float>;  ///< filtered parameter t
  std::vector<float>* LOC0_smt =
      new std::vector<float>;  ///< smoothed parameter local x
  std::vector<float>* LOC1_smt =
      new std::vector<float>;  ///< smoothed parameter local y
  std::vector<float>* PHI_smt =
      new std::vector<float>;  ///< smoothed parameter phi
  std::vector<float>* THETA_smt =
      new std::vector<float>;  ///< smoothed parameter theta
  std::vector<float>* QOP_smt =
      new std::vector<float>;  ///< smoothed parameter q/p
  std::vector<float>* T_smt = new std::vector<float>;  ///< smoothed parameter t

  std::vector<float>* res_LOC0_prt =
      new std::vector<float>;  ///< residual of predicted parameter local x
  std::vector<float>* res_LOC1_prt =
      new std::vector<float>;  ///< residual of predicted parameter local y
  std::vector<float>* res_PHI_prt =
      new std::vector<float>;  ///< residual of predicted parameter phi
  std::vector<float>* res_THETA_prt =
      new std::vector<float>;  ///< residual of predicted parameter theta
  std::vector<float>* res_QOP_prt =
      new std::vector<float>;  ///< residual of predicted parameter q/p
  std::vector<float>* res_T_prt =
      new std::vector<float>;  ///< residual of predicted parameter t
  std::vector<float>* res_LOC0_flt =
      new std::vector<float>;  ///< residual of filtered parameter local x
  std::vector<float>* res_LOC1_flt =
      new std::vector<float>;  ///< residual of filtered parameter local y
  std::vector<float>* res_PHI_flt =
      new std::vector<float>;  ///< residual of filtered parameter phi
  std::vector<float>* res_THETA_flt =
      new std::vector<float>;  ///< residual of filtered parameter theta
  std::vector<float>* res_QOP_flt =
      new std::vector<float>;  ///< residual of filtered parameter q/p
  std::vector<float>* res_T_flt =
      new std::vector<float>;  ///< residual of filtered parameter t
  std::vector<float>* res_LOC0_smt =
      new std::vector<float>;  ///< residual of smoothed parameter local x
  std::vector<float>* res_LOC1_smt =
      new std::vector<float>;  ///< residual of smoothed parameter local y
  std::vector<float>* res_PHI_smt =
      new std::vector<float>;  ///< residual of smoothed parameter phi
  std::vector<float>* res_THETA_smt =
      new std::vector<float>;  ///< residual of smoothed parameter theta
  std::vector<float>* res_QOP_smt =
      new std::vector<float>;  ///< residual of smoothed parameter q/p
  std::vector<float>* res_T_smt =
      new std::vector<float>;  ///< residual of smoothed parameter t

  std::vector<float>* pull_LOC0_prt =
      new std::vector<float>;  ///< pull of predicted parameter local x
  std::vector<float>* pull_LOC1_prt =
      new std::vector<float>;  ///< pull of predicted parameter local y
  std::vector<float>* pull_PHI_prt =
      new std::vector<float>;  ///< pull of predicted parameter phi
  std::vector<float>* pull_THETA_prt =
      new std::vector<float>;  ///< pull of predicted parameter theta
  std::vector<float>* pull_QOP_prt =
      new std::vector<float>;  ///< pull of predicted parameter q/p
  std::vector<float>* pull_T_prt =
      new std::vector<float>;  ///< pull of predicted parameter t
  std::vector<float>* pull_LOC0_flt =
      new std::vector<float>;  ///< pull of filtered parameter local x
  std::vector<float>* pull_LOC1_flt =
      new std::vector<float>;  ///< pull of filtered parameter local y
  std::vector<float>* pull_PHI_flt =
      new std::vector<float>;  ///< pull of filtered parameter phi
  std::vector<float>* pull_THETA_flt =
      new std::vector<float>;  ///< pull of filtered parameter theta
  std::vector<float>* pull_QOP_flt =
      new std::vector<float>;  ///< pull of filtered parameter q/p
  std::vector<float>* pull_T_flt =
      new std::vector<float>;  ///< pull of filtered parameter t
  std::vector<float>* pull_LOC0_smt =
      new std::vector<float>;  ///< pull of smoothed parameter local x
  std::vector<float>* pull_LOC1_smt =
      new std::vector<float>;  ///< pull of smoothed parameter local y
  std::vector<float>* pull_PHI_smt =
      new std::vector<float>;  ///< pull of smoothed parameter phi
  std::vector<float>* pull_THETA_smt =
      new std::vector<float>;  ///< pull of smoothed parameter theta
  std::vector<float>* pull_QOP_smt =
      new std::vector<float>;  ///< pull of smoothed parameter q/p
  std::vector<float>* pull_T_smt =
      new std::vector<float>;  ///< pull of smoothed parameter t

  std::vector<float>* g_x_prt = new std::vector<float>;
  std::vector<float>* g_y_prt = new std::vector<float>;
  std::vector<float>* g_z_prt = new std::vector<float>;
  std::vector<float>* g_x_flt = new std::vector<float>;
  std::vector<float>* g_y_flt = new std::vector<float>;
  std::vector<float>* g_z_flt = new std::vector<float>;
  std::vector<float>* g_x_smt = new std::vector<float>;
  std::vector<float>* g_y_smt = new std::vector<float>;
  std::vector<float>* g_z_smt = new std::vector<float>;

  std::vector<int>* volume_id = new std::vector<int>;  ///< volume_id
  std::vector<int>* layer_id = new std::vector<int>;   ///< layer_id
  std::vector<int>* module_id = new std::vector<int>;  ///< module_id

  std::vector<bool>* predicted = new std::vector<bool>;  ///< prediction status
  std::vector<bool>* filtered = new std::vector<bool>;   ///< filtering status
  std::vector<bool>* smoothed = new std::vector<bool>;   ///< smoothing status

  unsigned int nStates, nMeasurements;
};


//// Struct used for reading seed info writen out by the TrackFinderPerformanceWriter
struct SeedReader : public TreeReader{
  SeedReader() = delete;

  SeedReader(TTree* tree_, bool sortEvents) : TreeReader(tree_) {
    tree->SetBranchAddress("event_id", &trkEventId);
    tree->SetBranchAddress("track_id", &trkTrackId);
    tree->SetBranchAddress("size", &trkNumHits);
    tree->SetBranchAddress("nparticles", &trkNumParticles);
    tree->SetBranchAddress("particle_id", &trkParticleId);
    tree->SetBranchAddress("particle_nhits_total", &trkParticleNumHitsTotal);
    tree->SetBranchAddress("particle_nhits_on_track", &trkParticleNumHitsOnTrack);
 
    // It's not necessary if you just need to read one file, but please do it to
    // synchronize events if multiple root files are read
    if (sortEvents) {
      entryNumbers.resize(tree->GetEntries());
      tree->Draw("event_id", "", "goff");
      // Sort to get the entry numbers of the ordered events
      TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(),
                  false);
    }
  }

  uint64_t trkEventId;
  uint64_t trkTrackId;
  // track content
  // number of hits on track
  uint32_t trkNumHits;
  // number of particles contained in the track
  uint32_t trkNumParticles;
  // track particle content; for each contributing particle, largest first
  std::vector<uint64_t>* trkParticleId = new std::vector<uint64_t>;
  // total number of hits generated by this particle
  std::vector<uint32_t>* trkParticleNumHitsTotal = new std::vector<uint32_t>;
  // number of hits within this track
  std::vector<uint32_t>* trkParticleNumHitsOnTrack = new std::vector<uint32_t>;
};


/// Struct used for reading track summary info written out by the
/// RootTrajectorySummaryWriter
///
struct TrackSummaryReader : public TreeReader {
  // Delete the default constructor
  TrackSummaryReader() = delete;

  // The constructor
  TrackSummaryReader(TTree* tree_, bool sortEvents, bool duplicateMDC = false) : TreeReader(tree_) {
    tree->SetBranchAddress("event_nr", &eventId);
    tree->SetBranchAddress("nStates", &nStates);
    tree->SetBranchAddress("nMeasurements", &nMeasurements);
    tree->SetBranchAddress("nOutliers", &nOutliers);
    tree->SetBranchAddress("nHoles", &nHoles);
    tree->SetBranchAddress("chi2Sum", &chi2Sum);
    tree->SetBranchAddress("measurementChi2", &measurementChi2);
    tree->SetBranchAddress("measurementIndex", &measurementIndex);
    if(duplicateMDC){
      tree->SetBranchAddress("twinsMeasurementIndex", &twinsMeasurementIndex);
    } 
    tree->SetBranchAddress("NDF", &NDF);
    tree->SetBranchAddress("measurementVolume", &measurementVolume);
    tree->SetBranchAddress("measurementLayer", &measurementLayer);
    tree->SetBranchAddress("outlierVolume", &outlierVolume);
    tree->SetBranchAddress("outlierLayer", &outlierLayer);
    tree->SetBranchAddress("nMajorityHits", &nMajorityHits);
    tree->SetBranchAddress("nSharedHits", &nSharedHits);
    tree->SetBranchAddress("majorityParticleId", &majorityParticleId);
    tree->SetBranchAddress("majorityParticlePdg", &majorityParticlePdg);

    tree->SetBranchAddress("hasFittedParams", &hasFittedParams);

    tree->SetBranchAddress("t_theta", &t_theta);
    tree->SetBranchAddress("t_phi", &t_phi);
    tree->SetBranchAddress("t_eta", &t_eta);
    tree->SetBranchAddress("t_p", &t_p);
    tree->SetBranchAddress("t_pT", &t_pT);
    tree->SetBranchAddress("t_d0", &t_d0);
    tree->SetBranchAddress("t_z0", &t_z0);
    tree->SetBranchAddress("t_charge",&t_charge);
    tree->SetBranchAddress("t_time",&t_time);

    tree->SetBranchAddress("eLOC0_fit", &eLOC0_fit);
    tree->SetBranchAddress("eLOC1_fit", &eLOC1_fit);
    tree->SetBranchAddress("ePHI_fit", &ePHI_fit);
    tree->SetBranchAddress("eTHETA_fit", &eTHETA_fit);
    tree->SetBranchAddress("eQOP_fit", &eQOP_fit);
    tree->SetBranchAddress("eT_fit", &eT_fit);
    tree->SetBranchAddress("res_eLOC0_fit", &res_eLOC0_fit);
    tree->SetBranchAddress("res_eLOC1_fit", &res_eLOC1_fit);
    tree->SetBranchAddress("res_ePHI_fit", &res_ePHI_fit);
    tree->SetBranchAddress("res_eTHETA_fit", &res_eTHETA_fit);
    tree->SetBranchAddress("res_eQOP_fit", &res_eQOP_fit);
    tree->SetBranchAddress("res_eT_fit", &res_eT_fit);
    tree->SetBranchAddress("pull_eLOC0_fit", &pull_eLOC0_fit);
    tree->SetBranchAddress("pull_eLOC1_fit", &pull_eLOC1_fit);
    tree->SetBranchAddress("pull_ePHI_fit", &pull_ePHI_fit);
    tree->SetBranchAddress("pull_eTHETA_fit", &pull_eTHETA_fit);
    tree->SetBranchAddress("pull_eQOP_fit", &pull_eQOP_fit);
    tree->SetBranchAddress("pull_eT_fit", &pull_eT_fit);
    tree->SetBranchAddress("err_eLOC0_fit", &err_eLOC0_fit);
    tree->SetBranchAddress("err_eLOC1_fit", &err_eLOC1_fit);
    tree->SetBranchAddress("err_ePHI_fit", &err_ePHI_fit);
    tree->SetBranchAddress("err_eTHETA_fit", &err_eTHETA_fit);
    tree->SetBranchAddress("err_eQOP_fit", &err_eQOP_fit);
    tree->SetBranchAddress("err_eT_fit", &err_eT_fit);

    // It's not necessary if you just need to read one file, but please do it to
    // synchronize events if multiple root files are read
    if (sortEvents) {
      entryNumbers.resize(tree->GetEntries());
      tree->Draw("event_nr", "", "goff");
      // Sort to get the entry numbers of the ordered events
      TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(),
                  false);
    }
  }


  std::vector<RecoTrackInfo> getPreSelectedTracks (size_t ievent, unsigned int nMeasurementsMin = 5, unsigned int nOutliersMax = 1000, unsigned int nHolesMax = 1000, double absCosTheta=0.94 , bool duplicateMDC = false) {
    std::vector<RecoTrackInfo> preSelectedRecoTracks;
    getEntry(ievent);
    for (size_t itrack = 0; itrack <  nStates->size(); ++itrack) {

        bool hasFittedParameters_ = hasFittedParams->at(itrack);
        auto nMeasurements_ = nMeasurements->at(itrack);
        auto nOutliers_ = nOutliers->at(itrack);
        auto nHoles_ = nHoles->at(itrack);
        auto theta = eTHETA_fit->at(itrack);
        auto qop = eQOP_fit->at(itrack);
        auto pt = std::abs(1 / qop) * std::sin(theta);
        auto t_pt = t_pT->at(itrack);
        auto t_p_ = t_p->at(itrack);
        auto t_Eta = t_eta->at(itrack);
        auto eta = std::atanh(std::cos(theta));
        auto nMajorityHits_ = nMajorityHits->at(itrack);
        auto nSharedHits_ = nSharedHits->at(itrack);
        std::vector<double> res_fit;
        std::vector<double> pull_fit;
        res_fit.push_back(res_eLOC0_fit->at(itrack));
        res_fit.push_back(res_eLOC1_fit->at(itrack));
        res_fit.push_back(res_ePHI_fit->at(itrack));
        res_fit.push_back(res_eTHETA_fit->at(itrack));
        res_fit.push_back(res_eQOP_fit->at(itrack));
        res_fit.push_back(res_eT_fit->at(itrack));
        res_fit.push_back(pt -  t_pt);
        res_fit.push_back( std::abs(1 / qop) -  t_p_);

        pull_fit.push_back(pull_eLOC0_fit->at(itrack));
        pull_fit.push_back(pull_eLOC1_fit->at(itrack));
        pull_fit.push_back(pull_ePHI_fit->at(itrack));
        pull_fit.push_back(pull_eTHETA_fit->at(itrack));
        pull_fit.push_back(pull_eQOP_fit->at(itrack));
        pull_fit.push_back(pull_eT_fit->at(itrack));


        auto majorityParticleId_ = majorityParticleId->at(itrack);
        auto majorityParticlePdg_ = majorityParticlePdg->at(itrack);
        
	auto measurementIndex_double = measurementIndex->at(itrack);
	std::vector<uint64_t> measurementIndex_;
        measurementIndex_.reserve(measurementIndex_double.size());
	std::vector<uint64_t> twinsMeasurementIndex_;
        twinsMeasurementIndex_.reserve(measurementIndex_double.size());
	
        for(auto & mid: measurementIndex_double){
          measurementIndex_.push_back(uint64_t(mid));
          if(not duplicateMDC){ 
	    twinsMeasurementIndex_.push_back(uint64_t(mid));
	  }
	}
        
	if(duplicateMDC){
          auto twinsMeasurementIndex_double = twinsMeasurementIndex->at(itrack);
          for(auto & mid: twinsMeasurementIndex_double){
	     twinsMeasurementIndex_.push_back(uint64_t(mid));
          }
        } 

        // Select the track, e.g. you might also want to add cuts on the
        // nOutliers, nHoles
        if ((!hasFittedParameters_) or nMeasurements_ < nMeasurementsMin or
            nOutliers_ > nOutliersMax  or nHoles_ > nHolesMax or std::abs(std::cos(theta))>absCosTheta) {
          continue;
        }

        preSelectedRecoTracks.push_back(
            {eta, pt, qop, nMajorityHits_, nMeasurements_, nOutliers_, nSharedHits_, majorityParticleId_, majorityParticlePdg_, t_pt, t_Eta, measurementIndex_, twinsMeasurementIndex_, res_fit, pull_fit});
    }
    return preSelectedRecoTracks; 	
  }


  // The variables
  uint32_t eventId;
  std::vector<unsigned int>* nStates = new std::vector<unsigned int>;
  std::vector<unsigned int>* nMeasurements = new std::vector<unsigned int>;
  std::vector<unsigned int>* nOutliers = new std::vector<unsigned int>;
  std::vector<unsigned int>* nHoles = new std::vector<unsigned int>;
  std::vector<unsigned int>* nSharedHits = new std::vector<unsigned int>;
  std::vector<float>* chi2Sum = new std::vector<float>;
  std::vector<unsigned int>* NDF = new std::vector<unsigned int>;
  std::vector<std::vector<double>>* measurementChi2 =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* measurementIndex =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* twinsMeasurementIndex =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* outlierChi2 =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* measurementVolume =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* measurementLayer =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* outlierVolume =
      new std::vector<std::vector<double>>;
  std::vector<std::vector<double>>* outlierLayer =
      new std::vector<std::vector<double>>;
  std::vector<unsigned int>* nMajorityHits = new std::vector<unsigned int>;
  std::vector<uint64_t>* majorityParticleId = new std::vector<uint64_t>;
  std::vector<int>* majorityParticlePdg = new std::vector<int>;

  std::vector<bool>* hasFittedParams = new std::vector<bool>;

  // True parameters
  std::vector<float>* t_d0 = new std::vector<float>;
  std::vector<float>* t_z0 = new std::vector<float>;
  std::vector<float>* t_phi = new std::vector<float>;
  std::vector<float>* t_theta = new std::vector<float>;
  std::vector<float>* t_eta = new std::vector<float>;
  std::vector<float>* t_p = new std::vector<float>;
  std::vector<float>* t_pT = new std::vector<float>;
  std::vector<float>* t_time = new std::vector<float>;
  std::vector<int>* t_charge = new std::vector<int>;

  // Estimated parameters
  std::vector<float>* eLOC0_fit = new std::vector<float>;
  std::vector<float>* eLOC1_fit = new std::vector<float>;
  std::vector<float>* ePHI_fit = new std::vector<float>;
  std::vector<float>* eTHETA_fit = new std::vector<float>;
  std::vector<float>* eQOP_fit = new std::vector<float>;
  std::vector<float>* eT_fit = new std::vector<float>;

  std::vector<float>* err_eLOC0_fit = new std::vector<float>;
  std::vector<float>* err_eLOC1_fit = new std::vector<float>;
  std::vector<float>* err_ePHI_fit = new std::vector<float>;
  std::vector<float>* err_eTHETA_fit = new std::vector<float>;
  std::vector<float>* err_eQOP_fit = new std::vector<float>;
  std::vector<float>* err_eT_fit = new std::vector<float>;

    // The residual of fitted parameters
  std::vector<float>* res_eLOC0_fit = new std::vector<float>;
  std::vector<float>* res_eLOC1_fit = new std::vector<float>;
  std::vector<float>* res_ePHI_fit = new std::vector<float>;
  std::vector<float>* res_eTHETA_fit = new std::vector<float>;
  std::vector<float>* res_eQOP_fit = new std::vector<float>;
  std::vector<float>* res_eT_fit = new std::vector<float>;

  std::vector<float>* pull_eLOC0_fit = new std::vector<float>;
  std::vector<float>* pull_eLOC1_fit = new std::vector<float>;
  std::vector<float>* pull_ePHI_fit = new std::vector<float>;
  std::vector<float>* pull_eTHETA_fit = new std::vector<float>;
  std::vector<float>* pull_eQOP_fit = new std::vector<float>;
  std::vector<float>* pull_eT_fit = new std::vector<float>;

};

/// Struct used for reading particles written out by the
/// TrackFinderPerformanceWriter
///
struct ParticleReader : public TreeReader {
  // Delete the default constructor
  ParticleReader() = delete;

  // The constructor
  ParticleReader(TTree* tree_, bool sortEvents) : TreeReader(tree_) {
    tree->SetBranchAddress("event_id", &eventId);
    tree->SetBranchAddress("particle_id", &particleId);
    tree->SetBranchAddress("particle_type", &particleType);
    tree->SetBranchAddress("vx", &vx);
    tree->SetBranchAddress("vy", &vy);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("vt", &vt);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("m", &m);
    tree->SetBranchAddress("q", &q);
    tree->SetBranchAddress("nhits", &nHits);
    tree->SetBranchAddress("ntracks", &nTracks);
    tree->SetBranchAddress("ntracks_majority", &nTracksMajority);

    // It's not necessary if you just need to read one file, but please do it to
    // synchronize events if multiple root files are read
    if (sortEvents) {
      entryNumbers.resize(tree->GetEntries());
      tree->Draw("event_id", "", "goff");
      // Sort to get the entry numbers of the ordered events
      TMath::Sort(tree->GetEntries(), tree->GetV1(), entryNumbers.data(),
                  false);
    }
  }

  // Get all the particles with requested event id
  std::vector<ParticleInfo> getParticles(const uint32_t& eventNumber) const {
    // Find the start entry and the batch size for this event
    std::string eventNumberStr = std::to_string(eventNumber);
    std::string findStartEntry = "event_id<" + eventNumberStr;
    std::string findParticlesSize = "event_id==" + eventNumberStr;
    size_t startEntry = tree->GetEntries(findStartEntry.c_str());
    size_t nParticles = tree->GetEntries(findParticlesSize.c_str());
    std::vector<ParticleInfo> particles;
    if (nParticles == 0) {
       //std::cout<<"WARNING: Event " << eventNumber  << " has no particles "<< std::endl;      
       return particles; 
      //throw std::invalid_argument(
      //    "No particles found. Please check the input file.");
    }
    particles.reserve(nParticles);
    for (unsigned int i = 0; i < nParticles; ++i) {
      getEntry(startEntry + i);
      auto pt = std::hypot(px, py);
      auto p = std::hypot(pt, pz);
      auto eta = std::atanh(pz / p * 1.);
      particles.push_back({particleId, particleType, eta, p, pt, q, nHits});
    }
    return particles;
  }

  // The variables
  ULong64_t eventId = 0;
  ULong64_t particleId = 0;
  Int_t particleType = 0;
  float vx = 0, vy = 0, vz = 0;
  float vt = 0;
  float px = 0, py = 0, pz = 0;
  float m = 0;
  float q = 0;
  UShort_t nHits = 0;
  UShort_t nTracks = 0;
  UShort_t nTracksMajority = 0;
};



std::vector<std::pair<const RecoTrackInfo*, const RecoTrackInfo*>> removeDuplicateTracks(const std::vector<RecoTrackInfo>&  preSelectedRecoTracks, size_t nSharedHitsMax, size_t numMeasurements_cut , bool forMDC = false){      
      
      std::vector<std::pair<const RecoTrackInfo*, const RecoTrackInfo*>> vecPairTrackMastertrack;
      for(const auto & track: preSelectedRecoTracks){
        vecPairTrackMastertrack.push_back(std::make_pair(&track, nullptr));
      }

      for(auto & [track_l, mastertrack_l]: vecPairTrackMastertrack){
        if(mastertrack_l){
          continue;
        }
        uint64_t numSharedHits_l = track_l->nSharedHits;
        if(!forMDC and !numSharedHits_l){
          continue;
        }
        const std::vector<uint64_t>& vecHitId_l = track_l->measurementIndex;
        const std::vector<uint64_t>& vecHitIdTwins_l = track_l->twinsMeasurementIndex;
        for(const auto & [track_r, mastertrack_r]: vecPairTrackMastertrack){
          if(track_l == track_r){
            continue;
          }
          if(mastertrack_r){
            continue;
          }
          uint64_t numSharedHits_r = track_r->nSharedHits;
          if(!forMDC and !numSharedHits_r){
            continue;
          }
          uint64_t numMeasurements_r = track_r->nMeasurements;
          if(numMeasurements_r < numMeasurements_cut){
            continue;
          }
          const std::vector<uint64_t>& vecHitId_r = track_r->measurementIndex;
          const std::vector<uint64_t>& vecHitIdTwins_r = track_r->twinsMeasurementIndex;
          uint64_t numSharedHits_pair = 0;
	  //for(auto hitId_l : vecHitId_l){
	  for(int ihitId_l=0; ihitId_l<vecHitId_l.size(); ++ihitId_l){
	    auto hitId_l = vecHitId_l[ihitId_l];
	    auto hitIdTwins_l = vecHitIdTwins_l[ihitId_l];
            if(forMDC){
	        for(int ihitId_r=0; ihitId_r<vecHitId_r.size(); ++ihitId_r){
	          auto hitId_r = vecHitId_r[ihitId_r];
	          auto hitIdTwins_r = vecHitIdTwins_r[ihitId_r];
		  if(hitId_l== hitId_r or hitIdTwins_l== hitId_r){
		  //if(hitIdTwins_l== hitId_r){
                     numSharedHits_pair  ++;
		     //std::cout<<"hitId_l = " << hitId_l <<", hitId_r = " << hitId_r <<", hitIdTwins_l = " <<  hitIdTwins_l << std::endl; 
		  }
                }
	    } else { 
	      if(std::find(vecHitId_r.begin(),vecHitId_r.end(), hitId_l)!=vecHitId_r.end()){
                numSharedHits_pair  ++;
              }
            }
          }
          
	  uint64_t numMeasurements_l = track_l->nMeasurements;
          if(numSharedHits_pair > nSharedHitsMax && numMeasurements_l <= numMeasurements_r ){
            mastertrack_l = track_r;
            //std::cout<<"nSharedHits " << numSharedHits_pair<<", removed id "<< track_l->majorityParticleId <<" master id" << mastertrack_l->majorityParticleId<<std::endl;
          }
        }
     }
     //std::cout<<"Returned vecPairTrackMastertrack.size() " << vecPairTrackMastertrack.size() << std::endl;
     return  vecPairTrackMastertrack; 
}
