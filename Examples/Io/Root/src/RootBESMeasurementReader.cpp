// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootBESMeasurementReader.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

ActsExamples::RootBESMeasurementReader::RootBESMeasurementReader(
    const ActsExamples::RootBESMeasurementReader::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_logger{Acts::getDefaultLogger(name(), level)},
      m_cfg(config),
      m_events(0),
      m_mcParticleTreeReader(nullptr),
      m_mcPixHitTreeReader(nullptr),
      m_recMdcHitTreeReader(nullptr) {
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument(
        "Missing simulated particles output collection");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits output collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  TFile* file = TFile::Open(m_cfg.filePath.c_str(), "READ");

  m_mcParticleTreeReader =
      new TTreeReader(m_cfg.mcParticleTreeName.c_str(), file);
  m_mcPixHitTreeReader = new TTreeReader(m_cfg.mcPixHitTreeName.c_str(), file);
  m_recMdcHitTreeReader =
      new TTreeReader(m_cfg.recMdcHitTreeName.c_str(), file);

  // Set the branches

  particlePDG =
      new TTreeReaderArray<int>(*m_mcParticleTreeReader, "particleId");
  particleIndex =
      new TTreeReaderArray<int>(*m_mcParticleTreeReader, "particleIndex");
  particleVertexX = new TTreeReaderArray<Acts::ActsScalar>(
      *m_mcParticleTreeReader, "xInitialPosition");
  particleVertexY = new TTreeReaderArray<Acts::ActsScalar>(
      *m_mcParticleTreeReader, "yInitialPosition");
  particleVertexZ = new TTreeReaderArray<Acts::ActsScalar>(
      *m_mcParticleTreeReader, "zInitialPosition");
  particleMomentumX =
      new TTreeReaderArray<double>(*m_mcParticleTreeReader, "xInitialMomentum");
  particleMomentumY =
      new TTreeReaderArray<double>(*m_mcParticleTreeReader, "yInitialMomentum");
  particleMomentumZ =
      new TTreeReaderArray<double>(*m_mcParticleTreeReader, "zInitialMomentum");
//  particleTrackID =
//      new TTreeReaderArray<int>(*m_mcParticleTreeReader, "particleIndex");

  PIXparticleIndex =
      new TTreeReaderArray<uint>(*m_mcPixHitTreeReader, "particleIndex");
  PIXparticlePDG =
      new TTreeReaderArray<int>(*m_mcPixHitTreeReader, "particleId");
  PIXpositionX = new TTreeReaderArray<Acts::ActsScalar>(*m_mcPixHitTreeReader,
                                                        "xPosition");
  PIXpositionY = new TTreeReaderArray<Acts::ActsScalar>(*m_mcPixHitTreeReader,
                                                        "yPosition");
  PIXpositionZ = new TTreeReaderArray<Acts::ActsScalar>(*m_mcPixHitTreeReader,
                                                        "zPosition");
  PIXmomentumX = new TTreeReaderArray<Acts::ActsScalar>(*m_mcPixHitTreeReader,
                                                        "xMomentum");
  PIXmomentumY = new TTreeReaderArray<Acts::ActsScalar>(*m_mcPixHitTreeReader,
                                                        "yMomentum");
  PIXmomentumZ = new TTreeReaderArray<Acts::ActsScalar>(*m_mcPixHitTreeReader,
                                                        "zMomentum");

  MDCcellID = new TTreeReaderArray<double>(*m_recMdcHitTreeReader, "wire");
  MDClayerID = new TTreeReaderArray<double>(*m_recMdcHitTreeReader, "layer");
  MDCparticleIndex =
      new TTreeReaderArray<int>(*m_recMdcHitTreeReader, "particleIndex");
  MDCparticlePDG =
      new TTreeReaderArray<int>(*m_recMdcHitTreeReader, "particleId");
  MDCwireR = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
                                                        "r");
  MDCwirePhi = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
                                                        "phi");
 // MDCpositionX = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
 //                                                       "xPosition");
 // MDCpositionY = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
 //                                                       "yPosition");
 // MDCpositionZ = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
 //                                                       "zPosition");
 // MDCmomentumX = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
 //                                                       "xMomentum");
 // MDCmomentumY = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
 //                                                       "yMomentum");
 // MDCmomentumZ = new TTreeReaderArray<Acts::ActsScalar>(*m_recMdcHitTreeReader,
 //                                                       "zMomentum");
  MDCdriftDistance = new TTreeReaderArray<Acts::ActsScalar>(
      *m_recMdcHitTreeReader, "driftDist");
  MDCdriftDistanceError = new TTreeReaderArray<Acts::ActsScalar>(
      *m_recMdcHitTreeReader, "driftDistError");

  ACTS_DEBUG("Adding File " << m_cfg.filePath << " to tree '"
                            << m_cfg.mcParticleTreeName << "'.");

  m_events = m_mcParticleTreeReader->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

std::pair<size_t, size_t>
ActsExamples::RootBESMeasurementReader::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::RootBESMeasurementReader::~RootBESMeasurementReader() {
  delete particlePDG;
  delete particleIndex;
  delete particleVertexX;
  delete particleVertexY;
  delete particleVertexZ;
  delete particleMomentumX;
  delete particleMomentumY;
  delete particleMomentumZ;

  delete PIXparticleIndex;
  delete PIXparticlePDG;
  delete PIXpositionX;
  delete PIXpositionY;
  delete PIXpositionZ;
  delete PIXmomentumX;
  delete PIXmomentumY;
  delete PIXmomentumZ;

  delete MDCcellID;
  delete MDClayerID;
  delete MDCparticleIndex;
  delete MDCparticlePDG;
  delete MDCwireR;
  delete MDCwirePhi;
 // delete MDCpositionX;
 // delete MDCpositionY;
 // delete MDCpositionZ;
 // delete MDCmomentumX;
 // delete MDCmomentumY;
 // delete MDCmomentumZ;
  delete MDCdriftDistance;
  delete MDCdriftDistanceError;
}

ActsExamples::ProcessCode ActsExamples::RootBESMeasurementReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read BES tracker measurements.");

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(context);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  // Prepare output containers
  // need list here for stable addresses
  SimParticleContainer particles;
  SimParticleContainer::sequence_type unordered_particles;

  // 1) The simHits have a valid geometry id, but might have meaningless
  // particle index 2) For MDC, the simHit is actually not a real sim hit, it's
  // actually the same as the measurement
  SimHitContainer simHits;
  SimHitContainer::sequence_type unordered_hits;

  std::list<IndexSourceLink> sourceLinkStorage;
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurements;
  IndexMultimap<ActsFatras::Barcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;

  // read in the fitted track parameters and particles
  if (m_recMdcHitTreeReader != nullptr && m_mcPixHitTreeReader != nullptr &&
      m_mcParticleTreeReader != nullptr && context.eventNumber < m_events) {
    // The index of the sim hit among all the sim hits from this particle
    std::map<int, int> particleHitIdx;
    std::map<int, int> particlePIXHitIdx;

    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);

    int nParticles = 0;
    int nPIXHits = 0;
    int nMDCHits = 0;

    // now read
    if (m_mcPixHitTreeReader->Next()) {
      std::cout << "Reading mcPixHit in event " << m_evtCounter++ << std::endl;

      // Reading PIX hits
      for (size_t i = 0; i < PIXpositionX->GetSize(); ++i) {
        // int parentID = (*PIXparentID)[i];
        // if(parentID!=0) continue;
        nPIXHits++;

        Acts::Vector3 pos((*PIXpositionX)[i], (*PIXpositionY)[i],
                          (*PIXpositionZ)[i]);
        Acts::Vector3 mom((*PIXmomentumX)[i]/1000, (*PIXmomentumY)[i]/1000,
                          (*PIXmomentumZ)[i]/1000);

        int layerID = 0;
        Acts::GeometryIdentifier geoId = Acts::GeometryIdentifier()
                                             .setVolume(m_volumeIDs[0])
                                             .setLayer(2 * (layerID + 1))
                                             .setSensitive(1);
        std::cout<<"Pix hit at " << geoId << std::endl;
        std::cout<<"r = " << std::hypot((*PIXpositionX)[i],
         (*PIXpositionY)[i]) <<", phi = " << std::atan2((*PIXpositionY)[i],
         (*PIXpositionX)[i]) << std::endl;
        const Acts::Surface* surfacePtr =
            m_cfg.trackingGeometry->findSurface(geoId);
        if(surfacePtr==nullptr){
          std::cout<<"surface with geoId " << geoId << " is not found "<< std::endl;
	}	
	auto intersection = surfacePtr->intersect(context.geoContext, pos,
                                                  mom.normalized(), true);
        Acts::Vector3 posUpdated = intersection.intersection.position;

        auto cylinderSurface =
            dynamic_cast<const Acts::CylinderSurface*>(surfacePtr);
        if (cylinderSurface == nullptr) {
          std::cout << "Cast PIX surface to cylinder surface failed "
                    << std::endl;
        } else {
          auto bounds = cylinderSurface->bounds();
          auto values = bounds.values();
          // std::cout<<"PIX surface r = "<< values[0] << std::endl;
        }

        int particleInd = (*PIXparticleIndex)[i];
        int particlePdg = (*PIXparticlePDG)[i];
        double mass = 0;
        double charge = 0;
        if (m_particleMassCharge.find(particlePdg) !=
                                       m_particleMassCharge.end()) {
          mass = m_particleMassCharge[particlePdg][0];
          charge = m_particleMassCharge[particlePdg][1];
       // } else {
       //   std::cout<<"WARNING: the particle has pdgId " << particlePdg << std::endl;
	}

        ActsFatras::Hit::Vector4 pos4{
            posUpdated.x() * Acts::UnitConstants::mm,
            posUpdated.y() * Acts::UnitConstants::mm,
            posUpdated.z() * Acts::UnitConstants::mm,
            0 * Acts::UnitConstants::ns,
        };
        auto energy = std::sqrt(mom.x() * mom.x() + mom.y() * mom.y() +
                                mom.z() * mom.z() + mass * mass);
        ActsFatras::Hit::Vector4 mom4{
            mom.x() / 1000 * Acts::UnitConstants::GeV,
            mom.y() / 1000 * Acts::UnitConstants::GeV,
            mom.z() / 1000 * Acts::UnitConstants::GeV,
            energy / 1000. * Acts::UnitConstants::GeV,
        };
        ActsFatras::Hit::Vector4 delta4{
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
        };

        // The particleIndex of the hit could be meaningless, i.e. the hit might
        // a noise hit
        ActsFatras::Hit hit(geoId, particleInd, pos4, mom4, mom4 + delta4,
                            particleHitIdx[particleInd]);
        unordered_hits.push_back(std::move(hit));

        particleHitIdx[particleInd]++;
        particlePIXHitIdx[particleInd]++;
      }
    }

    if (m_recMdcHitTreeReader->Next()) {
      std::cout << "Reading recMdcHit in event " << m_evtCounter << std::endl;

      // Reading MDC hits
      for (size_t i = 0; i < MDCcellID->GetSize(); ++i) {
        // int parentID = (*MDCparentID)[i];
        // if(parentID!=0) continue;
        nMDCHits++;


	//This is totally wrong because we don't have the truth hit info in MDC. So don't use the truth hit info (except its corresponding particleId)
	double x = (*MDCwireR)[i]*std::cos((*MDCwirePhi)[i]);
///////////////////////////////////////////////////////////////////////////////////////////
	double y = (*MDCwireR)[i]*std::sin((*MDCwirePhi)[i]);
        Acts::Vector3 pos(x, (*MDCdriftDistance)[i]*10, (*MDCdriftDistanceError)[i]*10);
        Acts::Vector3 mom(1, 1, 1);
///////////////////////////////////////////////////////////////////////////////////////////


        // double dang = 2*M_PI/m_MDCnCells[(*MDClayerID)[i]];
        int layerID = (*MDClayerID)[i] >=36 ? (37 + ((*MDClayerID)[i]-36)*2) : (*MDClayerID)[i];
        //int layerID = (*MDClayerID)[i];
	int cellID = m_MDCnCells[(*MDClayerID)[i]] - (*MDCcellID)[i];
	
        Acts::GeometryIdentifier geoId = Acts::GeometryIdentifier()
                                             .setVolume(m_volumeIDs[1])
                                             .setLayer((layerID + 1) * 2)
                                             .setSensitive(cellID);
        const Acts::Surface* surfacePtr =
            m_cfg.trackingGeometry->findSurface(geoId);
     
        std::cout<<"MDC hit at " << geoId << std::endl;
        
	if(surfacePtr==nullptr){
          std::cout<<"layerID = "<< layerID <<", cellID = "<< cellID << std::endl;
	}

     //   auto intersection = surfacePtr->intersect(context.geoContext, pos,
     //                                             mom.normalized(), true);
     //   Acts::Vector3 posUpdated = intersection.intersection.position;
     //   auto lpResult = surfacePtr->globalToLocal(context.geoContext,
     //                                             posUpdated, mom.normalized());
     //   if (not lpResult.ok()) {
     //     ACTS_FATAL("Global to local transformation did not succeed.");
     //     return ProcessCode::ABORT;
     //   }

        int particleInd = (*MDCparticleIndex)[i];
        int particlePdg = (*MDCparticlePDG)[i];
        double mass = 0;
        double charge = 0;
        if (m_particleMassCharge.find(particlePdg) !=
                                       m_particleMassCharge.end()) {
          mass = m_particleMassCharge[particlePdg][0];
          charge = m_particleMassCharge[particlePdg][1];
        //} else {
        //  std::cout<<"WARNING: the particle has pdgId " << particlePdg << std::endl;
	}

        // This is totally wrong
        ActsFatras::Hit::Vector4 pos4{
            pos.x() * Acts::UnitConstants::mm,
            pos.y() * Acts::UnitConstants::mm,
            pos.z() * Acts::UnitConstants::mm,
            0 * Acts::UnitConstants::ns,
        };
        auto energy = std::sqrt(mom.x() * mom.x() + mom.y() * mom.y() +
                                mom.z() * mom.z() + mass * mass);
        // This is totally wrong
        ActsFatras::Hit::Vector4 mom4{
            mom.x() / 1000 * Acts::UnitConstants::GeV,
            mom.y() / 1000 * Acts::UnitConstants::GeV,
            mom.z() / 1000 * Acts::UnitConstants::GeV,
            energy / 1000. * Acts::UnitConstants::GeV,
        };
        ActsFatras::Hit::Vector4 delta4{
            0 * Acts::UnitConstants::GeV, 0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV, 0 * Acts::UnitConstants::GeV,
        };

        // The particleIndex of the hit could be meaningless, i.e. the hit might
        // a noise hit
        ActsFatras::Hit hit(geoId, particleInd, pos4, mom4, mom4 + delta4,
                            particleHitIdx[particleInd]);
        unordered_hits.push_back(std::move(hit));
        particleHitIdx[particleInd]++;
      }
    }

    simHits.insert(unordered_hits.begin(), unordered_hits.end());

    sourceLinks.reserve(simHits.size());
    measurements.reserve(simHits.size());
    measurementParticlesMap.reserve(simHits.size());
    measurementSimHitsMap.reserve(simHits.size());

    ACTS_DEBUG("Starting loop over modules ...");
    for (auto simHitsGroup : groupByModule(simHits)) {
      // Manual pair unpacking instead of using
      //   auto [moduleGeoId, moduleSimHits] : ...
      // otherwise clang on macos complains that it is unable to capture the
      // local binding in the lambda used for visiting the smearer below.
      Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
      const auto& moduleSimHits = simHitsGroup.second;

      const Acts::Surface* surfacePtr =
          m_cfg.trackingGeometry->findSurface(moduleGeoId);

      if (surfacePtr == nullptr) {
        // this is either an invalid geometry id or a misconfigured smearer
        // setup; both cases can not be handled and should be fatal.
        ACTS_ERROR("Could not find surface " << moduleGeoId
                                             << " for configured smearer");
        return ProcessCode::ABORT;
      }

      for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
        const auto& simHit = *h;
        const auto simHitIdx = simHits.index_of(h);
        auto pos = simHit.position();
        //auto dir = simHit.unitDirection();
        auto particleInd = simHit.particleId();

        Index measurementIdx = measurements.size();
        sourceLinkStorage.emplace_back(moduleGeoId, measurementIdx);
        IndexSourceLink& sourceLink = sourceLinkStorage.back();
        sourceLinks.insert(sourceLinks.end(), sourceLink);

        if (surfacePtr->type() == Acts::Surface::SurfaceType::Cylinder) {
          std::array<Acts::BoundIndices, 2> indices = {Acts::eBoundLoc0,
                                                       Acts::eBoundLoc1};

          Acts::ActsVector<2> par{m_PIXRadius[moduleGeoId.layer() / 2 - 1] *
                                      (Acts::VectorHelpers::phi(pos) +
                                       m_pixSmear[0] * stdNormal(rng)),
                                       //0),
                                  pos.z() + m_pixSmear[1] * stdNormal(rng)};
                                  //pos.z()};
          Acts::ActsSymMatrix<2> cov = Acts::ActsSymMatrix<2>::Identity();
          cov(0, 0) = m_pixSmear[0] * m_pixSmear[0];
          cov(1, 1) = m_pixSmear[1] * m_pixSmear[1];

          measurements.emplace_back(Acts::Measurement<Acts::BoundIndices, 2>(
              sourceLink, indices, par, cov));
        } else if (surfacePtr->type() == Acts::Surface::SurfaceType::Straw) {
          std::array<Acts::BoundIndices, 1> indices = {Acts::eBoundLoc0};
          auto driftDistance = pos[1]; //We store the drift distance in the posY above
          auto driftDistanceError = pos[2]; //We store the drift distance in the posZ above

          //double sigma = 0.125;
          // std::uniform_real_distribution<double> uniform(0,1);
          // if(uniform(rng)<0.782){
          //   sigma= 0.09;
          // }  else {
          //   sigma = 0.24;
          // }

          Acts::ActsVector<1> par{driftDistance};
          Acts::ActsSymMatrix<1> cov = Acts::ActsSymMatrix<1>::Identity();
          cov(0, 0) = driftDistanceError*driftDistanceError;
          measurements.emplace_back(Acts::Measurement<Acts::BoundIndices, 1>(
              sourceLink, indices, par, cov));

        } else {
          ACTS_ERROR("The surface type must be Cylinder or Straw");
          return ProcessCode::ABORT;
        }

        // The measurement might have meaningless particle index
        measurementParticlesMap.emplace_hint(measurementParticlesMap.end(),
                                             measurementIdx, particleInd);
        measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                           measurementIdx, simHitIdx);
      }
    }

    //    for (size_t i = 0; i < measurements.size(); i++) {
    //      const auto sourceLink_ = sourceLinks.nth(i);
    //      auto geoId_ = sourceLink_->get().geometryId();
    //    }

    if (m_mcParticleTreeReader->Next()) {
      std::cout << "Reading mcParticle in event " << m_evtCounter
                << std::endl;

      // Reading truth particles
      for (size_t i = 0; i < particlePDG->GetSize(); ++i) {
        if (particlePIXHitIdx[(*particleIndex)[i]] < 1) {
          continue;
        };
        nParticles++;

        Acts::Vector3 pos((*particleVertexX)[i]*10, (*particleVertexY)[i]*10,
                          (*particleVertexZ)[i]*10);
        Acts::Vector3 mom((*particleMomentumX)[i], (*particleMomentumY)[i],
                          (*particleMomentumZ)[i]);
       double charge = m_particleMassCharge[(*particlePDG)[i]][1]; 
       std::cout<<"particle charge = " << charge << std::endl; 
       // double charge = 0;
       // if ((*particlePDG)[i] == 13 or (*particlePDG)[i] == 11 or
       //     (*particlePDG)[i] == -211 or (*particlePDG)[i] == -2212) {
       //   charge = -1;
       // } else if ((*particlePDG)[i] == -13 or (*particlePDG)[i] == -11 or
       //            (*particlePDG)[i] == 211 or (*particlePDG)[i] == 2212) {
       //   charge = 1;
       // }
        ActsFatras::Particle particle(
            ActsFatras::Barcode((*particleIndex)[i]), Acts::PdgParticle((*particlePDG)[i]),
            charge * Acts::UnitConstants::e,
            m_particleMassCharge[(*particlePDG)[i]][0] *
                Acts::UnitConstants::GeV);
        // 1 means "Undefined"
        particle.setProcess(static_cast<ActsFatras::ProcessType>(1));
        particle.setPosition4(pos.x() * Acts::UnitConstants::mm,
                              pos.y() * Acts::UnitConstants::mm,
                              pos.z() * Acts::UnitConstants::mm,
                              0 * Acts::UnitConstants::ns);
        //// Only used for direction; normalization/units do not matter
        particle.setDirection(mom.x(), mom.y(), mom.z());  // in GeV
        particle.setAbsoluteMomentum(std::sqrt(mom.x() * mom.x() +
                                               mom.y() * mom.y() +
                                               mom.z() * mom.z()) *
                                     Acts::UnitConstants::GeV);
        unordered_particles.push_back(std::move(particle));
      }
    }

    for (const auto& [particleInd, nHits] : particleHitIdx) {
      std::cout << "particle " << particleInd << " has " << nHits << " nHits"
                << ", nPIXHits = " << particlePIXHitIdx[particleInd]
                << std::endl;
    }

    std::cout << "nPIXHits = " << nPIXHits << ", nMDCHits = " << nMDCHits
              << ", nParticles (nPIXHits>0) = " << nParticles << std::endl;

    // Write the collections to the EventStore
    context.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
    context.eventStore.add(m_cfg.outputSourceLinks + "__storage",
                           std::move(sourceLinkStorage));
    context.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
    context.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                           std::move(measurementParticlesMap));
    context.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                           std::move(measurementSimHitsMap));

    particles.insert(unordered_particles.begin(), unordered_particles.end());
    context.eventStore.add(m_cfg.outputParticles, std::move(particles));

    context.eventStore.add(m_cfg.outputSimHits, std::move(simHits));
  } else {
    ACTS_WARNING("Could not read McParticleHitCol in event.");
  }

  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
