// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <mutex>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

class TTreeReader;
// class TTreeReaderArray;

namespace ActsExamples {

/// @class RootBESMeasurementReader
///
/// @brief Reads in TrackParameter information from a root file
/// and fills it into a Acts::BoundTrackParameter format
class RootBESMeasurementReader : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Output simulated (truth) hits collection.
    std::string outputSimHits;
    /// Which particle collection to read into.
    std::string outputParticles;
    /// Output source links collection.
    std::string outputSourceLinks = "sourcelinks_";
    /// Output measurements collection.
    std::string outputMeasurements = "measurements_";
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap = "measurement_particles_map";
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap = "measurement_simhits_map";
    /// Tracking geometry required to get approximated sim hits on surface.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// Random numbers service.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;

    bool readRecMdcHit = false;
    bool runBesUpgrade = false;

    std::string mcParticleTreeName =
        "McParticleCol";  ///< name of the input tree for truth particles
    std::string mcPixHitTreeName =
        "PixMcHitCol";  ///< name of the input tree for truth pixel hits
    std::string mcMdcHitTreeName =
        "MdcMcHitCol";  ///< name of the input tree for truth mdc hits
    std::string recMdcHitTreeName =
        "MdcRecHitCol";    ///< name of the input tree for truth pixel hits
    std::string filePath;  ///< The name of the input file

    /// Whether the events are ordered or not
    bool orderedEvents = true;
  };

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootBESMeasurementReader(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootBESMeasurementReader();

  /// Framework name() method
  std::string name() const final override { return "RootBESMeasurementReader"; }

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(
      const ActsExamples::AlgorithmContext& context) final override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  size_t m_events = 0;

  std::map<int, std::array<double, 2>> m_particleMassCharge = {
      {11, {0.000511, -1}},    {-11, {0.000511, 1}},   {13, {0.105658, -1}},
      {-13, {0.105658, 1}},    {211, {0.139570, 1}},   {-211, {0.139570, -1}},
      {321, {0.493677, 1}},    {-321, {0.493677, -1}}, {2212, {0.938272, 1}},
      {-2212, {0.938272, -1}},
  };

  std::array<int, 43> m_MDCnCells = {
      40,  44,  48,  56,  64,  72,  80,  80,  76,  76,  88,  88,  100, 100, 112,
      112, 128, 128, 140, 140, 160, 160, 160, 160, 176, 176, 176, 176, 208, 208,
      208, 208, 240, 240, 240, 240, 256, 256, 256, 256, 288, 288, 288};

  std::array<int, 2> m_volumeIDs = {6, 8};
  std::array<Acts::ActsScalar, 1> m_PIXRadius = {35};
  std::array<double, 2> m_pixSmear = {0.00866,
                                      0.05774};

  size_t m_evtCounter = 0;

  /// The input tree name
  TTreeReader* m_mcParticleTreeReader = nullptr;
  TTreeReader* m_mcPixHitTreeReader = nullptr;
  TTreeReader* m_mcMdcHitTreeReader = nullptr;
  TTreeReader* m_recMdcHitTreeReader = nullptr;

  TTreeReaderArray<int>* particlePDG = nullptr;
  TTreeReaderArray<int>* particleIndex = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleVertexX = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleVertexY = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleVertexZ = nullptr;
  TTreeReaderArray<double>* particleMomentumX = nullptr;
  TTreeReaderArray<double>* particleMomentumY = nullptr;
  TTreeReaderArray<double>* particleMomentumZ = nullptr;

  // Pix MC
  TTreeReaderArray<uint>* PIXparticleIndex = nullptr;
  TTreeReaderArray<int>* PIXparticlePDG = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* PIXpositionX = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* PIXpositionY = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* PIXpositionZ = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* PIXmomentumX = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* PIXmomentumY = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* PIXmomentumZ = nullptr;

  // Mdc MC
  TTreeReaderArray<uint>* mcMDCcellID = nullptr;
  TTreeReaderArray<uint>* mcMDClayerID = nullptr;
  TTreeReaderArray<uint>* mcMDCparticleIndex = nullptr;
  TTreeReaderArray<int>* mcMDCparticlePDG = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCpositionX = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCpositionY = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCpositionZ = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCmomentumX = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCmomentumY = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCmomentumZ = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* mcMDCdriftDistance = nullptr;
  TTreeReaderArray<int>* mcMDCdriftSign = nullptr;
  TTreeReaderArray<int>* mcMDCinCellStatus = nullptr;

  // Mdc REC
  TTreeReaderArray<double>* recMDCcellID = nullptr;
  TTreeReaderArray<double>* recMDClayerID = nullptr;
  TTreeReaderArray<int>* recMDCparticleIndex = nullptr;
  TTreeReaderArray<int>* recMDCparticlePDG = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* recMDCwireR = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* recMDCwirePhi = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* recMDCdriftDistance = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* recMDCdriftDistanceError = nullptr;
};

}  // namespace ActsExamples
