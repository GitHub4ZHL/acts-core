// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#ifdef ACTS_PLUGIN_ONNX
#include "Acts/Plugins/Onnx/MLTrackClassifier.hpp"
#endif
#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>

#include <memory>

#include <boost/filesystem.hpp>

#include "RecInput.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace boost::filesystem;
using namespace std::placeholders;

int runRecCKFTracks(int argc, char* argv[],
                    std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  detector->addOptions(desc);
  Options::addBFieldOptions(desc);
  Options::addTrackFindingOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readBField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto hitSmearingCfg = setupSimHitSmearing(
      vm, sequencer, rnd, trackingGeometry, simHitReaderCfg.outputSimHits);

  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 500_MeV;
  particleSelectorCfg.nHitsMin = 9;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Run the particle smearing
  auto particleSmearingCfg =
      setupParticleSmearing(vm, sequencer, rnd, inputParticles);

  // Setup the track finding algorithm with CKF
  // It takes all the source links created from truth hit smearing, seeds from
  // truth particle smearing and source link selection config
  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
  trackFindingCfg.inputMeasurements = hitSmearingCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters =
      particleSmearingCfg.outputTrackParameters;
  trackFindingCfg.outputTrajectories = "trajectories";
  trackFindingCfg.findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFindingAlgorithm>(trackFindingCfg, logLevel));

  // write track states from CKF
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. Thsi could be avoided when a seperate track
  // selection algorithm is used.
  trackStatesWriter.inputParticles = particleReader.outputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      hitSmearingCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.outputDir = outputDir;
  trackStatesWriter.outputFilename = "trackstates_ckf.root";
  trackStatesWriter.outputTreename = "trackstates_ckf";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));

  // write track parameters from CKF
  RootTrajectoryParametersWriter::Config trackParamsWriter;
  trackParamsWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. Thsi could be avoided when a seperate track
  // selection algorithm is used.
  trackParamsWriter.inputParticles = particleReader.outputParticles;
  trackParamsWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackParamsWriter.outputDir = outputDir;
  trackParamsWriter.outputFilename = "trackparams_ckf.root";
  trackParamsWriter.outputTreename = "trackparams_ckf";
  sequencer.addWriter(std::make_shared<RootTrajectoryParametersWriter>(
      trackParamsWriter, logLevel));

  // Write CKF performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = inputParticles;
  perfWriterCfg.inputTrajectories = trackFindingCfg.outputTrajectories;
  perfWriterCfg.inputParticles = inputParticles;
  perfWriterCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  perfWriterCfg.outputDir = outputDir;
#ifdef ACTS_PLUGIN_ONNX
  // Onnx plugin related options
  // Path to default demo ML model for track classification
  path currentFilePath(__FILE__);
  path parentPath = currentFilePath.parent_path();
  path demoModelPath =
      canonical(parentPath / "MLAmbiguityResolutionDemo.onnx").native();
  // Threshold probability for neural network to classify track as duplicate
  double decisionThreshProb = 0.5;
  // Initialize OnnxRuntime plugin
  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "MLTrackClassifier");
  Acts::MLTrackClassifier neuralNetworkClassifier(env, demoModelPath.c_str());
  perfWriterCfg.duplicatedPredictor =
      std::bind(&Acts::MLTrackClassifier::isDuplicate, &neuralNetworkClassifier,
                std::placeholders::_1, decisionThreshProb);
#endif
  sequencer.addWriter(
      std::make_shared<CKFPerformanceWriter>(perfWriterCfg, logLevel));

  return sequencer.run();
}
