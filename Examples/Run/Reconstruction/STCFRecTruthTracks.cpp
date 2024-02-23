// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Detector/TGeoDetectorWithOptions.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootSTCFMeasurementReader.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/SpacePointMakerOptions.hpp"
#include "ActsExamples/Options/TrackFittingOptions.hpp"
#include "ActsExamples/Options/TruthSeedSelectorOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // auto detector = std::make_shared<ActsExamples::TGeoDetector>();
  const std::shared_ptr<ActsExamples::IBaseDetector>& detector =
      std::make_shared<ActsExamples::TGeoDetectorWithOptions>();

  using boost::program_options::value;

  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addFittingOptions(desc);
  Options::addDigitizationOptions(desc);
  Options::addParticleSmearingOptions(desc);
  Options::addSpacePointMakerOptions(desc);
  Options::addTruthSeedSelectorOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputfiles = vm["input-files"].template as<std::vector<std::string>>();
  if (inputfiles.size() > 1) {
    throw std::invalid_argument("Sorry, only one input file is possible here");
  } else if (inputfiles.size() < 1) {
    throw std::invalid_argument("Sorry, please provide one input file");
  }
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<const ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vm);

  // Read truth hits from CSV files
  RootSTCFMeasurementReader::Config STCFMeasurementReaderCfg;
  STCFMeasurementReaderCfg.filePath = inputDir + "/" + inputfiles.front();
  STCFMeasurementReaderCfg.outputSimHits = "hits";
  STCFMeasurementReaderCfg.outputParticles = "particles";
  STCFMeasurementReaderCfg.trackingGeometry = trackingGeometry;
  STCFMeasurementReaderCfg.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<RootSTCFMeasurementReader>(
      STCFMeasurementReaderCfg, logLevel));

  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = STCFMeasurementReaderCfg.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 30._MeV;
  particleSelectorCfg.nHitsMin = 5;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Run the particle smearing
  auto particleSmearingCfg =
      setupParticleSmearing(vm, sequencer, rnd, inputParticles);

  // The fitter needs the measurements (proto tracks) and initial
  // track states (proto states). The elements in both collections
  // must match and must be created from the same input particles.
  // Create truth tracks
  TruthTrackFinder::Config trackFinderCfg;
  trackFinderCfg.inputParticles = inputParticles;
  trackFinderCfg.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  trackFinderCfg.inputMeasurementSimHitsMap =
      STCFMeasurementReaderCfg.outputMeasurementSimHitsMap;
  trackFinderCfg.inputSimulatedHits = STCFMeasurementReaderCfg.outputSimHits;
  trackFinderCfg.removeHitsFromLoops = true;
  trackFinderCfg.outputProtoTracks = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));

  // setup the fitter
  const double reverseFilteringMomThreshold = 0.0;
  TrackFittingAlgorithm::Config fitter;
  fitter.inputMeasurements = STCFMeasurementReaderCfg.outputMeasurements;
  fitter.inputSourceLinks = STCFMeasurementReaderCfg.outputSourceLinks;
  fitter.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  fitter.inputInitialTrackParameters =
      particleSmearingCfg.outputTrackParameters;
  fitter.outputTracks = "tracks";
  fitter.pickTrack = vm["fit-pick-track"].as<int>();
  fitter.fit = makeKalmanFitterFunction(
      trackingGeometry, magneticField,
      vm["fit-multiple-scattering-correction"].as<bool>(),
      vm["fit-energy-loss-correction"].as<bool>(), reverseFilteringMomThreshold,
      Acts::FreeToBoundCorrection(
          vm["fit-ftob-nonlinear-correction"].as<bool>()));
  fitter.calibrator = std::make_shared<PassThroughCalibrator>();
  sequencer.addAlgorithm(
      std::make_shared<TrackFittingAlgorithm>(fitter, logLevel));

  // write track states from fitting
  RootTrackStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTracks = fitter.outputTracks;
  trackStatesWriter.inputParticles = inputParticles;
  trackStatesWriter.inputSimHits = STCFMeasurementReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      STCFMeasurementReaderCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.filePath = outputDir + "/trackstates_fitter.root";
  sequencer.addWriter(
      std::make_shared<RootTrackStatesWriter>(trackStatesWriter, logLevel));

  // write track summary from CKF
  RootTrackSummaryWriter::Config trackSummaryWriter;
  trackSummaryWriter.inputTracks = fitter.outputTracks;
  trackSummaryWriter.inputParticles = inputParticles;
  trackSummaryWriter.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  trackSummaryWriter.filePath = outputDir + "/tracksummary_fitter.root";
  sequencer.addWriter(
      std::make_shared<RootTrackSummaryWriter>(trackSummaryWriter, logLevel));

  // write reconstruction performance data
  TrackFinderPerformanceWriter::Config perfFinder;
  perfFinder.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  perfFinder.inputParticles = inputParticles;
  perfFinder.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  perfFinder.filePath = outputDir + "/performance_track_finder.root";
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(perfFinder, logLevel));

  TrackFitterPerformanceWriter::Config perfFitter;
  perfFitter.inputTracks = fitter.outputTracks;
  perfFitter.inputParticles = inputParticles;
  perfFitter.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  perfFitter.filePath = outputDir + "/performance_track_fitter.root";
  sequencer.addWriter(
      std::make_shared<TrackFitterPerformanceWriter>(perfFitter, logLevel));

  return sequencer.run();
}
