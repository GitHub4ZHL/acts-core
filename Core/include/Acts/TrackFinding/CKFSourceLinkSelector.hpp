// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <limits>
#include <map>

namespace Acts {

/// Selection cuts for associating calibrated measurements with predicted track
/// parameters on a surface.
///
/// The default configuration only takes the best matching calibrated
/// measurement without a cut on the local chi2.
struct SourceLinkSelectorCuts {
  /// Maximum local chi2 contribution.
  double chi2CutOff = std::numeric_limits<double>::max();
  /// Maximum number of associated measurements on a single surface.
  size_t numSourcelinksCutOff = 1;
};

/// @brief Calibrated measurement selection struct selecting those calibrated
/// measurements compatible with the given track parameter against provided
/// criteria on one surface
///
/// The selection criteria could be allowed maximum chi2
/// and allowed maximum number of measurements on one surface
///
/// If there is no compatible measurement, the measurement with the mininum
/// chi2 will be selected and the status will be tagged as an outlier
///
struct CKFSourceLinkSelector {
 public:
  /// Geometry-dependent cut configuration.
  ///
  /// Different components on the geometry can require different cut settings.
  /// The configuration must either contain explicit settings for all geometry
  /// components that are used or contain a global default.
  using Config = Acts::GeometryHierarchyMap<SourceLinkSelectorCuts>;

  /// @brief Default constructor
  CKFSourceLinkSelector() = default;
  /// @brief Constructor with config and (non-owning) logger
  ///
  /// @param config a config instance
  /// @param logger a logger instance
  CKFSourceLinkSelector(Config cfg) : m_config(std::move(cfg)) {}

  /// @brief Operater that select the calibrated measurements compatible with
  /// the given track parameter on a surface
  ///
  /// @tparam fittable_measurement_t The type of calibrated measurement
  ///
  /// @param predictedParams The predicted track parameter on a surface
  /// @param measurements The pool of calibrated measurements
  /// @param measChi2 The container for index and chi2 of intermediate
  /// measurement candidates
  /// @param measCandidateIndices The container for index of final measurement
  /// candidates
  /// @param isOutlier The indicator for outlier or not
  ///
  template <typename fittable_measurement_t>
  Result<void> operator()(
      const BoundTrackParameters& predictedParams,
      const std::vector<fittable_measurement_t>& measurements,
      std::vector<std::pair<size_t, double>>& measChi2,
      std::vector<size_t>& measCandidateIndices, bool& isOutlier,
      LoggerWrapper logger) const {
    ACTS_VERBOSE("Invoked CKFSourceLinkSelector");

    // Return error if no measurement
    if (measurements.empty()) {
      return CombinatorialKalmanFilterError::SourcelinkSelectionFailed;
    }

    // Get geoID of this surface
    auto surface = &predictedParams.referenceSurface();
    auto geoID = surface->geometryId();

    // Find the appropriate cuts
    auto cuts = m_config.find(geoID);
    if (cuts == m_config.end()) {
      // for now we consider missing cuts an unrecoverable error
      // TODO consider other options e.g. do not add measurements at all (not
      // even as outliers)
      return CombinatorialKalmanFilterError::SourcelinkSelectionFailed;
    }
    const auto chi2CutOff = cuts->chi2CutOff;
    const auto numSourcelinksCutOff = cuts->numSourcelinksCutOff;
    ACTS_VERBOSE("Allowed maximum chi2: " << chi2CutOff);
    ACTS_VERBOSE(
        "Allowed maximum number of measurements: " << numSourcelinksCutOff);

    measChi2.resize(measurements.size());
    double minChi2 = std::numeric_limits<double>::max();
    size_t minIndex = 0;
    size_t index = 0;
    size_t nInitialCandidates = 0;
    // Loop over all measurements to select the compatible measurements
    for (const auto& measurement : measurements) {
      std::visit(
          [&](const auto& calibrated) {
            // The measurement surface should be the same as parameter surface
            assert(&calibrated.referenceObject() == surface);

            // Take the projector (measurement mapping function)
            const auto& H = calibrated.projector();
            // Take the parameter covariance
            const auto& predictedCovariance = *predictedParams.covariance();
            // Get the residual
            const auto& residual =
                calibrated.residual(predictedParams.parameters());
            // Get the chi2
            double chi2 = (residual.transpose() *
                           ((calibrated.covariance() +
                             H * predictedCovariance * H.transpose()))
                               .inverse() *
                           residual)
                              .eval()(0, 0);

            ACTS_VERBOSE("Chi2: " << chi2);
            // Push the measurement index and chi2 if satisfying the criteria
            if (chi2 < chi2CutOff) {
              measChi2.at(nInitialCandidates) = {index, chi2};
              nInitialCandidates++;
            }
            // Search for the measurement with the min chi2
            if (chi2 < minChi2) {
              minChi2 = chi2;
              minIndex = index;
            }
          },
          measurement);
      index++;
    }

    // Get the number of measurement candidates with provided constraint
    // considered
    size_t nFinalCandidates =
        std::min(nInitialCandidates, numSourcelinksCutOff);

    // If there is no selected measurement, return the measurement with the best
    // chi2 and tag it as an outlier
    if (nFinalCandidates == 0) {
      measCandidateIndices.resize(1);
      measCandidateIndices.at(0) = minIndex;
      ACTS_DEBUG("No measurement candidate. Return an outlier measurement.");
      isOutlier = true;
      return Result<void>::success();
    }

    ACTS_VERBOSE("Number of measurement candidates: " << nFinalCandidates);
    measCandidateIndices.resize(nFinalCandidates);
    // Sort the initial measurement candidates based on chi2 in ascending order
    std::sort(measChi2.begin(), measChi2.begin() + nInitialCandidates,
              [](const std::pair<size_t, double>& lchi2,
                 const std::pair<size_t, double>& rchi2) {
                return lchi2.second < rchi2.second;
              });
    // Get only allowed number of measurement candidates, i.e. nFinalCandidates,
    // from the front and reset the values in the container
    size_t nRecorded = 0;
    for (const auto& [id, chi2] : measChi2) {
      if (nRecorded >= nFinalCandidates) {
        break;
      }
      measCandidateIndices.at(nRecorded) = id;
      nRecorded++;
    }
    isOutlier = false;
    return Result<void>::success();
  }

  /// The config
  Config m_config;
};

}  // namespace Acts
