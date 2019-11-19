// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
/// Enum which determines the stage of outlier search
enum class OutlierSearchStage { Filtering = 0, Smoothing = 1 };

using OutlierFinder =
    std::function<bool(double, const Surface*, OutlierSearchStage)>;

/// @brief void Measurement calibrator and converter
struct VoidKalmanComponents {
  /// @brief Public call mimicking a calibrator
  ///
  /// @tparam measurement_t Type of the measurement
  /// @tparam parameter_t Type of the parameters for calibration
  ///
  /// @param m Measurement to be moved through
  /// @param pars Parameters to be used for calibration
  ///
  /// @return void-calibrated measurement
  template <typename measurement_t, typename parameters_t>
  Result<measurement_t> operator()(measurement_t m,
                                   const parameters_t& /*pars*/) const {
    return m;
  }

  /// @brief void measurement converter only moves the
  /// the measurement through for further processing
  ///
  /// @tparam measurement_container_t Type of the measurement
  ///
  /// @param ms Measurements to be moved through
  ///
  /// @return moved measurements
  template <typename measurements_t>
  Result<measurements_t> operator()(measurements_t ms) const {
    return std::move(ms);
  }
};

/// @brief Void measurement calibrator for filtering
struct VoidMeasurementCalibrator {
  /// Main calibration call. In this implementation, it will dereference the
  /// given source link and expect it to result in something convertible to
  /// @c FittableMeasurement<source_link_t>.
  /// @tparam source_link_t Source link type which identifier the uncalibrated
  /// measurement
  /// @tparam parameters_t Parameters type (unused)
  /// @param sl Source link to turn into a measurement
  /// @param pars The parameters to calibrate with (unused)
  /// @note If the deref operator on @c source_link_t returns a reference, this
  /// will copy it before returning. If it is already returned by-value (for
  /// instance for a newly created measurement instance), return value
  /// optimizitaion should auto-move the result.
  /// @note This will not make the "calibrated" measurement point to the
  /// uncalibrated measurement via sourcelink, it's just a copy.
  template <typename source_link_t, typename parameters_t>
  FittableMeasurement<source_link_t> operator()(
      const source_link_t& sl, const parameters_t& /*pars*/) const {
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does fulfill SourceLinkConcept.");
    static_assert(
        concept ::converts_to<FittableMeasurement<source_link_t>,
                              concept ::detail_slc::dereferenceable_t,
                              source_link_t>,
        "For DefaultMeasurementCalibrator, source link needs to implement "
        "dereference operator");

    return *sl;
  }
};

/// @brief void Kalman updater
struct VoidKalmanUpdater {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam measurement_t Type of the measurement to be used
  /// @tpredicted_state_t Type of the (bound) predicted state
  ///
  /// @param m The measurement
  /// @param predicted The predicted parameters
  ///
  /// @return The copied predicted parameters
  template <typename track_state_t, typename predicted_state_t>
  auto operator()(track_state_t& /*m*/,
                  const predicted_state_t& predicted) const {
    return &(predicted.parameters);
  }
};

/// @brief void Kalman smoother
struct VoidKalmanSmoother {
  /// @brief Public call mimicking an updater
  ///
  /// @tparam track_states_t Type of the track states
  ///
  /// @param states The track states to be smoothed
  ///
  /// @return The resulting
  template <typename parameters_t, typename track_states_t>
  const parameters_t* operator()(track_states_t& /*states*/) const {
    return nullptr;
  }
};

/// @brief void outlier rejector
struct VoidOutlierFinder {
  /// @brief Public call mimicking an outlier rejector
  ///
  /// @tparam chi2 The chisq from fitting
  ///
  /// @param surface The surface of the measurement (in case the outlier
  /// criteria is detector-specific)
  /// @param searchStage The outlier search stage
  ///
  /// @return The resulting
  bool operator()(double /*chi2*/, const Surface* /*surface*/,
                  OutlierSearchStage /*searchStage*/) const {
    return false;
  }
};

}  // namespace Acts
