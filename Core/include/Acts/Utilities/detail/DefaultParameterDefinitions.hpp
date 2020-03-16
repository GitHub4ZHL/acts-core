// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <cmath>

// Acts includes
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterTypes.hpp"

namespace Acts {

/// Components of a bound track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum BoundParametersIndices : unsigned int {
  eBoundLoc0 = 0,
  eBoundLoc1 = 1,
  eBoundPhi = 2,
  eBoundTheta = 3,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  eBoundSignedInverseP = 4,
  eBoundTime = 5,
  // Last uninitialized value contains the total number of components
  eBoundParametersSize,
  // The following alias without prefix are
  // Generic spatial coordinates on the local surface
  eLOC_0 = eBoundLoc0,
  eLOC_1 = eBoundLoc1,
  // Spatial coordinates on a disk in polar coordinates
  eLOC_R = eLOC_0,
  eLOC_PHI = eLOC_1,
  // Spatial coordinates on a disk in Cartesian coordinates
  eLOC_X = eLOC_0,
  eLOC_Y = eLOC_1,
  // Spatial coordinates on a cylinder
  eLOC_RPHI = eLOC_0,
  eLOC_Z = eLOC_1,
  // Distance-of-closest approach on a virtual perigee surface
  eLOC_D0 = eLOC_0,
  eLOC_Z0 = eLOC_1,
  ePHI = eBoundPhi,
  eTHETA = eBoundTheta,
  eQOP = eBoundSignedInverseP,
  eT = eBoundTime,
  BoundParsDim = eBoundParametersSize,
};

/// Underlying fundamental scalar type for bound track parameters.
using BoundParametersScalar = double;

/// Components of a free track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum FreeParametersIndices : unsigned int {
  // Cartesian position
  // The spatial position components must be stored as one continous block.
  eFreePos0 = 0u,
  eFreePos1 = eFreePos0 + 1u,
  eFreePos2 = eFreePos0 + 2u,
  // Time
  eFreeTime = 3u,
  // Cartesian (unit) direction
  // The direction components must be stored as one continous block.
  eFreeDir0 = 4u,
  eFreeDir1 = eFreeDir0 + 1u,
  eFreeDir2 = eFreeDir0 + 2u,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  eFreeSignedInverseP = 7u,
  // Last uninitialized value contains the total number of components
  eFreeParametersSize,
  // For backward compatibility
  FreeParsDim = eFreeParametersSize,
};

/// Underlying fundamental scalar type for free track parameters.
using FreeParametersScalar = double;

/// Components of a space point vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum SpacePointIndices : unsigned int {
  // For spatial points
  // The spatial position components must be stored as one continous block.
  eSpX = 0u,
  eSpY = eSpX + 1u,
  eSpZ = eSpX + 2u,
  eSpT = 3u,
  // Alias to allow clear code when accessing momentum vectors
  eSpPx = eSpX,
  eSpPy = eSpY,
  eSpPz = eSpZ,
  eSpE = eSpT,
  // Last uninitialized value contains the total number of components
  eSpacePointSize,
  // for backward compatibility
  SpacePointDim = eSpacePointSize,
};

/// Underlying fundamental scalar type for space points.
using SpacePointScalar = double;

using ParDef = BoundParametersIndices;
using ParID_t = BoundParametersIndices;
using ParValue_t = double;

///
/// Type namings with bound parameters
///

/// Vector of bound parameters
using BoundVector = ActsVector<ParValue_t, BoundParsDim>;
/// Row vector of bound parameters
using BoundRowVector = ActsRowVector<ParValue_t, BoundParsDim>;
/// Matrix of bound-to-bound parameters
using BoundMatrix = ActsMatrix<ParValue_t, BoundParsDim, BoundParsDim>;
/// Symmetrical matrix of bound-to-bound parameters
using BoundSymMatrix = ActsSymMatrix<ParValue_t, BoundParsDim>;

///
/// Type naming with free parameters
///

/// Vector of free track parameters
using FreeVector = ActsVector<ParValue_t, FreeParsDim>;
/// Matrix of free-to-free parameters
using FreeMatrix = ActsMatrix<ParValue_t, FreeParsDim, FreeParsDim>;
/// Symmetrical matrix of free-to-free parameters
using FreeSymMatrix = ActsSymMatrix<ParValue_t, FreeParsDim>;

///
/// Type namings with bound & free parameters
///

/// Matrix of bound-to-free parameters
using BoundToFreeMatrix = ActsMatrix<ParValue_t, FreeParsDim, BoundParsDim>;
/// Matrix of free-to-bound parameters
using FreeToBoundMatrix = ActsMatrix<ParValue_t, BoundParsDim, FreeParsDim>;

///
/// Type namings with space points
///

/// Vector with space point parameters
using SpacePointVector = ActsVector<ParValue_t, SpacePointDim>;
/// Symmetrical matrix of space point-to-space point
using SpacePointSymMatrix = ActsSymMatrix<ParValue_t, SpacePointDim>;

///
/// Type namings with space points & bound parameters
///

/// Matrix of space point-to-bound parameters
using SpacePointToBoundMatrix =
    ActsMatrix<ParValue_t, BoundParsDim, SpacePointDim>;
/// Matrix with bound parameters-to-space point
using BoundToSpacePointMatrix =
    ActsMatrix<ParValue_t, SpacePointDim, BoundParsDim>;

template <ParID_t>
struct par_type;

template <ParID_t par>
using par_type_t = typename par_type<par>::type;

template <>
struct par_type<ParDef::eLOC_0> {
  using type = local_parameter;
};

template <>
struct par_type<ParDef::eLOC_1> {
  using type = local_parameter;
};

template <>
struct par_type<ParDef::ePHI> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eTHETA> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eQOP> {
  using type = unbound_parameter;
};

template <>
struct par_type<ParDef::eT> {
  using type = unbound_parameter;
};
}  // namespace Acts