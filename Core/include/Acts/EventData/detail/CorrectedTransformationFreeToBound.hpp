// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
//#include <unsupported/Eigen/MatrixFunctions>

#include <type_traits>

namespace Acts {
namespace detail {

/// @brief Parameters sampling based on covariance matrix sqrt root
///
struct CorrectedFreeToBoundTransformer {
  /// The parameter to tune the weight
  //ActsScalar kappa = 4;
  ActsScalar alpha = 0.1; // 0 < alpha <=1
  ActsScalar beta=2;


  /// Get the non-linearity corrected bound parameters and its covariance
  std::optional<std::tuple<BoundVector, BoundSymMatrix, BoundToFreeMatrix>>
  operator()(const FreeVector& freeParams, const FreeSymMatrix& freeCovariance,
             const FreeVector& freeDerivatives, const Surface& surface,
             const GeometryContext& geoContext,
             NavigationDirection navDir = forward) {
    // std::cout << "Invoking CorrectedFreeToBoundTransformer operator()"
    //          << std::endl;
    // Get the incidence angle
    Vector3 dir = freeParams.segment<3>(eFreeDir0);
    Vector3 normal = surface.normal(geoContext);
    ActsScalar incidenceAng = dir.dot(normal);
    //std::cout<<"Incidence angle = \n " <<  incidenceAng << std::endl;
    if(abs(incidenceAng)<0.1){
      return std::nullopt;
    }

    size_t sampleSize = 2 * (eFreeSize) + 1;
    std::vector<std::tuple<FreeVector, ActsScalar, ActsScalar>> sampledFreeParams;
    sampledFreeParams.reserve(sampleSize);

    // Initialize the covariance sqrt root matrix
    FreeSymMatrix covSqrt = FreeSymMatrix::Zero();
    // ActsDynamicMatrix covSqrt = ActsDynamicMatrix::Zero(6,6);
    //   std::cout << "freeCovariance =  \n" << freeCovariance << std::endl;
    // Get the covariance sqrt root matrix
    Eigen::JacobiSVD<FreeSymMatrix> svd(
        freeCovariance, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // U*S*V^-1
    auto S = svd.singularValues();
    auto U = svd.matrixU();
    // std::cout << "Its singular values are:" << std::endl << S << std::endl;
    /*
    auto V = svd.matrixV();
        if (U != V) {
          std::cout << "U != V " << std::endl;
        }
        std::cout << "S = \n " << S << std::endl;
        std::cout << "U = \n " << U << std::endl;
        std::cout << "V = \n " << V << std::endl;
    */
    FreeMatrix D = FreeMatrix::Zero();
    // ActsDynamicMatrix D =  ActsDynamicMatrix::Zero(6,6);
    for (unsigned i = 0; i < eFreeSize; ++i) {
      if (S(i) == 0) {
        // std::cout<<"Zero singular values found"<< std::endl;
        continue;
      }
      D(i, i) = std::sqrt(S(i));
    }

    FreeMatrix UP = FreeMatrix::Zero();
    // ActsDynamicMatrix UP =  ActsDynamicMatrix::Zero(6,6);
    for (unsigned i = 0; i < eFreeSize; ++i) {
      for (unsigned j = 0; j < eFreeSize; ++j) {
        UP(i, j) = U(i, j);
      }
    }
    // covSqrt = UP * D * UP.transpose();
    covSqrt = UP * D;

    // auto test = covSqrt * covSqrt;
    //  std::cout << "covSqrt = \n" << covSqrt << std::endl;
    double kappa = alpha*alpha*eFreeSize;
    double skappa = std::sqrt(kappa);

    // Sample the free parameters
    // 1. the baseline parameter
    sampledFreeParams.push_back({freeParams, (kappa - eFreeSize) / kappa, (kappa-eFreeSize)/kappa+(1.0-alpha*alpha+beta)});
    // 2. the shifted parameters
    for (unsigned i = 0; i < eFreeSize; ++i) {
      ActsScalar kappaSqrt = std::sqrt(kappa);
      //     std::cout << "delta =  \n " << covSqrt.col(i) * kappaSqrt << std::endl;
      sampledFreeParams.push_back(
          {freeParams + covSqrt.col(i) * kappaSqrt, 0.5 / kappa, 0.5/kappa});
      sampledFreeParams.push_back(
          {freeParams - covSqrt.col(i) * kappaSqrt, 0.5 / kappa, 0.5/kappa});
    }

    //    std::cout<<"sampledFreeParams.size() " << sampledFreeParams.size() <<
    //    std::endl;
    // Initialize the mean of the bound parameters
    BoundVector bpMean = BoundVector::Zero();
    // Initialize the bound covariance
    BoundSymMatrix bv = BoundSymMatrix::Zero();
    // Initialize the correlation matrix
    BoundToFreeMatrix cross = BoundToFreeMatrix::Zero();

    // The transformed bound parameters
    std::vector<std::pair<BoundVector, ActsScalar>> transformedBoundParams;
    // Loop over the sample points of the free parameters to get the weighted
    // bound parameters
    for (const auto& [params, mweight, cweight] : sampledFreeParams) {
      // Calculate the derivative of path length at the final surface or the
      // point-of-closest approach w.r.t. free parameters
      FreeToPathMatrix freeToPath =
          surface.freeToPathDerivative(geoContext, freeParams);

      auto correctedFreeParams = params;
      // Reintersect to get the corrected free params
      auto intersection =
          surface.intersect(geoContext, params.segment<3>(eFreePos0),
                            navDir * params.segment<3>(eFreeDir0), true);
      // The new position and time?
      correctedFreeParams.segment<3>(eFreePos0) =
          intersection.intersection.position;
      // std::cout << "intersection.pathLength = "
      //          << intersection.intersection.pathLength << std::endl;

      // if(abs( intersection.intersection.pathLength) < 1e-20) {
      //  continue;
      //}

      // should we use pathlength * freeDerivatives?

      // auto correctedFreeParams =
      //    params + (FreeSymMatrix::Identity() + freeDerivatives * freeToPath)
      //    *
      //                 (params - freeParams);

      // Transform the free to bound
      auto result = detail::transformFreeToBoundParameters(correctedFreeParams,
                                                           surface, geoContext);
      if (not result.ok()) {
        std::cout << "transformation failed!" << std::endl;
        continue;
      }

      auto bp = result.value();
      //      std::cout << "Transformation succeeded with transformed bp = \n" << bp << std::endl;
      transformedBoundParams.push_back({bp, cweight});

      bpMean = bpMean + mweight * bp;
    }

    //    std::cout << "transformedBoundParams.size() "
    //              << transformedBoundParams.size() << std::endl;
    if (transformedBoundParams.empty()) {
      return std::nullopt;
    }

    // Get the weighted bound covariance
    for (unsigned int isample = 0; isample < sampleSize; ++isample) {
      FreeVector fSigma = std::get<0>(sampledFreeParams[isample]) - freeParams;
      BoundVector bSigma = transformedBoundParams[isample].first - bpMean;

      // std::cout << "weight " << sampledFreeParams[isample].second
      //         << " for sigma = \n"
      //        << bSigma << std::endl;
      bv = bv + transformedBoundParams[isample].second * bSigma * bSigma.transpose();
      cross = cross +
              transformedBoundParams[isample].second * fSigma * bSigma.transpose();
    }

    return std::make_tuple(bpMean, bv, cross);
  }
};

}  // namespace detail
}  // namespace Acts
