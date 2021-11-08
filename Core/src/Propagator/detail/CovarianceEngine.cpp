// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/EventData/detail/CorrectedTransformationBoundToFree.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {
namespace {
/// Some type defs
using Jacobian = BoundMatrix;

using BoundState = std::tuple<BoundTrackParameters, Jacobian, double,
                              BoundVector, BoundSymMatrix, Jacobian>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to bound parameters at the final surface
///
/// @note Modifications of the jacobian related to the
/// projection onto a surface is considered. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] boundToFreeJacobian The projection jacobian from start local
/// to start free parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] freeDerivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] fullTransportJacobian The full jacobian from start local to
/// bound parameters at the final surface
/// @param [in, out] startBoundToFinalFreeJacobian The jacobian from start local to
/// free parameters at the final surface
/// @param [in] surface The final surface onto which the projection should be
/// performed
void boundToBoundJacobian(const GeometryContext& geoContext,
                          const FreeVector& freeParameters,
                          const BoundToFreeMatrix& boundToFreeJacobian,
                          const FreeMatrix& freeTransportJacobian,
                          const FreeVector& freeDerivatives,
                          BoundMatrix& fullTransportJacobian,
                          BoundToFreeMatrix& startBoundToFinalFreeJacobian,
                          FreeMatrix& finalFreeCorrection,
                          const Surface& surface) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  startBoundToFinalFreeJacobian = freeTransportJacobian * boundToFreeJacobian;
  finalFreeCorrection = FreeMatrix::Identity() + freeDerivatives * freeToPath;
  fullTransportJacobian =
      freeToBoundJacobian * finalFreeCorrection * startBoundToFinalFreeJacobian;
  //  std::cout<<"startBoundToFinalFreeJacobian =\n" << startBoundToFinalFreeJacobian << std::endl;
  //  std::cout<<"finalFreeCorrection =\n" << FreeMatrix::Identity() + freeDerivatives * freeToPath << std::endl;
  //  std::cout<<"freeToBoundJacobian =\n" << freeToBoundJacobian << std::endl;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to final curvilinear parameters
///
/// @note Modifications of the jacobian related to the
/// projection onto a curvilinear surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This is an approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] boundToFreeJacobian The projection jacobian from local start
/// to global final parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] freeDerivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] jacFull The full jacobian from start local to curvilinear
/// parameters
///
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
void boundToCurvilinearJacobian(const Vector3& direction,
                                const BoundToFreeMatrix& boundToFreeJacobian,
                                const FreeMatrix& freeTransportJacobian,
                                const FreeVector& freeDerivatives,
                                BoundMatrix& fullTransportJacobian) {
  // Calculate the derivative of path length at the the curvilinear surface
  // w.r.t. free parameters
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * direction;
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian = freeToCurvilinearJacobian(direction);
  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in] geoContext The geometry context
/// @param [in, out] freeTransportJacobian The transport jacobian from start
/// free to final free parameters
/// @param [in, out] freeToPathDerivatives Path length derivatives of the free,
/// nominal parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] surface The reference surface of the local parametrisation
Result<void> reinitializeJacobians(const GeometryContext& geoContext,
                                   FreeMatrix& freeTransportJacobian,
                                   FreeVector& freeToPathDerivatives,
                                   BoundToFreeMatrix& boundToFreeJacobian,
                                   const FreeVector& freeParameters,
                                   const Surface& surface) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  //  freeToPathDerivatives = FreeVector::Zero();

  // Get the local position
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (not lpResult.ok()) {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("CovarianceEngine", Logging::INFO));
    ACTS_FATAL(
        "Inconsistency in global to local transformation during propagation.")
  }
  // Transform from free to bound parameters
  Result<BoundVector> boundParameters = detail::transformFreeToBoundParameters(
      freeParameters, surface, geoContext);
  if (!boundParameters.ok()) {
    return boundParameters.error();
  }
  // Reset the jacobian from local to global
  boundToFreeJacobian =
      surface.boundToFreeJacobian(geoContext, *boundParameters);
  return Result<void>::success();
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in, out] freeTransportJacobian The transport jacobian from start
/// free to final free parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
void reinitializeJacobians(FreeMatrix& freeTransportJacobian,
                           FreeVector& freeToPathDerivatives,
                           BoundToFreeMatrix& boundToFreeJacobian,
                           const Vector3& direction) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();
  boundToFreeJacobian = BoundToFreeMatrix::Zero();

  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;

  boundToFreeJacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  boundToFreeJacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  boundToFreeJacobian(eFreePos1, eBoundLoc0) = cosPhi;
  boundToFreeJacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  boundToFreeJacobian(eFreePos2, eBoundLoc1) = sinTheta;
  boundToFreeJacobian(eFreeTime, eBoundTime) = 1;
  boundToFreeJacobian(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  boundToFreeJacobian(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  boundToFreeJacobian(eFreeDir2, eBoundTheta) = -sinTheta;
  boundToFreeJacobian(eFreeQOverP, eBoundQOverP) = 1;
}
}  // namespace

namespace detail {

Result<BoundState> boundState(
    const GeometryContext& geoContext, BoundSymMatrix& boundMatrix,
    FreeSymMatrix& freeMatrix, FreeToBoundMatrix& localToGlobalCorrelation,
    BoundMatrix& jacobian, BoundMatrix& correctedJacobian,
    BoundToFreeMatrix& startBoundToFinalFreeJacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& boundToFreeJacobian, FreeVector& parameters,
    bool covTransport, double accumulatedPath, const Surface& surface,
    bool localToGlobalCorrection, bool globalToLocalCorrection) {
  // Covariance transport
  std::optional<BoundSymMatrix> boundCov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    correctedJacobian = BoundMatrix::Identity();
    startBoundToFinalFreeJacobian = BoundToFreeMatrix::Identity();
    // Calculate the jacobian and transport the boundMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToBound(
        geoContext, boundMatrix, freeMatrix, localToGlobalCorrelation, jacobian,
        correctedJacobian, startBoundToFinalFreeJacobian, transportJacobian,
        derivatives, boundToFreeJacobian, parameters, surface,
        localToGlobalCorrection, globalToLocalCorrection);
    // The boundToFreeJacobian and derivatives are reinitialized
  }
  if (boundMatrix != BoundSymMatrix::Zero()) {
    boundCov = boundMatrix;
  } else {
    std::cout << "WARNING: boundMatrix will be zero" << std::endl;
  }

  // Create the bound parameters
  Result<BoundVector> bv =
      detail::transformFreeToBoundParameters(parameters, surface, geoContext);
  if (!bv.ok()) {
    return bv.error();
  }

  /*
  // Calculate corrected bound vec, bound cov and jacobian (they are meaning
  // only when globalToLocalCorrection is required) The impact of
  // localToGlobalCorrection is already contained in the bp, bv and jacobian
  /////////////////////////////////////////////////////////////////////////
  BoundVector correctedBp = BoundVector::Zero();
  BoundSymMatrix correctedBv = BoundSymMatrix::Zero();
  BoundMatrix correctedJacobian = BoundMatrix::Zero();

  // The correctedJacobian has different meaning for globalToLocalCorrection is
  // true or false
  // - if localToGlobal is true, the correctedJacobian includes the bound
  // covariance at the final surface
  // - if localToGlobal is false, the correctedJacobian is a real jacobian
  if (globalToLocalCorrection) {
    std::cout << "Without correction, bound vec = \n"
              << bv.value() << "\n bound cov = \n"
              << boundCov.value() << "\n jacobian = \n"
              << jacobian << std::endl;
    // Note the derivative must not be reinitialized
    std::cout << "derivatives are \n" << derivatives << std::endl;

    auto transformer = detail::CorrectedFreeToBoundTransformer();
    auto correctedRes =
        transformer(parameters, freeMatrix, derivatives, surface, geoContext);
    if (correctedRes.has_value()) {
      auto correctedValue = correctedRes.value();
      correctedBp = std::get<0>(correctedValue);
      correctedBv = std::get<1>(correctedValue);
      BoundToFreeMatrix fCross = std::get<2>(correctedValue);

      //std::cout << "freeMatrix.inverse() = \n"
      //          << freeMatrix.inverse() << std::endl;
      correctedJacobian = fCross.transpose() *
                          (freeMatrix.inverse()).transpose() *
                          startBoundToFinalFreeJacobian;

      // This term below should be similar to freeToBoundJacobian *
      // finalFreeCorrection in the else, but actually they are quite different
      // So this option does not work well
      //std::cout<<"fCross.transpose() * (freeMatrix.inverse()).transpose() = \n" << fCross.transpose() * (freeMatrix.inverse()).transpose() << std::endl;

      std::cout << "With correction, bound vec = \n"
                << correctedBp << "\n bound cov = \n"
                << correctedBv << " \n cross term = \n"
                << correctedJacobian << std::endl;
    }
  } else {
    // No correction for final free->bound jacobian
    FreeToBoundMatrix freeToBoundJacobian =
        surface.freeToBoundJacobian(geoContext, parameters);

    const FreeToPathMatrix freeToPath =
        surface.freeToPathDerivative(geoContext, parameters);
    // The derivatives are already initialized to zero
    FreeMatrix finalFreeCorrection =
        (FreeMatrix::Identity() + derivatives * freeToPath);

    std::cout << "freeToBoundJacobian * finalFreeCorrection =\n"
              << freeToBoundJacobian * finalFreeCorrection << std::endl;
    //std::cout<<"startBoundToFinalFreeJacobian =\n" << startBoundToFinalFreeJacobian << std::endl;
    //std::cout<<"finalFreeCorrection =\n" << finalFreeCorrection << std::endl;
    //std::cout<<"freeToBoundJacobian =\n" << freeToBoundJacobian << std::endl;
    correctedJacobian = freeToBoundJacobian * finalFreeCorrection *
                        startBoundToFinalFreeJacobian;
  }

  /////////////////////////////////////////////////////////////////////////

  // The correctedJacobian should be equal to jacobian in case no gtol and ltog
  // correction
  std::cout << "The jacobian without correction is \n" << jacobian << std::endl;
  std::cout << "The jacobian with correction is \n"
            << correctedJacobian << std::endl;
*/

  // BoundVector correctedBp = BoundVector::Zero();
  // BoundSymMatrix correctedBv = BoundSymMatrix::Zero();
  //  BoundMatrix cJacobian = BoundMatrix::Zero();

  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), *bv, std::move(boundCov)),
      jacobian, accumulatedPath, *bv, *boundCov, correctedJacobian);
}

CurvilinearState curvilinearState(BoundSymMatrix& boundMatrix,
                                  BoundMatrix& jacobian,
                                  FreeMatrix& transportJacobian,
                                  FreeVector& derivatives,
                                  BoundToFreeMatrix& boundToFreeJacobian,
                                  const FreeVector& parameters,
                                  bool covTransport, double accumulatedPath) {
  const Vector3& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSymMatrix> boundCov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the boundMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToCurvilinear(boundMatrix, jacobian, transportJacobian,
                                     derivatives, boundToFreeJacobian,
                                     direction);
  }
  if (boundMatrix != BoundSymMatrix::Zero()) {
    boundCov = boundMatrix;
  }

  // Create the curvilinear parameters
  Vector4 pos4 = Vector4::Zero();
  pos4[ePos0] = parameters[eFreePos0];
  pos4[ePos1] = parameters[eFreePos1];
  pos4[ePos2] = parameters[eFreePos2];
  pos4[eTime] = parameters[eFreeTime];
  CurvilinearTrackParameters curvilinearParams(
      pos4, direction, parameters[eFreeQOverP], std::move(boundCov));
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParams), jacobian,
                         accumulatedPath);
}

// After the transportCovarianceToBound, insures the free parameters and bound
// covariance is reliable but the free covariance might not be correct
void transportCovarianceToBound(
    const GeometryContext& geoContext, BoundSymMatrix& boundCovariance,
    FreeSymMatrix& freeCovariance, FreeToBoundMatrix& localToGlobalCorrelation,
    BoundMatrix& fullTransportJacobian, BoundMatrix& correctedJacobian,
    BoundToFreeMatrix& startBoundToFinalFreeJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, FreeVector& freeParameters,
    const Surface& surface, bool localToGlobalCorrection,
    bool globalToLocalCorrection) {

  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  FreeMatrix finalFreeCorrection = FreeMatrix::Zero();
  // @note this can be removed?
  boundToBoundJacobian(geoContext, freeParameters, boundToFreeJacobian,
                       freeTransportJacobian, freeDerivatives,
                       fullTransportJacobian, startBoundToFinalFreeJacobian,
                       finalFreeCorrection, surface);

  // const FreeToPathMatrix freeToPath =
  //    surface.freeToPathDerivative(geoContext, freeParameters);
  // FreeMatrix finalFreeCorrection = (FreeMatrix::Identity() + freeDerivatives
  // * freeToPath);

  // 1. localToGlobalCorrection determines if the starting freeCov is trustible
  // or not.
  // - If localToGlobalCorrection is true, then use the starting freeCov
  // (which is more accurate);
  // - Otherwise, use the starting boundCov
  //
  // 2. globalToLocalCorrection determins if the final free covariance considers
  // the path correction or not, and whether to do the correction
  // - If globalToLocalCorrection is true, then don't
  // consider the path correction at this point
  if (localToGlobalCorrection) {
    // Should this be calculated using corrected parameters if
    // globalToLocalCorrection is true?
    FreeToBoundMatrix freeToBoundJacobian =
        surface.freeToBoundJacobian(geoContext, freeParameters);

    // Reset (the meaning has been changed)
    startBoundToFinalFreeJacobian =
        freeTransportJacobian * localToGlobalCorrelation.transpose();

    if (globalToLocalCorrection) {
      // The free cov does not include the path correction now
      freeCovariance = freeTransportJacobian * freeCovariance *
                       freeTransportJacobian.transpose();

      auto transformer = detail::CorrectedFreeToBoundTransformer();
      auto correctedRes = transformer(freeParameters, freeCovariance,
                                      freeDerivatives, surface, geoContext);
      if (correctedRes.has_value()) {
        auto correctedValue = correctedRes.value();
        BoundVector boundParams = std::get<0>(correctedValue);
        // 1. Update the free parameters with the corrected bound parameters
        freeParameters = detail::transformBoundToFreeParameters(
            surface, geoContext, boundParams);

        // 2. Update the bound covariance
        boundCovariance = std::get<1>(correctedValue);

        // 3. Get the corrected jacobian (which is not exactly a jacobian)
        BoundToFreeMatrix fCross = std::get<2>(correctedValue);
        correctedJacobian = fCross.transpose() *
                            (freeCovariance.inverse()).transpose() *
                            startBoundToFinalFreeJacobian;
      } else {
        //std::cout << "WARNING: globalToLocal correction failed" << std::endl;
        // Correction failed. No correction
        boundCovariance = freeToBoundJacobian * finalFreeCorrection *
                          freeCovariance * finalFreeCorrection.transpose() *
                          freeToBoundJacobian.transpose();
        correctedJacobian = fullTransportJacobian;
      }
    } else {
      //std::cout << "No globalToLocalCorrection " << std::endl;
      freeCovariance = finalFreeCorrection * freeTransportJacobian *
                       freeCovariance * freeTransportJacobian.transpose() *
                       finalFreeCorrection.transpose();
      boundCovariance = freeToBoundJacobian * freeCovariance *
                        freeToBoundJacobian.transpose();

      // Get the corrected jacobian (which is not exactly a jacobian)
      correctedJacobian = freeToBoundJacobian * finalFreeCorrection *
                          startBoundToFinalFreeJacobian;
    }
  } else {
    //std::cout << "No localToGlobalCorrection " << std::endl;
    // startBoundToFinalFreeJacobian is already transfortJacobian *
    // J^start(B->F)

    if (globalToLocalCorrection) {
      //std::cout << "globalToLocalCorrection " << std::endl;
      // The free cov does not include the path correction now
      freeCovariance = startBoundToFinalFreeJacobian * boundCovariance *
                       startBoundToFinalFreeJacobian.transpose();

      auto transformer = detail::CorrectedFreeToBoundTransformer();
      auto correctedRes = transformer(freeParameters, freeCovariance,
                                      freeDerivatives, surface, geoContext);
      if (correctedRes.has_value()) {
        auto correctedValue = correctedRes.value();
        BoundVector boundParams = std::get<0>(correctedValue);
        // 1. Update the free parameters
        freeParameters = detail::transformBoundToFreeParameters(
            surface, geoContext, boundParams);

        // 2. Update the bound covariance
        boundCovariance = std::get<1>(correctedValue);

        // 3. Get the corrected jacobian
        BoundToFreeMatrix fCross = std::get<2>(correctedValue);
        correctedJacobian = fCross.transpose() *
                            (freeCovariance.inverse()).transpose() *
                            startBoundToFinalFreeJacobian;
      } else {
        //std::cout << "WARNING: globalToLocal correction failed" << std::endl;
        correctedJacobian = fullTransportJacobian;
        boundCovariance = fullTransportJacobian * boundCovariance *
                          fullTransportJacobian.transpose();
      }
    } else {
      //std::cout << "No globalToLocalCorrection " << std::endl;
      // No need to update the freeCovariance since it's not used?
      freeCovariance = finalFreeCorrection * startBoundToFinalFreeJacobian *
                       boundCovariance *
                       startBoundToFinalFreeJacobian.transpose() *
                       finalFreeCorrection.transpose();
      boundCovariance = fullTransportJacobian * boundCovariance *
                        fullTransportJacobian.transpose();
      correctedJacobian = fullTransportJacobian;
    }
  }

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, freeTransportJacobian, freeDerivatives,
                        boundToFreeJacobian, freeParameters, surface);

  if (localToGlobalCorrection) {
    // update the localToGlobalCorrection to that for the final surface
    // Is this necessary?
    localToGlobalCorrelation =
        boundCovariance * boundToFreeJacobian.transpose();
  }
}

void transportCovarianceToCurvilinear(BoundSymMatrix& boundCovariance,
                                      BoundMatrix& fullTransportJacobian,
                                      FreeMatrix& freeTransportJacobian,
                                      FreeVector& freeDerivatives,
                                      BoundToFreeMatrix& boundToFreeJacobian,
                                      const Vector3& direction) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  boundToCurvilinearJacobian(direction, boundToFreeJacobian,
                             freeTransportJacobian, freeDerivatives,
                             fullTransportJacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  boundCovariance = fullTransportJacobian * boundCovariance *
                    fullTransportJacobian.transpose();

  // Reinitialize jacobian components:
  // ->The free transportJacobian is reinitialized to Identity
  // ->The path derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeTransportJacobian, freeDerivatives,
                        boundToFreeJacobian, direction);
}

}  // namespace detail
}  // namespace Acts
