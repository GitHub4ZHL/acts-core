// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"

#include "Acts/Material/Interactions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {
namespace detail {
void PointwiseMaterialInteraction::evaluatePointwiseMaterialInteraction(
    bool multipleScattering, bool energyLoss) {
  if (energyLoss) {
     Eloss = computeEnergyLossBethe(slab, pdg, mass, qOverP, q);
    
    //Eloss = computeEnergyLossLandau(slab, pdg, mass, qOverP, q);
 
//    double p_ = 1./std::abs(qOverP); 
//    if(p_<0.08){
//       Eloss = computeEnergyLossBethe(slab, pdg, mass, qOverP, q);
//    } else if (p_<0.14){
//       Eloss = Eloss*0.8;
//    }


  }
  // Compute contributions from interactions
  if (performCovarianceTransport) {
    covarianceContributions(multipleScattering, energyLoss);
  }
}

void PointwiseMaterialInteraction::covarianceContributions(
    bool multipleScattering, bool energyLoss) {
  // Compute contributions from interactions
  if (multipleScattering) {
    // TODO use momentum before or after energy loss in backward mode?
    const auto theta0 =
        computeMultipleScatteringTheta0(slab, pdg, mass, qOverP, q);
    // sigmaPhi = theta0 / sin(theta)
    const auto sigmaPhi = theta0 * (dir.norm() / VectorHelpers::perp(dir));
    variancePhi = sigmaPhi * sigmaPhi;
    // sigmaTheta = theta0
    varianceTheta = theta0 * theta0;
  }
  // TODO just ionisation loss or full energy loss?
  if (energyLoss) {
    const auto sigmaQoverP =
        computeEnergyLossLandauSigmaQOverP(slab, pdg, mass, qOverP, q);
    varianceQoverP = sigmaQoverP * sigmaQoverP;
  }
}

double PointwiseMaterialInteraction::updateVariance(
    double variance, double change, NoiseUpdateMode updateMode) const {
  // Add/Subtract the change
  // Protect the variance against becoming negative
  //return std::max(0., variance + std::copysign(change, updateMode));
  return std::max(0., variance + change);
}
}  // namespace detail
}  // end of namespace Acts
