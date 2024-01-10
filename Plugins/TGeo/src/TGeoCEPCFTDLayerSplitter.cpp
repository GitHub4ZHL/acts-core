// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoCEPCFTDLayerSplitter.hpp"

#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "TGeoArb8.h"
#include "TGeoBBox.h"

Acts::TGeoCEPCFTDLayerSplitter::TGeoCEPCFTDLayerSplitter(
    const TGeoCEPCFTDLayerSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoCEPCFTDLayerSplitter::split(
    const GeometryContext& /*gctx*/,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  // Thickness
  auto tgIdentifier = tgde->identifier();
  const TGeoNode& tgNode = tgde->tgeoNode();
  std::cout << "TGeoCEPCFTDLayerSplitter collecting sibling detElement "
            << tgNode.GetName() << std::endl;

  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
      tgDetectorElements = {};

  TGeoTrap* trap = dynamic_cast<TGeoTrap*>(tgNode.GetVolume()->GetShape());
  const TGeoMatrix* mymatrix = tgNode.GetMatrix();
  TGeoHMatrix mytransform = TGeoCombiTrans(*mymatrix);
  const Double_t* mytranslation = mytransform.GetTranslation();
  ActsScalar zPos = mytranslation[2];
  std::cout << tgNode.GetName() << " has zPos =  " << zPos << std::endl;

  // Get its mother volumes
  auto siblings = tgNode.GetMotherVolume()->GetNodes();

  TIter iObj(siblings);
  while (TObject* obj = iObj()) {
    TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
    if (node == nullptr) {
      throw std::runtime_error("The node is invalid");
    }

    const TGeoMatrix* siblingmatrix = node->GetMatrix();

    TGeoHMatrix siblingtransform = TGeoCombiTrans(*siblingmatrix);
    const Double_t* siblingtranslation = siblingtransform.GetTranslation();
    ActsScalar siblingzPos = siblingtranslation[2];

    if (std::abs(siblingzPos - zPos) > 1) {
      continue;
    }

    std::string nodeName = node->GetName();
    std::cout << tgNode.GetName() << " found siblings with name " << nodeName
              << std::endl;

    std::vector<std::string> component0Names = {
        "component0_0", "component0_00x1", "component0_00x2", "component0_00x3",
        "component0_00x4"};
    std::vector<std::string> component1Names = {
        "component1_1", "component1_10x1", "component1_10x2", "component1_10x3",
        "component1_10x4"};

    TGeoNode* component0 = nullptr;
    TGeoNode* component1 = nullptr;
    for (const auto& name : component0Names) {
      component0 = node->GetVolume()->FindNode(name.c_str());
      if (component0) {
        break;
      }
    }
    for (const auto& name : component1Names) {
      component1 = node->GetVolume()->FindNode(name.c_str());
      if (component1) {
        break;
      }
    }
    if (component0 == nullptr) {
      throw std::runtime_error(
          "The node doesn't have daughter named with component0_0");
    }
    if (component1 == nullptr) {
      throw std::runtime_error(
          "The node doesn't have daughter named with component1_1");
    }

    //std::vector<TGeoNode*> dNodes = {component0, component1};
    std::vector<TGeoNode*> dNodes = {component0};

    for (const auto& dNode : dNodes) {
      const TGeoMatrix* sensorMatrix = dNode->GetMatrix();

      TGeoHMatrix transform =
          TGeoCombiTrans(*siblingmatrix) * TGeoCombiTrans(*sensorMatrix);

      // Get the placement and orientation in respect to its mother
      const Double_t* rotation = transform.GetRotationMatrix();
      const Double_t* translation = transform.GetTranslation();
      // Create a eigen transform
      Vector3 t(translation[0] * m_cfg.unitScalar,
                translation[1] * m_cfg.unitScalar,
                translation[2] * m_cfg.unitScalar);
      // std::cout<<"supposed translation " << t << std::endl;
      Vector3 cx(rotation[0], rotation[3], rotation[6]);
      Vector3 cy(rotation[1], rotation[4], rotation[7]);
      Vector3 cz(rotation[2], rotation[5], rotation[8]);
      auto etrf = TGeoPrimitivesHelper::makeTransform(cx, cy, cz, t);

      // Get bounds
      TGeoTrap* siblingTrap =
          dynamic_cast<TGeoTrap*>(node->GetVolume()->GetShape());
      if (siblingTrap == nullptr) {
        throw std::runtime_error("The node should be a trapezoid");
      }

      ActsScalar bl1 = siblingTrap->GetBl1() * m_cfg.unitScalar;
      ActsScalar H1 = siblingTrap->GetH1() * m_cfg.unitScalar;
      ActsScalar tl1 = siblingTrap->GetTl1() * m_cfg.unitScalar;
      ActsScalar dz = siblingTrap->GetDZ() * m_cfg.unitScalar;

      ActsScalar midY = std::hypot(t.x(), t.y());
      ActsScalar minR = (midY - H1) / std::cos(M_PI * 2 / 16 / 2);
      ActsScalar maxR = (midY + H1) / std::cos(M_PI * 2 / 16 / 2);
      ActsScalar halfXminR = minR * std::sin(M_PI * 2 / 16 / 2);
      ActsScalar halfXmaxR = maxR * std::sin(M_PI * 2 / 16 / 2);
      //std::cout << "componet has H1: " << H1 << ", bl1: " << bl1 << ", tl1 "
      //          << tl1 << ", thickness: " << dz;
      //std::cout << ", minR: " << minR << ", maxR: " << maxR
      //          << ", midY: " << midY << ", and z " << t.z() << std::endl;

      auto bounds = std::make_shared<Acts::DiscTrapezoidBounds>(
          halfXminR, halfXmaxR, minR, maxR);

      auto tgDetectorElement = std::make_shared<Acts::TGeoDetectorElement>(
          tgIdentifier, *node, etrf, bounds, dz);

      // test
      // auto surface = &(tgDetectorElement->surface());
      // auto center = surface->center(GeometryContext());
      // std::cout<<"center r = "<< hypot(center.x(), center.y()) <<" and center
      // z " << center.z() << std::endl;

      tgDetectorElements.push_back(tgDetectorElement);
    }
  }  // loop over all siblings

  std::cout << "Found " << tgDetectorElements.size()
            << " FTD detElements on this layer" << std::endl;

  if (not tgDetectorElements.empty()) {
    return tgDetectorElements;
  }

  return {tgde};
}
