// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoDriftChamberLayerSplitter.hpp"

#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "TGeoTube.h"

Acts::TGeoDriftChamberLayerSplitter::TGeoDriftChamberLayerSplitter(
    const TGeoDriftChamberLayerSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoDriftChamberLayerSplitter::split(
    const GeometryContext& /*gctx*/,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  // Thickness
  auto tgIdentifier = tgde->identifier();
  const TGeoNode& tgNode = tgde->tgeoNode();
  ACTS_DEBUG("TGeoDriftChamberLayerSplitter splitting detElement "
             << tgNode.GetName());

  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
      tgDetectorElements = {};

  // The objects (no concept of cell for CEPC DC) on this layer
  auto daugthers = tgNode.GetVolume()->GetNodes();
  // The relative position of layer w.r.t. chamber 
  const TGeoMatrix* lmatrix = tgNode.GetMatrix();
  //std::cout<<"layer matrix print " << std::endl;
  //lmatrix->Print();
  
  TGeoTube* layer = dynamic_cast<TGeoTube*>(tgNode.GetVolume()->GetShape());
  ActsScalar parameters[5];
  layer->GetBoundingCylinder(parameters);
  ActsScalar minR = layer->GetRmin() * m_cfg.unitScalar;
  ActsScalar maxR = layer->GetRmax() * m_cfg.unitScalar;

  ActsScalar thickness = maxR-minR;

  std::cout<<"minR = "<< minR << ", maxR " << maxR << std::endl;

  int nSenseWires = 0;
  TIter iObj(daugthers);
  while (TObject* obj = iObj()) {
    TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
    if (node != nullptr) {
      std::string nodeName = node->GetName(); 
      if (nodeName.find("SignalWire_") == std::string::npos) {
        continue;
      }
      const TGeoMatrix* smatrix = node->GetMatrix();
      // Find the inner W wire
      TGeoNode* WNode = node->GetVolume()->FindNode("SignalWireW_0");
      if (WNode != nullptr) {
        TGeoTube* wire = dynamic_cast<TGeoTube*>(WNode->GetVolume()->GetShape());
        if (wire == nullptr) {
          ACTS_WARNING(
              "cast bad (signal wire is always a tube. This is not supposed to "
              "happen)");
        } else {
          const TGeoMatrix* wmatrix = WNode->GetMatrix();
          ActsScalar halfZ = wire->GetDz() * m_cfg.unitScalar;
          ACTS_DEBUG("half length of the sense wire: " << halfZ);
          //std::cout<<"half length of the sense wire: " << halfZ << std::endl;;

          //std::cout<<"wire matrix print " << std::endl;
          //wmatrix->Print();
          TGeoHMatrix transform =
              TGeoCombiTrans(*lmatrix) * TGeoCombiTrans(*smatrix) * TGeoCombiTrans(*wmatrix);

          // Get the placement and orientation in respect to its mother
          const Double_t* rotation = transform.GetRotationMatrix();
          const Double_t* translation = transform.GetTranslation();
          // Create a eigen transform
          Vector3 t(translation[0] * m_cfg.unitScalar,
                    translation[1] * m_cfg.unitScalar,
                    translation[2] * m_cfg.unitScalar);
          //std::cout<<"supposed translation " << t << std::endl;	
          Vector3 cx(rotation[0], rotation[3], rotation[6]);
          Vector3 cy(rotation[1], rotation[4], rotation[7]);
          Vector3 cz(rotation[2], rotation[5], rotation[8]);
          auto etrf = TGeoPrimitivesHelper::makeTransform(cx, cy, cz, t);
          //std::cout<<"cx = "<< cx.transpose() << std::endl;
          //std::cout<<"cy = "<< cy.transpose() << std::endl;
          //std::cout<<"cz = "<< cz.transpose() << std::endl;
          if(rotation[8]!=1){
            throw std::runtime_error("Only axial wire for CEPC DC!");
	  }

          // make a lineBounds
          //auto tgWire = std::make_shared<Acts::LineBounds>(distanceToOuterCorner, halfZ);
          auto tgWire = std::make_shared<Acts::LineBounds>(thickness/2, halfZ);

          // Create a new detector element per split
          auto tgDetectorElement =
              std::make_shared<Acts::TGeoDetectorElement>(
                  tgIdentifier, *node, etrf, tgWire, thickness);

          tgDetectorElements.push_back(tgDetectorElement);
        
          nSenseWires++;
        }
      }
    }
  }

  ACTS_DEBUG("Found " << nSenseWires << " sense wires on this layer");
  std::cout<<"Found " << nSenseWires << " sense wires on this layer" << std::endl;

  if (not tgDetectorElements.empty()) {
    return tgDetectorElements;
  }

  return {tgde};
}
