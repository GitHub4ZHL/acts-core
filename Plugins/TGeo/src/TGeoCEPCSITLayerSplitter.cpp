// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoCEPCSITLayerSplitter.hpp"

#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "TGeoBBox.h"

Acts::TGeoCEPCSITLayerSplitter::TGeoCEPCSITLayerSplitter(
    const TGeoCEPCSITLayerSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoCEPCSITLayerSplitter::split(
    const GeometryContext& /*gctx*/,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  // Thickness
  auto tgIdentifier = tgde->identifier();
  const TGeoNode& tgNode = tgde->tgeoNode();
  std::cout << "TGeoCEPCSITLayerSplitter splitting detElement "
            << tgNode.GetName() << std::endl;

  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
      tgDetectorElements = {};

  // The ladders of the 'layer_**' on this layer
  auto daugthers = tgNode.GetVolume()->GetNodes();
  // The relative position of layer w.r.t. volume 
  const TGeoMatrix* lmatrix = tgNode.GetMatrix();

  TIter iObj(daugthers);
  while (TObject* obj = iObj()) {
    TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
    if (node != nullptr) {
      std::string nodeName = node->GetName();
      // The relative position of module w.r.t. layer
      const TGeoMatrix* nmatrix = node->GetMatrix();
      if (nodeName.find("_LadderLogical_") == std::string::npos) {
        continue; 
      }
      std::cout << tgNode.GetName() << " has ladder: " << nodeName
         << std::endl;

      auto senEnvelopeLogicals = node->GetVolume()->GetNodes();
      TIter jObj(senEnvelopeLogicals);
      int nSenEnvelopeLogical=0; 
      while (TObject* dobj = jObj()) {
        TGeoNode* dNode = dynamic_cast<TGeoNode*>(dobj);
        if (dNode == nullptr){
          throw std::runtime_error("null ptr in SIT/SET splitter");
        }
       
        std::string dNodeName = dNode->GetName();
        const TGeoMatrix* dmatrix = dNode->GetMatrix();
        if(dNodeName.find("SenEnvelopeLogical_") == std::string::npos){
          continue;    
        }	    

        auto senLogicals = dNode->GetVolume()->GetNodes();
        TIter kObj(senLogicals);
        while (TObject* ddobj = kObj()) {
          TGeoNode* ddNode = dynamic_cast<TGeoNode*>(ddobj);
          if (ddNode == nullptr){
            throw std::runtime_error("null ptr in SIT/SET splitter");
          }
          const TGeoMatrix* ddmatrix = ddNode->GetMatrix();

          TGeoBBox* sensor =
          dynamic_cast<TGeoBBox*>(ddNode->GetVolume()->GetShape());
      
          if (sensor == nullptr) {
            throw std::runtime_error("cast bad (SIT/SET sensor is always a box. This should not happen)");
          } 
          
          ActsScalar dx = sensor->GetDX()*m_cfg.unitScalar;
          ActsScalar dy = sensor->GetDY()*m_cfg.unitScalar;
          ActsScalar dz = sensor->GetDZ()*m_cfg.unitScalar;


          TGeoHMatrix transform =
              TGeoCombiTrans(*lmatrix) * TGeoCombiTrans(*nmatrix) * TGeoCombiTrans(*dmatrix) * TGeoCombiTrans(*ddmatrix);

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

          auto tgRectangle =
          std::make_shared<Acts::RectangleBounds>(dx, dy);

          auto tgDetectorElement = std::make_shared<Acts::TGeoDetectorElement>(
              tgIdentifier, *node, etrf, tgRectangle, dz);


          tgDetectorElements.push_back(tgDetectorElement);
	} // loop over the sensors  
      } // loop over the objects of this ladder
    } // if the ladder has node 
  } // loop over ladders

  std::cout << "Found " << tgDetectorElements.size()
            << " SIT/SET detElements on this layer" << std::endl;

  if (not tgDetectorElements.empty()) {
    return tgDetectorElements;
  }

  return {tgde};
}
