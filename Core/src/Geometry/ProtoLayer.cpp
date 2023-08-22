// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ProtoLayer.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>

using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace Acts {

ProtoLayer::ProtoLayer(const GeometryContext& gctx,
                       const std::vector<const Surface*>& surfaces)
    : m_surfaces(surfaces) {
  measure(gctx, surfaces);
}

ProtoLayer::ProtoLayer(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<const Surface>>& surfaces)
    : m_surfaces(unpack_shared_vector(surfaces)) {
  measure(gctx, m_surfaces);
}

double ProtoLayer::min(BinningValue bval, bool addenv) const {
  if (addenv) {
    //std::cout << "extent.min" << extent.min(bval) << std::endl;
    return extent.min(bval) - envelope[bval].first;
  }
  return extent.min(bval);
}

double ProtoLayer::max(BinningValue bval, bool addenv) const {
  if (addenv) {
    //std::cout << "extent.max" << extent.max(bval) << std::endl;
    return extent.max(bval) + envelope[bval].second;
  }
  return extent.max(bval);
}

double ProtoLayer::medium(BinningValue bval, bool addenv) const {
  if(not m_surfaces.empty()){ 
    if (bval == binR and m_surfaces[0]->type() == Surface::Straw) {
      //std::cout<<"ProtoLayer::medium" << std::endl;      
      return (m_rMin + m_rMax) / 2;
    }
  }
  return 0.5 * (min(bval, addenv) + max(bval, addenv));
}

double ProtoLayer::range(BinningValue bval, bool addenv) const {
  if(not m_surfaces.empty()){ 
    if (bval == binR and m_surfaces[0]->type() == Surface::Straw) {
      //std::cout<<"ProtoLayer::range" << std::endl;      
      auto detElement = m_surfaces[0]->associatedDetectorElement(); 
      //return detElement->thickness()*0.2;
      //return 8;
      //return 2;
      return 1;
    } else if (bval == binR and m_surfaces[0]->type() == Surface::Cylinder){
      //return 6;
      //return 2;
      return 1;
    }
  }
  return std::abs(max(bval, addenv) - min(bval, addenv));
}

std::ostream& ProtoLayer::toStream(std::ostream& sl) const {
  sl << "ProtoLayer with dimensions (min/max)" << std::endl;
  extent.toStream(sl);
  return sl;
}

void ProtoLayer::measure(const GeometryContext& gctx,
                         const std::vector<const Surface*>& surfaces) {
  for (const auto& sf : surfaces) {
    Vector3 center = sf->center(gctx);
    ActsScalar radius = std::hypot(center.x(), center.y());
    if (radius < m_rMin) {
      m_rMin = radius;
    }
    if (radius > m_rMax) {
      m_rMax = radius;
    }
    auto sfPolyhedron = sf->polyhedronRepresentation(gctx, 1);
    const DetectorElementBase* element = sf->associatedDetectorElement();
    if (element != nullptr) {
      // Take the thickness in account if necessary
      double thickness = element->thickness();
      // We need a translation along and opposite half thickness
      Vector3 sfNormal = sf->normal(gctx, sf->center(gctx));
      std::vector<double> deltaT = {-0.5 * thickness, 0.5 * thickness};
      for (const auto& dT : deltaT) {
        Transform3 dtransform = Transform3::Identity();
        dtransform.pretranslate(dT * sfNormal);
        extent.extend(sfPolyhedron.extent(dtransform));
      }
      continue;
    }
    extent.extend(sfPolyhedron.extent());
  }
}

void ProtoLayer::add(const GeometryContext& gctx, const Surface& surface) {
  m_surfaces.push_back(&surface);
  measure(gctx, m_surfaces);
}

}  // namespace Acts
