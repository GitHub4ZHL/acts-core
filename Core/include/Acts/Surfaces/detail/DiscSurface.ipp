// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Vector2D DiscSurface::localPolarToCartesian(
    const Vector2D& lpolar) const {
  return Vector2D(lpolar[eBoundLoc0] * cos(lpolar[eBoundLoc1]),
                  lpolar[eBoundLoc0] * sin(lpolar[eBoundLoc1]));
}

inline Vector2D DiscSurface::localCartesianToPolar(
    const Vector2D& lcart) const {
  return Vector2D(sqrt(lcart[eBoundLoc0] * lcart[eBoundLoc0] +
                       lcart[eBoundLoc1] * lcart[eBoundLoc1]),
                  atan2(lcart[eBoundLoc1], lcart[eBoundLoc0]));
}

inline SurfaceIntersection DiscSurface::intersect(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // Get the contextual transform
  auto gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction);
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status != Intersection3D::Status::unreachable and bcheck and
      m_bounds != nullptr) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    const Vector3D vecLocal(intersection.position - tMatrix.block<3, 1>(0, 3));
    const Vector2D lcartesian =
        tMatrix.block<3, 2>(0, 0).transpose() * vecLocal;
    if (bcheck.type() == BoundaryCheck::Type::eAbsolute and
        m_bounds->coversFullAzimuth()) {
      double tolerance = s_onSurfaceTolerance + bcheck.tolerance()[eBoundLoc0];
      if (not m_bounds->insideRadialBounds(VectorHelpers::perp(lcartesian),
                                           tolerance)) {
        intersection.status = Intersection3D::Status::missed;
      }
    } else if (not insideBounds(localCartesianToPolar(lcartesian), bcheck)) {
      intersection.status = Intersection3D::Status::missed;
    }
  }
  return {intersection, this};
}

inline LocalCartesianToBoundLocalMatrix
DiscSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3D& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3D localPos = sTransform.inverse() * position;
  const double lr = perp(localPos);
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  LocalCartesianToBoundLocalMatrix loc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  loc3DToLocBound << lcphi, lsphi, 0, -lsphi / lr, lcphi / lr, 0;

  return loc3DToLocBound;
}

inline BoundLocalToLocalCartesianMatrix
DiscSurface::boundLocalToLocalCartesianDerivative(
    const GeometryContext& gctx, const Vector3D& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3D localPos = sTransform.inverse() * position;
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  BoundLocalToLocalCartesianMatrix locBoundToLoc3D =
      BoundLocalToLocalCartesianMatrix::Zero();
  locBoundToLoc3D << lcphi, -localPos.y(), lsphi, localPos.x(), 0, 0;

  return locBoundToLoc3D;
}

inline Vector3D DiscSurface::normal(const GeometryContext& gctx,
                                    const Vector2D& /*unused*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline Vector3D DiscSurface::binningPosition(const GeometryContext& gctx,
                                             BinningValue bValue) const {
  if (bValue == binR) {
    double r = m_bounds->binningValueR();
    double phi = m_bounds->binningValuePhi();
    return Vector3D(r * cos(phi), r * sin(phi), center(gctx).z());
  }
  return center(gctx);
}

inline double DiscSurface::binningPositionValue(const GeometryContext& gctx,
                                                BinningValue bValue) const {
  // only modify binR
  if (bValue == binR) {
    return VectorHelpers::perp(center(gctx)) + m_bounds->binningValueR();
  }
  return GeometryObject::binningPositionValue(gctx, bValue);
}

inline double DiscSurface::pathCorrection(const GeometryContext& gctx,
                                          const Vector3D& position,
                                          const Vector3D& direction) const {
  /// we can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, position).dot(direction));
}
