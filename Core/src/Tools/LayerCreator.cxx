///////////////////////////////////////////////////////////////////
// LayerCreator.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Core module
#include "ACTS/Utilities/Definitions.h"
// GeometryTools
#include "ACTS/Tools/LayerCreator.h"
#include "ACTS/Utilities/OverlapDescriptor.h"
#include "ACTS/Layers/CylinderLayer.h"
#include "ACTS/Layers/DiscLayer.h"
#include "ACTS/Detector/GenericOverlapDescriptor.h"
#include "ACTS/Surfaces/CylinderBounds.h"
#include "ACTS/Surfaces/RadialBounds.h"
#include "ACTS/Surfaces/PlanarBounds.h"

// constructor
Acts::LayerCreator::LayerCreator(const Acts::LayerCreator::Config& lcConfig) :
  m_config()
{
    setConfiguration(lcConfig);
}

void Acts::LayerCreator::setConfiguration(const Acts::LayerCreator::Config& lcConfig)
{
    // @TODO check consistency
    // copy the configuration
    m_config = lcConfig;
}

Acts::LayerPtr Acts::LayerCreator::cylinderLayer(const std::vector<const Surface*>& surfaces,
                                                 double envelopeR, double envelopeZ,
                                                 size_t binsPhi, size_t binsZ) const
{
    // loop over the surfaces and estimate 
    double minR = 10e10;
    double maxR = -10e10;
    double minZ = 10e10;
    double maxZ = -10e10;
    double minPhi = 10;
    double maxPhi = -10;
    // make a surface loop and get the extends
    for (auto& surface : surfaces)
        moduleExtend(*surface, minR, maxR, minPhi, maxPhi, minZ, maxZ);
    
    // harmonize the phi boundaries @TODO - allow for sectorally filled arrrays later
    double phiStep = (maxPhi-minPhi)/(binsPhi-1);
    minPhi -= 0.5*phiStep;
    maxPhi += 0.5*phiStep;    
    // remaining layer parameters
    double layerR = 0.5*(minR+maxR);
    double layerHalfZ = maxZ;
    double layerThickness = (maxR-minR)+2*envelopeR;    
    
    // adjust the layer radius 
    layerR = 0.5*(minR+maxR);
    //MSG_VERBOSE("Creating a cylindrical Layer:");
    //MSG_VERBOSE(" - with layer R    = " << layerR);
    //MSG_VERBOSE(" - from R min/max  = " << minR            << " / "  << maxR);
    //MSG_VERBOSE(" - with z min/max  = " << -layerHalfZ     << " / " << layerHalfZ);
    //MSG_VERBOSE(" - and phi min/max = " <<  minPhi         << " / " << maxPhi);
    //MSG_VERBOSE(" - # of modules    = " << surfaces.size() << " ordered in ( " << binsPhi << " x " << binsZ << ")");
    
    // create the surface array
    SurfaceArray* sArray = m_config.surfaceArrayCreator->surfaceArrayOnCylinder(surfaces, layerR, minPhi, maxPhi, layerHalfZ, binsPhi, binsZ);

    // create the layer and push it back
    std::shared_ptr<const CylinderBounds> cBounds(new CylinderBounds(layerR, layerHalfZ+envelopeZ));
    
    // create the layer
    LayerPtr cLayer = CylinderLayer::create(nullptr, cBounds, sArray, layerThickness, new GenericOverlapDescriptor(), nullptr, active);
    
    // now return
    return cLayer;
} 

Acts::LayerPtr Acts::LayerCreator::discLayer(const std::vector<const Surface*>& surfaces,
                                             double envelopeMinR, double envelopeMaxR, double envelopeZ,
                                             size_t binsR, size_t binsPhi,
                                             const std::vector<double>& rBoundaries) const
{
    // loop over the surfaces and estimate 
    double minR = 10e10;
    double maxR = 0.;
    double minZ = 10e10;
    double maxZ = -10e10;
    double minPhi = 10;
    double maxPhi = -10;
    
    // make a surface loop and get the extends
    for (auto& surface : surfaces)
        moduleExtend(*surface, minR, maxR, minPhi, maxPhi, minZ, maxZ);
    // harmonize the phi boundaries @TODO - allow for sectorally filled arrrays later
    double phiStep = (maxPhi-minPhi)/(binsPhi-1);
    minPhi -= 0.5*phiStep;
    maxPhi += 0.5*phiStep;    
    // layer parametres
    double layerZ         = 0.5*(minZ+maxZ);
    double layerThickness = (maxZ-minZ)+2*envelopeZ;
    
    // create the surface array
    SurfaceArray* sArray = m_config.surfaceArrayCreator->surfaceArrayOnDisc(surfaces, minR, maxR, minPhi, maxPhi, binsR, binsPhi, rBoundaries);
    
    // create the share disc bounds
    std::shared_ptr<const DiscBounds> dBounds(new RadialBounds(minR-envelopeMinR,maxR+envelopeMaxR));
    
    // create the layer transforms
    Transform3D* transform = new Transform3D(Transform3D::Identity());
    transform->translation() = Vector3D(0.,0.,layerZ);

    // create the layers
    LayerPtr dLayer = DiscLayer::create(std::shared_ptr<Transform3D>(transform), 
                                        dBounds,
                                        sArray,
                                        layerThickness,
                                        new GenericOverlapDescriptor(),
                                        nullptr,
                                        active);
    // return the layer                                    
    return dLayer;
}

Acts::LayerPtr Acts::LayerCreator::planeLayer(const std::vector<const Surface*>& /**surfaces*/,
                                              double /**envelopeXY*/, double /**envelopeZ*/,
                                              size_t /**binsX*/, size_t /**binsY*/) const
{
    //@TODO implement
    return nullptr;
}

void Acts::LayerCreator::moduleExtend(const Surface& sf, 
                                      double& minR, double& maxR, 
                                      double& minPhi, double& maxPhi, 
                                      double& minZ, double& maxZ) const
{
    // get the associated detector element
    const DetectorElementBase* element = sf.associatedDetectorElement();
    if (element){
        // get the thickness
        double thickness = element->thickness();
        // check the shape
        const PlanarBounds* pBounds = dynamic_cast<const PlanarBounds*>(&(sf.bounds()));
        if (pBounds){
            // phi is always from the center for planar surfaces
            takeSmallerBigger(minPhi,maxPhi,sf.center().phi());
            // get the vertices
            std::vector< Vector2D > vertices = pBounds->vertices();
            size_t nVertices = vertices.size();
            // loop over the two sides of the module 
            for (int side = 0; side < 2; ++side){
                // loop over the vertex combinations 
                for (size_t iv = 0; iv < nVertices; ++iv){
                     size_t ivp = iv ? iv-1 : nVertices-1;
                     // thickness
                     double locz = side ? 0.5*thickness : -0.5*thickness;
                     // p1 & p2 vectors 
                     Vector3D p2(sf.transform()*Vector3D(vertices[iv].x(),vertices[iv].y(),locz));
                     Vector3D p1(sf.transform()*Vector3D(vertices[ivp].x(),vertices[ivp].y(),locz));
                     // let's get
                     takeSmallerBigger(minZ,maxZ,p2.z());
                     takeBigger(maxR,p2.perp());
                     takeSmaller(minR,radialDistance(p1,p2));
                }
            }
        } else 
            //MSG_WARNING("Not implemented yet for Non-Planar bounds")
            ;
    }
}       

double Acts::LayerCreator::radialDistance(const Vector3D& pos1, const Vector3D& pos2) const 
{
    
    // following nominclature found in header file and doxygen documentation
    // line one is the straight track
    const Vector3D&  ma  = pos1;
    const Vector3D   ea  = (pos2-pos1).unit();
    // line two is the line surface
    Vector3D mb(0.,0.,0);
    Vector3D eb(0.,0.,1.);
    // now go ahead and solve for the closest approach
    Vector3D  mab(mb - ma);
    double eaTeb = ea.dot(eb);
    double denom = 1 - eaTeb*eaTeb;
    if (fabs(denom)>10e-7){
       double lambda0 = (mab.dot(ea) - mab.dot(eb)*eaTeb)/denom;
       // evaluate validaty in terms of bounds
       if (lambda0 < 1. && lambda0 > 0.) return (ma+lambda0*ea).perp();
       return lambda0 < 0. ? pos1.perp() : pos2.perp();
    }
    return 10e101;    
}
