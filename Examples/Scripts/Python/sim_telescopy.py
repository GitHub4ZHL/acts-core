#!/usr/bin/env python3
import pathlib, acts, acts.examples
from pathlib import Path
from typing import Optional, Union
from acts import UnitConstants as u

from acts.examples import (
    FixedMultiplicityGenerator,
)

from acts.examples.geant4 import TelescopeG4DetectorConstructionFactory 
from acts.examples import TelescopeDetector


from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
    addGeant4
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg, 
    addCKFTracks,
    TrackSelectorConfig,
    CkfConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)


teleG4Config=TelescopeDetector.Config();
teleG4Config.bounds=[200, 200]
teleG4Config.positions=[30, 60, 90, 100, 120, 150, 180]
teleG4Config.thickness=[0.08, 0.08, 0.08, 0.001, 0.08, 0.08, 0.08]
teleG4Config.stereos=[0, 0, 0, 0, 0, 0, 0]
teleG4Config.binValue=0



u = acts.UnitConstants
detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
    bounds=[200, 200],
    positions=[30, 60, 90, 100, 120, 150, 180],
    thickness=[0.08, 0.08, 0.08, 0.001, 0.08, 0.08, 0.08],
    stereos=[0, 0, 0, 0, 0, 0, 0],
    binValue=0,
)


# field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
field = acts.ConstantBField(acts.Vector3(0.0 * u.T, 0, 2 * u.T))
outputDir = Path.cwd() / "result_telescope/multiplicity5"
if not outputDir.exists():
    outputDir.mkdir()

rnd = acts.examples.RandomNumbers(seed=43)
s = acts.examples.Sequencer(events=100, numThreads=1, outputDir=str(outputDir))

addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.2, 0.2, uniform=True),
    PhiConfig(-35*u.degree, 35 * u.degree), 
    #PhiConfig(0.0*u.degree, 0.0 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False), #The first num specifies the number of particles from each vertex
    multiplicity=5, #Fixed number of vertex.TODO: change it to number sampled using poission distribution
    rnd=rnd,
    vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(1*u.mm, 0.01*u.mm, 0.01*u.mm, 0.5*u.ns)),
)

'''
addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    preSelectParticles=ParticleSelectorConfig(
        rho=(0.0 * u.mm, 300.0 * u.mm),
        absZ=(0.0 * u.mm, 200.0 * u.mm),
        eta=(-0.8, 0.8),
        pt=(1 * u.GeV, None),
        removeNeutral=True,
    ),
    outputDirRoot=outputDir,
)
'''

addGeant4(
    s,
    detector=None,
    trackingGeometry=trackingGeometry,
    field=field,
    rnd=rnd,
    g4DetectorConstructionFactory=TelescopeG4DetectorConstructionFactory(teleG4Config),
    preSelectParticles=ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(150 * u.MeV, None),
            removeNeutral=True,
    ),
    outputDirRoot=outputDir,
    killVolume=trackingGeometry.worldVolume,
    killAfterTime=1000 * u.ns,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=Path("../../../Examples/Algorithms/Digitization/share/default-digi-config-telescope.json"),
    #digiConfigFile=Path("../../../Examples/Algorithms/Digitization/share/default-smearing-config-telescope.json"),
    #digiConfigFile=Path("../../../Examples/Algorithms/Digitization/share/default-smearing-config-telescope-with-time.json"),
    outputDirRoot=outputDir,
    rnd=rnd,
    logLevel=acts.logging.VERBOSE
)

"""
addSeeding(
    s,
    trackingGeometry,
    field,
    # TruthSeedRanges(pt=(500.0 * u.MeV, None), nHits=(9, None)),
    TruthSeedRanges(pt=(1 * u.GeV, None), eta=(-0.8, 0.8), nHits=(6, None)),

    # ParticleSmearingSigmas(pRel=0.01),  # only used by SeedingAlgorithm.TruthSmeared
    SeedFinderConfigArg(
        # r=(None, 200 * u.mm),  # rMin=default, 33mm
        r=(1, 141.5 * u.mm), # 100*sqrt(2)

        # deltaR=(1 * u.mm, 60 * u.mm),
        deltaR=(0.05 * u.mm, 100 * u.mm),

        collisionRegion=(-10 * u.mm, 10 * u.mm),
        z=(-200 * u.mm, 200 * u.mm), # half 200
        maxSeedsPerSpM=1,
        sigmaScattering=10,
        radLengthPerSeed=0.5,

        minPt=200 * u.MeV,

        cotThetaMax=1, #This is related with the particle eta
        impactMax=10 * u.mm,
    ),

    SeedFinderOptionsArg( 
        beamPos=(0,0), 
        bFieldInZ=2 * u.T, 
    ),

    initialVarInflation=(100,100,100,100,100,1),
    geoSelectionConfigFile=Path("/home/xiaocong/Software/Acts/acts/Examples/Algorithms/TrackFinding/share/geoSelection-telescope.json"),
    outputDirRoot=outputDir,
    rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
    #seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    logLevel=acts.logging.DEBUG
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    CkfConfig(
        chi2CutOff=15.0, 
        numMeasurementsCutOff=10,
    ),
    outputDirRoot=outputDir,
)


addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(maximumSharedHits=3),
    outputDirRoot=outputDir,
)

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Iterative,
    outputDirRoot=outputDir,
)
"""

s.run()
