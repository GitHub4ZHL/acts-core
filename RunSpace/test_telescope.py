#!/usr/bin/python
import pathlib, acts, acts.examples
from pathlib import Path
from typing import Optional, Union
from acts import UnitConstants as u
from acts.examples.geant4 import TelescopeG4DetectorConstructionFactory 
from acts.examples import TelescopeDetector

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addTelescopeSeeding,
    addCKFTracks,
    addTruthTrackingGsf,
    TrackSelectorConfig,
    CkfConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)

# Config
u = acts.UnitConstants
field = acts.ConstantBField(acts.Vector3(0, 0, 0)) # u.T


teleG4Config=TelescopeDetector.Config();
teleG4Config.bounds=[14.08, 28.16]
teleG4Config.positions=[30, 60, 90, 105, 120, 150, 180]
teleG4Config.stereos=[0, 0, 0, 0, 0, 0, 0]
teleG4Config.thickness = [80*u.um, 80*u.um, 80*u.um, 1*u.um, 80*u.um, 80*u.um, 80*u.um]
teleG4Config.binValue=0

detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
    bounds=[14.08, 28.16],
    positions=[30, 60, 90, 105, 120, 150, 180],
    stereos=[0, 0, 0, 0, 0, 0, 0],
    thickness=[80*u.um, 80*u.um, 80*u.um, 1*u.um, 80*u.um, 80*u.um, 80*u.um],
    binValue=0,
)

#chi2Cut=21.11
chi2Cut=18.42
nbranch=2
mul=1

digifile="../Examples/Algorithms/Digitization/share/default-digi-config-telescope.json"
#digifile="../Examples/Algorithms/Digitization/share/default-digi-config-telescope-time.json"
#digifile="../Examples/Algorithms/Digitization/share/default-digi-config-telescope-time2.json"


outputDir = Path.cwd() / f"result-without-time-mul{mul}-chi2Cut{chi2Cut}-nbranch{nbranch}-test"
#outputDir = Path.cwd() / f"result-with-time-mul{mul}-chi2Cut{chi2Cut}-nbranch{nbranch}-test"
#outputDir = Path.cwd() / f"result-with-time2-mul{mul}-chi2Cut{chi2Cut}-nbranch{nbranch}-test"

field = acts.ConstantBField(acts.Vector3(0, 0, 0)) # u.T
if not outputDir.exists():
    outputDir.mkdir()

rnd = acts.examples.RandomNumbers(seed=32)
s = acts.examples.Sequencer(events=10000, numThreads=1, outputDir=str(outputDir))

addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0113, 0.0113, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.649 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False),
    multiplicity=mul,
    #multiplicity=mul, # To modify the value of multiplicity
    rnd=rnd,
    vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0.*u.mm, 3.04*u.mm, 3.04*u.mm, 1.52*u.ns)),
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
        eta=(-3, 3),
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
            rho=(0.0, 30 * u.mm),
            absZ=(-30.0, 30.0 * u.m),
            eta=(-1.0, 1.0),
            pt=(1 * u.GeV, None),
            removeNeutral=True,
    ),
    outputDirRoot=outputDir,
    killVolume=trackingGeometry.worldVolume,
    killAfterTime=1000 * u.ns,
    #logLevel=acts.logging.VERBOSE,
    )


addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=Path(digifile),
    #digiConfigFile=Path("../Examples/Algorithms/Digitization/share/default-smearing-config-telescope.json"),
    outputDirRoot=outputDir,
    rnd=rnd,
    #logLevel=acts.logging.VERBOSE
)


addTelescopeSeeding(
    s,
    trackingGeometry,
    initialSigmas={1, 1, 1, 1, 1, 1},
    initialVarInflation={1, 1, 1, 1, 1, 1},
    outputDirRoot=outputDir,
    #logLevel=acts.logging.VERBOSE,
)

addCKFTracks(
    s, 
    trackingGeometry,
    field, 
    CkfConfig(
        chi2CutOff=chi2Cut,
        numMeasurementsCutOff=nbranch,
        #chi2CutOff=15.0,
        #numMeasurementsCutOff=10,
    ),
    outputDirRoot=outputDir,
    logLevel=acts.logging.DEBUG,
)

s.run()
