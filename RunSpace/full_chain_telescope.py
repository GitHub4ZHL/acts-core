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

teleG4Config=TelescopeDetector.Config();
teleG4Config.bounds=[14.08, 28.16]
teleG4Config.positions=[30, 60, 90, 105, 120, 150, 180]
teleG4Config.stereos=[0, 0, 0, 0, 0, 0, 0]
teleG4Config.thickness = [80*u.um, 80*u.um, 80*u.um, 1*u.um, 80*u.um, 80*u.um, 80*u.um]
teleG4Config.binValue=0

u = acts.UnitConstants

detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
    bounds=[14.08, 28.16],
    positions=[30, 60, 90, 105, 120, 150, 180],
    stereos=[0, 0, 0, 0, 0, 0, 0],
    thickness=[80*u.um, 80*u.um, 80*u.um, 1*u.um, 80*u.um, 80*u.um, 80*u.um],
    binValue=0,
)

# ParticleGun Configurationi
# Modify multiplicity
'''
mul = int(input("multiplicity = "))
if mul == 1:
    outputDir = Path.cwd() / "wot_multiplicity_1"
elif mul == 2:
    outputDir = Path.cwd() / "wot_multiplicity_2"
elif mul == 3:
    outputDir = Path.cwd() / "wot_multiplicity_3"
elif mul == 4:
    outputDir = Path.cwd() / "wot_multiplicity_4"
elif mul == 5:
    outputDir = Path.cwd() / "wot_multiplicity_5"
else:
    print("multiplicity error")
'''

# Modify position stddev
'''
stddev_p = float(input("value of Y-Z stddev"))
if stddev_p == 0:
    outputDir = Path.cwd() / "wot_pos_stddev_0"
elif stddev_p == 1:
    outputDir = Path.cwd() / "wot_pos_stddev_1"
elif stddev_p == 2:
    outputDir = Path.cwd() / "wot_pos_stddev_2"
elif stddev_p == 3:
    outputDir = Path.cwd() / "wot_pos_stddev_3"
elif stddev_p == 4:
    outputDir = Path.cwd() / "wot_pos_stddev_4"
elif stddev_p == 5:
    outputDir = Path.cwd() / "wot_pos_stddev_4"
else:
    print("pos_stddev error")
'''

# Modify time stddev
'''
stddev_t = float(input("value of time stddev"))
if stddev_t == 0:
    outputDir = Path.cwd() / "wot_time_stddev_0"
elif stddev_t == 1:
    outputDir = Path.cwd() / "wot_time_stddev_1"
elif stddev_t == 2:
    outputDir = Path.cwd() / "wot_time_stddev_2"
elif stddev_t == 3:
    outputDir = Path.cwd() / "wot_time_stddev_3"
else:
    print("time_stddev error")
'''

outputDir = Path.cwd() / "result-without-time"

field = acts.ConstantBField(acts.Vector3(0, 0, 0)) # u.T
if not outputDir.exists():
    outputDir.mkdir()

rnd = acts.examples.RandomNumbers(seed=42)
s = acts.examples.Sequencer(events=100, numThreads=1, outputDir=str(outputDir))

addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0113, 0.0113, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.649 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False),
    multiplicity=20,
    #multiplicity=mul, # To modify the value of multiplicity
    rnd=rnd,
    vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0.*u.mm, 3.04*u.mm, 3.04*u.mm, 1.52*u.ns)),
    #vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0*u.mm, stddev_p*u.mm, stddev_p*u.mm, 1.52*u.ns)),
    #vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0*u.mm, 5.*u.mm, 5.*u.um, stddev_t*u.ns)),
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
    logLevel=acts.logging.VERBOSE,
    #volumeMappings = ["Layer #0 Phys"],
    g4DetectorConstructionFactory=TelescopeG4DetectorConstructionFactory(teleG4Config),
    preSelectParticles=ParticleSelectorConfig(
            rho=(0.0, 300 * u.mm),
            absZ=(-200.0, 200.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(1 * u.GeV, None),
            removeNeutral=True,
    ),
    outputDirRoot=outputDir,
    #logLevel=acts.logging.VERBOSE,
    killVolume=trackingGeometry.worldVolume,
    killAfterTime=1000 * u.ns,
)


addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=Path("../Examples/Algorithms/Digitization/share/default-digi-config-telescope.json"),
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
        chi2CutOff=15.0,
        numMeasurementsCutOff=10,
    ),
    outputDirRoot=outputDir,
    #logLevel=acts.logging.DEBUG,
)

s.run()
