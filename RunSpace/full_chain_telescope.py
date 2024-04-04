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
rnd = acts.examples.RandomNumbers(seed=42)

teleG4Config=TelescopeDetector.Config();
teleG4Config.bounds=[14.08, 28.04]
teleG4Config.positions=[30, 60, 90, 105, 120, 150, 180]
teleG4Config.stereos=[0, 0, 0, 0, 0, 0, 0]
teleG4Config.thickness = [80*u.um, 80*u.um, 80*u.um, 1*u.um, 80*u.um, 80*u.um, 80*u.um]
teleG4Config.binValue=0

detector, trackingGeometry, decorators = acts.examples.TelescopeDetector.create(
    bounds=[14.08, 28.04],
    positions=[30, 60, 90, 105, 120, 150, 180],
    stereos=[0, 0, 0, 0, 0, 0, 0],
    thickness=[80*u.um, 80*u.um, 80*u.um, 1*u.um, 80*u.um, 80*u.um, 80*u.um],
    binValue=0,
)

max_multiplicity = int(input("Loop multiplicity from 1 to "))

# Begin multiplicity loop
for multiplicity in range(1, max_multiplicity+1):
    # Without time
    outputDir = Path.cwd() / f"wot_multiplicity_{multiplicity}"
    if not outputDir.exists():
        outputDir.mkdir()
    s = acts.examples.Sequencer(events=1000, numThreads=1, outputDir=str(outputDir))

    addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0113, 0.0113, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.649 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False),
    multiplicity=multiplicity,
    rnd=rnd,
    vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0.*u.mm, 3.04*u.mm, 3.04*u.mm, 1.52*u.ns)),
    #logLevel=acts.logging.VERBOSE,
    )

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
        chi2CutOff=50,
        numMeasurementsCutOff=1,
    ),
    outputDirRoot=outputDir,
    #logLevel=acts.logging.VERBOSE,
    )

    #''' Temporarily unable to run successfully
    addGSFTracks(
        s,
        trackingGeometry,
        field,
        outputDirRoot=outputDir,
        #logLevel=acts.logging.VERBOSE,
    )
    #'''
    continue
    s.run()
    print(f"Finished telescope simulation and reconstruction without time of multiplicity = {multiplicity}")
# End multiplicity loop

# Begin multiplicity loop
for multiplicity in range(1, max_multiplicity+1):
    # With time
    outputDir = Path.cwd() / f"wt_multiplicity_{multiplicity}"
    if not outputDir.exists():
        outputDir.mkdir()
    s = acts.examples.Sequencer(events=1000, numThreads=1, outputDir=str(outputDir))

    addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0113, 0.0113, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.649 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False),
    multiplicity=multiplicity,
    rnd=rnd,
    vtxGen=acts.examples.GaussianVertexGenerator(mean=acts.Vector4(0, 0, 0, 0), stddev=acts.Vector4(0.*u.mm, 3.04*u.mm, 3.04*u.mm, 1.52*u.ns)),
    #logLevel=acts.logging.VERBOSE,    
    )

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
    digiConfigFile=Path("../Examples/Algorithms/Digitization/share/default-digi-config-telescope-time.json"),
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
        chi2CutOff=50,
        numMeasurementsCutOff=1,
    ),
    outputDirRoot=outputDir,
    #logLevel=acts.logging.VERBOSE,
    )

    ''' Temporarily unable to run successfully
    addGSFTracks(
        s,
        trackingGeometry,
        field,
        outputDirRoot=outputDir,
        #logLevel=acts.logging.VERBOSE,
    )
    '''

    s.run()
    print(f"Finished telescope simulation and reconstruction with time of multiplicity = {multiplicity}")

# End multiplicity loop
print("All OK!")
