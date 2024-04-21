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

max_chi2 = int(input("Loop chi2 from 10 to "))

# Begin multiplicity loop
for chi2 in range(10, max_chi2+5, 5):
    # Without time
    outputDir = Path.cwd() / f"wt_1_chi2_{chi2}"
    if not outputDir.exists():
        outputDir.mkdir()
    s = acts.examples.Sequencer(events=10000, numThreads=1, outputDir=str(outputDir))

    addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0113, 0.0113, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.649 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False),
    multiplicity=1,
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
        chi2CutOff=chi2,
        numMeasurementsCutOff=1,
    ),
    outputDirRoot=outputDir,
    #logLevel=acts.logging.VERBOSE,
    )

    s.run()
# End multiplicity loop

# Begin multiplicity loop
for chi2 in range(10, max_chi2+5, 5):
    # With time
    outputDir = Path.cwd() / f"wt_15_chi2_{chi2}"
    if not outputDir.exists():
        outputDir.mkdir()
    s = acts.examples.Sequencer(events=10000, numThreads=1, outputDir=str(outputDir))

    addParticleGun(
    s,
    MomentumConfig(4 * u.GeV, 4 * u.GeV, transverse=True),
    EtaConfig(-0.0113, 0.0113, uniform=True),
    PhiConfig(-0.0 * u.degree, 0.649 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eElectron, randomizeCharge=False),
    multiplicity=15,
    rnd=rnd,
    vtxGen=acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0.*u.mm, 3.04*u.mm, 3.04*u.mm, 1.52*u.ns)),
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
        #http://courses.atlas.illinois.edu/spring2016/STAT/STAT200/pchisq.html 
        #p-Value here means the probability of rejecting a good hit 
        #13.82 and 16.27 correpsonds to p-Value of 0.001 for chisq with 2 and 3 degree of freedom, respectively
        #18.42 and 21.11 correpsonds to p-Value of 0.0001 for chisq with 2 and 3 degree of freedom, respectively
        #19.81 and 22.55 correpsonds to p-Value of 0.00005 for chisq with 2 and 3 degree of freedom, respectively
        chi2CutOff=chi2, 
        numMeasurementsCutOff=1,
    ),
    outputDirRoot=outputDir,
    #logLevel=acts.logging.VERBOSE,
    )

    s.run()

# End multiplicity loop
print("All OK!")
