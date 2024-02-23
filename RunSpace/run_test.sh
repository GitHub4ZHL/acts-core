
tgeoRootFile=STCF_tracker.root
tgeoConfigFile=tgeo_stcf_config.json
matInputFile=material_map_stcf.json
smearConfigFile=STCF-digi-smearing-config.json


# Check geometry building
../build/bin/ActsExampleGeometryTGeo \
	--geo-tgeo-filename=${tgeoRootFile} \
	--geo-tgeo-jsonconfig=${tgeoConfigFile} \
	--output-obj --output-dir=obj  \
	--geo-volume-loglevel=0

# Propagation
../build/bin/ActsExamplePropagationTGeo \
	-n 1000  \
	--prop-ntests=1 \
	--prop-eta-range=-1.73:1.73 --prop-pt-range=0.1:2 \
	--geo-tgeo-filename=${tgeoRootFile}  \
	--geo-tgeo-jsonconfig=${tgeoConfigFile} \
	--output-dir=data/propagation \
	--output-root \
	--bf-constant-tesla=0:0:1


# Fatras simulation
../build/bin/ActsExampleFatrasTGeo  \
	--events=1000 \
	--gen-eta=-1.73:1.73 \
	--gen-mom-transverse \
	--gen-mom-gev=0.1:2 \
	--geo-tgeo-filename=${tgeoRootFile}  \
	--geo-tgeo-jsonconfig=${tgeoConfigFile} \
        --mat-input-type file  \
	--mat-input-file=${matInputFile}  \
	--output-dir=data/fastsim  \
	--output-root \
        --output-csv  \
	--bf-constant-tesla=0:0:1 


# Truth fitting using Fatras fast sim
../build/bin/ActsExampleTruthTracksTGeo  \
        --input-dir=data/fastsim \
        --output-dir=data/reco_fastsim \
        --geo-tgeo-filename=${tgeoRootFile}  \
        --geo-tgeo-jsonconfig=${tgeoConfigFile} \
        --mat-input-type file  \
        --mat-input-file=${matInputFile}  \
        --digi-config-file=${smearConfigFile} \
        --bf-constant-tesla=0:0:1

# Truth fitting using Oscar full sim output, e.g. fullsim_pipijpsi.root 
../build/bin/ActsExampleTruthTracksSTCF  \
	--input-dir=data/fullsim \
        --input-files=fullsim_pipijpsi.root \
	--output-dir=data/reco_fullsim \
	--geo-tgeo-filename=${tgeoRootFile}  \
	--geo-tgeo-jsonconfig=${tgeoConfigFile} \
        --mat-input-type file  \
	--mat-input-file=${matInputFile}  \
        --bf-constant-tesla=0:0:1

