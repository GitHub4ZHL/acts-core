#!/bin/bash

detType="upgrade"
fieldType="nonUniformField"
#fieldType="uniformField"
tag=CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.simNonUniField.recUniField
#tag=CKF.smeared.MdcRecHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.simUniField.recUniField
#tag=CKF.smeared.MdcMcHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.simNonUniField.recUniField
#tag=CKF.smeared.MdcMcHits.PreditedDriftSign.chi2Cut30.maxPropSteps330.betheLoss.simUniField.recUniField


ckfPropSteps="330"
configFilePath=/home/xiaocong/Software/bes3Acts/acts/build/bin

momentums=("100" "125" "150" "175" "200" "250" "300" "350" "400" "450" "600" "800" "1000" "1200" "1400" "1600" "1800")
#momentums=("125")
absCosThetas=("0" "0.5" "0.8")
#absCosThetas=("0")
particles=("mu-" "pi-")
#particles=("pi-")
momentumTypes=("p" "pt")
#momentumTypes=("pt")

for ((m=0; m<${#momentumTypes[@]};++m)); do
  momentumType=${momentumTypes[m]}
  
  for ((i=0; i<${#momentums[@]};++i)); do
  for ((j=0; j<${#absCosThetas[@]};++j)); do
  for ((k=0; k<${#particles[@]};++k)); do
  
    momentum=${momentums[i]}
    absCosTheta=${absCosThetas[j]}
    particle=${particles[k]}
  
#    if [ "${momentumType}" = "pt" ]
#    then
#       absSinTheta=`echo "scale=3; sqrt(1-${absCosTheta}^2)" | bc`
#       realMomentum=`printf '%.2f\n' "$(echo "scale=2; ${momentum}/${absSinTheta}/1000" | bc)"`
#    else
#       realMomentum=`printf '%.2f\n' "$(echo "scale=2; ${momentum}/1000" | bc)"`
#    fi

    runBesUprade="" 
    if [ "${detType}" = "upgrade" ] 
    then
      runBesUprade="--run-bes-upgrade"
    fi

    
    #if [[ $particle == "proton" ]];then
    #   if [[ $momentum == "50" || $momentum == "75" || $momentum == "100" || $momentum == "125" || $momentum == "150" ]]; then 
    #     continue
    #   fi
    #fi
    
    #ckfchi2max=30
    #seedSigmaScattering="200"
    #ckfPropTolerance="0.0001"
    #ckfPropMass="139.57018"
    #if [[ $momentum == "50" || $momentum == "75" || $momentum == "100" ]];then
    #  seedSigmaScattering=2
    #fi
    #
    #if [[ $momentum == "50" ]]; then
    #  ckfPropTolerance="0.01"
    #elif [[ $momentum == "75" || $momentum == "100" || $momentum == "125" ]];then
    #  ckfPropTolerance="0.0008"
    #fi
    #
    #if [[ $particle == "mu-" ]];then
    #  ckfPropMass="105.6583755"
    #elif [[ $particle == "proton" ]];then
    #  ckfPropMass="938.272088"
    #fi
    
            #--ckf-reco-mdc \

    inputDir=/home/xiaocong/Software/bes3Acts/acts/RunSpace/samples/${detType}/${fieldType}/ana_momentum_angle_scan/${momentumType}/${particle}
    inputFile=single_${particle}_absCosThetaDeg_${absCosTheta}_${absCosTheta}_momentumMev_${momentum}_${momentum}.root
    outputDir=/home/xiaocong/Software/bes3Acts/acts/RunSpace/perf/${detType}/ana_momentum_angle_scan/${tag}/${momentumType}/${particle}/absCosThetaDeg_${absCosTheta}_${absCosTheta}_momentumMev_${momentum}_${momentum}/
    
    mkdir -p ${outputDir}
    
    # -j 1 \
 
    run="${configFilePath}/ActsExampleCKFTracksBES \
            --ckf-prop-steps 330 \
            ${runBesUprade} \
            --input-dir=${inputDir} \
            --input-files=${inputFile} \
            --output-dir=${outputDir} \
            --ckf-truth-smeared-seeds \
            --ckf-selection-nmax 1 \
            --ckf-selection-chi2max 30.0 \
            --ckf-initial-variance-inflation=100:100:100:100:100:1  \
            --mat-input-type file  \
            --mat-input-file=${configFilePath}/mat-bes_${detType}_alt5_manual_MDI2_MDCLayer100.json  \
            --geo-tgeo-filename=${configFilePath}/bes_pixel_mdc_Pip.root  \
            --geo-tgeo-jsonconfig=${configFilePath}/tgeo_bes_${detType}_config.json  \
            --bf-constant-tesla=0:0:-1"
  
    echo $run
    
    eval $run
    #echo ======================================== 
  
  done
  done
  done

done


