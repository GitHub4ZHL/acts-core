source /opt/ROOT/root_v6_26_06/INSTALL/bin/thisroot.sh 

source /home/xiaocong/Software/Oscar/geant4-v11.0.2/install/bin/geant4.sh
#source /home/xiaocong/Software/Oscar/geant4-v11.0.2.new.xerces/install/bin/geant4.sh

#geant
export LD_LIBRARY_PATH=/home/xiaocong/Software/Oscar/geant4-v11.0.2/install/lib64:$LD_LIBRARY_PATH
export G4ENSDFSTATEDATA=/home/xiaocong/Software/Oscar/geant4-v11.0.2/install/share/Geant4-11.0.2/data/G4ENSDFSTATE2.3
export G4LEVELGAMMADATA=/home/xiaocong/Software/Oscar/geant4-v11.0.2/install/share/Geant4-11.0.2/data/PhotonEvaporation5.7
export G4LEDATA=/home/xiaocong/Software/Oscar/geant4-v11.0.2/install/share/Geant4-11.0.2/data/G4EMLOW8.0
export G4SAIDXSDATA=/home/xiaocong/Software/Oscar/geant4-v11.0.2/install/share/Geant4-11.0.2/data/G4SAIDDATA2.0
#export G4NEUTRONXSDATA=/home/xiaocong/Software/Oscar/geant4-v11.0.2/install/share/Geant4-11.0.2/data/G4NEUTRONXS1.4


source /home/xiaocong/Software/Oscar/DD4hep-01-22/install/bin/thisdd4hep.sh

source  /home/xiaocong/Software/stcfActs/acts1/install/bin/this_acts.sh

source /home/xiaocong/Software/stcfActs/acts1/build/python/setup.sh

