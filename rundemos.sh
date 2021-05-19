#!/bin/bash
cd examples/geometry
name=(dna/CharltonDNA.txt dna/LinearDNA.txt dna/CircularPlasmid.txt\
      cells/EllipsoidCell.txt cells/FibroblastCell.txt cells/SphericalCell.txt\
      other/generateRandomCylinders.txt other/readBackRandomCylindres.txt)
for f in `seq 0 7`
  do
    topas ${name[$f]} > out; grep -c Execution out | awk '{if ($1 == 1)  print  "--success-- '${name[$f]}'";   else  print "--failed -- '${name[$f]}'"}'
done
rm -f *.xyz *.phsp *.header *.csv out
cd ../../

cd examples/processes
name=(ActiveChemistryDefault.txt ActiveChemistryExtended.txt ActiveChemistryRevised.txt\
      ActiveCustomizablePhysics.txt G4DNAModelPerRegion.txt RemoveChemicalSpeciesInVolume.txt)
for f in `seq 0 5`
  do
    topas ${name[$f]} > out; grep -c Execution out | awk '{if ($1 == 1)  print  "--success-- '${name[$f]}'";   else  print "--failed -- '${name[$f]}'"}'
done
rm -f out
cd ../../

cd examples/scorers
name=(DBSCAN.txt PDB4DNA.txt particleTuple.txt IonizationDetailInRandomCylinders.txt\
      simpleSSBandDSBInPlasmid.txt simpleSSBandDSBInPlasmidWithDBSCAN.txt NumberOfChemicalSpeciesAtTime.txt\
      GvalueG4DNADefault.txt GvalueRevisedPhysicsChemistry.txt)
for f in `seq 0 8`
  do
    topas ${name[$f]} > out; grep -c Execution out | awk '{if ($1 == 1)  print  "--success-- '${name[$f]}'";   else  print "--failed -- '${name[$f]}'"}'
done
rm -f *.root *.phsp *.header *csv out
cd ../../

