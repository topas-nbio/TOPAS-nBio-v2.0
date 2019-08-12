# General information
Run the demos with:
   source rundemos.csh

In some examples, a pause before quit is enabled. Then, just write `exit` at the terminal.

## Geometry examples
### DNA models
#### CharltonDNA.txt
 Simplified model of DNA. It is represented as a  cylinder with base pairs assambled 
 by two cylindrical segments of 0.34 nm length placed around a cylinder of 1 nm 
 diameter that represents the base. The DNA diameter is 2.3 nm.
 
 #### DNA.txt
 DNA model representing a whole nuclear DNA. The example is taken from 
 the Geant4 DNA example wholeNuclearDNA and implemented into TOPAS-nBio.

#### Plasmid.txt
 A circular plasmid of 2000 base pairs. The DNA is cylidrical with 2.3 nm diameter. 
 The bases are 0.34 nm length.


### Cell models 
#### EllipsoidCell.txt 
#### FibroblastCell.txt
#### SphericalCell.txt

### Sampling geometries
#### NikjooCylinders.txt

#### generateRandomCylinders.txt
Shows how to generate random cylinders enclosed by a user-defined TsSphere,
TsCylinder or TsBox. In an iterative process, the algorithm generates randomly 
a cylinder and checks with existing cylinders for overlaps. It rejects the new cylinder
if overlaps exists. The positions in nm and rotations in deg (6 columns) are stored 
in an ascii file. Accordingly, many warning messages about overlapping geometries are expected.

#### readBackRandomCylindres.txt
Shows how to read back the output ascii file from example 
`generateRandomCylinders.txt`. Users still need to set the shape and 
dimensions of the envelope geometry.

## Processes examples
### Physics
#### ActiveCustomizablePhysics.txt
Shows how to select the elastic/inelastic scattering models for e- and the elastic scattering
model for protons.

#### G4DNAModelPerRegion.txt
Shows how to set the Geant4DNA track-structure physics per region and condensed history elsewhere.

### Chemistry
#### ActiveChemistryDefault.txt
Shows how to activate the step-by-step transport of chemical species

#### ActiveChemistryRevised.txt
Shows how to set the chemistry parameters different than those provided by the Geant4-DNA default 
chemistry list.

#### ActiveChemistryExtended.txt
Shows how to customize the chemistry parameters using an extended list of reactions.

#### RemoveChemicalSpeciesInVolume.txt
Shows how to terminate the transport of chemical species in a user-defined volume with specific
material. 

## Scorer examples
#### DBSCAN.txt
Shows how to calcualte single and double strand breaks in a homogeneous geometry using the DBSCAN algorithm

#### DBSCAN_VRT.txt
Shows how to calcualte single and double strand breaks in a homogeneous geometry using the DBSCAN algorithm
with a variance reduction technique (particle splitting).

#### PDB4DNA.txt
Shows how to calculate single and double strand breaks using an overlaid geometry from a PDB file. The scoring
method is inherited from the PDB4DNA example of Geant4-DNA.

#### PDB4DNA_VRT.txt
Shows how to calculate single and double strand breaks using an overlaid geometry from a PDB file and 
variance reduction (particle splitting). The scoring method is inherited from the PDB4DNA example of Geant4-DNA.

#### NumberOfChemicalSpeciesAtTime.txt
Shows how to calculate the number of chemical species as a function of time. Time is discretized in a log-scale way.

#### GvalueG4DNADefault.txt
Shows how to calcualte the G-value (number of chemical species per 100 eV of energy deposit) 

#### GvalueRevisedPhysicsChemistry.txt
Shows how to calcualte the G-value (number of chemical species per 100 eV of energy deposit) using a revised set
of chemistry parameters.

#### particleTuple.txt
Shows how to get the physical chemical track information. For the chemical track, the information is retrived only at
TimeCut.

#### IonizationDetailInRandomCylinders.txt
Shows how to get the physical track information including only essential information.

#### simpleSSBandDSBInPlasmid.txt
Shows how to score SSB and DSB in a plasmid geometry using the simple definition of DSB (2 SSB separated within 10 bps).
The classification of the DSB is performed via a PDB4DNA-like algorithm, example `PDB4DNA.txt` 

#### simpleSSBandDSBInPlasmidWithDBSCAN.txt
Shows how to score SSB and DSB in a plasmid geometry using the simple definition of DSB (2 SSB separated within 10 bps).
The classification of the DSB is performed via the DBSCAN algorithm.


