# TOPAS-CellModels
1) Description:

Cell Models for TOPAS/Geant4 and the inclusion of nano particles in particle scattering simulations.

The C++ classes in this repository extend the functionality of the TOPAS (http://www.topasmc.org/) Monte-Carlo program, which is itself a wrapper of the Geant4 MCS Toolkit (http://geant4.org).


2) Installation:

Installation:

Navigate to the TOPAS extension directory:

  cd ~/topas_extensions/

Clone or download the sourcecode into your TOPAS extension directory:
 
  git clone https://github.com/BAMresearch/TOPAS-CellModels.git
 
Change to your topas directory:

  cd ~/Topas/

Install it:

  cmake ./ -DTOPAS_EXTENSIONS_DIR=~/topas_extensions/TOPAS-CellModels &&  make -j4


3) Description:

A simple spherical cell with nanoparticles can be generated in a fast manner.
The user has the option to include nanoparticles and different organelles in the cell, e.g. a nucleus, mitochondria, a cell membrane.
Details can be found in https://doi.org/10.1038/s41598-021-85964-2

4) Usage:

Examples can be found in the  "examples/" directory.


 
5) Bugs:

Please report bugs to hahn@physik.fu-berlin.de or on https://github.com/BAMresearch/TOPAS-CellModels


6) Literature:

If you use this extension please cite the following literature:

Hahn, M.B., Zutta Villate, J.M. "Combined cell and nanoparticle models for TOPAS to study radiation dose enhancement in cell organelles."
Sci Rep 11, 6721 (2021). 

https://doi.org/10.1038/s41598-021-85964-2


7) Etc:

Tags: topas topasmc topasmcs topas-mc topas-mcs topas-nbio mcs 
