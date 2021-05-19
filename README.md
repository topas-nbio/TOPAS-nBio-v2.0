# TOPAS-nBio
This is the TOPAS-nBio extension repository, a Monte Carlo simulation framework for (sub-) cellular radiobiology.

TOPAS-nBio is described here: https://topas-nbio.readthedocs.io/. 
This page includes a class documentation and the license.

TOPAS-nBio is an extension of TOPAS (TOol for PArticle Simulations), which can be obtained from www.topasmc.org. The TOPAS documentation can be found at https://topas.readthedocs.io/. 

The TOPAS-nBio package was described in Schuemann et al., Radiation Research, 2019, 191(2), p.125. This reference should be cited for all work using the TOPAS-nBio package.

We encourage contribution of user-developed extensions that fit within the TOPAS-nBio scheme. If you would like to contribute code, please contact the developers.


## General information

0) Pre-requisites:

   TOPAS installed with recommended `OS system, c++ and cmake versions`, see 
   topas https://topas.readthedocs.io/en/latest/getting-started/install.htm

1) We recommend having a global directory for extensions named topas_extension and move in TOPAS-nBio there 

   Linux:
        mkdir ~/topas_extensions

   Mac: 
        mkdir /Applications/topas_extensions

2) Unzip TOPAS-nBio directory in topas_extensions and navigate to the topas directory

   Linux:
        cd ~/topas
   Mac:
        cd /Applications/topas

3) Unzip the Geant4Headers.zip

   Linux:
        unzip -e Geant4Headers.zip
   Mac:
        unzip -e Geant4Headers.zip
        
4) Build the extensions

   Linux:
        cmake ./ -DTOPAS_EXTENSIONS_DIR=~/topas_extensions/TOPAS-nBio
        make -j4
   Mac:
        cmake ./ -DTOPAS_EXTENSIONS_DIR=/Applications/topas_extensions/TOPAS-nBio
        make -j4
 
5) Run the demos. For some demos, a pause before quit is enabled, then, write `exit` at the terminal prompt.

   Linux:
        source rundemos.csh

   Mac:
        source rundemos.csh

