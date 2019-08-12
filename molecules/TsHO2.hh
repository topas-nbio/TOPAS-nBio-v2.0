//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//
#ifndef TsHO2_h
#define TsHO2_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4MoleculeDefinition.hh"

// ######################################################################
// ###                         Peroxyde                               ###
// ######################################################################

class TsHO2 : public G4MoleculeDefinition
{
private:
      static /*G4ThreadLocal*/ TsHO2* theInstance;
      TsHO2() {}
      virtual ~TsHO2() {}
  
  public:
      static TsHO2* Definition();
  };
  
#endif

