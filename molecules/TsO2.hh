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
#ifndef TsO2_h
#define TsO2_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4MoleculeDefinition.hh"

// ######################################################################
// ###                         Peroxyde                               ###
// ######################################################################

class TsO2 : public G4MoleculeDefinition
{
private:
      static /*G4ThreadLocal*/ TsO2* theInstance;
      TsO2() {}
      virtual ~TsO2() {}
  
  public:
      static TsO2* Definition();
  };
  
#endif

