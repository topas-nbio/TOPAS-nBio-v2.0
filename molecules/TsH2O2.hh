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
#ifndef TsH2O2_h
#define TsH2O2_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4MoleculeDefinition.hh"

// ######################################################################
// ###                         Peroxyde                               ###
// ######################################################################

class TsH2O2 : public G4MoleculeDefinition
{
private:
      static /*G4ThreadLocal*/ TsH2O2* theInstance;
      TsH2O2() {}
      virtual ~TsH2O2() {}
  
  public:
      static TsH2O2* Definition();
  };
  
#endif

