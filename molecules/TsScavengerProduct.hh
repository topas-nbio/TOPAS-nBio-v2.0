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
#ifndef TsScavengerProduct_h
#define TsScavengerProduct_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4MoleculeDefinition.hh"

// ######################################################################
// ###                   ScavengerProduct (Generic)                   ###
// ######################################################################

class TsScavengerProduct : public G4MoleculeDefinition
{
private:
      static /*G4ThreadLocal*/ TsScavengerProduct* theInstance;
      TsScavengerProduct() {}
      virtual ~TsScavengerProduct() {}
  
  public:
      static TsScavengerProduct* Definition();
  };
  
#endif

