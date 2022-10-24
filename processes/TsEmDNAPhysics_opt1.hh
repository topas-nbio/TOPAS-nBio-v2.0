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

#ifndef TsEmDNAPhysics_opt1_h
#define TsEmDNAPhysics_opt1_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class TsParameterManager;

class TsEmDNAPhysics_opt1 : public G4VPhysicsConstructor
{
public:
    
    explicit TsEmDNAPhysics_opt1(G4int ver=1, const G4String& name="");
    TsEmDNAPhysics_opt1(TsParameterManager* pM);
    
    virtual ~TsEmDNAPhysics_opt1();
    
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
private:
    TsParameterManager* fPm;
    G4int  verbose;
};

#endif
