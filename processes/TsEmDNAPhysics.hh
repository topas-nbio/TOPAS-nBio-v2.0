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

#ifndef TsEmDNAPhysics_h
#define TsEmDNAPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class TsParameterManager;

class TsEmDNAPhysics : public G4VPhysicsConstructor
{
public:
    
    explicit TsEmDNAPhysics(G4int ver=1, const G4String& name="");
    TsEmDNAPhysics(TsParameterManager* pM);
    
    virtual ~TsEmDNAPhysics();
    
    virtual void ConstructParticle();
    virtual void ConstructProcess();
    
private:
    G4String fName;

    TsParameterManager* fPm;
    G4int  verbose;
};

#endif