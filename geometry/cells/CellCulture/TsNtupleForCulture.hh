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

#ifndef TsNtupleForCulture_hh
#define TsNtupleForCulture_hh

#include "TsVNtupleScorer.hh"

#include "G4VProcess.hh"

class TsNtupleForCulture : public TsVNtupleScorer
{
public:
    TsNtupleForCulture(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsNtupleForCulture();

    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    

    
protected:
    // Output variables
    G4float fPosX;
    G4float fPosY;
    G4float fPosZ;
    G4float fEnergy;
    G4float fEnergyDep;
    G4int fParticleType;
    G4int fProcess;
    G4int fTrackID;
    G4int fRunID;
    G4int fEventID;
    G4String fVolName;
    G4int fNucleusID;

};
#endif
