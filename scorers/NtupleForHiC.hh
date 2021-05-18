/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/


#ifndef NtupleForHiC_hh
#define NtupleForHiC_hh

#include "TsVNtupleScorer.hh"

#include "G4VProcess.hh"


#include "G4PSDirectionFlag.hh"

#include "TsVGeometryComponent.hh"
//#include "TsGeometryManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

class NtupleForHiC : public TsVNtupleScorer
{
public:
    NtupleForHiC(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM,
                    TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName,
                    G4bool isSubScorer);

    virtual ~NtupleForHiC();

    G4bool ProcessHits(G4Step*,G4TouchableHistory*);

    void UserHookForEndOfEvent();

    void UserHookForEndOfRun();

protected:
    // Output variables
    G4double fPosX=-1000.0;
    G4double fPosY=-1000.0;
    G4double fPosZ=-1000.0;

    G4int fCause=-1;
    G4int fEventID=-1;
    G4int fRunID=-1;
    G4int fStrand=1;
    G4int fChromosome=-1;
    G4int fBeadID=1;
    
    //Scoring Parameters
    G4double SensitiveFraction=0.141328949;
    G4double EMin=5.0;
    G4double EMax=37.5;
    
    



private:

};
#endif
