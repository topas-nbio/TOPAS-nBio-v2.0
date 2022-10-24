// Scorer for TsNtupleForCulture
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
//Tuple scorer for scoring events deposited in the Culture organelles

#include "TsNtupleForCulture.hh"

#include "G4SystemOfUnits.hh"

#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

TsNtupleForCulture::TsNtupleForCulture(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //SetScorer();
    fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
    fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
    fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
    fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
    fNtuple->RegisterColumnF(&fEnergyDep, "Energy Deposited", "MeV");
    fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
    fNtuple->RegisterColumnI(&fTrackID, "Track ID");
    fNtuple->RegisterColumnI(&fRunID, "Run ID");
    fNtuple->RegisterColumnI(&fEventID, "Event ID");
    fNtuple->RegisterColumnS(&fVolName, "Volume Name");
    fNtuple->RegisterColumnI(&fNucleusID, "Nucleus ID");

}


TsNtupleForCulture::~TsNtupleForCulture() {;}


G4bool TsNtupleForCulture::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    ResolveSolid(aStep);
    
    G4double flagEnergyDep      = aStep->GetTotalEnergyDeposit();
    
    G4StepPoint* theStepPoint=aStep->GetPreStepPoint();
    
    //Find volume name
    G4Track* aTrack = aStep->GetTrack();
    G4String volumeName = aTrack->GetVolume()->GetName();
    
    G4TouchableHistory* theTouchable =  (G4TouchableHistory*)(aTrack->GetTouchable());
    
    G4int motherCopyNo = theTouchable->GetReplicaNumber(1);
    
    //Get position
    G4ThreeVector pos = theStepPoint->GetPosition();
    
    //Score events that deposit energy in the Culture organelles:
    if ((flagEnergyDep > 0) && (volumeName != "MyCulture")){
        
        //Get position
        fPosX = pos.x();
        fPosY = pos.y();
        fPosZ = pos.z();
        
        //Get particle Energy
        fEnergy = theStepPoint->GetKineticEnergy();

        //Get Edep
        fEnergyDep = flagEnergyDep;
        
        //Get particle type
        fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

        // Get IDs
        fTrackID = aStep->GetTrack()->GetTrackID();
        fRunID   = GetRunID();
        fEventID = GetEventID();
        
        //Get volume Name
        fVolName = volumeName;
        
        //Get the cell ID of the nucleus
        fNucleusID = motherCopyNo;

        fNtuple->Fill();

        return true;
    }
    return false;   
}
