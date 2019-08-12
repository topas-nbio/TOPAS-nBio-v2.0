// Scorer for TsNtupleForNeuron
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


#include "TsNtupleForNeuron.hh"

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

TsNtupleForNeuron::TsNtupleForNeuron(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //SetSurfaceScorer();
    
    fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
    fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
    fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
    fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
    fNtuple->RegisterColumnF(&fEnergyDep, "Energy Deposited", "MeV");
    fNtuple->RegisterColumnS(&fVolume, "Volume Name");
    fNtuple->RegisterColumnI(&fEventID, "Event ID");
    fNtuple->RegisterColumnI(&fTrackID, "Track ID");
    
}


TsNtupleForNeuron::~TsNtupleForNeuron() {;}


G4bool TsNtupleForNeuron::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    ResolveSolid(aStep);
    
    fEnergyDep      = aStep->GetTotalEnergyDeposit();
    
    G4StepPoint* theStepPoint=aStep->GetTrack()->GetStep()->GetPreStepPoint();
    G4TouchableHandle theTouchable = theStepPoint->GetTouchableHandle();
    G4String volumeName = theTouchable->GetVolume()->GetName();
    
    G4int VolNum = 0;
    
    if (volumeName == "Neuron/Soma") {VolNum = 1;}
    if (volumeName == "Neuron/SomaUnion") {VolNum = 2;}
    if (volumeName == "Neuron/Axon") {VolNum = 3;}
    if (volumeName == "Neuron/BDend") {VolNum = 4;}
    if (volumeName == "Neuron/ADend") {VolNum = 5;}
    
    
    
    if ((fEnergyDep > 0) && (VolNum != 0)){
        
        G4cout << volumeName << G4endl;
        
        //Get position
        G4ThreeVector pos = theStepPoint->GetPosition();
        fPosX = pos.x();
        fPosY = pos.y();
        fPosZ = pos.z();

        //Get kinetic energy
        fEnergy = theStepPoint->GetKineticEnergy();
        
        //Get energy deposition
        fEnergyDep = aStep->GetTrack()->GetStep()->GetTotalEnergyDeposit();
        
        //Get particle type
        //fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
        
        //Get Volume name
        fVolume = volumeName;
        
        // Get IDs
        fTrackID = aStep->GetTrack()->GetTrackID();
        fEventID = GetEventID();
        
        fNtuple->Fill();
        return true;
    }
    return false;   
}
