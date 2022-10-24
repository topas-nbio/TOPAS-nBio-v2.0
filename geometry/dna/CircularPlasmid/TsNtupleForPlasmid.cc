// Scorer for TsNtupleForPlasmid
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

//An example ntuple scorer for the DNA circular plasmid geometry.
//File may be edited to score different quantities.

#include "TsNtupleForPlasmid.hh"

#include "G4PSDirectionFlag.hh"

#include "TsVGeometryComponent.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

TsNtupleForPlasmid::TsNtupleForPlasmid(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //Define ntuple columns

    fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
    fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
    fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
    fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
    fNtuple->RegisterColumnI(&fVolume, "Volume 1=BP, 2=Backbone1, 3=Backbone2");
    fNtuple->RegisterColumnI(&fRepNum, "BasePair Replication Number");
    fNtuple->RegisterColumnF(&fEnergyDep, "Energy Deposited", "MeV");
    fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
    fNtuple->RegisterColumnI(&fProcess, "Process Type");
    fNtuple->RegisterColumnI(&fEventID, "Event ID");
    fNtuple->RegisterColumnI(&fTrackID, "Track ID");

}


TsNtupleForPlasmid::~TsNtupleForPlasmid() {;}


G4bool TsNtupleForPlasmid::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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
    
    G4int VolNum;

    if (volumeName == "MyDNA/BasePair") VolNum = 1;
    if (volumeName == "MyDNA/Backbone1") VolNum = 2;
    if (volumeName == "MyDNA/Backbone2") VolNum = 3;
    

    if ((fEnergyDep > 0) && ((VolNum == 1) || (VolNum == 2) || (VolNum == 3))) {
    
        //Get Volume name
        fVolume = VolNum;
        
        //The basepair number
        fRepNum = theTouchable->GetReplicaNumber(1);
        G4cout << theTouchable->GetReplicaNumber(1) << G4endl;
    
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
        fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();


        G4String processName = aStep->GetTrack()->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

        //Define process type:
        if (processName=="msc") flagProcess =1;
        if (processName=="e-_G4DNAElastic")		flagProcess =2;
        if (processName=="e-_G4DNAExcitation")		flagProcess =3;
        if (processName=="e-_G4DNAIonisation")		flagProcess =4;
        if (processName=="e-_G4DNAAttachment")		flagProcess =5;
        if (processName=="e-_G4DNAVibExcitation")	flagProcess =6;
        if (processName=="eCapture")			flagProcess =7;

        if (processName=="proton_G4DNAExcitation")	flagProcess =8;
        if (processName=="proton_G4DNAIonisation")	flagProcess =9;
        if (processName=="proton_G4DNAChargeDecrease")	flagProcess =10;

        if (processName=="hydrogen_G4DNAExcitation")	 flagProcess =11;
        if (processName=="hydrogen_G4DNAIonisation")	 flagProcess =12;
        if (processName=="hydrogen_G4DNAChargeIncrease")flagProcess =13;

        if (processName=="alpha_G4DNAExcitation")	flagProcess =14;
        if (processName=="alpha_G4DNAIonisation")	flagProcess =15;
        if (processName=="alpha_G4DNAChargeDecrease")	flagProcess =16;

        if (processName=="alpha+_G4DNAExcitation")	flagProcess =17;
        if (processName=="alpha+_G4DNAIonisation")	flagProcess =18;
        if (processName=="alpha+_G4DNAChargeDecrease")	flagProcess =19;
        if (processName=="alpha+_G4DNAChargeIncrease")	flagProcess =20;

        if (processName=="helium_G4DNAExcitation")	flagProcess =21;
        if (processName=="helium_G4DNAIonisation")	flagProcess =22;
        if (processName=="helium_G4DNAChargeIncrease")	flagProcess =23;

        if (processName=="hIoni")	flagProcess =24;
        if (processName=="eIoni")	flagProcess =25;

        fProcess = flagProcess;

        // Get IDs
        fTrackID = aStep->GetTrack()->GetTrackID();
        fRunID   = GetRunID();
        fEventID = GetEventID();

        fNtuple->Fill();
        return true;
    }
    return false;
}
