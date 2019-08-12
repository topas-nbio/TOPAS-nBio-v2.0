// Scorer for Tuple
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

#include "TsScoreMoleculeTuple.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Molecule.hh"

TsScoreMoleculeTuple::TsScoreMoleculeTuple(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fIncludeKineticEnergy(false), fIncludeEventID(false),  fIncludeTrackID(false), fIncludeParentID(false), fIncludeStepNumber(false),
fIncludeParticleName(false), fIncludeProcessName(false), fIncludeVolumeName(false), fIncludeVolumeCopyNumber(false),
fIncludeGlobalTime(false), fIncludeEnergyDeposited(false), fIncludeVertex(false)
{
	SetUnit("");
	
	fNtuple->RegisterColumnI(&fMoleculeID, "MoleculeID or ParticlePDG");
	fNtuple->RegisterColumnF(&fPosX, "Position X", "um");
	fNtuple->RegisterColumnF(&fPosY, "Position Y", "um");
	fNtuple->RegisterColumnF(&fPosZ, "Position Z", "um");
	
	fTimeCut = 1.0 * us;
	if ( fPm->ParameterExists(GetFullParmName("TimeCut") ) )
		fTimeCut = fPm->GetDoubleParameter(GetFullParmName("TimeCut"),"Time");
	
	fIncludeChemistry = true;
	fIncludePhysics = true;
	if ( fPm->ParameterExists(GetFullParmName("IncludeChemicalTrack")) )
		fIncludeChemistry = fPm->GetBooleanParameter(GetFullParmName("IncludeChemicalTrack"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludePhysicalTrack")) )
		fIncludePhysics = fPm->GetBooleanParameter(GetFullParmName("IncludePhysicalTrack"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeKineticEnergy")) )
		fIncludeKineticEnergy = fPm->GetBooleanParameter(GetFullParmName("IncludeKineticEnergy"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeEventID")) )
		fIncludeEventID = fPm->GetBooleanParameter(GetFullParmName("IncludeEventID"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeTrackID")) )
		fIncludeTrackID = fPm->GetBooleanParameter(GetFullParmName("IncludeTrackID"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeParentID")) )
		fIncludeParentID = fPm->GetBooleanParameter(GetFullParmName("IncludeParentID"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeStepNumber")) )
		fIncludeStepNumber = fPm->GetBooleanParameter(GetFullParmName("IncludeStepNumber"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeParticleName")) )
		fIncludeParticleName = fPm->GetBooleanParameter(GetFullParmName("IncludeParticleName"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludePhysicalProcessName")) )
		fIncludeProcessName = fPm->GetBooleanParameter(GetFullParmName("IncludePhysicalProcessName"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeVolumeName")) )
		fIncludeVolumeName = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeName"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeVolumeCopyNumber")) )
		fIncludeVolumeCopyNumber = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeCopyNumber"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeGlobalTime")) )
		fIncludeGlobalTime = fPm->GetBooleanParameter(GetFullParmName("IncludeGlobalTime"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeEnergyDeposited")) )
		fIncludeEnergyDeposited = fPm->GetBooleanParameter(GetFullParmName("IncludeEnergyDeposited"));
	
	if ( fPm->ParameterExists(GetFullParmName("IncludeVertexPosition")) )
		fIncludeVertex = fPm->GetBooleanParameter(GetFullParmName("IncludeVertexPosition"));
	
	if ( fIncludeEventID )
		fNtuple->RegisterColumnI(&fEvt,      "EventID");
	
	if ( fIncludeTrackID )
		fNtuple->RegisterColumnI(&fTrackID,  "Track ID");
	
	if ( fIncludeStepNumber)
		fNtuple->RegisterColumnI(&fStepNumber, "Step number");
	
	if ( fIncludeParticleName )
		fNtuple->RegisterColumnS(&fParticleName, "Particle name");
	
	if ( fIncludeProcessName )
		fNtuple->RegisterColumnS(&fProcessName, "Process name");
		
	if ( fIncludeVolumeName )
		fNtuple->RegisterColumnS(&fVolumeName, "Volume name");
	
	if ( fIncludeVolumeCopyNumber )
		fNtuple->RegisterColumnI(&fVolumeCopyNumber, "Volume copy number");
	
	if ( fIncludeParentID ) {
		fNtuple->RegisterColumnI(&fParentAID, "ParentA ID");
		fNtuple->RegisterColumnI(&fParentBID, "ParentB ID");
	}
	
	if ( fIncludeVertex ) {
		fNtuple->RegisterColumnF(&fVertexPositionX, "Vertex position x", "um");
		fNtuple->RegisterColumnF(&fVertexPositionY, "Vertex position y", "um");
		fNtuple->RegisterColumnF(&fVertexPositionZ, "Vertex position z", "um");
	}
	
	if ( fIncludeGlobalTime )
		fNtuple->RegisterColumnF(&fTime, "Global time", "ps");
	
	if ( fIncludeEnergyDeposited )
		fNtuple->RegisterColumnF(&fEnergyDeposited, "Energy deposited", "keV");
	
	if ( fIncludeKineticEnergy )
		fNtuple->RegisterColumnF(&fKineticEnergy, "Kinetice energy", "keV");
	
}


TsScoreMoleculeTuple::~TsScoreMoleculeTuple() {;}


G4bool TsScoreMoleculeTuple::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	G4Track* aTrack = aStep->GetTrack();
	
	if ( fIncludeChemistry && aTrack->GetTrackID() < 0 ) {
		fTime = aStep->GetPreStepPoint()->GetGlobalTime();
		
		if ( fTime >= fTimeCut ) {
			G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
			fVolumeName = touchable->GetVolume()->GetName();
			fVolumeCopyNumber  = touchable->GetVolume()->GetCopyNo();
			
			fEvt = GetEventID();
			fKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
			G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
			fPosX = pos.x();
			fPosY = pos.y();
			fPosZ = pos.z();
			if ( fIncludeVertex ) {
				G4ThreeVector vpos = aTrack->GetVertexPosition();
				fVertexPositionX = vpos.x();
				fVertexPositionY = vpos.y();
				fVertexPositionZ = vpos.z();
			}
			
			fParentAID = -1;
			fParentBID = -1;
			fTrackID = aTrack->GetTrackID();
			
			fParticleName = GetMolecule(aTrack)->GetName();
			fMoleculeID = GetMolecule(aTrack)->GetMoleculeID();
			
			GetMolecule(aTrack)->GetParentID(fParentAID, fParentBID);
			fEnergyDeposited = 0.0;
			fProcessName = "none";
			fStepNumber = aTrack->GetCurrentStepNumber();
			
			fNtuple->Fill();
			
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return true;
		}
	} else if ( fIncludePhysics && aStep->GetTotalEnergyDeposit() > 0 ){
		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		fVolumeName = touchable->GetVolume()->GetName();
		fVolumeCopyNumber  = touchable->GetVolume()->GetCopyNo();

		fEvt = GetEventID();
		G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
		fPosX = pos.x();
		fPosY = pos.y();
		fPosZ = pos.z();
		if ( fIncludeVertex ) {
			G4ThreeVector vpos = aTrack->GetVertexPosition();
			fVertexPositionX = vpos.x();
			fVertexPositionY = vpos.y();
			fVertexPositionZ = vpos.z();
		}
		fParentAID = aTrack->GetParentID();
		fParentBID = -1;
		fTrackID = aTrack->GetTrackID();
		fEnergyDeposited = aStep->GetTotalEnergyDeposit();
		fStepNumber = aTrack->GetCurrentStepNumber();
		fKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
		fMoleculeID = aTrack->GetParticleDefinition()->GetPDGEncoding();
		fParticleName = aTrack->GetParticleDefinition()->GetParticleName();
		fProcessName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
		
		fNtuple->Fill();
		
		return true;
	}
	
	return false;
}

