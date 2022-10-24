// Scorer for TsIRTInterTrack
//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreWithIRTMultipleTracks.hh"
#include "TsIRT.hh"
#include "TsIRTConfiguration.hh"

#include "G4ITTrackHolder.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Scheduler.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"
#include "G4Electron_aq.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

TsScoreWithIRTMultipleTracks::TsScoreWithIRTMultipleTracks(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fEnergyLossKill(0), fName(scorerName)
{
	SetUnit("");
	
	fIRT = new TsIRT(fPm, fName);
	
	fNtuple->RegisterColumnD(&fGValue, "GValue: number of molecules per 100 eV of energy deposit", "");
	fNtuple->RegisterColumnD(&fGValueError, "GValue statistical error", "");
	fNtuple->RegisterColumnD(&fTime,    "Time", "picosecond");
	fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
	
	fEnergyLossKill = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"),"Energy");
	fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"),"Energy");
	fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");
	
	fTCut = 1.0 * ps;
	if ( fPm->ParameterExists(GetFullParmName("CutoffTime")) )
		fTCut = fPm->GetDoubleParameter(GetFullParmName("CutoffTime"), "Time");
	
	fNumberOfSteps = 0;
	fMaximumNumberOfSteps = 0;
	fUseMaximumNumberOfSteps = false;
	if ( fPm->ParameterExists(GetFullParmName("MaximumNumberOfSteps")) ) {
		fMaximumNumberOfSteps = fPm->GetIntegerParameter(GetFullParmName("MaximumNumberOfSteps"));
		fUseMaximumNumberOfSteps = true;
	}
	
	fUseMultipleTracks = false;
	fNumberOfTracksPerEvent = 0;
	fNumberOfMultipleTracks = 1;
	if (fPm->ParameterExists(GetFullParmName("UseMultipleTracks")) &&
		fPm->GetBooleanParameter(GetFullParmName("UseMultipleTracks"))) {
		fUseMultipleTracks = true;
		fVTimeDelay = fPm->GetDoubleVector(GetFullParmName("TimeDelay"),"Time");
		fNumberOfMultipleTracks = fPm->GetVectorLength(GetFullParmName("TimeDelay"));
		
		fSpatialOffsetX  = fPm->GetDoubleVector(GetFullParmName("SpatialOffsetX"), "Length");
		fSpatialOffsetY  = fPm->GetDoubleVector(GetFullParmName("SpatialOffsetY"), "Length");
		fSpatialOffsetZ  = fPm->GetDoubleVector(GetFullParmName("SpatialOffsetZ"), "Length");
	}

	for ( int i = 0; i < fNumberOfMultipleTracks; i++ )
		fVEnergyDeposit.push_back(0.0);
	
	fNbOfScoredEvents = 0;
	fEnergyLoss = 0.0;
	fTotalTrackLength = 0.0;
	fEnergyDepositPlusEnergyKinetic = 0.0;
	fTheTotalEdep = 0.0;
}


TsScoreWithIRTMultipleTracks::~TsScoreWithIRTMultipleTracks()
{;}


G4bool TsScoreWithIRTMultipleTracks::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	if ( -1 < aStep->GetTrack()->GetTrackID() ) {
		
		if ( aStep->GetTrack()->GetParentID() == 0 ) {
			
			if ( !fUseMaximumNumberOfSteps ) {
				fTotalTrackLength += aStep->GetStepLength();
				
				G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
				G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;
				fEnergyLoss += eLoss;
				
				fLET += aStep->GetTotalEnergyDeposit();
				
				fEnergyDepositPlusEnergyKinetic += aStep->GetTotalEnergyDeposit();
				
				if ( fEnergyLoss > fEnergyLossAbort ) {
					G4RunManager::GetRunManager()->AbortEvent();
					return true;
				}
				if ( fEnergyLoss >= fEnergyLossKill ) {
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					return true;
				}
				if ( fTotalTrackLength >= fMaximumTrackLength ) {
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					return true;
				}
			} else {
				fNumberOfSteps++;
				if ( fMaximumNumberOfSteps < fNumberOfSteps ) {
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				}
			}
		} else if ( aStep->GetTrack()->GetParentID() == 1 && aStep->GetTrack()->GetCurrentStepNumber() == 1 ) {
			G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();
			if ( ekin < 100.0 * eV ) {
				fEnergyDepositPlusEnergyKinetic += ekin;
				fLET += ekin;
			}
		}
		
		G4double edep = aStep->GetTotalEnergyDeposit();
		if ( edep > 0 ) {
			fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
			if ( fUseMultipleTracks )
				fVEnergyDeposit[fNumberOfTracksPerEvent] += edep * aStep->GetPreStepPoint()->GetWeight();
			return true;
		}
	} else {
		G4double time = aStep->GetTrack()->GetGlobalTime();
		if ( time >= fTCut ) {
			time = 1.0 * ps;
			G4ThreeVector p = aStep->GetPreStepPoint()->GetPosition();
			const G4String name = GetMolecule(aStep->GetTrack())->GetName();
			if ( fUseMultipleTracks ) {
				time += fVTimeDelay[fNumberOfTracksPerEvent];
				fIRT->AddMolecule(aStep, time, 0, G4ThreeVector(fSpatialOffsetX[fNumberOfTracksPerEvent],
																fSpatialOffsetY[fNumberOfTracksPerEvent],
																fSpatialOffsetZ[fNumberOfTracksPerEvent]));
			} else {
				fIRT->AddMolecule(aStep, time, 0, G4ThreeVector());
			}
			
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
		return true;
	}
	return false;
}


void TsScoreWithIRTMultipleTracks::UserHookForEndOfEvent() {
	fNumberOfTracksPerEvent++;
	fTheTotalEdep += fEnergyDepositPerEvent;
	
	if(!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()
	   && (fNumberOfTracksPerEvent == fNumberOfMultipleTracks) ) {
		
		G4cout << " - For event " << GetEventID() << ": " << G4endl;
		G4cout << " --- IRT start for event " << GetEventID() << G4endl;
		fIRT->runIRT();
		std::map<G4String, std::map<G4double, G4int>> irt = fIRT->GetGValues();
		G4cout << " --- IRT ends for event " << GetEventID() << G4endl;
		for ( auto& nameTimeAndGvalue : irt ) {
			G4String name = nameTimeAndGvalue.first;
			for ( auto& timeAndGvalue : (nameTimeAndGvalue.second) ) {
				G4double time = timeAndGvalue.first;
				G4double gvalue = timeAndGvalue.second;
				
				//if ( fUseMultipleTracks ) {
				//	// Assuming 2 tracks:
				//	if (      time >= fVTimeDelay[0] && time < fVTimeDelay[1] )
				//		gvalue *= 100/(fVEnergyDeposit[0]/eV);
				//	else
				//		gvalue *= 100/((fVEnergyDeposit[0]+fVEnergyDeposit[1])/eV);
				//} else {
				//	gvalue *= 100/(fEnergyDepositPerEvent/eV);
				//}
				fGValuePerSpeciePerTime[name][time] += gvalue;
				fGValuePerSpeciePerTime2[name][time] += gvalue*gvalue;
			}
		}
		
		fNbOfScoredEvents++;
		irt.clear();
		
		fEnergyDepositPlusEnergyKinetic /= fTotalTrackLength;
		fLET /= fTotalTrackLength;
		G4cout << " ---- Edep1 : " << fVEnergyDeposit[0]/eV << " eV at " << fVTimeDelay[0]/ps << " ps " << G4endl;
		G4cout << " ---- Edep2 : " << fVEnergyDeposit[1]/eV << " eV at " << fVTimeDelay[1]/ps << " ps " << G4endl;
		G4cout << " ---- T1-T2 : " << fEnergyDepositPlusEnergyKinetic / (keV/um) << " keV/um ---- ";
		G4cout << " ---- LETEdep: " << fLET/(keV/um) << " keV/um." << G4endl;
		
		fNumberOfTracksPerEvent = 0;
		
		fXMin = 1e9*nm;
		fYMin = 1e9*nm;
		fZMin = 1e9*nm;
		
		fXMax = 0e0*nm;
		fYMax = 0e0*nm;
		fZMax = 0e0*nm;
		
		fTheTotalEdep = 0;
		for ( int i = 0; i < fNumberOfMultipleTracks; i++ )
			fVEnergyDeposit[i] = 0.0;
		
		fSpecies.clear();
	}
	
	fEnergyDepositPerEvent = 0.0;
	fEnergyLoss = 0.0;
	fTotalTrackLength = 0.0;
	fNumberOfSteps = 0;
	fEnergyDepositPlusEnergyKinetic = 0.0;
	fLET = 0.0;
}


void TsScoreWithIRTMultipleTracks::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
	
	TsScoreWithIRTMultipleTracks* workerGvalueScorer = dynamic_cast<TsScoreWithIRTMultipleTracks*>(workerScorer);
	
	fNbOfScoredEvents += workerGvalueScorer->fNbOfScoredEvents;
	
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
	std::map<G4String, std::map<G4double, G4double> >::iterator mIter;
	
	for ( wIter = workerGvalueScorer->fGValuePerSpeciePerTime.begin(); wIter != workerGvalueScorer->fGValuePerSpeciePerTime.end(); wIter++) {
		mIter = fGValuePerSpeciePerTime.find( wIter->first );
		if ( mIter == fGValuePerSpeciePerTime.end() ) {
			fGValuePerSpeciePerTime.insert(std::pair<G4String, std::map<G4double, G4double> > ( wIter->first, wIter->second));
		} else {
			std::map<G4double, G4double>::iterator witer;
			std::map<G4double, G4double>::iterator miter;
			for ( witer = (wIter->second).begin(); witer != (wIter->second).end(); witer++ ) {
				miter = (mIter->second).find( witer->first );
				miter->second += witer->second;
			}
		}
	}
	
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter2;
	std::map<G4String, std::map<G4double, G4double> >::iterator mIter2;
	
	for ( wIter2 = workerGvalueScorer->fGValuePerSpeciePerTime2.begin(); wIter2 != workerGvalueScorer->fGValuePerSpeciePerTime2.end(); wIter2++) {
		mIter2 = fGValuePerSpeciePerTime2.find( wIter2->first );
		if ( mIter2 == fGValuePerSpeciePerTime2.end() ) {
			fGValuePerSpeciePerTime2.insert(std::pair<G4String, std::map<G4double, G4double> > ( wIter2->first, wIter2->second));
		} else {
			std::map<G4double, G4double>::iterator witer2;
			std::map<G4double, G4double>::iterator miter2;
			for ( witer2 = (wIter2->second).begin(); witer2 != (wIter2->second).end(); witer2++ ) {
				miter2 = (mIter2->second).find( witer2->first );
				miter2->second += witer2->second;
			}
		}
	}
	
	workerGvalueScorer->fGValuePerSpeciePerTime.clear();
	workerGvalueScorer->fGValuePerSpeciePerTime2.clear();
	workerGvalueScorer->fNbOfScoredEvents = 0;
	workerGvalueScorer->fEnergyDepositPerEvent = 0.0;
	workerGvalueScorer->fSpecies.clear();
}


void TsScoreWithIRTMultipleTracks::Output() {
	
	if ( fNbOfScoredEvents == 0 )
		return;
	
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter2;
	std::map<G4double, G4double>::iterator iter;
	std::map<G4double, G4double>::iterator iter2;
	
	for ( wIter = fGValuePerSpeciePerTime.begin(); wIter != fGValuePerSpeciePerTime.end(); wIter++ ) {
		wIter2 = fGValuePerSpeciePerTime2.find(wIter->first);
		
		for ( iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++ ) {
			iter2 = (wIter2->second).find( iter->first );
			
			fGValue = iter->second/fNbOfScoredEvents;
			if ( fNbOfScoredEvents > 1 ) {
				fGValueError = sqrt( (1.0/(fNbOfScoredEvents-1)) * ( (iter2->second)/fNbOfScoredEvents - fGValue*fGValue));
			} else {
				fGValueError = 1.0;
			}
			fTime = iter->first;
			
			fMoleculeName = wIter->first;
			fNtuple->Fill();
		}
	}
	
	fNtuple->Write();
}


void TsScoreWithIRTMultipleTracks::Clear() {
	fGValuePerSpeciePerTime.clear();
	fNbOfScoredEvents = 0;
	UpdateFileNameForUpcomingRun();
}


