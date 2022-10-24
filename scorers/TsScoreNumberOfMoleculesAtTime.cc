// Scorer for TsSBSMoleculesAtTime
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

#include "TsScoreNumberOfMoleculesAtTime.hh"

#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4MoleculeCounter.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4RunManager.hh"

TsScoreNumberOfMoleculesAtTime::TsScoreNumberOfMoleculesAtTime(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM)
{
    G4MoleculeCounter::Instance()->Use();
    G4DNAChemistryManager::Instance()->ResetCounterWhenRunEnds(false);

    fNtuple->RegisterColumnI(&fNumberOfMolecules, "Total number of molecules at time");
    fNtuple->RegisterColumnF(&fMeanNumberOfMolecules, "Mean number of molecules per event at time", "");
    fNtuple->RegisterColumnF(&fTime,    "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");

	G4int tBins = 10;
	if ( fPm->ParameterExists(GetFullParmName("TimeBins")))
		tBins = fPm->GetIntegerParameter(GetFullParmName("TimeBins"));

	G4double tBinMin = 1.0 * ps;
	if ( fPm->ParameterExists(GetFullParmName("TimeBinMin")))
		tBinMin = fPm->GetDoubleParameter(GetFullParmName("TimeBinMin"),"Time");

	G4double tBinMax = 1.0 * us;
	if ( fPm->ParameterExists(GetFullParmName("TimeBinMax")))
		tBinMax = fPm->GetDoubleParameter(GetFullParmName("TimeBinMax"),"Time");

	G4bool tBinLog = true;
	if ( fPm->ParameterExists(GetFullParmName("TimeBinLog")))
		tBinLog = fPm->GetBooleanParameter(GetFullParmName("TimeBinLog"));
	
	if ( !tBinLog ) {
		for ( int i = 0; i < tBins; i++ )
			fTimeToRecord.push_back(tBinMin + i * (tBinMax - tBinMin)/tBins);
		
	} else {
		if ( tBinMin < 1.0 * ps ) {
			G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
			G4cerr << "The scorer named: " << GetName() << " has TimeBinMin lower than 1 ps" << G4endl;
			fPm->AbortSession(1);
		}
		
		G4double logXMin = std::log10(tBinMin);
		G4double logXMax = std::log10(tBinMax);
		for ( int i = 0; i < tBins; i++ )
			fTimeToRecord.push_back( std::pow(10, logXMin + i * (logXMax - logXMin)/(tBins-1)) );
		
	}

    fEnergyLossKill = 0.0;
    if ( fPm->ParameterExists(GetFullParmName("KillPrimaryIfEnergyLossExceeds")) )
        fEnergyLossKill = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"),"Energy");

    fEnergyLossAbort = 0.0;
    if ( fPm->ParameterExists(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds") ) )
        fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"), "Energy");

    fMaximumTrackLength = 0.0;
    if ( fPm->ParameterExists(GetFullParmName("KillPrimaryBasedOnTrackLength")) )
        fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");

    fNbOfScoredEvents = 0;
}


TsScoreNumberOfMoleculesAtTime::~TsScoreNumberOfMoleculesAtTime()
{;}


G4bool TsScoreNumberOfMoleculesAtTime::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    if ( -1 < aStep->GetTrack()->GetTrackID() ) { // physical tracks
        if ( aStep->GetTrack()->GetParentID() == 0 ) {
            G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
            G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;
            fEnergyLoss += eLoss;
            fTotalTrackLength += aStep->GetStepLength();

            if ( 0 < fEnergyLossAbort && fEnergyLoss > fEnergyLossAbort ) {
                std::cout << " -- Aborting event " << GetEventID() << std::endl;
                G4RunManager::GetRunManager()->AbortEvent();
                return true;
            }

            if ( 0 < fEnergyLossKill && fEnergyLoss >= fEnergyLossKill ) {
                std::cout << " -- Killing primary based on maximum accumulated energ loss." << std::endl;
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return true;
            }

            if ( 0 < fMaximumTrackLength && fTotalTrackLength >= fMaximumTrackLength ) {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                std::cout << " -- Killing primary based on maximum accumulated track-length." << std::endl;
                return true;
            }
        }

        if ( aStep->GetTotalEnergyDeposit() > 0 ) 
            return true;
       
    }

    return false;
}


void TsScoreNumberOfMoleculesAtTime::AccumulateEvent() {

    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) { 
        G4MoleculeCounter::RecordedMolecules species = G4MoleculeCounter::Instance()->GetRecordedMolecules();

        if ( species.get() == 0 ) return;

        auto molecule_it = species->begin();

        while ( molecule_it != species->end() ) {
            for ( int iTime = 0; iTime < (int)fTimeToRecord.size(); iTime++ ) {
                G4int moleculeAtSpecificTime = G4MoleculeCounter::Instance()->GetNMoleculesAtTime((*molecule_it), fTimeToRecord[iTime]);
                G4String moleculeID = (*molecule_it)->GetName();
                fMoleculesPerTime[moleculeID][fTimeToRecord[iTime]] += moleculeAtSpecificTime;
            }
            molecule_it++;
        }

        fNbOfScoredEvents++;
    }
    G4MoleculeCounter::Instance()->ResetCounter();
    fEnergyLoss = 0.0;
    fTotalTrackLength = 0.0;
}


void TsScoreNumberOfMoleculesAtTime::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

    TsScoreNumberOfMoleculesAtTime* workerMoleculesAtTimeScorer = dynamic_cast<TsScoreNumberOfMoleculesAtTime*>(workerScorer);

    fNbOfScoredEvents += workerMoleculesAtTimeScorer->fNbOfScoredEvents;

    std::map<G4String, std::map<G4double, G4int> >::iterator wIter;
    std::map<G4String, std::map<G4double, G4int> >::iterator mIter;

    for ( wIter = workerMoleculesAtTimeScorer->fMoleculesPerTime.begin(); wIter != workerMoleculesAtTimeScorer->fMoleculesPerTime.end(); wIter++) {
        mIter = fMoleculesPerTime.find( wIter->first );
        if ( mIter == fMoleculesPerTime.end() ) {
            fMoleculesPerTime.insert(std::pair<G4String, std::map<G4double, G4int> > ( wIter->first, wIter->second));
        } else {
            std::map<G4double, G4int>::iterator witer;
            std::map<G4double, G4int>::iterator miter;
            for ( witer = (wIter->second).begin(); witer != (wIter->second).end(); witer++ ) {
                miter = (mIter->second).find( witer->first );
                miter->second += witer->second;
            }
        }
    }

    workerMoleculesAtTimeScorer->fMoleculesPerTime.clear();
    workerMoleculesAtTimeScorer->fNbOfScoredEvents = 0;

}


void TsScoreNumberOfMoleculesAtTime::Output() {

    if ( fNbOfScoredEvents == 0 )
        return;

    std::map<G4String, std::map<G4double, G4int> >::iterator wIter;
    std::map<G4double, G4int>::iterator iter;

    for ( wIter = fMoleculesPerTime.begin(); wIter != fMoleculesPerTime.end(); wIter++ ) {
        for ( iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++ ) {
            fNumberOfMolecules = iter->second;
            fMeanNumberOfMolecules = fNumberOfMolecules/fNbOfScoredEvents;
            fTime = iter->first;
            fMoleculeName = wIter->first;
            fNtuple->Fill();
        }
    }

    fNtuple->Write();

}


void TsScoreNumberOfMoleculesAtTime::Clear() {
    fMoleculesPerTime.clear();
    fNbOfScoredEvents = 0;
    UpdateFileNameForUpcomingRun();
}

