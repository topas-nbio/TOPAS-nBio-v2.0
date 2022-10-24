// Scorer for simplessbanddsbwithdbscan
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

#include "TsScoreSimpleSSBandDSBWithDBSCAN.hh"
#include "TsTrackInformation.hh"

#include "ClusteringAlgo.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

#include <map>

TsScoreSimpleSSBandDSBWithDBSCAN::TsScoreSimpleSSBandDSBWithDBSCAN(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
																   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("");
	
	// Default parameters for DBSCAN.
	G4double epsilon = 3.2 * nm;
	G4double probability = 0.16;
	G4double minEnergy = 5.0 * eV;
	G4double maxEnergy = 37.5 * eV;
	G4int minNumOfElements = 2;
	fNbOfAlgo = 1;
	
	G4bool sampleEnergyFromRampDistribution = true;
	
	if ( fPm->ParameterExists(GetFullParmName("MaximumDistanceForDefiningDSB")) )
		epsilon = fPm->GetDoubleParameter(GetFullParmName("MaximumDistanceForDefiningDSB"),"Length");
	if ( fPm->ParameterExists(GetFullParmName("ScoreHitsWithProbability")) )
		probability = fPm->GetUnitlessParameter(GetFullParmName("ScoreHitsWithProbability"));
	if ( fPm->ParameterExists(GetFullParmName("EnergyThresholdForHavingSSB")) )
		minEnergy = fPm->GetDoubleParameter(GetFullParmName("EnergyThresholdForHavingSSB"), "Energy");
	if ( fPm->ParameterExists(GetFullParmName("EnergyLimitForHavingSSB")) )
		maxEnergy = fPm->GetDoubleParameter(GetFullParmName("EnergyLimitForHavingSSB"), "Energy");
	if ( fPm->ParameterExists(GetFullParmName("MinimumNumberOfSSBtoFormDSB")) )
		minNumOfElements = fPm->GetIntegerParameter(GetFullParmName("MinimumNumberOfSSBtoFormDSB"));
	if ( fPm->ParameterExists(GetFullParmName("SampleEnergyFromRampDistribution")))
		sampleEnergyFromRampDistribution = fPm->GetBooleanParameter(GetFullParmName("SampleEnergyFromRampDistribution"));
	
	// This is for variance reduction.
	if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
		fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
	
	for ( int i = 0; i < fNbOfAlgo; i++ ) // Not binning in energy or time are valid
		fVClustering.push_back( new ClusteringAlgo(epsilon, minNumOfElements, probability,
												   minEnergy, maxEnergy, sampleEnergyFromRampDistribution) );
	
	fBasePairDepth = 0;
	if (fPm->ParameterExists(GetFullParmName("BasePairPositionAtGeometricHierarchy")))
		fBasePairDepth = fPm->GetIntegerParameter(GetFullParmName("BasePairPositionAtGeometricHierarchy"));
	
	G4String strand1MaterialName = "G4_WATER";
	if (fPm->ParameterExists(GetFullParmName("Strand1MaterialName")))
		strand1MaterialName = fPm->GetStringParameter(GetFullParmName("Strand1MaterialName"));

	G4String strand2MaterialName = "G4_WATER";
	if (fPm->ParameterExists(GetFullParmName("Strand2MaterialName")))
		strand2MaterialName = fPm->GetStringParameter(GetFullParmName("Strand2MaterialName"));
	
	fStrand1Material = GetMaterial(strand1MaterialName);
	fStrand2Material = GetMaterial(strand2MaterialName);
	
	fNtuple->RegisterColumnI(&fEventID, "Event number");
	fNtuple->RegisterColumnI(&fSSB,     "Single strand breaks");
	fNtuple->RegisterColumnI(&fDSB,     "Double strand breaks");
	fNtuple->RegisterColumnI(&fCSB,     "Complex strand breaks");
	fNtuple->RegisterColumnI(&fSizeX,   "Clusters sizes");
	fNtuple->RegisterColumnI(&fSizeY,   "Clusters size weights");
	
	SuppressStandardOutputHandling();
}


TsScoreSimpleSSBandDSBWithDBSCAN::~TsScoreSimpleSSBandDSBWithDBSCAN() {
	for ( int i = 0; i < fNbOfAlgo; i++ )
		delete fVClustering[i];
}


G4bool TsScoreSimpleSSBandDSBWithDBSCAN::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	G4double edep = aStep->GetTotalEnergyDeposit();
	
	if ( edep > 0. ) {
		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
		
		if ( material == fStrand1Material || material == fStrand2Material ) {
			
			//G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
			//G4int numStrand  = touchable->GetVolume(fBasePairDepth)->GetCopyNo();
			//G4int parentDepth = touchable->GetVolume(fBasePairDepth+1)->GetCopyNo();
			
			G4Track* aTrack = aStep->GetTrack();
			G4int index = 1;
			
			if ( 1 < fNbOfAlgo ) {
				TsTrackInformation* trackInformation = (TsTrackInformation*)aTrack->GetUserInformation();
				index = trackInformation->GetSplitTrackID();
			}
			
			G4int sideOfStrand = material == fStrand1Material ? 1 : 2;
			
			if ( 2 < index ) {
				fVClustering[index-3]->RegisterDamage(preStep->GetPosition(), edep, sideOfStrand);
			} else {
				for ( int i = 0; i < fNbOfAlgo; i++ ) {
					fVClustering[i]->RegisterDamage(preStep->GetPosition(), edep, sideOfStrand);
				}
			}
			return true;
		}
	}
	return false;
}


void TsScoreSimpleSSBandDSBWithDBSCAN::UserHookForEndOfEvent() {
	
	fEventID = GetEventID();
	for ( int i = 0; i < fNbOfAlgo; i++ ) {
		std::map<G4int, G4int> sizeDistribution = fVClustering[i]->RunClustering();
		fSSB = fVClustering[i]->GetSSB();
		fCSB = fVClustering[i]->GetComplexSSB();
		fDSB = fVClustering[i]->GetDSB();
		
		while ( !sizeDistribution.empty() ) {
			fSizeX = sizeDistribution.begin()->first;
			fSizeY = sizeDistribution.begin()->second;
			sizeDistribution.erase(sizeDistribution.begin());
		}
		if ( fSSB > 0 || fDSB > 0 )
			fNtuple->Fill();
		fVClustering[i]->Purge();
	}
}


