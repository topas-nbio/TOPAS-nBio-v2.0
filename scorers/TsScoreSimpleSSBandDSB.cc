// Scorer for simplessbanddsb
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

#include "TsScoreSimpleSSBandDSB.hh"
#include "TsTrackInformation.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <map>

TsScoreSimpleSSBandDSB::TsScoreSimpleSSBandDSB(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
											   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)

{
	SetUnit("");
	
	// Default parameters
	fThresDistForDSB = 10;
	fThresEdepForSSB = 17.5 * eV;
	
	fNbOfAlgo = 1;
	
	if ( fPm->ParameterExists(GetFullParmName("BasePairDistanceForDefiningDSB")) )
		fThresDistForDSB = fPm->GetIntegerParameter(GetFullParmName("BasePairDistanceForDefiningDSB"));
	if ( fPm->ParameterExists(GetFullParmName("EnergyThresholdForHavingSSB")) )
		fThresEdepForSSB = fPm->GetDoubleParameter(GetFullParmName("EnergyThresholdForHavingSSB"), "Energy");
	
	fBasePairDepth = 0;
	if ( fPm->ParameterExists(GetFullParmName("BasePairPositionAtGeometricHierarchy")))
		fBasePairDepth = fPm->GetIntegerParameter(GetFullParmName("BasePairPositionAtGeometricHierarchy"));
	
	G4String strand1MaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("Strand1MaterialName")))
		strand1MaterialName = fPm->GetStringParameter(GetFullParmName("Strand1MaterialName"));

	G4String strand2MaterialName = "G4_WATER";
	if ( fPm->ParameterExists(GetFullParmName("Strand2MaterialName")))
		strand2MaterialName = fPm->GetStringParameter(GetFullParmName("Strand2MaterialName"));
	
	fStrand1Material = GetMaterial(strand1MaterialName);
	fStrand2Material = GetMaterial(strand2MaterialName);
	
	// This is for variance reduction
	if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
		fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
	
	//fThresEdepForSSB /= eV;
	
	fNtuple->RegisterColumnI(&fEventID, "Event number");
	fNtuple->RegisterColumnI(&fDNAParent, "DNA parent geometry");
	fNtuple->RegisterColumnI(&fSSB,       "Single strand breaks");
	fNtuple->RegisterColumnI(&fDSB,       "Double strand breaks");
	
	SuppressStandardOutputHandling();
	
}


TsScoreSimpleSSBandDSB::~TsScoreSimpleSSBandDSB() {
}


G4bool TsScoreSimpleSSBandDSB::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	G4double edep = aStep->GetTotalEnergyDeposit(); ///eV;
	
	if ( edep > 0 ) {
		G4StepPoint* preStep = aStep->GetPreStepPoint();
		G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
		
		if ( material == fStrand1Material || material == fStrand2Material ) {
			G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
			G4int numStrand  = touchable->GetVolume(fBasePairDepth)->GetCopyNo();
			G4int parentDepth = touchable->GetVolume(fBasePairDepth+1)->GetCopyNo();
			
			G4int index = 1;
			if ( 1 < fNbOfAlgo ) {
				TsTrackInformation* trackInformation = (TsTrackInformation*)aStep->GetTrack()->GetUserInformation();
				index = trackInformation->GetSplitTrackID();
			}
			
			if ( 2 < index ) {
				if ( material == fStrand1Material )
					fGenVEdepStrand1[parentDepth][index-3][numStrand] += edep;
				else
					fGenVEdepStrand2[parentDepth][index-3][numStrand] += edep;
			} else {
				if ( material == fStrand1Material )
					for ( int i = 0; i < fNbOfAlgo; i++ )
						fGenVEdepStrand1[parentDepth][i][numStrand] += edep;
				else
					for ( int i = 0; i < fNbOfAlgo; i++ )
						fGenVEdepStrand2[parentDepth][i][numStrand] += edep;
			}
			
			return true;
		}
	}
	return false;
}


void TsScoreSimpleSSBandDSB::UserHookForEndOfEvent() {
	fEventID = GetEventID();
	for ( auto& energyAtStrands : fGenVEdepStrand1 ) {
		G4int parentDepth = energyAtStrands.first;
		fVEdepStrand1 = fGenVEdepStrand1[parentDepth];
		fVEdepStrand2 = fGenVEdepStrand2[parentDepth];
		fDNAParent = parentDepth;
		
		for ( int i = 0; i < fNbOfAlgo; i++ ) {
			G4int sb[2] = {0, 0};
			ComputeStrandBreaks(sb, i);
			fSSB = sb[0];
			fDSB = sb[1];
			if ( fSSB > 0 || fDSB > 0 )
				fNtuple->Fill();
		}
		fVEdepStrand1.erase(fVEdepStrand1.begin(), fVEdepStrand1.end());
		fVEdepStrand2.erase(fVEdepStrand2.begin(), fVEdepStrand2.end());
	}
	fGenVEdepStrand1.erase(fGenVEdepStrand1.begin(), fGenVEdepStrand1.end());
	fGenVEdepStrand2.erase(fGenVEdepStrand2.begin(), fGenVEdepStrand2.end());
}


// This class was taken from Geant4/examples/extended/medical/dna/pdb4dna

void TsScoreSimpleSSBandDSB::ComputeStrandBreaks(G4int* sb, G4int cluster)
{
	// sb quantities
	//
	G4int ssb1=0;
	G4int ssb2=0;
	G4int dsb=0;
	
	// nucleotide id and energy deposit for each strand
	G4int nucl1;
	G4int nucl2;
	G4double edep1;
	G4double edep2;
	
	//Read strand1
	//
	while ( !fVEdepStrand1[cluster].empty() )
	{
		nucl1 = fVEdepStrand1[cluster].begin()->first;
		edep1 = fVEdepStrand1[cluster].begin()->second;
		fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
		
		// SSB in strand1
		//
		if ( edep1 >= fThresEdepForSSB )
		{
			ssb1++;
		}
		
		// Look at strand2
		//
		if ( !fVEdepStrand2[cluster].empty() )
		{
			do
			{
				nucl2 = fVEdepStrand2[cluster].begin()->first;
				edep2 = fVEdepStrand2[cluster].begin()->second;
				if ( edep2 >= fThresEdepForSSB )
				{
					ssb2++;
				}
				fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
			} while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fVEdepStrand2[cluster].empty()) );
			
			// no dsb
			//
			if ( nucl2-nucl1 > fThresDistForDSB )
			{
				fVEdepStrand2[cluster][nucl2]=edep2;
				if ( edep2 >= fThresEdepForSSB )
				{
					ssb2--;
				}
			}
			
			// one dsb
			//
			if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
			{
				if ( ( edep2 >= fThresEdepForSSB ) &&
					( edep1 >= fThresEdepForSSB ) )
				{
					ssb1--;
					ssb2--;
					dsb++;
				}
			}
		}
	}
	
	// End with not processed data
	//
	while ( !fVEdepStrand1[cluster].empty() )
	{
		nucl1 = fVEdepStrand1[cluster].begin()->first;
		edep1 = fVEdepStrand1[cluster].begin()->second;
		if ( edep1 >= fThresEdepForSSB )
		{
			ssb1++;
		}
		fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
	}
	
	while ( !fVEdepStrand2[cluster].empty() )
	{
		nucl2 = fVEdepStrand2[cluster].begin()->first;
		edep2 = fVEdepStrand2[cluster].begin()->second;
		if ( edep2 >= fThresEdepForSSB )
		{
			ssb2++;
		}
		fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
	}
	sb[0]=ssb1+ssb2;
	sb[1]=dsb;
}


