// Scorer for dbscan
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

#include "TsScoreDBSCAN.hh"
#include "TsTrackInformation.hh"

#include "ClusteringAlgo.hh"
#include "G4SystemOfUnits.hh"

#include <map>

TsScoreDBSCAN::TsScoreDBSCAN(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
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
    
    if ( fPm->ParameterExists(GetFullParmName("MinimumDistanceForDSB")) )
        epsilon = fPm->GetDoubleParameter(GetFullParmName("MinimumDistanceForDSB"),"Length");
    if ( fPm->ParameterExists(GetFullParmName("SampleHitsWithProbability")) )
        probability = fPm->GetUnitlessParameter(GetFullParmName("SampleHitsWithProbability"));
    if ( fPm->ParameterExists(GetFullParmName("LowerEnergyForSamplingSSB")) )
        minEnergy = fPm->GetDoubleParameter(GetFullParmName("LowerEnergyForSamplingSSB"), "Energy");
    if ( fPm->ParameterExists(GetFullParmName("UpperEnergyForSamplingSSB")) )
        maxEnergy = fPm->GetDoubleParameter(GetFullParmName("UpperEnergyForSamplingSSB"), "Energy");
    if ( fPm->ParameterExists(GetFullParmName("MinimumNumberofSSBtoFormDSB")) )
        minNumOfElements = fPm->GetIntegerParameter(GetFullParmName("MinimumNumberofSSBtoFormDSB"));
    
    // This is for variance reduction.
    if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
        fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
    
    for ( int i = 0; i < fNbOfAlgo; i++ ) // Not binning in energy or time are valid
        fVClustering.push_back( new ClusteringAlgo(epsilon, minNumOfElements, probability,
                                                   minEnergy, maxEnergy, sampleEnergyFromRampDistribution) );
    
    fNtuple->RegisterColumnI(&fEventID, "Event number");
    fNtuple->RegisterColumnI(&fSSB,     "Single strand breaks");
    fNtuple->RegisterColumnI(&fDSB,     "Double strand breaks");
    fNtuple->RegisterColumnI(&fCSB,     "Complex strand breaks");
    fNtuple->RegisterColumnI(&fSizeX,   "Clusters sizes");
    fNtuple->RegisterColumnI(&fSizeY,   "Clusters size weights");
    
    SuppressStandardOutputHandling();
}


TsScoreDBSCAN::~TsScoreDBSCAN() {
    for ( int i = 0; i < fNbOfAlgo; i++ )
        delete fVClustering[i];
}


G4bool TsScoreDBSCAN::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    G4double edep = aStep->GetTotalEnergyDeposit();
    
    if ( edep > 0. ) {
        
        G4Track* aTrack = aStep->GetTrack();
        G4int index = 1;
        
        if ( 1 < fNbOfAlgo ) {
            TsTrackInformation* trackInformation = (TsTrackInformation*)aTrack->GetUserInformation();
            index = trackInformation->GetSplitTrackID();
        }
        
        if ( 2 < index ) {
            fVClustering[index-3]->RegisterDamage(aStep->GetPreStepPoint()->GetPosition(), edep, -1);
        } else {
            for ( int i = 0; i < fNbOfAlgo; i++ ) {
                fVClustering[i]->RegisterDamage(aStep->GetPreStepPoint()->GetPosition(), edep, -1);
            }
        }
        return true;
    }
    return false;
}


void TsScoreDBSCAN::UserHookForEndOfEvent() {

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
        
        fNtuple->Fill();
        fVClustering[i]->Purge();
    }
}
