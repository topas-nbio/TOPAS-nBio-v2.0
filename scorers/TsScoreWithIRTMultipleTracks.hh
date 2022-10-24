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

#ifndef TsScoreWithIRTMultipleTracks_hh
#define TsScoreWithIRTMultipleTracks_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRT;

class TsScoreWithIRTMultipleTracks : public TsVNtupleScorer
{
public:
    TsScoreWithIRTMultipleTracks(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                      G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsScoreWithIRTMultipleTracks();
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
 
    void UserHookForEndOfEvent();   

protected:
    
    void Output();
    void Clear();
    
    // Output variables
    G4double fGValue;
    G4double fGValueError;
    G4double fTime;
    G4String fMoleculeName;
	
    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime;
    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime2;
   
private:
    TsParameterManager* fPm;
	TsIRT* fIRT;

	std::vector<TsIRTConfiguration::TsMolecule> fSpecies;
	
    G4double fEnergyDepositPerEvent;
    G4int fNbOfScoredEvents;
    G4double fEnergyLossKill;
    G4double fEnergyLossAbort;
    G4double fEnergyLoss;
	G4String fName;
	
    G4double fTCut;
    G4double fTotalTrackLength;
    G4double fMaximumTrackLength;
    G4double fXMin;
    G4double fYMin;
    G4double fZMin;
    G4double fXMax;
    G4double fYMax;
    G4double fZMax;
    
    G4int  fMaximumNumberOfSteps;
    G4bool fUseMaximumNumberOfSteps;
    G4int  fNumberOfSteps;
    
    G4double fEnergyDepositPlusEnergyKinetic;
    G4double fLET;
	
	G4bool fUseMultipleTracks;
	G4int fNumberOfMultipleTracks;
	G4int fNumberOfTracksPerEvent;
	
	std::vector<G4double> fVEnergyDeposit;
	G4double* fVTimeDelay;
	G4double* fSpatialOffsetX;
	G4double* fSpatialOffsetY;
	G4double* fSpatialOffsetZ;
	G4double fTheTotalEdep;
};

#endif

