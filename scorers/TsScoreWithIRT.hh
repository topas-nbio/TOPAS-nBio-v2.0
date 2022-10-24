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

#ifndef TsScoreWithIRT_hh
#define TsScoreWithIRT_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRT;

class TsScoreWithIRT : public TsVNtupleScorer
{
public:
    TsScoreWithIRT(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                      G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsScoreWithIRT();
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
 
    void UserHookForEndOfEvent();
	
	virtual void UserHookForPreTimeStepAction();
	//virtual void UserHookForBeginOfChemicalTrack();

	//virtual void FillIRT(const G4Track*);
	
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

    std::map<G4int, std::map<G4double, G4double> > fDeltaGPerReactionPerTime;
    std::map<G4int, std::map<G4double, G4double> > fDeltaGPerReactionPerTime2;
   
private:
    TsParameterManager* fPm;
	TsIRT* fIRT;

	std::vector<TsIRTConfiguration::TsMolecule> fSpecies;

    G4Timer fIRTTimer[1];
    G4double fExecutionTime;
    G4double fExecutionTimeStdv;
	
    G4double fEnergyDepositPerEvent;
    G4int fNbOfScoredEvents;
    G4double fEnergyLossKill;
    G4double fEnergyLossAbort;
    G4double fEnergyLoss;
	G4String fName;
	
    G4double fTCut;
    G4double fTotalTrackLength;
    G4double fMaximumTrackLength;
    
    G4int  fMaximumNumberOfSteps;
    G4bool fUseMaximumNumberOfSteps;
    G4int  fNumberOfSteps;
    
    G4double fEnergyDepositPlusEnergyKinetic;
    G4double fLET;
		
	std::vector<G4double> fVEnergyDeposit;

    G4bool   fReportDelta;
    G4String fOutputFile;
};

#endif

