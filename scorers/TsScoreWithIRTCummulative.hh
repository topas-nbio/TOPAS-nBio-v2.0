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

#ifndef TsScoreWithIRTCummulative_hh
#define TsScoreWithIRTCummulative_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRT;
class TsIRTUtils;

class TsScoreWithIRTCummulative : public TsVNtupleScorer
{
public:
	TsScoreWithIRTCummulative(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreWithIRTCummulative();
	
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	
	void UserHookForEndOfEvent();
    
    virtual void UserHookForPreTimeStepAction();
		
protected:
	
	// Output variables
	G4double fGValue;
	G4double fTime;
	G4String fMoleculeName;
	G4double fEnergy;
	
private:
	G4bool Inside(G4ThreeVector);
	
	TsParameterManager* fPm;
	TsIRT* fIRT;
	TsIRTUtils* fUtils;
	
	std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerEvent;
	std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerEventEverywhere;
	std::vector<G4double> fRandomTimes;
	
	G4double fEnergyDepositPerEvent;
	G4double fEnergyDepositPerEventEverywhere;

	G4double fTotalDose;
	G4double fPrescribedDose;
	G4int fNbOfScoredEvents;
	G4int fNbOfScoredEventsEverywhere;

	G4double fTCut;
	G4String fName;
	
	G4double fTimeMean;
	G4double fTimeStdv;
	G4double fTimeFWHM;
	G4String fTimeDistribution;
	G4int fTimeDistributionType;
    
    G4int fNumberOfPulses;
    G4double fPulsesTimeDelay;
    G4double fDosePerPulse;
    G4double fPulseTimeShift;
	
	std::vector<G4double> fStepTimes;
	std::vector<G4double> fVEnergyDeposit;
	std::vector<G4double> fVEnergyDepositEverywhere;

	G4int fEventID;
	G4int fOldEvent;
	G4double fShiftTime;
	G4double fMinShiftTime;
	
	G4double* fTimeValues;
	G4double* fTimeWeights;
	std::vector<G4double> fTimeTops;
	G4int fNbOfTimes;
	std::ofstream fTimeOutFile;
	G4bool fLowLimitTo1ps;
	
	G4bool fTestIsInside;
    G4String fSensitiveVolume;
};

#endif


