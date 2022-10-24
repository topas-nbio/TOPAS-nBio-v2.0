#ifndef TsIRT_hh
#define TsIRT_hh

#include "TsIRTConfiguration.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"

#include <vector>
#include <map>

class TsParameterManager;
class TsIRTUtils;
class TsIRT {
	
public:
	TsIRT(TsParameterManager*, G4String);
	~TsIRT();
	
	void runIRT();
	
	void AddMolecule(TsIRTConfiguration::TsMolecule);
	void AddMolecule(G4Step*, G4double, G4int, G4ThreeVector);
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(const G4Track*, G4double, G4int, G4ThreeVector);
	
	void Clean();
	
	G4bool Inside(G4ThreeVector);
	
	inline std::vector<G4double> GetStepTimes() {return fStepTimes; };
	inline TsIRTUtils* GetUtils() { return fUtils; };
	
	std::map<G4String, std::map<G4double, G4int>> GetGValues();
	std::map<G4int, std::map<G4double, G4int>> GetDeltaGValues();
	std::map<G4String, std::map<G4double, G4int>> GetGValuesInVolume();
	std::map<G4int, std::pair<G4int,G4int>> GetReactedDNA();
	
	std::pair<G4String, G4String> GetReactants(G4int);
	
	std::vector<G4String> GetProducts(G4int);
	
private:
	void initializeScorers();
	void FindBinIndexes(G4ThreeVector thisPos, G4double rcutOff);
	void contactReactions(G4int i);
	G4String GetFullParmName(G4String name);
	
private:
	
	TsParameterManager* fPm;
	TsIRTConfiguration* fReactionConf;
	TsIRTUtils* fUtils;
	
	G4String fName;
	
	std::vector<TsIRTConfiguration::TsMolecule> fChemicalSpecies;
	std::vector<TsIRTConfiguration::TsMolecule> fVChemicalSpecies;
	std::vector<std::pair<G4double,G4int>> fVTimes;
	
	std::vector<G4double> fStepTimes;
	std::map<G4int, G4String> fMolecules;
	std::map<G4String, G4int> fMoleculesIDs;
	std::map<G4String, std::map<G4double, G4int>> fGValues;
	std::map<G4int, std::map<G4double, G4int>> fDeltaGValues;
	std::map<G4String, std::map<G4double, G4int>> fGValuesInVolume;
	std::map<G4int, std::pair<G4int,G4int>> fReactedDNA;
	
	G4String fSpinBehavior;
	G4bool fUseSpinScaled;
	G4bool fHighTimeScavenger;
	G4bool fReportDelta;
	G4double fRCutOff;
	G4double fBinWidth;
	G4double fXMin;
	G4double fXMax;
	G4double fYMin;
	G4double fYMax;
	G4double fZMin;
	G4double fZMax;
	G4int fNx;
	G4int fNy;
	G4int fNz;
	G4int fVerbosity;
	G4double fDx;
	G4double fDy;
	G4double fDz;
	
	G4int fxiniIndex;
	G4int fxendIndex;
	G4int fyiniIndex;
	G4int fyendIndex;
	G4int fziniIndex;
	G4int fzendIndex;
	
	G4bool fTestForContactReactions;
	G4bool fSortByTime;
	G4int fGlobalIndex;
	G4bool fTestIsInside;
	G4bool fScorersInitialized;
	G4bool fReactedByContact;
	
	std::map<G4int, std::map<G4int, G4int>> fTheGvalue;
	std::map<G4int, std::map<G4int, G4int>> fTheGvalueInVolume;
	std::map<G4int,std::map<G4int,std::map<G4int,std::vector<G4int>>>> fSpaceBinned;
	std::vector<G4bool> fUsed;
	
};
#endif

