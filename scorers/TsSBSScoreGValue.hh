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

#ifndef TsSBSScoreGValue_hh
#define TsSBSScoreGValue_hh

#include "TsVNtupleScorer.hh"
#include "Randomize.hh"

#include <stdint.h>

class G4MolecularConfiguration;

class TsSBSScoreGValue : public TsVNtupleScorer
{
public:
    TsSBSScoreGValue(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsSBSScoreGValue();
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
    void UserHookForEndOfEvent();
   
protected:
//    void AccumulateEvent();
    
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
    
    G4double fEnergyDepositPerEvent;
    G4double* fTimeToRecord;
    G4int fNbTimeToRecord;
    G4int fNbOfScoredEvents;
    G4double fEnergyLossKill;
    G4double fEnergyLossAbort;
    G4double fEnergyLoss;
    G4double fMaximumTrackLength;
    G4double fTotalTrackLength;
    G4int fNbOfMoleculesToScavenge;
    G4int* fMoleculeIDToScavenge;
    G4double* fScavengingCapacity;

    std::map<G4double, G4double> fScavengerProducts;
};

#endif

