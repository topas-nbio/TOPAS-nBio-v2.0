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

#ifndef TsScoreDBSCAN_hh
#define TsScoreDBSCAN_hh

#include "TsVNtupleScorer.hh"

#include <map>

class ClusteringAlgo;

class TsScoreDBSCAN : public TsVNtupleScorer
{
public:
    TsScoreDBSCAN(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsScoreDBSCAN();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    void UserHookForEndOfEvent();
    
private:
    std::vector<ClusteringAlgo*> fVClustering;
    G4int fNbOfAlgo;
    
    G4int fEventID;
    G4int fSSB;
    G4int fDSB;
    G4int fCSB;
    G4int fSizeX;
    G4int fSizeY;
    
};
#endif
