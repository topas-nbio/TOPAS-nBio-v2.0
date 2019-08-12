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

#ifndef TsScoreSimpleSSBandDSB_hh
#define TsScoreSimpleSSBandDSB_hh

#include "TsVNtupleScorer.hh"

#include <map>


class G4Material;

class TsScoreSimpleSSBandDSB : public TsVNtupleScorer
{
public:
    TsScoreSimpleSSBandDSB(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsScoreSimpleSSBandDSB();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    void UserHookForEndOfEvent();
    
private:
    void ComputeStrandBreaks(G4int*, G4int);
    
private:
    G4Material* fStrand1Material;
    G4Material* fStrand2Material;
    
    G4int fThresDistForDSB;
    G4double fThresEdepForSSB;
    
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand1;
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand2;
    
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand1;
    std::map<G4int, std::map<G4int, std::map<G4int, G4double> > > fGenVEdepStrand2;
    
    G4int fNbOfAlgo;
    G4int fEventID;
    G4int fDNAParent; // if any
    G4int fSSB;
    G4int fDSB;

    G4int fBasePairDepth;
    
    G4String fStrand1VolumeName;
    G4String fStrand2VolumeName;
    
    
};
#endif
