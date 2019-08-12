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

#ifndef TsScorePDB4DNA_hh
#define TsScorePDB4DNA_hh

#include "TsVNtupleScorer.hh"
#include "PDBlib.hh"
#include <map>

class ClusteringAlgo;

class TsScorePDB4DNA : public TsVNtupleScorer
{
public:
    TsScorePDB4DNA(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    
    virtual ~TsScorePDB4DNA();
    
    G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
    void UserHookForEndOfEvent();
    
private:
    void ComputeStrandBreaks(G4int*, G4int);
    
private:
    G4String fPDBFileName;
    PDBlib fPDBlib;
    Molecule* fMoleculeList;
    Barycenter* fBarycenterList;
    
    G4int fThresDistForDSB;
    G4double fThresEdepForSSB;
    
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand1;
    std::map<G4int, std::map<G4int, G4double> > fVEdepStrand2;
    
    G4int fNbOfAlgo;
    G4int fEventID;
    G4int fSSB;
    G4int fDSB;
    
};
#endif
