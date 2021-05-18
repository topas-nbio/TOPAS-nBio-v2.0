/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/


#ifndef NtupleForHiC_EndOfSession_hh
#define NtupleForHiC_EndOfSession_hh

#include "TsParameterManager.hh"

#include "NtupleForHiC_HitPoint.hh"
#include "NtupleForHiC_Clusters.hh"
#include "NtupleForHiC_ClusterAlgo.hh"
#include "ChromosomeParser.hh"
#include "TsDamagePhaseSpaceStore.hh"

//#include <map.hh>
//#include <vector.hh>


class NtupleForHiC_EndOfSession
{
public:
    NtupleForHiC_EndOfSession(TsParameterManager* pM);
    virtual ~NtupleForHiC_EndOfSession();
    
    
    
private:
    
    G4bool Start(TsParameterManager* pM);
    std::map<int,std::vector<NtupleForHiC_HitPoint*>> ReadOutput(G4String filename);
    std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> StartClustering(std::map<G4int,std::vector<NtupleForHiC_HitPoint*>> ExposureHits);
    void AddChromosomeFractions(G4String HiCVerticesFile,
                                std::map<G4int,G4int> &ChromosomeSize,
                                std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> &ExposureClusters,
                                G4ThreeVector &NucDim);
    G4bool SetChromosomeArm(G4int ChromID,
                            G4double FractionAlong);
    TsDamagePhaseSpaceStore* FillPhaseSpaceStore(std::map<G4int,G4int> ChromosomeSize,
                                                 std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> &ExposureClusters,
                                                 G4int &DamageEntries);
    void FillPhaseSpaceHeader(std::map<G4int,G4int> ChromosomeSize,
                              TsDamagePhaseSpaceStore*PSStore,
                              TsParameterManager* pM,
                              G4ThreeVector NucDim,
                              G4int DamageEntries);
    
    void WriteSDD(G4String Filename,
                  TsDamagePhaseSpaceStore*PSStore);
    
    G4bool WriteSSB=false;
    G4bool OverwriteSDD=true;
    G4bool WriteSDDHeader=true;
    

};
#endif
