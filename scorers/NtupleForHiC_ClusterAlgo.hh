/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#ifndef NtupleForHiC_ClusterAlgo_h
#define NtupleForHiC_ClusterAlgo_h 1

#include <vector>
#include <map>
#include <G4ThreeVector.hh>
#include "NtupleForHiC_HitPoint.hh"
#include "NtupleForHiC_Clusters.hh"



class NtupleForHiC_ClusterAlgo{
public:
    NtupleForHiC_ClusterAlgo();
    ~NtupleForHiC_ClusterAlgo();
    
   void RunDBSCAN(std::vector<NtupleForHiC_HitPoint*> &Hits);
    
    std::map<G4int,NtupleForHiC_Clusters*> GetClusterMap(){return clustermap;}
    
private:
    int minPts = 2;
    double minSep = 3.2*nm;
    
    void RegionQuery(std::vector<NtupleForHiC_HitPoint*> &Hits,
                     G4int PtNum,
                     G4int &NumlocalPts,
                     std::vector<G4int> &LocalPtsID);
    
    void ClusterCheck(NtupleForHiC_HitPoint*Pt1,
                      NtupleForHiC_HitPoint*Pt2);
    
    void ExpandCluster(G4int &maxClustID,
                       std::vector<NtupleForHiC_HitPoint*> &Hits,
                       std::vector<G4int> &LocalPtsID,
                       G4int NumlocalPts);
    
    int MergeClusters(std::vector<G4int> &ClusterIDs);
    
    void FindBoundaryPts(std::vector<NtupleForHiC_HitPoint*> &Hits,
                         G4int currPt,
                         G4int clusterID,
                         std::vector<G4int>&TempIDs);
    
    
    void MergeClustersCentre();
    
    vector<NtupleForHiC_Clusters*> listClust;
    
    void DeleteClusters();
        
    std::map<G4int,NtupleForHiC_Clusters*> clustermap;
    
    
};


#endif
