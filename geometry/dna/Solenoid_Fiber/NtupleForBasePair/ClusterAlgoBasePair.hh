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
//  ClusterAlgo.hh
//  Clustering
//
//  Created by Nick Henthorn on 12/10/2016.

#ifndef ClusterAlgoBasePair_h
#define ClusterAlgoBasePair_h 1

#include <vector>
#include <G4ThreeVector.hh>
#include "HitPoint.hh"
#include "Clusters.hh"
#include <map>


class ClusterAlgoBasePair{
public:
    ClusterAlgoBasePair();
    ~ClusterAlgoBasePair();

    void RunDBSCAN(vector <HitPoint*> &Hits);
    vector < vector <int> > GetICSD(){return ICSD;}

    void ClusterBases(vector<HitPoint*> &BaseHits);

    void PotentialClusters(vector<HitPoint*> &BackHits, vector<HitPoint*> &BaseHits);

private:
    int minPts = 2;
    int BaseSep = 3;
    int minSep = 3.4*nm;
    int minSepBP = 10;

    vector < vector <int> > ICSD;

    void RegionQuery(vector<HitPoint*> &Hits, int PtNum, int &NumlocalPts, vector<int> &LocalPtsID);

    void ClusterCheck(HitPoint*Pt1, HitPoint*Pt2);

    void ExpandCluster(int &maxClustID, vector<HitPoint*> &Hits, vector<int> &LocalPtsID, int NumlocalPts);

    int MergeClusters(vector <int> &ClusterIDs);

    void FindBoundaryPts(vector<HitPoint*> &Hits, int currPt, int clusterID, vector<int>&TempIDs);

    void ClusterICSD(vector<HitPoint*> &Hits);

    void MergeClustersCentre();

    Clusters * clust;
    vector <Clusters*> listClust;

    void DeleteClusters();

    //clusterID and pointer to cluster
    map<int, Clusters*> clustermap;

    void RenumberClusters();




};


#endif /* ClusterAlgo_h */
