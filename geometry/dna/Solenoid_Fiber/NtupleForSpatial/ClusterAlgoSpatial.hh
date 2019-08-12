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

#ifndef ClusterAlgoSpatial_h
#define ClusterAlgoSpatial_h 1

#include <vector>
#include <G4ThreeVector.hh>
#include "HitPoint.hh"
#include "Clusters.hh"
#include <map>

class ClusterAlgoSpatial{
public:
    ClusterAlgoSpatial();
    ~ClusterAlgoSpatial();

    void RunDBSCAN(vector <HitPoint*> &Hits);


    map<int, Clusters*> GetClusterMap(){return clustermap;}


private:
    int minPts = 2;
    double minSep = 3.32; //in nm

    vector < vector <int> > ICSD;

    void RegionQuery(vector<HitPoint*> &Hits, int PtNum, int &NumlocalPts, vector<int> &LocalPtsID);

    void ClusterCheck(HitPoint*Pt1, HitPoint*Pt2);

    void ExpandCluster(int &maxClustID, vector<HitPoint*> &Hits, vector<int> &LocalPtsID, int NumlocalPts);

    int MergeClusters(vector <int> &ClusterIDs);

    void FindBoundaryPts(vector<HitPoint*> &Hits, int currPt, int clusterID, vector<int>&TempIDs);

    void MergeClustersCentre();

    Clusters * clust;
    vector <Clusters*> listClust;

    void DeleteClusters();

    //clusterID and pointer to cluster
    map<int, Clusters*> clustermap;


};


#endif /* ClusterAlgo_h */
