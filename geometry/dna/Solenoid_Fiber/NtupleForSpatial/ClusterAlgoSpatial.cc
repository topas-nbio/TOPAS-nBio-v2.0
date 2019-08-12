// Extra Class for use by NtupleForSpatial
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
//  ClusterAlgo.cc
//  Clustering
//
//  Created by Nick Henthorn on 12/10/2016.

#include "ClusterAlgoSpatial.hh"
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <unistd.h>

ClusterAlgoSpatial::ClusterAlgoSpatial(){}
ClusterAlgoSpatial::~ClusterAlgoSpatial()
{
    //delete the clusters -- if using getclusters then a new map must be set up (dont pass by reference)
    DeleteClusters();
}


void ClusterAlgoSpatial::RunDBSCAN(vector <HitPoint*> &Hits)
{
    int size=Hits.size();
    //    if (!Parameters::Instance()->GetRunEventByEvent()){
    //        //cout<<endl<<"Running clustering on "<<size<<" Hits"<<endl;
    //    }
    int maxClustID=0;

    for (int i=0; i<size; i++){
        if (!Hits[i]->isVisited()){
            Hits[i]->SetIsVisited(true);
            int NumlocalPts=1;
            vector <int> LocalPtsID; //a list of the IDs for all points within local region

            RegionQuery(Hits, i, NumlocalPts, LocalPtsID); //determine nb pts within min sep on opp strand

            if (NumlocalPts>=minPts){
                ExpandCluster(maxClustID, Hits, LocalPtsID, NumlocalPts);
            }

        }
    }


    //check if any clusters can be merged
    MergeClustersCentre();
}


//Determine the number of points within minSep of point PtNum that are on opposite strand
void ClusterAlgoSpatial::RegionQuery(vector<HitPoint*> &Hits, int PtNum, int &NumlocalPts,
                                    vector<int> &LocalPtsID)
{
    int size = Hits.size();
    G4ThreeVector Pos1 = Hits[PtNum]->GetPosition();
    LocalPtsID.push_back(PtNum);

    for (int i=0; i<size; i++){
        if (i!=PtNum){
            double sep = (Pos1-Hits[i]->GetPosition()).mag();
            //if within min separation, and opposite strand, and in same chrom territory
            //i.e. cannot cluster with itself - if isolated Pt then NumlocalPts=0
            //also cannot cluster with a point in a different chromosome
            if (sep<=minSep && Hits[PtNum]->GetStrand()!= Hits[i]->GetStrand() && Hits[PtNum]->GetChromID()==Hits[i]->GetChromID()){
                NumlocalPts++;
                LocalPtsID.push_back(i);
            }
        }

    }
}

//Add all points to the same cluster
void ClusterAlgoSpatial::ExpandCluster(int &maxClustID, vector <HitPoint*> &Hits, vector <int> &LocalPtsID,
                                      int NumlocalPts)
{
    int clusterID = 0;

    vector <int> tempIDS;
    //Check if any points are already in a cluster and set all points as core
    //If multiple clusterIDs found, then return lowest clusterID and delete all other clusters
    for (int pt=0; pt<NumlocalPts; pt++){
        int pointID=LocalPtsID[pt];
        if (Hits[pointID]->GetClusterID()!=0){ //if in cluster
            clusterID=Hits[pointID]->GetClusterID();
            int tempIDSsize=tempIDS.size();
            bool found=false;
            for (int i=0;i<tempIDSsize;i++){
                if (tempIDS[i]==clusterID){found=true; break;}
            }
            if (!found){tempIDS.push_back(clusterID);}
        }
    }

    //if more than 1 cluster found in region - merge
    if (tempIDS.size()>1){
        clusterID=MergeClusters(tempIDS);
        //Add any local hits not in a cluster to the merged cluster
        for (int pt=0; pt<NumlocalPts; pt++){
            if(Hits[LocalPtsID[pt]]->GetClusterID()==0){
                clustermap[clusterID]->AddPoint(Hits[LocalPtsID[pt]]);
                clustermap[clusterID]->SetChromID(Hits[LocalPtsID[pt]]->GetChromID());
            }
        }
    }

    //otherwise add pts to cluster found or create new cluster
    else{
        //If no points in a cluster, create a new cluster
        if (clusterID==0){
            maxClustID++;
            clusterID=maxClustID;
            clust=new Clusters;
            clust->SetClusterID(clusterID);

            clustermap[clusterID]=clust;
        }

        //Find relevant cluster based on clusterID
        if (clustermap.find(clusterID)!=clustermap.end()){
            clust=clustermap[clusterID];
        }else{
            cout<<"Err in expand, cant find cluster: "<<clusterID<<endl;
        }

        vector<int> TempIDs; //a list of cluster IDs that must be merged, found when looking for boundary Pts
        for (int pt=0; pt<NumlocalPts; pt++){
            int pointID=LocalPtsID[pt];
            Hits[pointID]->SetIsVisited(true);

            clust->AddPoint(Hits[pointID]);
            clustermap[clusterID]->SetChromID(Hits[pointID]->GetChromID());

            //Find boundary pts -- Not Necessary
            FindBoundaryPts(Hits, pointID, clusterID, TempIDs);
        }

        if (TempIDs.size()>1){MergeClusters(TempIDs);}
    }
}


//if multiple points in the region all belong to different clusters -- slow?
int ClusterAlgoSpatial::MergeClusters(vector <int> &clusterIDs)
{
    int nClust = clusterIDs.size();
    //    cout<<"merging clusters: ";
    //    for (int i=0; i<nClust; i++){
    //        if (i!=nClust-1){cout<<clusterIDs[i]<<", ";}
    //        else{cout<<clusterIDs[i]<<endl;}
    //    }

    //find lowest cluster ID
    int lowestClusterID=numeric_limits<int>::max();
    for (int ID=0; ID<nClust; ID++){
        if (clusterIDs[ID]<lowestClusterID){lowestClusterID=clusterIDs[ID];}
    }

    //delete other clusters
    for (int i=0; i<nClust; i++){
        if (clusterIDs[i]!=lowestClusterID){
            if (clustermap.find(clusterIDs[i])!=clustermap.end()){
                //move points to lowest cluster
                vector <HitPoint*>Hits = clustermap[clusterIDs[i]]->GetHits();
                int HitSize=Hits.size();
                for (int j=0;j<HitSize;j++){
                    clustermap[lowestClusterID]->AddPoint(Hits[j]);

                    //not necessary to set chrom id of cluster again
                    //but will check that you're not changing the chrom ID incorrectly
                    clustermap[lowestClusterID]->SetChromID(Hits[j]->GetChromID());

                }
                delete clustermap[clusterIDs[i]];
                clustermap[clusterIDs[i]] = NULL;

                clustermap.erase(clusterIDs[i]);
            }
            else{
                cout<<"Err in merge, cant find cluster: "<<clusterIDs[i]<<endl;
            }
        }

    }
    //cout<<"lowest: "<<lowestClusterID<<endl;
    return lowestClusterID;
}


//Take a core point and find boundary points
void ClusterAlgoSpatial::FindBoundaryPts(vector <HitPoint*> &Hits, int currPt, int clusterID, vector<int>&TempIDs)
{
    int size = Hits.size();
    G4ThreeVector pos1=Hits[currPt]->GetPosition();
    G4int strand1=Hits[currPt]->GetStrand();

    //Find relevant cluster based on clusterID
    if (clustermap.find(clusterID)!=clustermap.end()){
        clust=clustermap[clusterID];
    }
    else{
        cout<<"Err in boundary, cant find cluster: "<<clusterID<<endl;
    }

    for (int pt=0; pt<size; pt++){
        double sep = (pos1-Hits[pt]->GetPosition()).mag();
        if (sep<=minSep && Hits[pt]->GetStrand()!=strand1 && Hits[pt]->GetChromID()==Hits[currPt]->GetChromID()){
            Hits[pt]->SetIsVisited(true);
            int tempID = Hits[pt]->GetClusterID();

            clust->AddPoint(Hits[pt]);
            clust->SetChromID(Hits[pt]->GetChromID());

            int TempIDsSize=TempIDs.size();
            bool found=false;
            for (int i=0;i<TempIDsSize;i++){
                if (TempIDs[i]==tempID){found=true; break;}
            }
            if (!found && tempID!=0){TempIDs.push_back(tempID);} // a list of unique cluster IDs that need to be merged

        }
    }
}

//merge clusters if centres are within min sep
void ClusterAlgoSpatial::MergeClustersCentre()
{
    vector<int> deletedClust;

    map<int, Clusters*>::iterator data1, data2;

    for (data1=clustermap.begin();data1 != clustermap.end(); data1++){
        G4ThreeVector centre1 = data1->second->GetCentre();

        for (data2=data1;data2 != clustermap.end(); data2++){
            //if not currently the same cluster but from the same chromosome territory
            if (data1->second->GetClusterID()!=data2->second->GetClusterID() &&
                data1->second->GetChromID()==data2->second->GetChromID()){
                double sep = (centre1-data2->second->GetCentre()).mag();



                if (sep<=minSep){


                    if (data1->second->GetClusterID()>data2->second->GetClusterID()){
                        vector<HitPoint*> Hits = data1->second->GetHits();
                        int HitsSize=Hits.size();
                        for (int k=0; k<HitsSize; k++){
                            data2->second->AddPoint(Hits[k]);
                        }
                        deletedClust.push_back(data1->first);
                    }
                    else{
                        vector<HitPoint*> Hits = data2->second->GetHits();
                        int HitsSize=Hits.size();
                        for (int k=0; k<HitsSize; k++){
                            data1->second->AddPoint(Hits[k]);
                        }
                        deletedClust.push_back(data2->first);
                    }
                }
            }
        }
    }

    int mergesize=deletedClust.size();
    for (int i=0; i<mergesize; i++){
        if (clustermap.find(deletedClust[i])!=clustermap.end()){
            delete clustermap[deletedClust[i]];
            clustermap[deletedClust[i]]=NULL;
            clustermap.erase(deletedClust[i]);
        }
    }

    //renumber clusters -- avoid missing cluster ID
    //if (mergesize>0){
    int ID=1;
    for (auto data:clustermap){
        data.second->SetClusterID(ID);
        ID++;
    }
    //  }

}


void ClusterAlgoSpatial::DeleteClusters()
{
    for (auto data:clustermap){
        delete data.second;
        data.second=NULL;
    }

    clustermap.clear();
}
