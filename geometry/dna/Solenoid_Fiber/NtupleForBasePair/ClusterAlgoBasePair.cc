// Extra Class for use by NtupleForBasePair
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

#include "ClusterAlgoBasePair.hh"

ClusterAlgoBasePair::ClusterAlgoBasePair(){}
ClusterAlgoBasePair::~ClusterAlgoBasePair(){
    //delete the clusters
    DeleteClusters();
}


void ClusterAlgoBasePair::RunDBSCAN(vector <HitPoint*> &Hits)
{
    int size=Hits.size();

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

    //check if any clusters can be merged -- not needed for DNA volumes
    //MergeClustersCentre();

    //renumber the clusters (avoid gaps due to merging)
    RenumberClusters();


    //generate cluster size distribution
    ClusterICSD(Hits);
}


//Determine the number of points within minSep of point PtNum that are on opposite strand
void ClusterAlgoBasePair::RegionQuery(vector<HitPoint*> &Hits, int PtNum, int &NumlocalPts,
                         vector<int> &LocalPtsID)
{

    int size = Hits.size();
    int BP1 = Hits[PtNum]->GetBP();
    LocalPtsID.push_back(PtNum);

    for (int i=0; i<size; i++){
        if (i!=PtNum){
            int BP2=Hits[i]->GetBP();
            int nBPGap=abs(BP1-BP2);
            //bool hasClustered=false;
            bool DiffStrand=false;
            if (Hits[PtNum]->GetStrand()!= Hits[i]->GetStrand()){
                DiffStrand=true;
            }

            //double sep = (Pos1-Hits[i]->GetPosition()).mag();
            //if within min separation and opposite strand
            //i.e. cannot cluster with itself - if isolated Pt then NumlocalPts=0
            if (nBPGap<=minSepBP && DiffStrand){
                NumlocalPts++;
                LocalPtsID.push_back(i);
            //    hasClustered=true;
            }


        }
    }
}

//Add all points to the same cluster
void ClusterAlgoBasePair::ExpandCluster(int &maxClustID, vector <HitPoint*> &Hits, vector <int> &LocalPtsID,
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

            //Find boundary pts -- Not Necessary
            FindBoundaryPts(Hits, pointID, clusterID, TempIDs);
        }
        if (TempIDs.size()>1){MergeClusters(TempIDs);}
    }


}


//if multiple points in the region all belong to different clusters -- slow?
int ClusterAlgoBasePair::MergeClusters(vector <int> &clusterIDs)
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
void ClusterAlgoBasePair::FindBoundaryPts(vector <HitPoint*> &Hits, int currPt, int clusterID, vector<int>&TempIDs)
{

    int size = Hits.size();
    //G4ThreeVector pos1=Hits[currPt]->GetPosition();
    int BP1 = Hits[currPt]->GetBP();
    G4int strand1=Hits[currPt]->GetStrand();

    //Find relevant cluster based on clusterID
    if (clustermap.find(clusterID)!=clustermap.end()){
        clust=clustermap[clusterID];
    }
    else{
        cout<<"Err in boundary, cant find cluster: "<<clusterID<<endl;
    }

    for (int pt=0; pt<size; pt++){
         int nBPGap=abs(BP1-Hits[pt]->GetBP());
        //double sep = (pos1-Hits[pt]->GetPosition()).mag();
        if (nBPGap<=minSepBP && Hits[pt]->GetStrand()!=strand1){
            Hits[pt]->SetIsVisited(true);
            int tempID = Hits[pt]->GetClusterID();

            clust->AddPoint(Hits[pt]);

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
void ClusterAlgoBasePair::MergeClustersCentre()
{
    vector<int> deletedClust;

    map<int, Clusters*>::iterator data1, data2;

    for (data1=clustermap.begin();data1 != clustermap.end(); data1++){
        G4ThreeVector centre1 = data1->second->GetCentre();

        for (data2=data1;data2 != clustermap.end(); data2++){
            if (data1->second->GetClusterID()!=data2->second->GetClusterID()){
                double sep = (centre1-data2->second->GetCentre()).mag();
                if (sep<=minSepBP){
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
    RenumberClusters();

}

//calculate the ionisation cluster size distribution
void ClusterAlgoBasePair::ClusterICSD(vector<HitPoint*> &Hits)
{


    int maxClustSize = 0;
    int nSSB=0;
    vector <int> ClusterSize;

    for (auto data:Hits){
        if(data->GetClusterID()==0){nSSB++;}
    }

    for (auto data:clustermap){
        if (data.second->GetNbPts()>maxClustSize){
            maxClustSize=data.second->GetNbPts();
        }
    }



//    G4cout<<"NClust "<<clustermap.size()<<" "<<nClusters<<G4endl;
//    for (auto data:clustermap){
//        cout<<"ID "<<data.second->GetClusterID()<<endl;
//        cout<<"Size "<<data.second->GetNbPts()<<endl;
//        cout<<"with Hits "<<endl;
//
//        vector<HitPoint*> test = data.second->GetHits();
//        for(int i=0;i<test.size();i++){
//            cout<<"        "<<test[i]->GetHitID()<<endl;
//        }
//        cout<<endl;
//    }


    vector<int> row(2);
    row[0]=0;
    row[1]=0;
    ICSD.push_back(row);
    row[0]=1;
    row[1]=nSSB;
    ICSD.push_back(row);

    for (int cluSize=2; cluSize<maxClustSize+1; cluSize++){
        int f=0;
        for (auto data:clustermap){
            if (data.second->GetNbPts() == cluSize){
                f++;
            }
        }
        row[0]=cluSize;
        row[1]=f;

        ICSD.push_back(row);
    }
}

void ClusterAlgoBasePair::DeleteClusters()
{
    for (auto data:clustermap){
        delete data.second;
        data.second=NULL;
    }

    clustermap.clear();
}

//add bases to clusters if within 3bp -> 3bp defined in parameters.hh
//3bp from the extreme backbones involved in the DSB
//example: (X=included damage, O=excluded damage)

//  -----------X---X---------------X-------------
//  O----------------X------X---------X-O--------
//  -----X--O--------X-X--------X----------------
//  --------X-------------X----------------------

void ClusterAlgoBasePair::ClusterBases(vector<HitPoint *> &BaseHits)
{
    int nBaseHits=BaseHits.size();
    for (int i=0;i<nBaseHits;i++){
        int BaseNum = BaseHits[i]->GetBP();
        int BaseStrand = BaseHits[i]->GetStrand();
        bool AddBase=true;
        for (auto clusters:clustermap){
            vector<HitPoint*> HitsInClust = clusters.second->GetHits();
            int nInClust=HitsInClust.size();
            int maxBP=0;
            int minBP=std::numeric_limits<int>::max();
            //bool hasClustered=false;
            for (int j=0;j<nInClust;j++){
                int BackNum = HitsInClust[j]->GetBP();
                int BackStrand = HitsInClust[j]->GetStrand();
                if (BackNum>maxBP){maxBP=BackNum;}
                if (BackNum<minBP){minBP=BackNum;}

                if ((BaseNum==BackNum && BaseStrand==BackStrand)){
                    AddBase=false;
                }

            }

            //add the base if within 3bp of the extremities and not part of the damage backbone
            if (BaseNum<maxBP+BaseSep && BaseNum>minBP-BaseSep && AddBase){
                clusters.second->AddPoint(BaseHits[i]);
            //    hasClustered=true;
            }

        }
    }
}


//If there is a damaged base on opp strand and <10bp from a damaged backbone
//form a cluster (with -ve cluster ID)
//These are potential isolated damages than can lead to DSB
//i.e. Base Excision Repair
void ClusterAlgoBasePair::PotentialClusters(vector<HitPoint *> &BackHits, vector<HitPoint *> &BaseHits)
{
    //set up vectors for isolated backbone and base damage
    vector<HitPoint*> IsoBack;
    vector<HitPoint*> IsoBase;

    int nBack=BackHits.size();
    for (int i=0;i<nBack;i++){
        if (BackHits[i]->GetClusterID()==0){
            IsoBack.push_back(BackHits[i]);
        }
    }

    int nBase=BaseHits.size();
    for (int i=0;i<nBase;i++){
        if (BaseHits[i]->GetClusterID()==0){
            IsoBase.push_back(BaseHits[i]);
        }
    }

    int nIsoBack=IsoBack.size();
    int nIsoBase=IsoBase.size();

    int runningClustID=0;

    for (int i=0;i<nIsoBack;i++){
        int backStrand=IsoBack[i]->GetStrand();
        int backBP=IsoBack[i]->GetBP();
        for (int j=0;j<nIsoBase;j++){
            int baseStrand=IsoBase[j]->GetStrand();
            int baseBP=IsoBase[j]->GetBP();

            //potential cluster between base and backbone
            if (abs(backBP-baseBP)<=minSepBP &&
                backStrand!=baseStrand){

                //are either of the volumes in a cluster?
                int ClusterID=0;
                int ClustBack=IsoBack[i]->GetClusterID();
                int ClustBase=IsoBase[j]->GetClusterID();

                //both in clust
                if (ClustBack!=0 && ClustBase!=0){
                    int ClustToRemove=0;

                    if (ClustBack<ClustBase){
                        ClusterID=ClustBack;
                        ClustToRemove=ClustBase;
                    }
                    else {
                        ClusterID=ClustBase;
                        ClustToRemove=ClustBack;
                    }

                    //get all the points from the cluster to remove
                    //move all points from higher cluster into lower cluster
                    vector<HitPoint*> MoveHits = clustermap[ClustToRemove]->GetHits();
                    int nHits=MoveHits.size();
                    for (int k=0;k<nHits;k++){
                        clustermap[ClusterID]->AddPoint(MoveHits[k]);
                    }

                    //add the points to the cluster
                    clustermap[ClusterID]->AddPoint(IsoBack[i]);
                    clustermap[ClusterID]->AddPoint(IsoBase[j]);
                }

                //neither in clust -- make new
                if (ClustBack==0 && ClustBase==0){
                    runningClustID--;
                    ClusterID=runningClustID;

                    clust=new Clusters;
                    clust->SetClusterID(ClusterID);
                    clust->AddPoint(IsoBack[i]);
                    clust->AddPoint(IsoBase[j]);
                    clustermap[ClusterID]=clust;
                }

                //one volume in a cluster
                if ((ClustBack!=0 && ClustBase==0) ||
                    (ClustBack==0 && ClustBase!=0)){
                    if (ClustBack<ClustBase){ClusterID=ClustBack;}
                    else {ClusterID=ClustBase;}

                    //add the points to the cluster
                    clustermap[ClusterID]->AddPoint(IsoBack[i]);
                    clustermap[ClusterID]->AddPoint(IsoBase[j]);
                }
            }
        }
    }

}



void ClusterAlgoBasePair::RenumberClusters()
{
    int ID=1;
    for (auto data:clustermap){
        data.second->SetClusterID(ID);
        ID++;
    }
}
