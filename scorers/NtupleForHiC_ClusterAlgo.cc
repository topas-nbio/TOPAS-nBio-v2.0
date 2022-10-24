// Extra Class for use by NtupleForHiC_EndOfSession
/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#include "NtupleForHiC_ClusterAlgo.hh"
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <unistd.h>

NtupleForHiC_ClusterAlgo::NtupleForHiC_ClusterAlgo(){}
NtupleForHiC_ClusterAlgo::~NtupleForHiC_ClusterAlgo()
{
//    DeleteClusters();
}


void NtupleForHiC_ClusterAlgo::RunDBSCAN(std::vector<NtupleForHiC_HitPoint*> &Hits)
{
    G4int size=Hits.size();
    G4int maxClustID=0;
    
    
    
    for (G4int i=0; i<size; i++){
        if (!Hits[i]->GetIsVisited()){
            Hits[i]->SetIsVisited(true);
            G4int NumlocalPts=1;
            
            //a list of the IDs for all points within local region
            std::vector<G4int> LocalPtsID;
            
            //determine nb pts within min sep on opp strand
            RegionQuery(Hits,
                        i,
                        NumlocalPts,
                        LocalPtsID);
            
            if (NumlocalPts>=minPts){
                ExpandCluster(maxClustID,
                              Hits,
                              LocalPtsID,
                              NumlocalPts);
            }
        }
    }
    
    
    //check if any clusters can be merged
    MergeClustersCentre();

    
    //Hits not in a cluster go into cluster ID 0
    NtupleForHiC_Clusters*clust=new NtupleForHiC_Clusters();
    clust->SetIsSSB(true);
    clustermap[0]=clust;
    for (unsigned int i=0;i<Hits.size();i++){
        if (Hits[i]->GetClusterID()<=0){
            clustermap[0]->AddPoint(Hits[i]);
        }
    }
    
}


//Determine the number of points within minSep of point PtNum that are on opposite strand
void NtupleForHiC_ClusterAlgo::RegionQuery(std::vector<NtupleForHiC_HitPoint*> &Hits,
                                           G4int PtNum,
                                           G4int &NumlocalPts,
                                           std::vector<G4int> &LocalPtsID)
{
    G4int size = Hits.size();
    G4ThreeVector Pos1 = Hits[PtNum]->GetPosition();
    LocalPtsID.push_back(PtNum);
    
    for (G4int i=0; i<size; i++){
        if (i!=PtNum){
            G4double sep = (Pos1-Hits[i]->GetPosition()).mag();
            
            //if within min separation, and opposite strand, and in same chrom territory
            //i.e. cannot cluster with itself - if isolated Pt then NumlocalPts=0
            //also cannot cluster with a point in a different chromosome
            if (sep<=minSep &&
                Hits[PtNum]->GetStrand()!= Hits[i]->GetStrand() &&
                Hits[PtNum]->GetChromID()==Hits[i]->GetChromID()){
                
                NumlocalPts++;
                LocalPtsID.push_back(i);
            }
        }

    }
}

//Add all points to the same cluster
void NtupleForHiC_ClusterAlgo::ExpandCluster(G4int &maxClustID,
                                             std::vector<NtupleForHiC_HitPoint*> &Hits,
                                             std::vector<G4int> &LocalPtsID,
                                             G4int NumlocalPts)
{
    G4int clusterID = 0;

    std::vector<G4int> tempIDS;
    //Check if any points are already in a cluster and set all points as core
    //If multiple clusterIDs found, then return lowest clusterID and delete all other clusters
    for (G4int pt=0; pt<NumlocalPts; pt++){
        G4int pointID=LocalPtsID[pt];
        if (Hits[pointID]->GetClusterID()!=0){ //if in cluster
            clusterID=Hits[pointID]->GetClusterID();
            G4int tempIDSsize=tempIDS.size();
            G4bool found=false;
            for (G4int i=0;i<tempIDSsize;i++){
                if (tempIDS[i]==clusterID){found=true; break;}
            }
            if (!found){tempIDS.push_back(clusterID);}
        }
    }
    
    //if more than 1 cluster found in region - merge
    if (tempIDS.size()>1){
        clusterID=MergeClusters(tempIDS);
        //Add any local hits not in a cluster to the merged cluster
        for (G4int pt=0; pt<NumlocalPts; pt++){
            if(Hits[LocalPtsID[pt]]->GetClusterID()==0){
                clustermap[clusterID]->AddPoint(Hits[LocalPtsID[pt]]);
                clustermap[clusterID]->SetChromID(Hits[LocalPtsID[pt]]->GetChromID());
                clustermap[clusterID]->SetChromatid(Hits[LocalPtsID[pt]]->GetChromatid());
                clustermap[clusterID]->SetChromCopy(Hits[LocalPtsID[pt]]->GetChromCopy());
                clustermap[clusterID]->SetIsDSB(true);
            }
        }
    }
    
    //otherwise add pts to cluster found or create new cluster
    else{
        NtupleForHiC_Clusters*clust;
        //If no points in a cluster, create a new cluster
        if (clusterID==0){
            maxClustID++;
            clusterID=maxClustID;
            clust=new NtupleForHiC_Clusters;
            clust->SetClusterID(clusterID);
            clust->SetIsDSB(true);
            
            clustermap[clusterID]=clust;
        }
        
        //Find relevant cluster based on clusterID
        if (clustermap.find(clusterID)!=clustermap.end()){
            clust=clustermap[clusterID];
        }else{
            std::cout<<"Err in expand, cant find cluster: "<<clusterID<<std::endl;
        }
        
        std::vector<G4int> TempIDs; //a list of cluster IDs that must be merged, found when looking for boundary Pts
        for (G4int pt=0; pt<NumlocalPts; pt++){
            G4int pointID=LocalPtsID[pt];
            Hits[pointID]->SetIsVisited(true);
            
            clust->AddPoint(Hits[pointID]);
            clustermap[clusterID]->SetChromID(Hits[pointID]->GetChromID());
            clustermap[clusterID]->SetChromatid(Hits[pointID]->GetChromatid());
            clustermap[clusterID]->SetChromCopy(Hits[pointID]->GetChromCopy());
            
            //Find boundary pts -- Not Necessary
            FindBoundaryPts(Hits,
                            pointID,
                            clusterID,
                            TempIDs);
        }
     
        if (TempIDs.size()>1){MergeClusters(TempIDs);}
    }
}


//if multiple points in the region all belong to different clusters -- slow?
int NtupleForHiC_ClusterAlgo::MergeClusters(std::vector<G4int> &clusterIDs)
{
    G4int nClust = clusterIDs.size();
//    std::cout<<"merging clusters: ";
//    for (int i=0; i<nClust; i++){
//        if (i!=nClust-1){std::cout<<clusterIDs[i]<<", ";}
//        else{std::cout<<clusterIDs[i]<<std::endl;}
//    }
    
    //find lowest cluster ID
    int lowestClusterID=numeric_limits<G4int>::max();
    for (G4int ID=0; ID<nClust; ID++){
        if (clusterIDs[ID]<lowestClusterID){lowestClusterID=clusterIDs[ID];}
    }
    
    //delete other clusters
    for (G4int i=0; i<nClust; i++){
        if (clusterIDs[i]!=lowestClusterID){
            if (clustermap.find(clusterIDs[i])!=clustermap.end()){
                //move points to lowest cluster
                std::vector<NtupleForHiC_HitPoint*>Hits = clustermap[clusterIDs[i]]->GetHits();
                G4int HitSize=Hits.size();
                for (G4int j=0;j<HitSize;j++){
                    clustermap[lowestClusterID]->AddPoint(Hits[j]);
                    
                    //not necessary to set chrom id of cluster again
                    //but will check that you're not changing the chrom ID incorrectly
                    clustermap[lowestClusterID]->SetChromID(Hits[j]->GetChromID());
                    clustermap[lowestClusterID]->SetChromatid(Hits[j]->GetChromatid());
                    clustermap[lowestClusterID]->SetChromCopy(Hits[j]->GetChromCopy());
                    
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
    //std::cout<<"lowest: "<<lowestClusterID<<std::endl;
    return lowestClusterID;
}


//Take a core point and find boundary points
void NtupleForHiC_ClusterAlgo::FindBoundaryPts(std::vector<NtupleForHiC_HitPoint*> &Hits,
                                               G4int currPt,
                                               G4int clusterID,
                                               std::vector<G4int>&TempIDs)
{
    NtupleForHiC_Clusters*clust;
    G4int size = Hits.size();
    G4ThreeVector pos1=Hits[currPt]->GetPosition();
    G4int strand1=Hits[currPt]->GetStrand();
    
    //Find relevant cluster based on clusterID
    if (clustermap.find(clusterID)!=clustermap.end()){
        clust=clustermap[clusterID];
    }
    else{
        std::cout<<"Err in boundary, cant find cluster: "<<clusterID<<std::endl;
    }
    
    for (G4int pt=0; pt<size; pt++){
        G4double sep = (pos1-Hits[pt]->GetPosition()).mag();
        if (sep<=minSep &&
            Hits[pt]->GetStrand()!=strand1 &&
            Hits[pt]->GetChromID()==Hits[currPt]->GetChromID() &&
            Hits[pt]->GetChromCopy() == Hits[currPt]->GetChromCopy() &&
            Hits[pt]->GetChromatid() == Hits[currPt]->GetChromatid()){
            Hits[pt]->SetIsVisited(true);
            G4int tempID = Hits[pt]->GetClusterID();
            
            clust->AddPoint(Hits[pt]);
            clust->SetChromID(Hits[pt]->GetChromID());
            clust->SetChromatid(Hits[pt]->GetChromatid());
            clust->SetChromCopy(Hits[pt]->GetChromCopy());
            
            G4int TempIDsSize=TempIDs.size();
            G4bool found=false;
            for (G4int i=0;i<TempIDsSize;i++){
                if (TempIDs[i]==tempID){found=true; break;}
            }
            if (!found && tempID!=0){TempIDs.push_back(tempID);} // a list of unique cluster IDs that need to be merged
            
        }
    }
}

//merge clusters if centres are within min sep
void NtupleForHiC_ClusterAlgo::MergeClustersCentre()
{
    std::vector<G4int> deletedClust;
    
    std::map<G4int, NtupleForHiC_Clusters*>::iterator data1, data2;
   
    for (data1=clustermap.begin();data1 != clustermap.end(); data1++){
        G4ThreeVector centre1 = data1->second->GetCentre();
        
        for (data2=data1;data2 != clustermap.end(); data2++){
            //if not currently the same cluster but from the same chromosome territory
            if (data1->second->GetClusterID()!=data2->second->GetClusterID() &&
                data1->second->GetChromID()==data2->second->GetChromID()){
                G4double sep = (centre1-data2->second->GetCentre()).mag();
                
                
                
                if (sep<=minSep){
                    
                    
                    if (data1->second->GetClusterID()>data2->second->GetClusterID()){
                        vector<NtupleForHiC_HitPoint*> Hits = data1->second->GetHits();
                        G4int HitsSize=Hits.size();
                        for (G4int k=0; k<HitsSize; k++){
                            data2->second->AddPoint(Hits[k]);
                        }
                        deletedClust.push_back(data1->first);
                    }
                    else{
                        vector<NtupleForHiC_HitPoint*> Hits = data2->second->GetHits();
                        G4int HitsSize=Hits.size();
                        for (G4int k=0; k<HitsSize; k++){
                            data1->second->AddPoint(Hits[k]);
                        }
                        deletedClust.push_back(data2->first);
                    }
                }
            }
        }
    }
    
    G4int mergesize=deletedClust.size();
    for (G4int i=0; i<mergesize; i++){
        if (clustermap.find(deletedClust[i])!=clustermap.end()){
            delete clustermap[deletedClust[i]];
            clustermap[deletedClust[i]]=NULL;
            clustermap.erase(deletedClust[i]);
        }
    }
    
    //renumber clusters -- avoid missing cluster ID
    //if (mergesize>0){
        G4int ID=1;
        for (auto data:clustermap){
            data.second->SetClusterID(ID);
            ID++;
        }
  //  }

}


void NtupleForHiC_ClusterAlgo::DeleteClusters()
{
    for (auto data:clustermap){
        delete data.second;
        data.second=NULL;
    }
    
    clustermap.clear();
}





