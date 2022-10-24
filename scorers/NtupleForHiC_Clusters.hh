/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#ifndef NtupleForHiC_Clusters_h
#define NtupleForHiC_Clusters_h 1

#include "NtupleForHiC_HitPoint.hh"
#include <G4ThreeVector.hh>

class NtupleForHiC_Clusters{
public:
    NtupleForHiC_Clusters(){}
    ~NtupleForHiC_Clusters(){}
    
    void AddPoint(NtupleForHiC_HitPoint*Hit){
        bool present=false;
        for (int i=0; i<GetNbPts(); i++){
            if (Hits[i]==Hit){present=true;}
        }
        
        //always reassign cluster ID to this cluster ID (if present should already be clusternum)
        Hit->SetClusterID(clusterID);
        
        if(!present){
            // Hit->SetClusterID(clusterID);
            Hits.push_back(Hit);
            Centre.push_back(Hit->GetPosition());
            chromID=Hit->GetChromID();
            chromCopy=Hit->GetChromCopy();
            chromatid=Hit->GetChromatid();
            ClusterTimesInMinutes.push_back(Hit->GetHitTimeInMinutes());
        }
        
        if (Hit->GetIsBase() && !present){nBases++;}
        else if (Hit->GetIsBack() && !present){nBacks++;}
        
    }
    
    void SetCentre(G4ThreeVector pos){cent=pos;}
    G4ThreeVector GetCent(){return cent;}
    
    G4ThreeVector GetCentre(){
        G4int size = Centre.size();
        G4ThreeVector COM(0.,0.,0.);
        for (G4int i=0; i<size; i++){
            COM+=Centre[i];
        }
        return COM/(double)size;
    }
    
    G4ThreeVector GetMinPos(){
        G4int size = Centre.size();
        G4ThreeVector Min(0.,0.,0.);
        G4double MinDis = std::numeric_limits<G4double>::max();
        for (G4int i=0; i<size; i++){
            G4double test = sqrt((Centre[i].x()*Centre[i].x())
                                 +(Centre[i].y()*Centre[i].y())
                                 +(Centre[i].z()*Centre[i].z()));
            if (test<MinDis){Min=Centre[i];}
        }
        return Min;
    }
    G4ThreeVector GetMaxPos(){
        G4int size = Centre.size();
        G4ThreeVector Max(0.,0.,0.);
        G4double MaxDis = 0.0;
        for (G4int i=0; i<size; i++){
            G4double test = sqrt((Centre[i].x()*Centre[i].x())
                               +(Centre[i].y()*Centre[i].y())
                               +(Centre[i].z()*Centre[i].z()));
            if (test>MaxDis){Max=Centre[i];}
        }
        return Max;
    }
    
    std::vector<NtupleForHiC_HitPoint*> GetHits(){return Hits;}
    
    G4int GetNbPts(){return Hits.size();}
    
    void SetClusterID(G4int ID){
        clusterID=ID;
        for (G4int i=0; i<GetNbPts(); i++){
            Hits[i]->SetClusterID(clusterID);
        }
    }
    G4int GetClusterID(){return clusterID;}
    
    void SetnBases(G4int num){
        nBaseSet=true;
        nBases=num;
    }
    G4int GetnBases(){
        G4int n=nBases;
        if (!nBaseSet){
            n=0;
            G4int nHits=Hits.size();
            for (G4int i=0;i<nHits;i++){
                if (Hits[i]->GetIsBase()){
                    n++;
                }
            }
        }
        return n;
    }
    
    void SetnBacks(G4int num){
        nBackSet=true;
        nBacks=num;
    }
    G4int GetnBacks(){
        G4int n=nBacks;
        if (!nBackSet){
            n=0;
            G4int nHits=Hits.size();
            for (G4int i=0;i<nHits;i++){
                if (Hits[i]->GetIsBack()){
                    n++;
                }
            }
        }
        
        if (GetIsDSB()){
            return n-2;
        }
        return n;
    }
    
    void SetIsDSB(G4bool opt){isDSB=opt;}
    G4bool GetIsDSB(){return isDSB;}
    
    void SetIsSSB(G4bool opt){isSSB=opt;}
    G4bool GetIsSSB(){return isSSB;}
    
    void SetIsIsoBase(G4bool opt){isIsoBase=opt;}
    G4bool GetIsIsoBase(){return isIsoBase;}
    
    void SetChromID(G4int id){chromID=id;}
    G4int GetChromID(){return chromID;}
    
    void SetChromCopy(G4int val){chromCopy=val;}
    G4int GetChromCopy(){return chromCopy;}
    
    void SetChromatid(G4int val){chromatid=val;}
    G4int GetChromatid(){return chromatid;}
    
    void SetChromArm(G4bool opt){chromArm=opt;}
    G4bool GetChromArm(){return chromArm;}
    
    void SetHiCBeadCopy(G4int val){HiCBeadCopy=val;}
    G4int GetHicBeadCopy(){
//        for (int i=0;i<GetHits().size();i++){
//            HiCBeadCopy=GetHits()[i]->GetVolumeCopy();
//        }
        //take the last hit as the copy number
        //(doesnt account for hits in diff beads)
        HiCBeadCopy=GetHits()[GetHits().size()-1]->GetVolumeCopy();
        return HiCBeadCopy;
        
    }
    
    
    G4double GetClusterTimeInMinutes(){
        G4int size = ClusterTimesInMinutes.size();
        double Tim=0.0;
        for (G4int i=0; i<size; i++){
            Tim+=ClusterTimesInMinutes[i];
        }
        return Tim/(double)size;
    }
    
    G4int GetPrimaryID(){
        std::vector<G4int> PIDs;
        for (unsigned int i=0;i<GetHits().size();i++){
            G4bool present=false;
            for (unsigned int j=0;j<PIDs.size();j++){
                if (Hits[i]->GetEventID()==PIDs[j]){present=true;}
            }
            if (!present){
                PIDs.push_back(Hits[i]->GetEventID());
            }
        }
        
        PIDs[0]+=1;
        if (PIDs.size()>1){
            if (clusterID>0){
//                std::cout<<"DSB cluster made from multiple primaries"<<std::endl;
            }
            return PIDs[0]; //always return first primary
        } else {
            return PIDs[0];
        }
        
        
    }
    
    void SetChromosomePosition(G4double po){chromosomeposition=po;}
    G4double GetChromosomePosition(){return chromosomeposition;}
    
    vector<G4int> GetDamageCause(){
        vector<G4int> Causes;
        G4int DirectDam=0;
        G4int IndirectDam=0;
        G4int OtherDam=0;
        for (unsigned int i=0;i<GetHits().size();i++){
            if (GetHits()[i]->GetIsPhys()){DirectDam++;}
            if (GetHits()[i]->GetIsChem()){IndirectDam++;}
            if (GetHits()[i]->GetIsChargeMig()){OtherDam++;}
        }
        if (DirectDam>0 && IndirectDam>0){Causes.push_back(2);}
        else if (DirectDam>0 && IndirectDam==0){Causes.push_back(0);}
        else if (DirectDam==0 && IndirectDam>0){Causes.push_back(1);}
        else {Causes.push_back(3);}
        
        Causes.push_back(DirectDam);
        Causes.push_back(IndirectDam);
        
        return Causes;
    }
    
    
private:
    std::vector<G4ThreeVector> Centre;
    G4ThreeVector cent;
    
    std::vector<NtupleForHiC_HitPoint*> Hits;
    
    G4int clusterID=0;
    
    G4int nBases=0; //num bases
    G4int nBacks=0; //num backbones
    
    G4bool isDSB=false;
    G4bool isSSB=false;
    G4bool isIsoBase=false;
    
    G4int chromID=-1;
    G4int chromCopy=-1;
    G4int chromatid=-1;
    G4bool chromArm=0;
    G4int HiCBeadCopy=-1;
    
    G4double chromosomeposition=-1;
    
    std::vector<G4double> ClusterTimesInMinutes;
    
    G4bool nBackSet=false;
    G4bool nBaseSet=false;
    
    
};
#endif
