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
//  Clusters.hh
//  Clustering
//
//  Created by Nick Henthorn on 13/10/2016.

#ifndef Clusters_h
#define Clusters_h 1

#include "HitPoint.hh"
#include <G4ThreeVector.hh>


using namespace std;

class Clusters{
public:
    Clusters(){}
    ~Clusters(){}

    void AddPoint(HitPoint*Hit){
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
        }

        if (Hit->GetIsBase() && !present){nBases++;}
        else if (Hit->GetIsBack() && !present){nBacks++;}

    }

    void SetCentre(G4ThreeVector pos){cent=pos;}
    G4ThreeVector GetCent(){return cent;}

    G4ThreeVector GetCentre(){
        int size = Centre.size();
        G4ThreeVector COM(0.,0.,0.);
        for (int i=0; i<size; i++){
            COM+=Centre[i];
        }
        return COM/(double)size;
    }

    vector <HitPoint*> GetHits(){return Hits;}

    int GetNbPts(){return Hits.size();}

    void SetClusterID(int ID){
        clusterID=ID;
        for (int i=0; i<GetNbPts(); i++){
            Hits[i]->SetClusterID(clusterID);
        }
    }
    int GetClusterID(){return clusterID;}

    int GetnBases(){
        int n=0;
        int nHits=Hits.size();
        for (int i=0;i<nHits;i++){
            if (Hits[i]->GetIsBase()){
                n++;
            }
        }
        return n;
    }

    int GetnBacks(){
        int n=0;
        int nHits=Hits.size();
        for (int i=0;i<nHits;i++){
            if (Hits[i]->GetIsBack()){
                n++;
            }
        }
        return n;
    }

    void SetIsDSB(bool opt){isDSB=opt;}
    bool GetIsDSB(){return isDSB;}

    void SetIsSSB(bool opt){isSSB=opt;}
    bool GetIsSSB(){return isSSB;}

    void SetIsIsoBase(bool opt){isIsoBase=opt;}
    bool GetIsIsoBase(){return isIsoBase;}

    void SetChromID(int id){chromID=id;}
    int GetChromID(){return chromID;}

    void SetChromCopy(int val){chromCopy=val;}
    int GetChromCopy(){return chromCopy;}

    void SetChromatid(int val){chromatid=val;}
    int GetChromatid(){return chromatid;}



    //For Damage Spec Output
    void SetNBack(int num){nBack=num;}
    int GetNBack(){return nBack;}

    void SetNBase(int num){nBase=num;}
    int GetNBase(){return nBase;}

    void SetSpec(string spe){Spec=spe;}
    string GetSpec(){return Spec;}

    void SetCause(int ca){cause=ca;}
    int GetCause(){return cause;}


private:
    vector<G4ThreeVector> Centre;
    G4ThreeVector cent;

    vector<HitPoint*> Hits;

    int clusterID=0;

    int nBases=0; //num bases
    int nBacks=0; //num backbones

    bool isDSB=false;
    bool isSSB=false;
    bool isIsoBase=false;


    int chromID=-1;
    int chromCopy=-1;
    int chromatid=-1;



    //For Damage Spec Output
    int nBack=0;
    int nBase=0;
    string Spec;
    int cause=-1;


};
#endif /* Clusters_h */
