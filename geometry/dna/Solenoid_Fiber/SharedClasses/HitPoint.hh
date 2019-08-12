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
//  Points.h
//  Clustering
//
//  Created by Nick Henthorn on 12/10/2016.

#ifndef HitPoint_h
#define HitPoint_h 1

#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include "Randomize.hh"
#include <G4ThreeVector.hh>



using namespace std;

class HitPoint{
public:
    HitPoint(){};
    ~HitPoint(){};

    void SetPosition(G4ThreeVector p){Position=p;}
    G4ThreeVector GetPosition(){return Position;}

    void SetStrand(G4int stra){Strand=stra;}
    G4int GetStrand(){return Strand;}

    void SetEnergy(G4double ene){Energy=ene;}
    G4double GetEnergy(){return Energy;}

    void SetIsBack(G4bool opt){isBack=opt;}
    G4bool GetIsBack(){return isBack;}
    void SetIsBase(G4bool opt){isBase=opt;}
    G4bool GetIsBase(){return isBase;}

    void SetBP(G4int ba){BP=ba;}
    G4int GetBP(){return BP;}

    void SetClusterID(G4int CI){ClusterID=CI;}
    G4int GetClusterID(){return ClusterID;}

    void SetChecked(G4bool opt){IsCheckedForDuplicates = opt;}
    G4bool GetChecked(){return IsCheckedForDuplicates;}

    void SetIsIonis(bool opt){IsIonisation = opt;}
    G4bool GetIsIonis(){return IsIonisation;}

    void SetMarkedForDelete(bool opt){MarkedForDelete=opt;}
    bool GetMarkedForDelete(){return MarkedForDelete;}

    void SetIsChem(bool opt){isChem=opt;}
    bool GetIsChem(){return isChem;}

    void SetIsPhys(bool opt){isPhys=opt;}
    bool GetIsPhys(){return isPhys;}

    void SetIsVisited(bool opt){Visited=opt;}
    bool isVisited(){return Visited;}

    void SetChromID(int id){chromID=id;}
    int GetChromID(){return chromID;}

    void SetChromCopy(int val){chromCopy=val;}
    int GetChromCopy(){return chromCopy;}

    void SetChromatid(int val){chromatid=val;}
    int GetChromatid(){return chromatid;}

    void SetRunID(int RID){RunID=RID;}
    int GetRunID(){return RunID;}

    void SetEventID(int EID){EventID=EID;}
    int GetEventID(){return EventID;}

    //for the DSB Reporting
    void SetNBack(int num){nBack=num;}
    int GetNBack(){return nBack;}

    void SetNBase(int num){nBase=num;}
    int GetNBase(){return nBase;}

    void SetSpec(string spe){Spec=spe;}
    string GetSpec(){return Spec;}

    void SetCause(int ca){cause=ca;}
    int GetCause(){return cause;}


private:
    G4double Energy=0.0;
    G4ThreeVector Position=G4ThreeVector(0,0,0);
    G4int Strand=0;

    G4int ClusterID=0;

    int chromID=-1;
    int chromCopy=-1;
    int chromatid=-1;


    bool Visited=false;

    bool isBack=false;
    bool isBase=false;
    G4int BP=0;

    bool IsIonisation=false;

    bool IsCheckedForDuplicates=false;

    bool MarkedForDelete=false;

    bool isPhys=false;
    bool isChem=false;

    int RunID=-1;
    int EventID=-1;

    //for the DSB reporting
    int nBack=0;
    int nBase=0;
    string Spec;
    int cause=-1;
};


#endif /* Points_h */
