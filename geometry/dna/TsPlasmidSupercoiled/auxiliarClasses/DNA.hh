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
/*
 *  Developed by Nicholas Henthorn, PRECISE, University of Manchester
 *  nicholas.henthorn@manchester.ac.uk
 *  https://tinyurl.com/y7xuxw9d
 *  17/10/2018, LINK TO PUBLICATION
 */

#ifndef DNA_h
#define DNA_h
#include <vector>
#include <string>
#include <cmath>
#include <globals.hh>
#include "TsVGeometryComponent.hh"
#include <G4SystemOfUnits.hh>
#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

using namespace std;

class DNA{
public:
    DNA(){}
    ~DNA(){}
    
    void SetIsBase(bool opt){
        if (opt && isBack){
            G4cout<<"Error, trying to set backbone as base"<<G4endl;
        } else {
            isBase=opt;
        }
    }
    bool GetIsBase(){return isBase;}
    
    void SetIsBack(bool opt){
        if (opt && isBase){
            G4cout<<"Error, trying to set base as backbone"<<G4endl;
        } else {
            isBack=opt;
        }
    }
    bool GetIsBack(){return isBack;}
    
    void SetPos(G4ThreeVector p){pos=p;}
    G4ThreeVector GetPos(){return pos;}
    
    void SetRot(G4RotationMatrix*r){rot=r;}
    G4RotationMatrix* GetRot(){return rot;}
    
    void SetStrand(G4int str){strand=str;}
    G4int GetStrand(){return strand;}
    
    void SetBP(G4int b){BP=b;}
    G4int GetBP(){return BP;}
    
private:
    bool isBase=false;
    bool isBack=false;
    
    G4ThreeVector pos;
    G4RotationMatrix *rot;
    
    G4int strand=0;
    
    G4int BP=0;
};


#endif /* DNA_h */
