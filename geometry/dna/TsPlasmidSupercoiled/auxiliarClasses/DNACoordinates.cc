// Extra Class for use by TsPlasmid
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

#include <stdio.h>
#include <fstream>
#include "DNACoordinates.hh"

using namespace std;
using namespace CLHEP;

void DNACoordinates::Generate(vector<G4ThreeVector> &path,
                              vector<DNA*> &DNAPts,
                              G4bool BuildHalfCyl,
                              G4bool BuildQuartCyl,
                              G4bool BuildSphere,
							  G4bool segment)
{
    //split the path into 0.34nm steps
	vector<G4ThreeVector> newPath; // = path;
	if ( segment )
    	Segment(path, newPath);
	else
		newPath = path;
    
    if (BuildSphere){
        BuildSphereDNA(newPath, DNAPts);
    } else if (BuildHalfCyl || BuildQuartCyl){
        BuildDNA(newPath, DNAPts);
    } else {
        BuildDNA(newPath, DNAPts);
    }
}

void DNACoordinates::Segment(vector<G4ThreeVector> &path, vector<G4ThreeVector> &newPath)
{
    double rise=0.34*nm;
    int nPts=path.size();
    int counter=0;
    newPath.push_back(path[0]);
    for (int i=0;i<nPts-1;i++){
        G4ThreeVector vector = (path[i+1]-newPath[counter]);
        double length = vector.mag();
        int nDiv = length/rise;
        G4ThreeVector unit = vector.unit();
        unit *= rise;
        for (int j=0;j<nDiv;j++){
            G4ThreeVector nuwe = newPath[counter] + unit;
            newPath.push_back(nuwe);
            counter++;
        }
    }
}

void DNACoordinates::BuildDNA(vector<G4ThreeVector> &newPath, vector<DNA*> &DNAPts)
{
    //double helixRadius = 1.0;
    double rotPair = ((2.0*pi)/10.0);   //10bp per turn
    //double majorGroove = 2.2*nm;
    //double minorGroove = 1.2*nm;
    
    int nBP = newPath.size();
    //nBP=50;
    
    
    //write the coordiantes of the DNA elements
    bool writePositions=false;
    
    ofstream backpos1, backpos2, basepos1, basepos2;

    if (writePositions){
        backpos1.open("BackPos1.dat");
        backpos2.open("BackPos2.dat");
        basepos1.open("BasePos1.dat");
        basepos2.open("BasePos2.dat");
    }
    
    
    for (int bp=0; bp<nBP-1; bp++){
        //Position of base + back in xy
        //Definitely gives right handed coil (checked) -- left handed in -ve z?
        double angle1 = -(double)bp * rotPair;
        //double angle2 = angle1 + pi;

        
        //Rotation to point to next plane
        G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
        G4ThreeVector norm (0.,0.,1.); //the normal to the plane (G4 build planes facing -z)
        G4double DotProd = norm.dot(vecNext);
        G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
        G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm
        
        //set up new 3Vectors for rotated pos
        G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.);

        //Apply rotation
        G4RotationMatrix *rot1 = new G4RotationMatrix;
        rot1->rotate(AngBetween, cross);
        rot1->rotateZ(angle1);
        
        
        G4RotationMatrix *rot2 = new G4RotationMatrix;
        rot2->rotate(AngBetween, cross);
        rot2->rotateZ(angle1);
        
        //Translate
        back1+=newPath[bp];
        back2+=newPath[bp];
        base1+=newPath[bp];
        base2+=newPath[bp];
        
        int bpID=bp+1;
        
        //Add the DNA
        DNA * bac1 = new DNA;
        bac1->SetIsBack(true);
        bac1->SetIsBase(false);
        bac1->SetPos(back1);
        bac1->SetBP(bpID);
        bac1->SetStrand(1);
        bac1->SetRot(rot1);
        DNAPts.push_back(bac1);
        
        DNA * bac2 = new DNA;
        bac2->SetIsBack(true);
        bac2->SetIsBase(false);
        bac2->SetPos(back2);
        bac2->SetBP(bpID);
        bac2->SetStrand(2);
        bac2->SetRot(rot2);
        DNAPts.push_back(bac2);
        
        DNA * bas1 = new DNA;
        bas1->SetIsBase(true);
        bas1->SetIsBack(false);
        bas1->SetPos(base1);
        bas1->SetBP(bpID);
        bas1->SetStrand(1);
        bas1->SetRot(rot1);
        DNAPts.push_back(bas1);
        
        DNA * bas2 = new DNA;
        bas2->SetIsBase(true);
        bas2->SetIsBack(false);
        bas2->SetPos(base2);
        bas2->SetBP(bpID);
        bas2->SetStrand(2);
        bas2->SetRot(rot2);
        DNAPts.push_back(bas2);
        

        if (writePositions){
            basepos1
            <<base1.x()/nm<<" "<<base1.y()/nm<<" "<<base1.z()/nm<<endl;
            basepos2
            <<base2.x()/nm<<" "<<base2.y()/nm<<" "<<base2.z()/nm<<endl;
            backpos1
            <<back1.x()/nm<<" "<<back1.y()/nm<<" "<<back1.z()/nm<<endl;
            backpos2
            <<back2.x()/nm<<" "<<back2.y()/nm<<" "<<back2.z()/nm<<endl;
        }
        
    }

    if (writePositions){
        basepos1.close();
        basepos2.close();
        backpos1.close();
        backpos1.close();
    }
    
}

void DNACoordinates::BuildSphereDNA(vector<G4ThreeVector> &newPath, vector<DNA*> &DNAPts)
{
    double helixRadius = 1.15 * nm;
	double rBack = helixRadius - 0.29*nm;// helixRadius - 0.24*nm;
	double rBase = helixRadius - 2.0*0.29*nm - 0.3*nm;// helixRadius - 2.0 * 0.24*nm - 0.208*nm;  //rBack - 0.24*nm - 0.208*nm;
    double rotPair = ((2.0*pi)/10.0);   //10bp per turn
    //double majorGroove = 2.2*nm;
    //double minorGroove = 1.2*nm;
    //double RatioMajor = majorGroove/(minorGroove+majorGroove);
    
    int nBP = newPath.size();
    //nBP=200;
    
    
    //write the coordiantes of the DNA elements
    bool writePositions=true;
    
	ofstream plasmidPositions;
    
    if (writePositions){
        plasmidPositions.open("plasmid.xyz");
    }
    
    
    for (int bp=0; bp<nBP-1; bp++){
        //Position of base + back in xy
        //Definitely gives right handed coil (checked) -- left handed in -ve z?
        double angle1 = -(double)bp * rotPair;
        double angle2 = angle1+pi;// + (120.0*pi/180.0); //offset for strand2 (major and minor groove)
        
        G4ThreeVector back1temp = G4ThreeVector((rBack*cos(angle1)), (rBack*sin(angle1)), 0.0);
        G4ThreeVector back2temp = G4ThreeVector((rBack*cos(angle2)), (rBack*sin(angle2)), 0.0);
        G4ThreeVector base1temp = G4ThreeVector((rBase*cos(angle1)), (rBase*sin(angle1)), 0.0);
        G4ThreeVector base2temp = G4ThreeVector((rBase*cos(angle2)), (rBase*sin(angle2)), 0.0);
        
        //Rotation to point to next plane
        G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
        G4ThreeVector norm (0.,0.,1.); //the normal to the plane (G4 build planes facing -z)
        G4double DotProd = norm.dot(vecNext);
        G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
        G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm
        
        //set up new 3Vectors for rotated pos
        G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.);
        
        //Apply rotation
        G4RotationMatrix *rot = new G4RotationMatrix;
        rot->rotate(AngBetween, cross);
		//rot->rotateZ(angle1);
		
        ApplyRotation(back1, back1temp, rot);
        ApplyRotation(back2, back2temp, rot);
        ApplyRotation(base1, base1temp, rot);
        ApplyRotation(base2, base2temp, rot);
        delete rot;
        
        //Translate
        back1+=newPath[bp];
        back2+=newPath[bp];
        base1+=newPath[bp];
        base2+=newPath[bp];
        
        int bpID=bp+1;
        
        //Add the DNA
        DNA * bac1 = new DNA;
        bac1->SetIsBack(true);
        bac1->SetIsBase(false);
        bac1->SetPos(back1);
        bac1->SetBP(bpID);
        bac1->SetStrand(1);
        bac1->SetRot(0);
        DNAPts.push_back(bac1);
        
        DNA * bac2 = new DNA;
        bac2->SetIsBack(true);
        bac2->SetIsBase(false);
        bac2->SetPos(back2);
        bac2->SetBP(bpID);
        bac2->SetStrand(2);
        bac2->SetRot(0);
        DNAPts.push_back(bac2);
        
        DNA * bas1 = new DNA;
        bas1->SetIsBase(true);
        bas1->SetIsBack(false);
        bas1->SetPos(base1);
        bas1->SetBP(bpID);
        bas1->SetStrand(1);
        bas1->SetRot(0);
        DNAPts.push_back(bas1);
        
        DNA * bas2 = new DNA;
        bas2->SetIsBase(true);
        bas2->SetIsBack(false);
        bas2->SetPos(base2);
        bas2->SetBP(bpID);
        bas2->SetStrand(2);
        bas2->SetRot(0);
        DNAPts.push_back(bas2);
        
        
        if (writePositions){
            plasmidPositions<< bpID << " " << 102 << " " << base1.x()/nm<<" "<<base1.y()/nm<<" "<<base1.z()/nm<<endl;
            plasmidPositions<< -bpID << " " << 103 << " " << base2.x()/nm<<" "<<base2.y()/nm<<" "<<base2.z()/nm<<endl;
            plasmidPositions<< bpID << " " << 104 << " " << back1.x()/nm<<" "<<back1.y()/nm<<" "<<back1.z()/nm<<endl;
            plasmidPositions<< -bpID << " " << 104 << " " << back2.x()/nm<<" "<<back2.y()/nm<<" "<<back2.z()/nm<<endl;
        }
        
    }
    
    if (writePositions){
		plasmidPositions.close();
    }
    
}


void DNACoordinates::ApplyRotation(G4ThreeVector &rotated, G4ThreeVector &vector, G4RotationMatrix *rot)
{
    rotated[0] = vector[0]*rot->xx() + vector[1]*rot->yx() + vector[2]*rot->zx();
    rotated[1] = vector[0]*rot->xy() + vector[1]*rot->yy() + vector[2]*rot->zy();
    rotated[2] = vector[0]*rot->xz() + vector[1]*rot->yz() + vector[2]*rot->zz();
    
}

