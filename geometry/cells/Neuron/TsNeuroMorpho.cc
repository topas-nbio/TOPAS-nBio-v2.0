// Component for TsNeuroMorpho
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
//Interface to the www.neuromorpho.org database. Download input files directly from the website.

#include "TsNeuroMorpho.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4PVPlacement.hh"

#include "G4UnionSolid.hh"

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>

using namespace std;


TsNeuroMorpho::TsNeuroMorpho(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsNeuroMorpho::~TsNeuroMorpho()
{;}

void TsNeuroMorpho::ResolveParameters() {
    
    // User defined parameters.
    
    G4String name = GetFullParmName("NeuroMorphoFileName");
    if (!fPm->ParameterExists(name)) {
        G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
        G4cerr << "Parameter " << name << " has to be specified for this scorer." << G4endl;
        exit(1);
    }
    G4String NeuroMorphoFileName = fPm->GetStringParameter(name);
    
    HLX = fPm->GetDoubleParameter(GetFullParmName("HLX"), "Length");
    HLY = fPm->GetDoubleParameter(GetFullParmName("HLY"), "Length");
    HLZ = fPm->GetDoubleParameter(GetFullParmName("HLZ"), "Length");
    
    /****************  Read in the NeuroMorpho file: ******************************
     NeuroMorpho file consists of format SWC contains the following information:
     
     1. an integer number as compartment identifier
     2. type of neuronal compartment
        0 - undefined
        1 - soma
        2 - axon
        3 - basal dendrite
        4 - apical dendrite
        5 - custom (user-defined preferences)
        6 - unspecified neurites
        7 - glia processes
     3. x coordinate of the compartment
     4. y coordinate of the compartment
     5. z coordinate of the compartment
     6. radius of the compartment
     7. parent compartment
     
     The header of each file begins with "#" - skip these lines.
     Every compartment has only one parent and the parent compartment for the first point in each file
     is always -1 (if the file does not include the soma information then the originating point of the
     tree will be connected to a parent of -1). The index for parent compartments are always less than
     child compartments. Loops and unconnected branches are excluded. All trees should originate from
     the soma and have parent type if the file includes soma information. Soma can be a single point or
     more than one point. When the soma is encoded as one line in the SWC, it is interpreted as a
     "sphere". When it is encoded by more than 1 line, it could be a set of tapering cylinders
     (as in some pyramidal cells) or even a 2D projected contour("circumference").
     ***************************************************************************************************/
    
    const char* FileName = NeuroMorphoFileName;
    std::string line = "";
    
    ifstream f(FileName, ios::in);
    fPositions.reserve(3000);
    fPartID.reserve(3000);
    fRadius.reserve(3000);
    fID1.reserve(3000);
    fID2.reserve(3000);
    
    int partID;
    G4double x, y, z, radius;
    int ID1, ID2;
    
    if (f.is_open())
    {
        while (!f.eof())
        {
            getline(f,line);
            if (line[0] != '#')
            {
                std::istringstream iss(line);
                iss >> ID1 >> partID >> x >> y >> z >> radius >> ID2;
                fPositions.push_back(new G4ThreeVector(x*um, y*um, z*um));
                fPartID.push_back(partID);
                fRadius.push_back(radius * um);
                fID1.push_back(ID1);
                fID2.push_back(ID2);
            }
        }
    }
    else
    {
        G4cout << "ERROR: Unable to open file " << FileName << G4endl;
        exit(1);
    }
    f.close();
    
    TotalParts = fPartID.size();
}


G4VPhysicalVolume* TsNeuroMorpho::Construct()
{
	BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : Box
    //***********************************************************************
    
	G4Box* gWrapper = new G4Box("Wrapper", HLX, HLY, HLZ);
    
    fEnvelopeLog = CreateLogicalVolume(gWrapper);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //***********************************************************************
    //              Neuron Geometry
    //***********************************************************************
    
    //Pre-checks for the cell components:
    //Check number of soma components. This is usually 1, but in some cases can be more:
    
    SomaNumber = 0;
    
    for (G4int i = 0; i < TotalParts; i++){
        if (fPartID[i] == 1) {SomaNumber++;}
    }
    
    AxonNumber = 0;
    for (G4int i = 0; i < TotalParts; i++){
        if (fPartID[i] == 2) {AxonNumber++;}
    }
    
    BasalDendriteNumber = 0;
    for (G4int i = 0; i < TotalParts; i++){
        if (fPartID[i] == 3) {BasalDendriteNumber++;}
    }
    
    ApicalDendriteNumber = 0;
    for (G4int i = 0; i < TotalParts; i++){
        if (fPartID[i] == 4) {ApicalDendriteNumber++;}
    }
    
    G4cout << "The neuron model contains: " << G4endl;
    G4cout << SomaNumber << " Soma components" << G4endl;
    G4cout << AxonNumber << " Axon components" << G4endl;
    G4cout << BasalDendriteNumber << " Basal dendrite components" << G4endl;
    G4cout << ApicalDendriteNumber << " Apical dendrite components" << G4endl;
    
    //Cell components (read in from the neuroMorpho file)
    //Undefined components are not currently supported and will result in an error
    
    //**************************************
    //                  SOMA
    //**************************************
    
    G4RotationMatrix* rotSoma = new G4RotationMatrix();
    
    rotSoma->rotateX(0);
    rotSoma->rotateY(0);
    
    G4double somaRadius = fRadius[0];
    G4double x0 = fPositions[0]->x();
    G4double y0 = fPositions[0]->y();
    G4double z0 = fPositions[0]->z();
    
    //First line of file should define a soma component, if this is not the case,
    //the neuron is one of the SWC file special cases, which is not currently supported. Check this:
    if (fPartID[0] != 1){
        G4cout << "ERROR: Neuron does not contain a soma/cell body. These are special cases in the NeuroMorpho database, which we do not currently support." << G4endl;
        exit(1);
    }
    

    if (SomaNumber == 1){
        // Soma is a sphere:
        G4Orb* Soma = new G4Orb("Soma", somaRadius);
        G4LogicalVolume* SomaLog = CreateLogicalVolume("Soma", Soma);
        //G4VPhysicalVolume* SomaPhys = 
        CreatePhysicalVolume("Soma", SomaLog, rotSoma, fPositions[0], fEnvelopePhys);
        }
    else if (SomaNumber > 1){
        //Soma consists of more than one component - unionize the components (combine spheres):
        G4ThreeVector PosSoma(x0,y0,z0);
        
        G4Orb* Soma1 = new G4Orb("Soma1", fRadius[0]);
        G4Orb* Soma2 = new G4Orb("Soma2", fRadius[1]);
        G4UnionSolid* SomaUnion = new G4UnionSolid("SomaUnion", Soma1, Soma2, 0, PosSoma);
        
        //Used for overlap checking:
        fSomaRadii.reserve(10);
        fSomaRadii.push_back(fRadius[0]);
        fSomaRadii.push_back(fRadius[1]);
        
        for (G4int n = 2; n < SomaNumber; n++){
            
            G4double x0 = fPositions[n]->x();
            G4double y0 = fPositions[n]->y();
            G4double z0 = fPositions[n]->z();
            G4ThreeVector PosSoma(x0,y0,z0);
            
            G4Orb* Soma_n = new G4Orb("Soma", fRadius[n]);
            SomaUnion = new G4UnionSolid("SomaUnion", SomaUnion, Soma_n, 0, PosSoma);
            
            fSomaRadii.push_back(fRadius[n]);
            
            
        }
        
        G4LogicalVolume* SomaUnionLog = CreateLogicalVolume("SomaUnion", SomaUnion);
        //G4VPhysicalVolume* SomaUnionPhys = 
        CreatePhysicalVolume("SomaUnion", SomaUnionLog, rotSoma, fPositions[0], fEnvelopePhys);
    }

    //***************************************************
    // Axon and Dendtrites
    //***************************************************
    
    G4RotationMatrix* rotm = new G4RotationMatrix();
    
    rotm->rotateX(0);
    rotm->rotateY(0);
    
    G4int j = 0;
    G4int k = 0;
    G4int l = 0;
   
    G4int a = 0;
    G4int b = 0;
    
    for (G4int i = 0; i <= TotalParts-(SomaNumber); i++){
        
        //The first component is always part of the soma, in very rare cases,
        //a soma is not defined, we currently do not support these rare neuron types.
        if (fPartID[i] == 0){
            G4cerr << "Topas is exiting due to a serious error in reading in NeuroMorpho file." << G4endl;
            G4cerr << "SWC file contains an undefined component" << G4endl;
            exit(1);
            
        }
        
        //Treat glial proccesses the same as dendtrites
        if (fPartID[i] == 7){fPartID[i] = 2;}
        
        if ((fPartID[i]!=1)){
        
            //Find the components of the axons and dendtrites:
            a = fID2[i]-1;
            b = fID1[i]-1;
        
            //Calculate cylinder length using distance formula:
            rx = (fPositions[b]->x()-fPositions[a]->x());
            ry = (fPositions[b]->y()-fPositions[a]->y());
            rz = (fPositions[b]->z()-fPositions[a]->z());
            r = sqrt(rx*rx + ry*ry + rz*rz);
        
            //The position of the cylinder (midpoint between two points)
            mx = (fPositions[a]->x()+fPositions[b]->x())/2;
            my = (fPositions[a]->y()+fPositions[b]->y())/2;
            mz = (fPositions[a]->z()+fPositions[b]->z())/2;
        
            //Find the direction of the cylinder (unit vector b):
            bx = (fPositions[b]->x()-fPositions[a]->x())/r;
            by = (fPositions[b]->y()-fPositions[a]->y())/r;
            bz = (fPositions[b]->z()-fPositions[a]->z())/r;
        
            //Rotate the cylinder (0,0,1) to lay along the vector b.
            //Calculate the matrix elements using the Rodrigues' Rotation Formula:
        
            G4double p = bx*bx + by*by;
        
            //Matrix elements are:
            mxx = bz + (by*by*(1-bz))/p;
            mxy = -1*bx*by*(1-bz)/p;
            mxz = bx*sqrt(1-bz*bz)/sqrt(p);
            myx = -1*bx*by*(1-bz)/p;
            myy = bz+(bx*bx*(1-bz)/p);
            myz = by*sqrt(1-bz*bz)/sqrt(p);
            mzx = -1*bx*sqrt(1-bz*bz)/sqrt(p);
            mzy = -1*by*sqrt(1-bz*bz)/sqrt(p);
            mzz = bz;
        
            //Special cases where the formula fails:
            if ( (bx == 0) && (by == 0) && (bz == 1) ){
            
                mxx = 1;
                mxy = 0;
                mxz = 0;
                myx = 0;
                myy = 1;
                myz = 0;
                mzx = 0;
                mzy = 0;
                mzz = 1;
            }
        
            if ( (bx == 0) && (by == 0) && (bz == -1) ){
            
                mxx = 1;
                mxy = 0;
                mxz = 0;
                myx = 0;
                myy = -1;
                myz = 0;
                mzx = 0;
                mzy = 0;
                mzz = -1;
            }
        
            //Define matrix vectors for rotation matrix (colX, colY, colZ)
            G4ThreeVector colX = G4ThreeVector(mxx, myx, mzx);
            G4ThreeVector colY = G4ThreeVector(mxy, myy, mzy);
            G4ThreeVector colZ = G4ThreeVector(mxz, myz, mzz);
            G4RotationMatrix rotDend = G4RotationMatrix(colX, colY, colZ);
            G4RotationMatrix* rotDendInv = new G4RotationMatrix(rotDend.inverse());
            
            //Check if component is in the soma volume
            
            G4bool InsideSoma = false;
            
            if (SomaNumber == 1){
                if ((mx-x0)*(mx-x0) + (my-y0)*(my-y0) + (mz-z0)*(mz-z0) < somaRadius*somaRadius)
                {InsideSoma = true;}
            }
            
            if (SomaNumber > 1){
                for (unsigned int jp = 0; jp < fSomaRadii.size(); jp++){
                    if ((mx-x0)*(mx-x0) + (my-y0)*(my-y0) + (mz-z0)*(mz-z0) < somaRadius*somaRadius){InsideSoma = true;}
                }
            }
            
            //If no overlaps, place volumes
            if (InsideSoma != true){
    
                //Axon component (cylinder)
                if (fPartID[i] == 2){
                    
                    G4ThreeVector* CylPos = new G4ThreeVector(mx, my, mz);
                
                    G4Tubs* Axon = new G4Tubs("Axon", 0, fRadius[i], r/2, 0*deg, 360*deg);
                    G4LogicalVolume* AxonLog = CreateLogicalVolume("Axon", Axon);
                    //G4VPhysicalVolume* AxonPhys = 
                    CreatePhysicalVolume("Axon", j, false, AxonLog, rotDendInv, CylPos, fEnvelopePhys);
                    j++;
                }
          
        
                //Basal dendrite component (cylinder)
                if (fPartID[i] == 3){
    
                    G4ThreeVector* CylPos = new G4ThreeVector(mx, my, mz);

                    G4Tubs* BDend = new G4Tubs("BasalDendrite", 0, fRadius[i], r/2, 0*deg, 360*deg);
                    G4LogicalVolume* BDendLog = CreateLogicalVolume("BasalDendrite", BDend);
                    //G4VPhysicalVolume* BDendPhys = 
                    CreatePhysicalVolume("BasalDendrite", k, false, BDendLog, rotDendInv, CylPos, fEnvelopePhys);
                    k++;
                }
        
                //Apical dendrite component (cylinder)
                if (fPartID[i] == 4){
                  
                    G4ThreeVector* CylPos = new G4ThreeVector(mx, my, mz);
                  
                    G4Tubs* ADend = new G4Tubs("ApicalDendrite", 0, fRadius[i], r/2, 0*deg, 360*deg);
                    G4LogicalVolume* ADendLog = CreateLogicalVolume("ApicalDendrite", ADend);
                    //G4VPhysicalVolume* ADendPhys = 
                    CreatePhysicalVolume("ApicalDendrite", l, false, ADendLog, rotDendInv, CylPos, fEnvelopePhys);
                    l++;
                  
            
                    }
                }
            }
        }
 
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
    













        
