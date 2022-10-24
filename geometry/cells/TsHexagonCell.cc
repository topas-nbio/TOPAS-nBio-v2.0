// Component for TsHexagonCell
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
// Irregular-shaped cell (3D hexagon).
// User has the option of including organelles: nucleus and/or mitochondria.
//

#include "TsHexagonCell.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"


TsHexagonCell::TsHexagonCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{;}


TsHexagonCell::~TsHexagonCell()
{;}


G4VPhysicalVolume* TsHexagonCell::Construct()
{
	BeginConstruction();
    
    //***********************************************************************
    //              Cell Geometry : Hexagon cell
    //***********************************************************************
    
    //Base shape:
    std::vector<G4TwoVector> cellpoly(6);
    cellpoly[0] = G4TwoVector(15*micrometer, 15*micrometer);
    cellpoly[1] = G4TwoVector(-15*micrometer, 15*micrometer);
    cellpoly[2] = G4TwoVector(-20*micrometer, 0*micrometer);
    cellpoly[3] = G4TwoVector(-15*micrometer, -15*micrometer);
    cellpoly[4] = G4TwoVector(15*micrometer,-15*micrometer);
    cellpoly[5] = G4TwoVector(20*micrometer, 0*micrometer);
    
    // Z height
    G4double hz = 10*micrometer;  //half length along z
    
    G4double CellRadius = 10*micrometer;
    
    
	G4ExtrudedSolid* gFibroCell = new G4ExtrudedSolid(fName,
                                                      cellpoly,
                                                      hz,
                                                      G4TwoVector(0,0),
                                                      0.5,
                                                      G4TwoVector(0,0),
                                                      1.0); // shape along z
    
	fEnvelopeLog = CreateLogicalVolume(gFibroCell);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************
    
    //***************************
    // Subcomponent: Nucleus
    //***************************
    
    G4double NuclRadius = 0.0*um;
    G4String name = GetFullParmName("Nucleus/NucleusRadius");
    if (fPm->ParameterExists(name)) {
        NuclRadius = fPm->GetDoubleParameter(name, "Length");
        G4String subComponentName1 = "Nucleus";
        
        G4RotationMatrix* rotNuc = new G4RotationMatrix();
        
        rotNuc->rotateX(0);
        rotNuc->rotateY(0);
        
        G4double transNucX = 0 * um;
        G4double transNucY = 0 * um;
        G4double transNucZ = 0 * um;
        
        G4String name1 = GetFullParmName("Nucleus/translateX");
        if (fPm -> ParameterExists(name1)){
            transNucX = fPm->GetDoubleParameter(name1, "Length");
        }
        
        name1 = GetFullParmName("Nucleus/translateY");
        if (fPm -> ParameterExists(name1)){
            transNucY = fPm->GetDoubleParameter(name1, "Length");
        }
        
        name1 = GetFullParmName("Nucleus/translateZ");
        if (fPm -> ParameterExists(name1)){
            transNucZ = fPm->GetDoubleParameter(name1, "Length");
        }
        
        
        G4ThreeVector* NucPos = new G4ThreeVector(transNucX,transNucY,transNucZ);
        
        G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
        G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
        G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, rotNuc, NucPos, fEnvelopePhys);
        
        G4bool OverlapCheck = pNucleus->CheckOverlaps();
        
        if (OverlapCheck == true){
            G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
            G4cerr << "Nucleus overlaps with the cell." << G4endl;
            fPm->AbortSession(1);
        }
        
    }
    
    
    //*******************************
    // Subcomponent: Mitochondria
    //*******************************
    
    name = GetFullParmName("Mitochondria/NumberOfMitochondria");
    if (fPm->ParameterExists(name)) {
        
        //number of mitochondria
        const G4int NbOfMito  = fPm->GetIntegerParameter( GetFullParmName("Mitochondria/NumberOfMitochondria") );
        
        //Semi-axis lengths of the ellpsoid/mitochondria (default values if none are specified)
        G4double EllA = 0.5*micrometer;
        G4double EllB = 0.3*micrometer;
        G4double EllC = 0.9*micrometer;
        
        name=GetFullParmName("Mitochondria/a");
        if (fPm->ParameterExists(name)){EllA = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/a"), "Length" );}
        
        name=GetFullParmName("Mitochondria/b");
        if (fPm->ParameterExists(name)){EllB = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/b"), "Length" );}
        
        name=GetFullParmName("Mitochondria/c");
        if (fPm->ParameterExists(name)){EllC = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/c"), "Length" );}
        
        
        G4String subComponentName2 = "Mitochondria";
        G4Ellipsoid* gMito = new G4Ellipsoid("gMito", EllA, EllB, EllC);
        G4LogicalVolume* lMito = CreateLogicalVolume(subComponentName2, gMito);
        
        //Randomly distribute mitochondria throughout cell volume
        for (int j = 0; j < NbOfMito; j++){
            
            G4bool Overlap = true;
            while (Overlap == true){
                
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(CellRadius - NuclRadius);
                G4double phi = G4UniformRand()*2*pi;
                G4double psi = G4UniformRand()*2*pi;
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
                
                x = (NuclRadius + dr)* std::cos(u) * std::sin(v);
                y = (NuclRadius + dr)* std::sin(u) * std::sin(v);
                z = (NuclRadius + dr)* std::cos(v);
                
                G4ThreeVector* position = new G4ThreeVector(x,y,z);
                
                G4RotationMatrix* rotm = new G4RotationMatrix();
                
                rotm->rotateX(psi);
                rotm->rotateY(phi);
                
                G4VPhysicalVolume* pMito = CreatePhysicalVolume(subComponentName2, j, true, lMito, rotm, position, fEnvelopePhys);
                
                G4bool OverlapCheck = pMito->CheckOverlaps();
                
                if (OverlapCheck == false){break;}
                if (OverlapCheck == true){
                    G4PhysicalVolumeStore::DeRegister(pMito);
                    G4cout << "**** Finding new position for volume " << subComponentName2 << ":" << j <<  " ****" << G4endl;
                }
            }
        }
    }

    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
