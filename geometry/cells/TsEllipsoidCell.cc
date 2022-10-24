// Component for TsEllipsoidCell
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
// An ellipsoid cell.
// User has the option of including organelles: nucleus and/or mitochondria.

#include "TsEllipsoidCell.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsEllipsoidCell::TsEllipsoidCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsEllipsoidCell::~TsEllipsoidCell()
{;}

void TsEllipsoidCell::ResolveParameters() {
    
    G4String name = GetFullParmName("xSemiAxis");
    if (fPm->ParameterExists(name)) {
        //User specified cell size values
        xSA = fPm->GetDoubleParameter(GetFullParmName("xSemiAxis"), "Length");
        ySA = fPm->GetDoubleParameter(GetFullParmName("ySemiAxis"), "Length");
        zSA = fPm->GetDoubleParameter(GetFullParmName("zSemiAxis"), "Length");
        
    }else{
        //Set default values for cell dimensions
        xSA = 20 * um;
        ySA = 10 * um;
        zSA = 15 * um;
    }
    
}


G4VPhysicalVolume* TsEllipsoidCell::Construct()
{
	BeginConstruction();

    //***********************************************************************
    //              Envelope Geometry : ellipsoid cell
    //***********************************************************************
    
	G4Ellipsoid* envelopeSolid = new G4Ellipsoid(fName, xSA, ySA, zSA);
	fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************
    
    //Find biggest semi-axis length
    G4double CellRadius;
    if ((xSA <= ySA) && (xSA <= zSA)){CellRadius = xSA;}
    if ((ySA <= xSA) && (ySA <= zSA)){CellRadius = ySA;}
    if ((zSA <= xSA) && (zSA <= ySA)){CellRadius = zSA;}
    
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
            if (transNucX > xSA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                fPm->AbortSession(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/translateY");
        if (fPm -> ParameterExists(name1)){
            transNucY = fPm->GetDoubleParameter(name1, "Length");
            if (transNucY > ySA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                fPm->AbortSession(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/translateZ");
        if (fPm -> ParameterExists(name1)){
            transNucZ = fPm->GetDoubleParameter(name1, "Length");
            if (transNucZ > zSA) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                fPm->AbortSession(1);
            }
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
