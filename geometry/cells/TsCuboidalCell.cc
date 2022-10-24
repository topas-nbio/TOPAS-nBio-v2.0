// Component for TsCuboidalCell
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
// Cube or rectanglar cell.
// User has the option of including organelles: nucleus and/or mitochondria.

#include "TsCuboidalCell.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsCuboidalCell::TsCuboidalCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsCuboidalCell::~TsCuboidalCell()
{;}


void TsCuboidalCell::ResolveParameters() {
    HLX = fPm->GetDoubleParameter(GetFullParmName("Cell_HLX"), "Length");
    HLY = fPm->GetDoubleParameter(GetFullParmName("Cell_HLY"), "Length");
    HLZ = fPm->GetDoubleParameter(GetFullParmName("Cell_HLZ"), "Length");
}


G4VPhysicalVolume* TsCuboidalCell::Construct()
{
	BeginConstruction();

    // Used for defining the area that mitochondria are distributed in
    G4double CellRadius = HLX;
    
    //***********************************************************************
    //              Envelope Geometry : cuboidal/columnar cell
    //***********************************************************************
    
    G4Box* gCell = new G4Box(fName, HLX, HLY, HLZ);
    fEnvelopeLog = CreateLogicalVolume(gCell);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************
    
    //***************************
    // Subcomponent: Nucleus
    //***************************
    
    G4double NucleusRadius = 0.0*um;
    G4String name = GetFullParmName("Nucleus/NucleusRadius");
    if (fPm->ParameterExists(name)) {
        NucleusRadius = fPm->GetDoubleParameter(name, "Length");
        G4String subComponentName1 = "Nucleus";
        
        G4RotationMatrix* rotNuc = new G4RotationMatrix();
        
        rotNuc->rotateX(0);
        rotNuc->rotateY(0);
        
        G4double transNucX = 0 * um;
        G4double transNucY = 0 * um;
        G4double transNucZ = 0 * um;
        
        G4String name1 = GetFullParmName("Nucleus/transNucX");
        if (fPm -> ParameterExists(name1)){
            transNucX = fPm->GetDoubleParameter(name1, "Length");
        }
        
        name1 = GetFullParmName("Nucleus/transNucY");
        if (fPm -> ParameterExists(name1)){
            transNucY = fPm->GetDoubleParameter(name1, "Length");
        }
        
        name1 = GetFullParmName("Nucleus/transNucZ");
        if (fPm -> ParameterExists(name1)){
            transNucZ = fPm->GetDoubleParameter(name1, "Length");
        }
        
        
        G4ThreeVector* NucPos = new G4ThreeVector(transNucX,transNucY,transNucZ);
        
        G4Orb* gNucleus = new G4Orb("gNucleus", NucleusRadius);
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
        
        //Semi-axis lengths of the ellpsoid
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
        G4Ellipsoid* gMito = new G4Ellipsoid(subComponentName2, EllA, EllB, EllC);
        G4LogicalVolume* lMito = CreateLogicalVolume(subComponentName2, gMito);
        
        //Randomly distribute mitochondria throughout cell volume outside nucleus (default)
        for (int j = 0; j < NbOfMito; j++){
            
            G4bool Overlap = true;
            while (Overlap == true){
                
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(CellRadius - NucleusRadius);
                G4double phi = G4UniformRand()*2*pi;
                G4double psi = G4UniformRand()*2*pi;
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
                
                x = (NucleusRadius + dr)* std::cos(u) * std::sin(v);
                y = (NucleusRadius + dr)* std::sin(u) * std::sin(v);
                z = (NucleusRadius + dr)* std::cos(v);
                
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
