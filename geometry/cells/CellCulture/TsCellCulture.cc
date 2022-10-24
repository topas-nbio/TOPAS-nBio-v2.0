// Component for TsCellCulture
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


// A simple cell culture consisting of random spherical cells.


#include "TsCellCulture.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsCellCulture::TsCellCulture(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}

TsCellCulture::~TsCellCulture()
{;}

void TsCellCulture::ResolveParameters() {
    HLX = fPm->GetDoubleParameter(GetFullParmName("Container_HLX"), "Length");
    HLY = fPm->GetDoubleParameter(GetFullParmName("Container_HLY"), "Length");
    HLZ = fPm->GetDoubleParameter(GetFullParmName("Container_HLZ"), "Length");
    
    CellRadius = fPm->GetDoubleParameter(GetFullParmName("CellRadius"), "Length");
    NbOfCells  = fPm->GetIntegerParameter(GetFullParmName("NumberOfCells"));
    NuclRadius = fPm->GetDoubleParameter(GetFullParmName("NucleusRadius"), "Length");
}


G4VPhysicalVolume* TsCellCulture::Construct()
{
    BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : Rectanglar container
    //***********************************************************************
    
    G4Box* gBox = new G4Box(fName, HLX, HLY, HLZ);
    fEnvelopeLog = CreateLogicalVolume(gBox);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //***********************************************************************
    //              Cell geometry : spherical
    //***********************************************************************
    //Cell geometry
    G4Orb* gCell = new G4Orb("cell", CellRadius);
    G4LogicalVolume* lCell = CreateLogicalVolume(gCell);
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************

    // Nucleus
    G4String subComponentName1 = "Nucleus";
    G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);

    // Randomly place cells in the volume
    for (int j = 0; j < NbOfCells; j++){

        G4bool Overlap = true;
        while (Overlap == true){

            G4double phi = 0;
            G4double psi = 0;
            G4double x = 0.0;
            G4double y = 0.0;
            G4double z = 0.0;

            x = (2*G4UniformRand()-1)*(HLX-CellRadius) ;
            y = (2*G4UniformRand()-1)*(HLY-CellRadius) ;
            z = (2*G4UniformRand()-1)*(HLZ-CellRadius) ;

            G4ThreeVector* position = new G4ThreeVector(x,y,z);
            G4ThreeVector* posNucl = new G4ThreeVector(0*mm,0*mm,0*mm);
            
            G4RotationMatrix* rotm = new G4RotationMatrix();

            rotm->rotateX(psi);
            rotm->rotateY(phi);

            G4VPhysicalVolume* pCell = CreatePhysicalVolume("Cell", j, true, lCell, rotm, position, fEnvelopePhys);
            G4VPhysicalVolume* pNucleus = CreatePhysicalVolume("Nucleus", j, true, lNucleus, rotm, posNucl, pCell);
        
            G4bool OverlapCheck = pCell->CheckOverlaps();

            if (OverlapCheck == false){
                OverlapCheck = pNucleus->CheckOverlaps();
                if (OverlapCheck == false)
                    break;
            }
            if (OverlapCheck == true){
                pCell = NULL;
                pNucleus = NULL;
                G4cout << "**** Finding new position for volume Cell : " << j <<  " ****" << G4endl;
            }
        }
    }

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
