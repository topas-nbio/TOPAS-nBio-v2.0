// Component for TsLinearDNA
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
// A strand of DNA.

#include "TsLinearDNA.hh"

#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "G4VisAttributes.hh"

TsLinearDNA::TsLinearDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsLinearDNA::~TsLinearDNA()
{;}

void TsLinearDNA::ResolveParameters() {
    fNumberOfBasePairs = fPm->GetIntegerParameter(GetFullParmName("NumberOfBasePairs"));
    HL = (0.34 * nm * fNumberOfBasePairs)/2;
}



G4VPhysicalVolume* TsLinearDNA::Construct()
{
	BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : Cylinder containing strand
    //***********************************************************************
    
    G4Tubs* gCyl = new G4Tubs("gCyl",
                              0*nm,
                              1.5*nm,
                              HL,
                              0*deg,
                              360*deg);
    
    fEnvelopeLog = CreateLogicalVolume(gCyl);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //**************************************************************************
    //                 Subcomponent 1: Base pair
    //**************************************************************************
    
    //Base pair - a cylinder of radius 0.5 nm and length 0.34 nm.
    
    G4String Subcomponent1 = "BasePair";
    G4Tubs* gBp = new G4Tubs(Subcomponent1, 0, 0.5*nm, 0.17*nm, 0.0*deg, 360.0*deg);
    G4LogicalVolume* lBp = CreateLogicalVolume(Subcomponent1, gBp);
    
    //**************************************************************************
    //                 Subcomponent 2: Sugar phosphate backbone
    //**************************************************************************
    
    //Phosphodiester group - two sugars each consisting of quarter cylinders
    //The sugars are wrapped around the base pair
    
    G4String Subcomponent2 = "Backbone1";
    G4Tubs* gSugarPhosphate1 = new G4Tubs(Subcomponent2, 0.5*nm, 1.185*nm, 0.17*nm, 0*deg, 45*deg);
    G4LogicalVolume* lSugarPhosphate1 = CreateLogicalVolume(Subcomponent2, gSugarPhosphate1);
    
    G4String Subcomponent3 = "Backbone2";
    G4Tubs* gSugarPhosphate2 = new G4Tubs(Subcomponent3, 0.5*nm, 1.185*nm, 0.17*nm, 180*deg, 45*deg);
    G4LogicalVolume* lSugarPhosphate2 = CreateLogicalVolume(Subcomponent3, gSugarPhosphate2);
    
    // Rotation of strands around the base pair
    
    G4double x = 0.0;
    G4double y = 0.0;
    G4double z0 = -HL + 0.17*nm;
    
    for (G4int j = 0; j < fNumberOfBasePairs ; j++){
        
        G4double theta = 36*deg*j;
        G4double z = z0 + j*0.34*nm;
        
        G4ThreeVector* position = new G4ThreeVector(x, y, z);
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot ->rotateZ(theta);
        
        //G4VPhysicalVolume* pBP = 
        CreatePhysicalVolume(Subcomponent1, j, true, lBp, rot, position, fEnvelopePhys);
        //G4VPhysicalVolume* pSugar1 = 
        CreatePhysicalVolume(Subcomponent2, j, true, lSugarPhosphate1, rot, position, fEnvelopePhys);
        //G4VPhysicalVolume* pSugar2 = 
        CreatePhysicalVolume(Subcomponent3, j, true, lSugarPhosphate2, rot, position, fEnvelopePhys);
        
    }
    
    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
