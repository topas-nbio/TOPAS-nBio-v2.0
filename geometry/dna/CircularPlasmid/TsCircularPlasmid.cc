// Component for TsCircularPlasmid
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
// Geometry for a circular DNA plasmid, user specifies the number of base pairs in the ring.
// The plasmid consists of a ring (containing all geometric components) and
// boxes arranged within the ring (each containing the base pair component).
// Each DNA segment consists of a spherical base pair and the surrounding sugar backbone (2 quarter spheres)


#include "TsCircularPlasmid.hh"
#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcommand.hh"

#include "G4VisAttributes.hh"


TsCircularPlasmid::TsCircularPlasmid(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                     TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsCircularPlasmid::~TsCircularPlasmid()
{;}


void TsCircularPlasmid::ResolveParameters() {
	fNumberOfBasePairs = fPm->GetIntegerParameter(GetFullParmName("NumberOfBasePairs"));
	fRMin = (0.34 * nm) * (fNumberOfBasePairs) / twopi;
	fRMax = fRMin + 2.4 * nm;
}


G4VPhysicalVolume* TsCircularPlasmid::Construct()
{
    BeginConstruction();
    
    //****************************************************************************
    //                             Ring for plasmid (envelope)
    //****************************************************************************
	
	G4Tubs* envelope = new G4Tubs(fName, fRMin, fRMax, 1.2*nm, 0.0*deg, 360*deg);
    fEnvelopeLog = CreateLogicalVolume(envelope);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	
    //****************************************************************************
    //                             Boxes for base pairs
    //****************************************************************************

    G4Box* gBox = new G4Box("Base", 1.185*nm, 1.185*nm, 0.17*nm);
	G4LogicalVolume* lBox = CreateLogicalVolume(gBox);
    
    //**************************************************************************
    //                 Subcomponent 1: Base pair
    //**************************************************************************
    //Base pair - a cylinder of radius 0.5 nm and length 0.34 nm.
    G4String subComponent1 = "BasePair";
    G4Tubs* gBp1 = new G4Tubs(subComponent1, 0, 0.5*nm, 0.17*nm, 0.0*deg, 360.0*deg);
    G4LogicalVolume* lBp1 = CreateLogicalVolume(subComponent1, gBp1);
	
    //**************************************************************************
    //                 Subcomponent 1: Sugar phosphate
    //**************************************************************************
    //Phosphodiester group - two sugars each consisting of quarter cylinders
    //The sugars are wrapped around the base pair
    G4String subComponent2 = "Backbone1";
    G4Tubs* gSugarPhosphate1 = new G4Tubs(subComponent2, 0.5*nm, 1.185*nm, 0.17*nm, 0*deg, 90*deg);
    G4LogicalVolume* lSugarPhosphate1 = CreateLogicalVolume(subComponent2, gSugarPhosphate1);
	
    G4String subComponent3 = "Backbone2";
    G4Tubs* gSugarPhosphate2 = new G4Tubs(subComponent3, 0.5*nm, 1.185*nm, 0.17*nm, 180*deg, 90*deg);
    G4LogicalVolume* lSugarPhosphate2 = CreateLogicalVolume(subComponent3, gSugarPhosphate2);
	
    //*************************************************************************
    //              Place components in ring
    //*************************************************************************
	
    for (G4int j = 0; j < fNumberOfBasePairs ; j++) {
        G4double phi;
        phi = 360.*deg/fNumberOfBasePairs*j;
		
        G4RotationMatrix rotm  = G4RotationMatrix();
        rotm.rotateX(90*deg);
        rotm.rotateY(36*deg*j);
        rotm.rotateZ(phi);
        G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
        G4ThreeVector position = (fRMin+1.2*nm)*uz;
        G4Transform3D transform = G4Transform3D(rotm,position);
        
        //place boxes within the ring
        G4String aName = fName + "_" + G4UIcommand::ConvertToString(j);
		G4VPhysicalVolume* pBasePair = CreatePhysicalVolume("Base",
															j, true, lBox,
															transform,
															fEnvelopePhys);
		if ( j == 0 ) {
		  CreatePhysicalVolume("BasePair", lBp1, pBasePair);
		  CreatePhysicalVolume("Backbone1", lSugarPhosphate1, pBasePair);
		  CreatePhysicalVolume("Backbone2", lSugarPhosphate2, pBasePair);
		}
    }
    
    
    G4VisAttributes* BoxVis = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    BoxVis->SetVisibility(false);
    lBox->SetVisAttributes(BoxVis);
    
    InstantiateChildren(fEnvelopePhys);
    
    return fEnvelopePhys;
}



