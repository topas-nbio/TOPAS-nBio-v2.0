// Component for TsMembraneRing
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
// Lipid bilayer in a ring

#include "TsMembraneRing.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4Color.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "TsParameterManager.hh"
#include "G4VisAttributes.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"

TsMembraneRing::TsMembraneRing(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsMembraneRing::~TsMembraneRing()
{;}

void TsMembraneRing::ResolveParameters() {
    
    fNumberOfLipids = fPm->GetIntegerParameter(GetFullParmName("NumberOfLipids"));
    fLipidHeadRadius = fPm->GetDoubleParameter(GetFullParmName("LipidHeadRadius"), "Length");
    fLipidTailLength = fPm->GetDoubleParameter(GetFullParmName("LipidTailHalfLength"), "Length");
    fLipidTailRadius = fPm->GetDoubleParameter(GetFullParmName("LipidTailRadius"), "Length");
    
    fLipidLength = (fLipidHeadRadius * 4) + (fLipidTailLength * 2);
    
    fRMin = (fLipidHeadRadius * 2) * (fNumberOfLipids) / twopi;
    fRMax = fRMin + (fLipidHeadRadius * 2);
    
}


G4VPhysicalVolume* TsMembraneRing::Construct()
{
	BeginConstruction();
    
    //****************************************************************************
    //                             Ring (envelope)
    //****************************************************************************
    G4Tubs* envelope = new G4Tubs(fName, fRMin, fRMax, fLipidLength/2, 0.0*deg, 360*deg);
    fEnvelopeLog = CreateLogicalVolume(envelope);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //****************************************************************************
    //                             Boxes for lipids
    //****************************************************************************
    
    G4Box* gBox = new G4Box("Lipidbox", fLipidHeadRadius, fLipidLength/2, fLipidHeadRadius);
    G4LogicalVolume* lBox = CreateLogicalVolume(gBox);
    
    //**************************************************************************
    //                 Subcomponent 1: Lipid bilayer
    //**************************************************************************

    //Consists of two heads (two spheres) and a tail (cylinder)
    
    G4Orb* ghead1 = new G4Orb("LipidHead1", fLipidHeadRadius);
    G4LogicalVolume* lhead1 = CreateLogicalVolume("LipidHead1", ghead1);
    
    G4Orb* ghead2 = new G4Orb("LipidHead2", fLipidHeadRadius);
    G4LogicalVolume* lhead2 = CreateLogicalVolume("LipidHead2", ghead2);
    
    G4Tubs* gTail = new G4Tubs("LipidTail", 0*nm, fLipidTailRadius, fLipidTailLength, 0*deg, 360*deg);
    G4LogicalVolume* lTail = CreateLogicalVolume("LipidTail", gTail);
    
    //*************************************************************************
    //              Place components in ring
    //*************************************************************************
    
    //Postions of lipid components:
    G4RotationMatrix* rot1 = new G4RotationMatrix();
    G4RotationMatrix* rot2 = new G4RotationMatrix();
    G4ThreeVector* position_head1 = new G4ThreeVector(0,fLipidTailLength+fLipidHeadRadius,0);
    G4ThreeVector* position_head2 = new G4ThreeVector(0,-fLipidTailLength-fLipidHeadRadius,0);
    G4ThreeVector* position_tail = new G4ThreeVector(0,0,0);
    G4double angle = 90*deg;
    rot1->rotateX(0);
    rot2->rotateX(angle);
    
    //Create ring and place components within the ring:
    for (G4int j = 0; j < fNumberOfLipids; j++) {
        G4double phi;
        phi = 360.*deg/fNumberOfLipids*j;
        
        G4RotationMatrix rotm  = G4RotationMatrix();
        rotm.rotateX(90*deg);
        rotm.rotateY(0*deg*j);
        rotm.rotateZ(phi);
        G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
        G4ThreeVector position = (fRMin+1.2*nm)*uz;
        G4Transform3D transform = G4Transform3D(rotm,position);
        
        //place boxes within the ring
        G4String aName = fName + "_" + G4UIcommand::ConvertToString(j);
        G4VPhysicalVolume* pBox = CreatePhysicalVolume("Box",
                                                        j, true, lBox,
                                                        transform,
                                                        fEnvelopePhys);
        if ( j == 0 ) {
            CreatePhysicalVolume("LipidHead1", lhead1, rot1, position_head1, pBox);
            CreatePhysicalVolume("LipidHead2", lhead2, rot1, position_head2, pBox);
            CreatePhysicalVolume("LipidTail", lTail, rot2, position_tail, pBox);
        }
    }
    
    //Vis attributes
    
    G4VisAttributes* wrapVisAtt = new G4VisAttributes(G4Color(0.5,0.5,0.5));
    wrapVisAtt->SetVisibility(false);
    lBox->SetVisAttributes(wrapVisAtt);
    fEnvelopeLog->SetVisAttributes(wrapVisAtt);
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
