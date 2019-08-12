// Component for TsChannel
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

#include "TsChannel.hh"

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

TsChannel::TsChannel(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
     ResolveParameters();
}


TsChannel::~TsChannel()
{;}

void TsChannel::ResolveParameters() {
    
    fLipidHeadRadius = fPm->GetDoubleParameter(GetFullParmName("LipidHeadRadius"), "Length");
    fLipidTailLength = fPm->GetDoubleParameter(GetFullParmName("LipidTailHalfLength"), "Length");
    fLipidTailRadius = fPm->GetDoubleParameter(GetFullParmName("LipidTailRadius"), "Length");
    
    fLipidLength = (fLipidHeadRadius * 4) + (fLipidTailLength * 2);
    
    fNumberOfRows = fPm->GetIntegerParameter(GetFullParmName("NumberOfRows"));
    fNumberOfCols = fPm->GetIntegerParameter(GetFullParmName("NumberOfColumns"));
    
    fChannelRadius = fPm->GetDoubleParameter(GetFullParmName("ChannelRmax"), "Length");
    fChannelInnerRadius = fPm->GetDoubleParameter(GetFullParmName("ChannelRmin"), "Length");
    fChannelLength = fPm->GetDoubleParameter(GetFullParmName("ChannelHalfLength"), "Length");
    
}


G4VPhysicalVolume* TsChannel::Construct()
{
	BeginConstruction();
    
    //****************************************************************************
    //                             Box (envelope)
    //****************************************************************************
    G4double HX = fLipidHeadRadius*fNumberOfRows + fLipidHeadRadius;
    G4double HY = fChannelLength;
    G4double HZ = fLipidHeadRadius*fNumberOfCols + fLipidHeadRadius;
    
    G4Box* gWrap = new G4Box("Wrapper", HX, HY, HZ);
    fEnvelopeLog = CreateLogicalVolume(gWrap);
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
    //          Subcomponent 2: simple channel
    //*************************************************************************
    
    G4Tubs* gChannel = new G4Tubs("Channel", fChannelInnerRadius, fChannelRadius, fChannelLength, 0*deg, 360*deg);
    G4LogicalVolume* lChannel = CreateLogicalVolume("Channel", gChannel);
    
    G4ThreeVector* positionIon = new G4ThreeVector(0,0,0);
    G4RotationMatrix* rotIon = new G4RotationMatrix();
    rotIon->rotateX(90*deg);
    
    CreatePhysicalVolume("Channel", lChannel, rotIon, positionIon, fEnvelopePhys);
    

    //*************************************************************************
    //              Place components in box
    //*************************************************************************
    
    G4RotationMatrix* rot1 = new G4RotationMatrix();
    G4RotationMatrix* rot2 = new G4RotationMatrix();
    G4ThreeVector* position_head1 = new G4ThreeVector(0,fLipidTailLength+fLipidHeadRadius,0);
    G4ThreeVector* position_head2 = new G4ThreeVector(0,-fLipidTailLength-fLipidHeadRadius,0);
    G4ThreeVector* position_tail = new G4ThreeVector(0,0,0);
    G4double angle = 90*deg;
    rot1->rotateX(0);
    rot2->rotateX(angle);
    
    
   //*************************************************************************
   //Place in a layer
   //*************************************************************************
    
    G4double x_row = 0.0;
    G4double y_const = 0.0;
    G4double z_col = 0.0;
    
    G4double start_x = fLipidHeadRadius*fNumberOfRows;
    G4double start_z = fLipidHeadRadius*fNumberOfCols;
    
    G4int row, col;
    
    row = 1;
    col = 1;
    G4double shiftx = 2*fLipidHeadRadius;
    G4double shiftz = 2*fLipidHeadRadius;
    
    G4int TotalNumber = fNumberOfCols*fNumberOfRows+1;
    
    for (int i = 1; i < TotalNumber; i++){
        
        if (i % fNumberOfCols != 0){
            x_row = start_x - shiftx*row;
            z_col = start_z - shiftz*col;
            col++;
            
        }
        else{
            x_row = start_x - shiftx*row;
            z_col = start_z - shiftz*col;
            row++;
            col = 1;
        }
        
        
        G4ThreeVector* position = new G4ThreeVector(x_row, y_const, z_col);
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot ->rotateX(0);
        
        //Create a hole for the channel:
        
        G4double radius_o = std::sqrt(x_row*x_row + z_col*z_col) - 2*fLipidHeadRadius;
        
        if (radius_o > fChannelRadius){
        
            G4VPhysicalVolume* pBox = CreatePhysicalVolume("Box", i, true, lBox, rot, position, fEnvelopePhys);
            
            if ( i == 1 ) {
                CreatePhysicalVolume("LipidHead1", lhead1, rot1, position_head1, pBox);
                CreatePhysicalVolume("LipidHead2", lhead2, rot1, position_head2, pBox);
                CreatePhysicalVolume("LipidTail", lTail, rot2, position_tail, pBox);
            }
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
