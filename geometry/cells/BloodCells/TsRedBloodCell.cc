// Component for TsRedBloodCell
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
// A single red blood cell (RBC).
// Red blood cells do not contain a nucleus or organelles but store hemoglobin, an oxygen binding protein.

#include "TsRedBloodCell.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

TsRedBloodCell::TsRedBloodCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsRedBloodCell::~TsRedBloodCell()
{;}

void TsRedBloodCell::ResolveParameters() {
    RedBloodCellRadius = fPm->GetDoubleParameter(GetFullParmName("RBCRadius"), "Length");
    
    TorusWidth = fPm->GetDoubleParameter(GetFullParmName("RBCWidth"), "Length");
    
}


G4VPhysicalVolume* TsRedBloodCell::Construct()
{
	BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : spherical cell
    //***********************************************************************
    
    G4Torus* solidTorus = new G4Torus("torus", 0.0*um, TorusWidth, RedBloodCellRadius, 0*degree, 360*degree);
    G4Tubs* solidDisk = new G4Tubs("disk", 0.0*um, RedBloodCellRadius-TorusWidth, 0.5*um, 0*degree, 360*degree);
    
    
    G4RotationMatrix* myRotation = new G4RotationMatrix();
    myRotation->rotateX(0.*deg);
    myRotation->rotateY(0.*deg);
    myRotation->rotateZ(0.*rad);
    G4UnionSolid* gRBC = new G4UnionSolid(fName, solidTorus, solidDisk, myRotation, G4ThreeVector(0,0,0));
    
    fEnvelopeLog = CreateLogicalVolume(gRBC);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
