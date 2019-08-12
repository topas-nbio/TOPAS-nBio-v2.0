// Component for TsNeutrophil
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
// Neutrophil.
// Contains irregular nucleus with many lobes (2-5 lobes).

#include "TsNeutrophil.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4UnionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsNeutrophil::TsNeutrophil(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsNeutrophil::~TsNeutrophil()
{;}

void TsNeutrophil::ResolveParameters() {
    NeutrophilRadius = fPm->GetDoubleParameter(GetFullParmName("NeutrophilRadius"), "Length");
}


G4VPhysicalVolume* TsNeutrophil::Construct()
{
	BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : spherical cell
    //***********************************************************************
    
    G4Orb* gNeutrophil = new G4Orb(fName, NeutrophilRadius);
    fEnvelopeLog = CreateLogicalVolume(gNeutrophil);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    // Irregular nucleus with 3 to 5 lobes
    //***********************************************************************
    
    G4String subComponentName1 = "Nucleus";
    G4double Sphere1Radius = 1.0*um;
    G4double Sphere2Radius = 2.0*um;
    G4double Sphere3Radius = 0.5*um;
    G4double Sphere4Radius = 1.5*um;
    G4double Sphere5Radius = 1.4*um;
    
    G4Orb* solidSphere1 = new G4Orb("Sphere1", Sphere1Radius);
    G4Orb* solidSphere2 = new G4Orb("Sphere2", Sphere2Radius);
    G4Orb* solidSphere3 = new G4Orb("Sphere3", Sphere3Radius);
    G4Orb* solidSphere4 = new G4Orb("Sphere4", Sphere4Radius);
    G4Orb* solidSphere5 = new G4Orb("Sphere5", Sphere5Radius);
    
    G4RotationMatrix* nuclRotation = new G4RotationMatrix();
    nuclRotation->rotateX(0.*deg);
    nuclRotation->rotateY(0.*deg);
    nuclRotation->rotateZ(0.*rad);
    
    G4UnionSolid* gNucleus1 = new G4UnionSolid(subComponentName1, solidSphere1, solidSphere2, nuclRotation, G4ThreeVector(0*um,3.0*um,0));
    G4UnionSolid* gNucleus2 = new G4UnionSolid(subComponentName1, gNucleus1, solidSphere3, nuclRotation, G4ThreeVector(2.5*um,2*um,0));
    G4UnionSolid* gNucleus3 = new G4UnionSolid(subComponentName1, gNucleus2, solidSphere4, nuclRotation, G4ThreeVector(3.3*um,0.4*um,0.5*um));
    G4UnionSolid* gNucleus = new G4UnionSolid(subComponentName1, gNucleus3, solidSphere5, nuclRotation, G4ThreeVector(3.5*um,-2.0*um,0.2*um));
    
    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
    CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);
    

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
