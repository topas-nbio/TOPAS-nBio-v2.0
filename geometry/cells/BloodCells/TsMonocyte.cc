// Component for TsMonocyte
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
// White blood cell - Monocyte.
// Monocytes have a large kidney-shaped nucleus.

#include "TsMonocyte.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Torus.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4UnionSolid.hh"


TsMonocyte::TsMonocyte(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsMonocyte::~TsMonocyte()
{;}

void TsMonocyte::ResolveParameters() {
    MonocyteRadius = fPm->GetDoubleParameter(GetFullParmName("MonocyteRadius"), "Length");
}


G4VPhysicalVolume* TsMonocyte::Construct()
{
	BeginConstruction();    
    
    //***********************************************************************
    //              Envelope Geometry : spherical cell
    //***********************************************************************
    
    G4Orb* gMonocyte = new G4Orb(fName, MonocyteRadius);
    fEnvelopeLog = CreateLogicalVolume(gMonocyte);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    // Nucleus is kidney-shaped.
    //***********************************************************************
    
    G4String subComponentName1 = "Nucleus";
    G4double torMin = 0.0*um;
    G4double torMax = 4.0*um;
    G4double nuclRadius = 5.0*um;
    G4double pPhi = 0.0*deg;
    G4double pAng = 250*deg;
    
    G4Torus* gNucleus = new G4Torus(subComponentName1, torMin, torMax, nuclRadius, pPhi, pAng);
    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
    CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);
        

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
