// Component for TsLymphocyte
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
// Lymphocyte

#include "TsLymphocyte.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsLymphocyte::TsLymphocyte(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}

TsLymphocyte::~TsLymphocyte()
{;}

void TsLymphocyte::ResolveParameters() {
    LymphocyteRadius = fPm->GetDoubleParameter(GetFullParmName("LymphocyteRadius"), "Length");
}


G4VPhysicalVolume* TsLymphocyte::Construct()
{
	BeginConstruction();    
    
    //***********************************************************************
    //              Envelope Geometry : spherical cell
    //***********************************************************************
    
    G4Orb* gLymphocyte = new G4Orb(fName, LymphocyteRadius);
    fEnvelopeLog = CreateLogicalVolume(gLymphocyte);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    // Nucleus is usuall large (fills the cell) and may be indented.
    //***********************************************************************
   
    G4String subComponentName1 = "Nucleus";
    G4double NuclRadius = 4*um;
    
    G4String name = GetFullParmName("Nucleus/NucleusRadius");
    if (fPm->ParameterExists(name)) {
        NuclRadius = fPm->GetDoubleParameter(name, "Length");
    }
        
    G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
    CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
