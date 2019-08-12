// Component for TsCylindericalChromatin
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

#include "TsCylindericalChromatin.hh"

#include "TsParameterManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Tubs.hh"

TsCylindericalChromatin::TsCylindericalChromatin(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{;}


TsCylindericalChromatin::~TsCylindericalChromatin()
{;}


G4VPhysicalVolume* TsCylindericalChromatin::Construct()
{
	BeginConstruction();

    G4String name = GetFullParmName("ChromatinRadius");
    G4double ChromatinRadius;
    if (!fPm->ParameterExists(name)) {
        ChromatinRadius = 12.5*nm;
    }
    else{
       ChromatinRadius = fPm->GetDoubleParameter(name, "Length");
    }
    
    name = GetFullParmName("ChromatinHalfLength");
    G4double ChromatinLength;
    if (!fPm->ParameterExists(name)) {
        ChromatinLength = 12.5*nm;
    }
    else{
        ChromatinLength = fPm->GetDoubleParameter(name, "Length");
    }
    
	G4Tubs* envelopeSolid = new G4Tubs(fName, 0, ChromatinRadius, ChromatinLength,0.0*deg, 360.0*deg);
	fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
