// Component for TsCylindericalNucleosome
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

#include "TsCylindericalNucleosome.hh"

#include "TsParameterManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Tubs.hh"

TsCylindericalNucleosome::TsCylindericalNucleosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{;}


TsCylindericalNucleosome::~TsCylindericalNucleosome()
{;}


G4VPhysicalVolume* TsCylindericalNucleosome::Construct()
{
	BeginConstruction();

    G4String name = GetFullParmName("NucleosomeRadius");
    G4double NucleosomeRadius;
    if (!fPm->ParameterExists(name)) {
        NucleosomeRadius = 5*nm;
    }
    else{
       NucleosomeRadius = fPm->GetDoubleParameter(name, "Length");
    }
    
    name = GetFullParmName("NucleosomeHalfLength");
    G4double NucleosomeLength;
    if (!fPm->ParameterExists(name)) {
        NucleosomeLength = 2.5*nm;
    }
    else{
        NucleosomeLength = fPm->GetDoubleParameter(name, "Length");
    }
    
	G4Tubs* envelopeSolid = new G4Tubs(fName, 0, NucleosomeRadius, NucleosomeLength, 0.0*deg, 360.0*deg);
	fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
