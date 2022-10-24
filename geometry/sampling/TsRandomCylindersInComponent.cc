// Component for TsRandomCylindersInComponent
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

#include "TsRandomCylindersInComponent.hh"
#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <fstream>
#include <math.h>

TsRandomCylindersInComponent::TsRandomCylindersInComponent(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
														   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
	fCreate = true;
	if (fPm->ParameterExists(GetFullParmName("GenerateCylinders")))
		fCreate = fPm->GetBooleanParameter(GetFullParmName("GenerateCylinders"));
	
	fOutputFile = "RandomCylinders_positions_nm_rotations_deg.txt";
	if ( fPm->ParameterExists(GetFullParmName("OutputFile") ))
		fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));

	if ( !fCreate ) 
		fOutputFile = fPm->GetStringParameter(GetFullParmName("InputFile"));
	
	fNumberOfCylinders = 10000;
	if ( fPm->ParameterExists(GetFullParmName("NumberOfCylinders") ))
		fNumberOfCylinders = fPm->GetIntegerParameter(GetFullParmName("NumberOfCylinders"));
	
	fTRMax = 0.5*2.3*nm;
	fTHL   = 0.5*3.4*nm;
	if (fPm->ParameterExists(GetFullParmName("Cylinders/RMax")))
		fTRMax = fPm->GetDoubleParameter(GetFullParmName("Cylinders/RMax"), "Length");
	if ( fPm->ParameterExists(GetFullParmName("Cylinders/HL")))
		fTHL = fPm->GetDoubleParameter(GetFullParmName("Cylinders/HL"), "Length");
	
	fEnvelopeType = fPm->GetStringParameter(GetFullParmName("EnvelopeShape"));
	fEnvelopeType.toLower();
	if ( fEnvelopeType.contains("tssphere") ) {
		fRMax = fPm->GetDoubleParameter(GetFullParmName("RMax"), "Length");
	} else if (fEnvelopeType.contains("tsbox")) {
		fHLX = fPm->GetDoubleParameter(GetFullParmName("HLX"), "Length");
		fHLY = fPm->GetDoubleParameter(GetFullParmName("HLY"), "Length");
		fHLZ = fPm->GetDoubleParameter(GetFullParmName("HLZ"), "Length");
	} else if ( fEnvelopeType.contains("tscylinder")) {
		fRMax = fPm->GetDoubleParameter(GetFullParmName("RMax"), "Length");
		fHL = fPm->GetDoubleParameter(GetFullParmName("HL"), "Length");
	} else {
		G4cerr << "TOPAS is exiting due to an error in parameter " << GetFullParmName("EnvelopeShape") << G4endl;
		G4cerr << ". Envelope shape " << fEnvelopeType << " is not found or not supported." <<
		"Only TsSphere, TsCylinder (full filled, not RMin) and TsBox are supported." << G4endl;
		   fPm->AbortSession(1);
	}
	
	if ( fCreate && fPm->GetIntegerParameter("Ts/NumberOfThreads") > 1 ) {
		G4cerr << "TOPAS is exiting due to an error. Conflict between component TsRandomCylindersInComponent " << G4endl;
		G4cerr << "and Ts/NumberOfThreads. " << GetFullParmName("GenerateCylinders") << " is set to true. Then "
		<< "Ts/NumberOfThredas must be set to 1." << G4endl;
	    fPm->AbortSession(1);
	}
}


TsRandomCylindersInComponent::~TsRandomCylindersInComponent()
{;}


G4VPhysicalVolume* TsRandomCylindersInComponent::Construct()
{
	BeginConstruction();
	
	G4double extensionX, extensionY, extensionZ;
	if ( fEnvelopeType.contains("tssphere" ) ) {
		fEnvelope = new G4Orb(fName, fRMax);
		extensionX = 2.0 * fRMax;
		extensionY = extensionX;
		extensionZ = extensionY;
	} else if ( fEnvelopeType.contains("tscylinder")) {
		fEnvelope = new G4Tubs(fName, 0.0, fRMax, fHL, 0.0, 360.0*deg);
		extensionX = 2.0 * fRMax;
		extensionY = extensionX;
		extensionZ = 2.0 * fHL;
	} else {
		fEnvelope = new G4Box(fName, fHLX, fHLY, fHLZ);
		extensionX = 2.0 * fHLX;
		extensionY = 2.0 * fHLY;
		extensionZ = 2.0 * fHLZ;
	}

	
	fEnvelopeLog = CreateLogicalVolume(fEnvelope);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	
	G4String subComponentName = "Cylinders";
	G4Tubs* cylinder = new G4Tubs(subComponentName, 0, fTRMax, fTHL, 0*deg, 360*deg);
	G4LogicalVolume* cylinderLog = CreateLogicalVolume(subComponentName, cylinder);
	
	G4double x, y, z, u, v, w;

	if ( fCreate ) {
		std::ofstream outFile(fOutputFile);
		G4int n = 0;
		G4int volID = 1;
        	
		while( n < fNumberOfCylinders ) {
			x = G4RandFlat::shoot(-0.5*extensionX,0.5*extensionX);
			y = G4RandFlat::shoot(-0.5*extensionY,0.5*extensionY);
			z = G4RandFlat::shoot(-0.5*extensionZ,0.5*extensionZ);
			
			EInside test_status = fEnvelope->Inside(G4ThreeVector(x, y, z));

			if ( test_status == kInside ) {
				G4RotationMatrix* rot = new G4RotationMatrix();
				u = G4RandFlat::shoot(0., 180.0*deg);
				v = G4RandFlat::shoot(0., 180.0*deg);
				w = G4RandFlat::shoot(0., 180.0*deg);
				rot->rotateX(u);
				rot->rotateY(v);
				rot->rotateZ(w);
				G4ThreeVector* pos = new G4ThreeVector(x, y, z);
				G4VPhysicalVolume* phys = CreatePhysicalVolume(fName, volID, true, cylinderLog,
															   rot, pos, fEnvelopePhys);
				if ( phys->CheckOverlaps(1000, 0., false, 1) ) {
					delete phys;
				} else {
					outFile << x/nm << " " << y/nm << " " << z/nm << " "
					        << u/deg << " " << v/deg << " " << w/deg << G4endl;
					
					volID++;
					n++;
				}
			}
		}
	} else {
		G4int volID = 1;
		std::ifstream inFile(fOutputFile);
		while(1) {
			inFile >> x >> y >> z >> u >> v >> w;
			if ( !inFile.good() ) break;
			G4RotationMatrix* rot = new G4RotationMatrix();
			rot->rotateX(u*deg);
			rot->rotateY(v*deg);
			rot->rotateZ(w*deg);
			CreatePhysicalVolume(fName, volID, true, cylinderLog, rot, new G4ThreeVector(x*nm,y*nm,z*nm), fEnvelopePhys);
			volID++;
		}
		
	}
	
	InstantiateChildren(fEnvelopePhys);
	return fEnvelopePhys;
}

