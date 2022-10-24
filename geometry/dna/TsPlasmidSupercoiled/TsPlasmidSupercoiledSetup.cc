// Component for TsPlasmidSupercoiledSetup
//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#include "TsPlasmidSupercoiledSetup.hh"
#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Point3D.hh"
#include "G4Scene.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"

#include "G4NistManager.hh"
#include <sstream>

TsPlasmidSupercoiledSetup::TsPlasmidSupercoiledSetup(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
										   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsPlasmidSupercoiledSetup::~TsPlasmidSupercoiledSetup()
{
}


G4VPhysicalVolume* TsPlasmidSupercoiledSetup::Construct()
{
	BeginConstruction();
	
	G4double radius = fPm->GetDoubleParameter(GetFullParmName("R"),"Length");
	G4String fileName = fPm->GetStringParameter(GetFullParmName("InputFile"));
	G4String ofileName = fPm->GetStringParameter(GetFullParmName("OutputFile"));
	G4bool random = true;
	if ( fPm->ParameterExists(GetFullParmName("Randomize")))
		random = fPm->GetBooleanParameter(GetFullParmName("Randomize"));
	
	fSolid = new G4Orb(fName, radius);
	fEnvelopeLog = CreateLogicalVolume(fSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	
	fOffsetX = 0 * nm;
	fOffsetY = 0 * nm;
	fOffsetZ = 0 * nm;
	ParseFile(fileName);
	
	G4String plasmidName = "Plasmid";
	G4Box* box_solid = new G4Box(plasmidName, 0.5*(fXMax-fXMin)+0.5*3.4*nm,
								 0.5*(fYMax-fYMin)+0.5*3.4*nm,
								 0.5*(fZMax-fZMin)+0.5*3.4*nm);
	
	G4LogicalVolume* box_logic = CreateLogicalVolume(plasmidName, box_solid);
	std::vector<G4RotationMatrix*> rotations;
	std::vector<G4ThreeVector*> translations;
	G4String theName = ofileName + "_envelopes.xyz";
	std::ofstream envFile(theName);
	G4double x, y, z, u, v, w;
	G4int n = 0, volID = 0;
	
	G4int nPlasmids = 0;
	
	if ( random ) {
		nPlasmids = fPm->GetIntegerParameter(GetFullParmName("NumberOfPlasmids"));
		while( n < nPlasmids ) {
			x = G4RandFlat::shoot(-radius,radius);
			y = G4RandFlat::shoot(-radius,radius);
			z = G4RandFlat::shoot(-radius,radius);
			
			EInside test_status = fSolid->Inside(G4ThreeVector(x, y, z));
			
			if ( test_status == kInside ) {
				G4RotationMatrix* rot = new G4RotationMatrix();
				u = G4RandFlat::shoot(0., 180.0*deg);
				v = G4RandFlat::shoot(0., 180.0*deg);
				w = G4RandFlat::shoot(0., 180.0*deg);
				rot->rotateX(u);
				rot->rotateY(v);
				rot->rotateZ(w);
				G4ThreeVector* pos = new G4ThreeVector(x, y, z);
				G4VPhysicalVolume* phys = CreatePhysicalVolume(plasmidName, volID, true, box_logic,
											rot, pos, fEnvelopePhys);
				if ( phys->CheckOverlaps(1000, 0., false, 1) ) {
					delete phys;
				} else {
					rotations.push_back(rot);
					translations.push_back(pos);
					envFile << x/nm << " " << y/nm << " " << z/nm << " " << u/deg << " " << v/deg << " " << w/deg << G4endl;
					delete phys;
					volID++;
					n++;
				}
			}
		}
	} else {
		G4int nx = fPm->GetIntegerParameter(GetFullParmName("nx"));
		G4int ny = fPm->GetIntegerParameter(GetFullParmName("ny"));
		G4int nz = fPm->GetIntegerParameter(GetFullParmName("nz"));
		nPlasmids = nx*ny*nz;
		G4double dx = (fXMax-fXMin)+3.4*nm;
		G4double dy = (fYMax-fYMin)+3.4*nm;
		G4double dz = (fZMax-fZMin)+3.4*nm;
		G4double xmin = -0.5 * nx * dx + 0.5 * dx;
		G4double ymin = -0.5 * ny * dy + 0.5 * dy;
		G4double zmin = -0.5 * nz * dz + 0.5 * dz;
		
		for ( int i = 0; i < nx; i++ ) {
			x = xmin + i * dx;
			for ( int j = 0; j < ny; j++ ) {
				y = ymin + j * dy;
				for ( int k = 0; k < nz; k++ ) {
					z = zmin + k * dz;
					G4RotationMatrix* rot = new G4RotationMatrix();
					u = 0;
					v = 0;
					w = 0;
					G4ThreeVector* pos = new G4ThreeVector(x, y, z);
					rotations.push_back(rot);
					translations.push_back(pos);
					envFile << x/nm << " " << y/nm << " " << z/nm << " " << u/deg << " " << v/deg << " " << w/deg << G4endl;
				}
			}
		}
	}

	envFile.close();
	
	theName = ofileName + "_plasmids.xyz";	
 	std::ofstream plFile(theName);
	G4int currentCopyNumber = 0, lastCopyNumber = 0;
	size_t t = 0;
	for ( int p = 0; p < nPlasmids; p++ ) {
		G4Point3D* aPoint = new G4Point3D(translations[p]->x(), translations[p]->y(), translations[p]->z());
		G4RotationMatrix* aInvRot = new G4RotationMatrix( rotations[p]->inverse() );
		for ( t = 0; t < fDNASpecies.size(); t++ ) {
			G4Point3D newPoint = G4Translate3D(*aInvRot * G4Point3D(fDNASpecies[t].second.second.x(),
																	fDNASpecies[t].second.second.y(),
																	fDNASpecies[t].second.second.z())) * (*aPoint);
			
			currentCopyNumber = fDNASpecies[t].first < 0 ? fDNASpecies[t].first - lastCopyNumber : fDNASpecies[t].first + lastCopyNumber;
			plFile << currentCopyNumber << " " << fDNASpecies[t].second.first << " "
			       << newPoint.x()/nm << " " << newPoint.y()/nm << " " << newPoint.z()/nm << G4endl;
		}
		lastCopyNumber = currentCopyNumber; //fDNASpecies[t-1].first;
		if ( lastCopyNumber < 0 )
			lastCopyNumber *= -1;
		
		delete aInvRot;
	}
	plFile.close();
	
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}


void TsPlasmidSupercoiledSetup::ParseFile(G4String fileName)
{
	// Clear the containers
	
	G4bool addMoleFlag = true;
	
	// Setup the input stream
	std::ifstream file(fileName.c_str());
	
	// Check if the file was correctly opened
	if(!file.is_open())
	{
		G4String msg = fileName+" could not be opened";
		G4Exception("PhysGeoImport::ParseFile()", "Geo_InputFileNotOpened", FatalException, msg);
	}
	
	fXMin = 1*mm, fYMin = 1*mm, fZMin = 1*mm;
	fXMax=0.0, fYMax=0.0, fZMax=0.0;
	
	// Define the line string variable
	std::string line;
	// Read the file line per line
	while(std::getline(file, line) )
	{
		// Check the line to determine if it is empty
		if(line.empty()) continue; // skip the line if it is empty
		
		// Data string stream
		std::istringstream issLine(line);
		
		// String to determine the first letter/word
		std::string firstItem;
		
		// Put the first letter/word within the string
		issLine >> firstItem;
		
		// Check first letter to determine if the line is data or comment
		if(firstItem=="#") continue; // skip the line if it is comment
		
		// Use the file
		else if(firstItem=="_Name")
		{
			std::string name;
			issLine >> name;
			
			fGeoName = name;
		}
		else if(firstItem=="_Size")
		{
			G4double size;
			issLine >> size;
			size *= nm;
			
			fSize = size;
		}
		else if(firstItem == "_Version")
		{
			
		}
		else if(firstItem=="_Number")
		{
			
		}
		else if(firstItem=="_Radius")
		{
			std::string name;
			issLine >> name;
			
			G4double radius;
			issLine >> radius;
			radius *= nm;
			
			G4double waterRadius;
			issLine >> waterRadius;
			waterRadius *= nm;
			
		}
		else if(firstItem=="_pl")
		{
			std::string name;
			issLine >> name;
			
			std::string material;
			issLine >> material;
			
			G4int strand;
			issLine >> strand;
			
			G4int copyNumber;
			issLine >> copyNumber;
			
			G4double x;
			issLine >> x;
			x *= nm;
			
			G4double y;
			issLine >> y;
			y *= nm;
			
			G4double z;
			issLine >> z;
			z *= nm;
			
			if ( fXMin > x )
				fXMin = x;
			if ( fXMax < x )
				fXMax = x;
			if ( fYMin > y )
				fYMin = y;
			if ( fYMax < y )
				fYMax = y;
			if ( fZMin > z )
				fZMin = z;
			if ( fZMax < z )
				fZMax = z;
			
			G4int moleculeID = -1;
			
			if(addMoleFlag)
			{
				if(name != "phosphate1" && name != "phosphate2") //phosphate and histone do not react with molecules
				{
					if(name == "deoxyribose1" || name == "deoxyribose2") {
						name = "Deoxyribose";
						moleculeID = 104;
					}
					if(name == "base_adenine") {
						name = "Adenine";
						moleculeID = 102;
					}
					if(name == "base_thymine") {
						name = "Thymine";
						moleculeID = 103;
					}
					if(name == "base_cytosine") {
						name = "Cytosine";
						moleculeID = 101;
					}
					if(name == "base_guanine") {
						name = "Guanine";
						moleculeID = 100;
					}
					if(name == "histone") {
						name = "Histone";
						moleculeID = 99;
					}
					if ( moleculeID > 0 ) {
						if ( strand == 1 )
							fDNASpecies.push_back(std::make_pair(copyNumber, std::make_pair(moleculeID, G4ThreeVector(x, y, z))));
						else
							fDNASpecies.push_back(std::make_pair(-copyNumber, std::make_pair(moleculeID, G4ThreeVector(x, y, z))));
					}
				}
			}
		}
		else
		{
			// Geant4 exception
			G4String msg = firstItem+" is not defined in the parser. Check the input file: "+fileName+".";
			G4Exception("PhysGeoImport::ParseFile()", "Geo_WrongParse", FatalException, msg);
		}
	}
	fOffsetX = (fXMin + fXMax)*0.5;
	fOffsetY = (fYMin + fYMax)*0.5;
	fOffsetZ = (fZMin + fZMax)*0.5;
	file.close();
}
