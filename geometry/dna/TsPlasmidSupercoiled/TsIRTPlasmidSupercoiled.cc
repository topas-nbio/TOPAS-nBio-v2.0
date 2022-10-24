// Component for TsIRTPlasmidSupercoiled
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

#include "TsIRTPlasmidSupercoiled.hh"
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

TsIRTPlasmidSupercoiled::TsIRTPlasmidSupercoiled(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
										   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsIRTPlasmidSupercoiled::~TsIRTPlasmidSupercoiled()
{
}


G4VPhysicalVolume* TsIRTPlasmidSupercoiled::Construct()
{
	BeginConstruction();
	G4int nPlasdmids = fPm->GetIntegerParameter(GetFullParmName("NumberOfPlasmids"));
	G4double radius = fPm->GetDoubleParameter(GetFullParmName("R"),"Length");
	G4String fileName = fPm->GetStringParameter(GetFullParmName("InputFile"));
	
	fSolid = new G4Orb(fName, radius);
	fEnvelopeLog = CreateLogicalVolume(fSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	
	fOffsetX = 0 * nm;
	fOffsetY = 0 * nm;
	fOffsetZ = 0 * nm;
	ParseFile(fileName);
	
	std::vector<G4RotationMatrix*> rotations;
	std::vector<G4ThreeVector*> translations;
	G4String theName = fPm->GetStringParameter(GetFullParmName("EnvelopeFileName"));
 
	std::ifstream envFile(theName); //"envelopes.xyz");
	G4double x, y, z, u, v, w;
	G4int n = 0, volID = 0;
	while( n < nPlasdmids ) {
		envFile >> x >> y >> z >> u >> v >> w;
		if ( !envFile.good() )
			break;
		x *= nm;
		y *= nm;
		z *= nm;
		u *= deg;
		v *= deg;
		w *= deg;
		G4RotationMatrix* rot = new G4RotationMatrix();
		rot->rotateX(u);
		rot->rotateY(v);
		rot->rotateZ(w);
		G4ThreeVector* pos = new G4ThreeVector(x, y, z);
		
		rotations.push_back(rot);
		translations.push_back(pos);
		
		volID++;
		n++;
	}
	envFile.close();
	
	G4LogicalVolume* subComponent = 0;
	G4String plasmidName = "plasmid";
	for ( int p = 0; p < nPlasdmids; p++ ) {
		if ( p == 0 )
			subComponent = CreateLogicVolume(fileName, p, rotations[p], translations[p]);
		else
			CreatePhysicalVolume(plasmidName, p, true, subComponent, rotations[p], translations[p], fEnvelopePhys);
	}
	
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}


G4LogicalVolume* TsIRTPlasmidSupercoiled::CreateLogicVolume(G4String, G4int copy, G4RotationMatrix* rot, G4ThreeVector* trans)
{
	G4String materialBox = fParentComponent->GetResolvedMaterialName();
	G4String materialDNA = GetResolvedMaterialName();
	
	std::string boxNameSolid = fGeoName+"_solid";
	G4Box* box_solid = new G4Box(boxNameSolid, 0.5*(fXMax-fXMin)+0.5*3.4*nm,
								 0.5*(fYMax-fYMin)+0.5*3.4*nm,
								 0.5*(fZMax-fZMin)+0.5*3.4*nm);
	G4String boxNameLogic = fGeoName+"_logic";
	G4LogicalVolume* box_logic = CreateLogicalVolume(boxNameLogic, box_solid);
	
	G4VPhysicalVolume* box_phys = CreatePhysicalVolume(boxNameLogic, copy, true, box_logic, rot, trans, fEnvelopePhys);
	
	G4ThreeVector offsetPosition(fOffsetX, fOffsetY, fOffsetZ);
	
	G4VisAttributes* purines = fPm->GetColor("red");
	purines->SetForceSolid();
	G4VisAttributes* pyrimidines = fPm->GetColor("yellow");
	pyrimidines->SetForceSolid();
	G4VisAttributes* deoxyriboses = fPm->GetColor("magenta");
	deoxyriboses->SetForceSolid();
	G4VisAttributes* phosphates = fPm->GetColor("grass");
	phosphates->SetForceSolid();
	G4VisAttributes* histones = fPm->GetColor("blue");
	histones->SetForceSolid();
	G4VisAttributes* wshell = new G4VisAttributes(G4Color(0,0,1,0.1));
	wshell->SetForceSolid();
	RegisterVisAtt(wshell);
	
	G4Orb* purineAndPyrimidine = 0;
	G4Orb* phosphate = 0;
	G4Orb* deoxyribose = 0;
	G4Orb* wpurineAndPyrimidine = 0;
	G4Orb* wphosphate = 0;
	G4Orb* wdeoxyribose = 0;
	G4Orb* histone = 0;
	G4bool solidPurinePyrimidine = true;
	G4bool solidPhosphate = true;
	G4bool solidDeoxyribose = true;
	G4bool solidHistone = true;
	
	G4cout << "#### Creating plasmid base pairs: begin" << G4endl;
	std::sort(fMolecules.begin(), fMolecules.end() );
	for(int i=0, ie=fMolecules.size(); i<ie; ++i)
	{
		G4String name = fMolecules[i].fName;
		G4double radius = fMolecules[i].fRadius;
		G4double waterRadius = fMolecules[i].fRadiusWater;
		G4ThreeVector moleculePosition = fMolecules[i].fPosition - offsetPosition;
		G4int copyNum = fMolecules[i].fCopyNumber;
		
		if ( solidPhosphate || solidPurinePyrimidine || solidDeoxyribose || solidHistone ) {
			if ( name.contains("phosphate") && solidPhosphate ) {
				phosphate = new G4Orb("phosphate", radius);
				wphosphate = new G4Orb("wphosphate", waterRadius);
				solidPhosphate = false;
				G4cout << "####  Built phosphate" << G4endl;
			} else if ( name.contains("base") && solidPurinePyrimidine ) {
				purineAndPyrimidine = new G4Orb("purineandpyrimidine", radius);
				wpurineAndPyrimidine = new G4Orb("wpurineandpyrimidine", waterRadius);
				solidPurinePyrimidine = false;
				G4cout << "####  Built purines and pyrimidines" << G4endl;
			} else if ( name.contains("deoxyribose") && solidDeoxyribose ) {
				deoxyribose = new G4Orb("deoxyribose", radius);
				wdeoxyribose = new G4Orb("wdeoxyribose", waterRadius);
				solidDeoxyribose = false;
				G4cout << "####  Built deoxyriboses" << G4endl;
			} else if ( name.contains("histone") && solidHistone ) {
				histone = new G4Orb("histone", radius);
				solidHistone = false;
			}
		}
		// Water hydration shell volume part
		G4VPhysicalVolume* moleculeWater_phys = 0;
		// If water radius != 0 then we have a water hydration shell
		G4double tol = 0.0001;
		G4VSolid* moleculeWaterCut_solid = 0;
		G4LogicalVolume* moleculeWater_logic = 0;
		if(waterRadius > (0 + tol)*nm)
		{
			G4String nameWaterSolid = name+"_waterShell";
			if ( name.contains("phosphate") ) {
				moleculeWaterCut_solid = CreateCutSolid(wphosphate, fMolecules[i], fMolecules, false);
			} else if ( name.contains("base") ) {
				moleculeWaterCut_solid = CreateCutSolid(wpurineAndPyrimidine, fMolecules[i], fMolecules, false);
			} else if ( name.contains("deoxyribose") ) {
				moleculeWaterCut_solid = CreateCutSolid(wdeoxyribose, fMolecules[i], fMolecules, false);
			}
			
			//moleculeWater_logic = CreateLogicalVolume(nameWaterSolid, materialDNA, moleculeWaterCut_solid);
			moleculeWater_logic = CreateLogicalVolume(nameWaterSolid, moleculeWaterCut_solid);
			moleculeWater_logic->SetVisAttributes(wshell);
			moleculeWater_phys = CreatePhysicalVolume(nameWaterSolid, copyNum, true, moleculeWater_logic, 0, new G4ThreeVector(moleculePosition), box_phys);
		}
		
		// Dna volume part
		G4VSolid* moleculeCut_solid = 0;
		G4String nameSolid = fMolecules[i].fName+"_solid";
		
		G4String nameLogic = name;
		G4LogicalVolume* molecule_logic = 0;
		
		if ( name.contains("phosphate") ) {
			moleculeCut_solid = CreateCutSolid(phosphate, fMolecules[i], fMolecules, true);
		} else if ( name.contains("base") ) {
			moleculeCut_solid = CreateCutSolid(purineAndPyrimidine, fMolecules[i], fMolecules, true);
		} else if ( name.contains("deoxyribose") ) {
			moleculeCut_solid = CreateCutSolid(deoxyribose, fMolecules[i], fMolecules, true);
		} else if ( name.contains("histone") ) {
			moleculeCut_solid = CreateCutSolid(histone, fMolecules[i], fMolecules, true);
		}
		
		molecule_logic = CreateLogicalVolume(nameLogic, moleculeCut_solid);
		
		G4ThreeVector position(0);
		G4String namePhys = name;
		if(waterRadius > (0 + tol)*nm)
			CreatePhysicalVolume(namePhys, copyNum, true, molecule_logic, 0, new G4ThreeVector(position), moleculeWater_phys);
		else
			CreatePhysicalVolume(namePhys, copyNum, true, molecule_logic, 0, new G4ThreeVector(moleculePosition), box_phys);
	}
	G4cout << "#### Creating plasmid base pairs: end" << G4endl;
	// Clear the containers
	fMolecules.clear();
	fRadiusMap.clear();
	fWaterRadiusMap.clear();
	return box_logic;
}

void TsIRTPlasmidSupercoiled::ParseFile(G4String fileName)
{
	// Clear the containers
	fMolecules.clear();
	fRadiusMap.clear();
	fWaterRadiusMap.clear();
	
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
			
			fRadiusMap[name] = radius;
			fWaterRadiusMap[name] = waterRadius;
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
			
			TempMolecule molecule(name, copyNumber, G4ThreeVector(x, y, z), fRadiusMap[name], fWaterRadiusMap[name], material, strand);
			fMolecules.push_back(molecule);
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
					if ( moleculeID > 99 ) {
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


G4VSolid* TsIRTPlasmidSupercoiled::CreateCutSolid(G4Orb *solidOrbRef,
											   TempMolecule &molRef,
											   std::vector<TempMolecule> &molList,
											   G4bool in)
{
	// The idea behing this method is to cut overlap volumes by selecting one of them (the reference) and checking all the other volumes (the targets).
	// If a reference and a target volumes are close enough to overlap they will be cut.
	// The reference is already selected when we enter this method.
	
	// Use the tiny space to differentiate the frontiers (may not be necessary)
	G4double tinySpace = 0.001*nm;
	
	// Cutted solid to be returned
	G4SubtractionSolid* solidCut(NULL);
	
	// Some flags
	G4bool isCutted = false;
	G4bool isOurVol = false;
	
	// Radius of the molecule to cut
	G4double radiusRef;
	if(molRef.fRadiusWater==0)
		radiusRef = molRef.fRadius;
	else
		radiusRef = molRef.fRadiusWater;
	
	// Reference volume position
	G4ThreeVector posRef = molRef.fPosition;
	
	// Look the other volumes of the voxel
	// Loop on all the target volumes (other volumes with potential overlaps)
	for(int i=0, ie=molList.size(); i<ie; ++i)
	{
		G4ThreeVector posTar = molList[i].fPosition;
		
		G4double rTar = posRef.z();
		G4double zTar = posTar.z();
		
		if(zTar>rTar+20*nm)
		{
			break;
		}
		else if(zTar<rTar-20*nm)
		{
			continue;
		}
		
		
		// Retrieve current target sphere informations
		G4double radiusTar;
		if(molList[i].fRadiusWater==0) radiusTar = molList[i].fRadius;
		else radiusTar = molList[i].fRadiusWater;
		
		// Compute the distance reference-target
		G4double distance = std::abs( (posRef - posTar).getR() );
		
		// Use the distance to check if the current target is also the reference.
		// This can only happen once per loop.
		if(distance==0 && !isOurVol)
		{
			// Target volume is also reference volume.
			
			// Set the flag
			isOurVol = true;
			
			// Next iteration
			continue;
		}
		// If the condition is correct more than one time then there is a mistake somewhere.
		else if(distance == 0 && isOurVol)
		{
			G4cerr<<"********************* Fatal Error **************************"<<G4endl;
			G4cerr<<"DetectorConstruction::CreateCutSolid: Two volumes are placed at the same position."<<G4endl;
			exit(EXIT_FAILURE);
		}
		
		// If the volumes are differents then we want to know if they are
		// close enough to overlap and, thus, to intiate a cut.
		else if(distance <= radiusRef+radiusTar)
		{
			// Volumes are close enough, there will be a cut
			
			// Box used to cut
			G4Box* solidBox = new G4Box("solid_box_for_cut", 2*radiusTar, 2*radiusTar, 2*radiusTar);
			
			// This part is tricky.
			// The goal is to calculate the position of the intersection center
			
			// diff vector to from ref to tar
			G4ThreeVector diff = posTar - posRef;
			
			// Find the intersection point and add to it half the length of the box used to cut
			G4double d = (pow(radiusRef,2)-pow(radiusTar,2)+pow(distance,2) ) / (2*distance) + solidBox->GetZHalfLength() - tinySpace;
			
			// If we are in another volume we choose to double the tinySpace to differentiate without ambiguities the inner and outer volume frontiers.
			// (may not be necessary)
			if(in) d -= 2*tinySpace;
			
			// Position of the box in order to achieve the cut.
			// "* ( diff/diff.getR() )" is necessary to get a vector in the right direction as output
			G4ThreeVector pos = d *( diff/diff.getR() );
			
			// Get the rotation angles because the box used to cut needs to be rotated
			// to give the right "cut angle".
			G4double phi = std::acos(pos.getZ()/pos.getR());
			G4double theta = std::acos( pos.getX() / ( pos.getR()*std::cos(M_PI/2.-phi) ) );
			
			if(pos.getY()<0) theta = -theta;
			
			G4ThreeVector rotAxisForPhi(1*nm,0.,0.);
			rotAxisForPhi.rotateZ(theta+M_PI/2);
			
			// Create the rotation matrix
			G4RotationMatrix *rotMat = new G4RotationMatrix;
			rotMat->rotate(-phi, rotAxisForPhi);
			
			// Rotate it again
			G4ThreeVector rotZAxis(0.,0.,1*nm);
			rotMat->rotate(theta, rotZAxis);
			
			// If the volume is cutted for the first time
			if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, pos);
			
			// For the other times
			else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, pos);
			
			// Set the cut flag
			isCutted = true;
		}
	}
	
	//delete rotMat;
	if(isCutted) return solidCut;
	
	// Otherwise, we return the original volume
	else return solidOrbRef;
}




