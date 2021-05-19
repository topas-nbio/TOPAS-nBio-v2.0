// Component for TsNucleus
//
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			*
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	 *
// *																  *
// ********************************************************************
//
// Authors: Hongyu Zhu, Nicholas Henthorn, Alexander Klapproth, Jan Schuemann

// This class creates either single histones, chromatin fibers, voxels filled with DNA 
// or even a whole nucleus. The DNA is arranged in:
// 1) Basepairs wrapped around histones
// 2) Histones with linkers arranged to form chromatin fibers
// 3) Fibers arranged in a Hilbert filling curve to fill a voxel
// 4) Voxels arranged in using a Fractal walk pattern to fill a whole nucleus with DNA

// The fiber is based on an extension developed by Nicholas Henthorn
// A hydration shell was around the backbone volume and the fibers were used to fill the 
// nucleus based on work by Hongyu Zhu
// The fractal walk pattern was designed by Aimee McNamara
// Additional suggestions and edits were made by Alexander Klapproth

// This file includes the geometry for a linear segment of the chromatin Fiber
// The fiber is a solenoid conformation of cylinrical histones,
// wrapped by 1.65 left handed turns of the double helix
// The Histones are placed, recording position and rotation
// the histone position and rotation is passed to the DNA path and a general path is made
// the general path is segmented in 0.34 nm steps
// the DNA position and rotation is created and DNA placed


#include "TsNucleus.hh"
#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TwoVector.hh"

#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4CutTubs.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4UIcommand.hh"

#include <fstream>
#include <sstream>

using namespace std;

TsNucleus::TsNucleus(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{}


TsNucleus::~TsNucleus()
{}


G4VPhysicalVolume* TsNucleus::Construct()
{
	BeginConstruction();

	// User defined parameters.
	// Specify the radius of the nucleus
	fNucleusRadius = fPm->GetDoubleParameter(GetFullParmName("NucleusRadius"),"Length");

	// Specify how many chromatin fiber loops (Hilbert Curve) to put in a voxel
	fHilbertCurveLayer = fPm->GetIntegerParameter(GetFullParmName("HilbertCurveLayer"));
  
	// Specify how many voxels to put in X/Y/Z directions
	fHilbertCurve3DRepeat = fPm->GetIntegerParameter(GetFullParmName("HilbertCurve3DRepeat"));

	// Specify the file needed to build the hilbert filling curve
	G4String HilbertCurveFileName = fPm->GetStringParameter(GetFullParmName("HilbertCurveFileName"));

	fCheckOverlap = false;
	if (fPm->ParameterExists(GetFullParmName("fCheckOverlap")))
	  fCheckOverlap = fPm->GetBooleanParameter(GetFullParmName("fCheckOverlap"));

	fFillCylindersWithDNA = true;
	if (fPm->ParameterExists(GetFullParmName("FillCylindersWithDNA")))
	  fFillCylindersWithDNA = fPm->GetBooleanParameter(GetFullParmName("FillCylindersWithDNA"));

	fAddHistones = true;
	if (fPm->ParameterExists(GetFullParmName("AddHistones")))
		fAddHistones = fPm->GetBooleanParameter(GetFullParmName("AddHistones"));
	
	fOnlyBuildOneHistone = false;
	if (fPm->ParameterExists(GetFullParmName("OnlyBuildOneHistone")))
	  fOnlyBuildOneHistone = fPm->GetBooleanParameter(GetFullParmName("OnlyBuildOneHistone"));

	fShowNucleus = false;
	if (fPm->ParameterExists(GetFullParmName("ShowNucleusOutline")))
	  fShowNucleus = fPm->GetBooleanParameter(GetFullParmName("ShowNucleusOutline"));

	fShowDNAVoxels = false;
	if (fPm->ParameterExists(GetFullParmName("ShowDNAVoxels")))
		fShowDNAVoxels = fPm->GetBooleanParameter(GetFullParmName("ShowDNAVoxels"));

	fShowChromatinCylinders = true;
	if (fPm->ParameterExists(GetFullParmName("ShowChromatinCylinders")))
	  fShowChromatinCylinders = fPm->GetBooleanParameter(GetFullParmName("ShowChromatinCylinders"));

	//For chromatin fiber
	fFiberRadius = (37.0879/2.0) *nm;
	if (fPm->ParameterExists(GetFullParmName("FiberRadius")))
	  fFiberRadius = fPm->GetDoubleParameter(GetFullParmName("FiberRadius"),"Length");
	
	fFiberLength = 120 *nm;
	if (fPm->ParameterExists(GetFullParmName("FiberLength")))
	  fFiberLength = fPm->GetDoubleParameter(GetFullParmName("FiberLength"),"Length");
	
	fAddBases = true;
	if (fPm->ParameterExists(GetFullParmName("AddBases")))
	  fAddBases = fPm->GetBooleanParameter(GetFullParmName("AddBases"));
	
	fAddBackbones = true;
	if (fPm->ParameterExists(GetFullParmName("AddBackbone")))
		fAddBackbones = fPm->GetBooleanParameter(GetFullParmName("AddBackbone"));
	
	fAddHydrationShell = true;
	if (fPm->ParameterExists(GetFullParmName("AddHydrationShell")))
	  fAddHydrationShell = fPm->GetBooleanParameter(GetFullParmName("AddHydrationShell"));

	fHydrationShellThickness = 0.16*nm;
	if (fPm->ParameterExists(GetFullParmName("HydrationShellThickness")))
	  fHydrationShellThickness = fPm->GetDoubleParameter(GetFullParmName("fHydrationShellThickness"),"Length");
	
	fDNAModel = "QuarterCylinder";
	if (fPm->ParameterExists(GetFullParmName("DNAModel")))
	  fDNAModel = fPm->GetStringParameter(GetFullParmName("DNAModel"));
	if (!(fDNAModel=="HalfCylinder" || fDNAModel=="Sphere" || fDNAModel=="QuarterCylinder"))
		G4cerr << "Trying to build DNA with undefined value" << fDNAModel << G4endl;


	// Not using this feature anymore
	/*fScoreOnBases = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnBases")) && fPm->ParameterExists(GetFullParmName("AddBases")))
		fScoreOnBases = fPm->GetBooleanParameter(GetFullParmName("ScoreOnBases"));

	fScoreOnBackbones = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnBackbones")) && fPm->ParameterExists(GetFullParmName("AddBackbone")))
		fScoreOnBackbones = fPm->GetBooleanParameter(GetFullParmName("ScoreOnBackbones"));

	fScoreOnHydrationShell = false;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnHydrationShell")) && fPm->ParameterExists(GetFullParmName("AddHydrationShell")))
		fScoreOnHydrationShell = fPm->GetBooleanParameter(GetFullParmName("ScoreOnHydrationShell"));

	fScoreOnHistones = false;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnHistones")))
		fScoreOnHistones = fPm->GetBooleanParameter(GetFullParmName("ScoreOnHistones"));*/


	//****************************************************************************
	//					 Read  Hilbert space filling  data					//
	//****************************************************************************
	const char* filename = HilbertCurveFileName;
	std::ifstream Hfile(filename, std::ios::in);
	G4double x, y, z;
 
	if (Hfile.is_open()) {
		while (Hfile >> x >> y >> z){
			fiberPosX.push_back(x);
			fiberPosY.push_back(y);
			fiberPosZ.push_back(z);
		}
	}
	else {
		G4cout << "ERROR: Unable to open file " << HilbertCurveFileName << G4endl;
		exit(1);
	}

	G4int TotalPoints = fiberPosX.size();
	G4cout << "Number of points " << TotalPoints << G4endl;

	//****************************************************************************
	//							   Set geometry size						  //
	//****************************************************************************
	
	G4double FiberEnvelopeRadius = fFiberRadius;	//Fiber radius
	G4double FiberEnvelopeLength = fFiberLength;
	G4double HilbertPointDistance = 40*nm+FiberEnvelopeLength;

	G4int HilbertFold = 1;
	if(TotalPoints==8) HilbertFold = 1;
	else if(TotalPoints==64) HilbertFold = 3;
	else if(TotalPoints==512) HilbertFold = 7;
	else if(TotalPoints==4096) HilbertFold = 15;
	else {
		G4cerr << "Hilbert folding does not match number available folds." << G4endl;
		exit(1);
	}

	fVoxelLength = HilbertPointDistance*HilbertFold + fHilbertCurveLayer*FiberEnvelopeRadius*2 +0.001*um; //add tolerance 0.001 um
	G4double ParaContainerHalfSize = fVoxelLength*fHilbertCurve3DRepeat/2;

	G4cout<<"*********************************************************************************"<<G4endl;
	G4cout<<"The Hilbert curve data file name: " << HilbertCurveFileName <<G4endl;
	G4cout<<"The nuclues was devided into " << fHilbertCurve3DRepeat << " subdomains(voxels) on X/Y/Z diection respectively."<<G4endl;
	G4cout<<"Voxel size :"<< fVoxelLength/um << " um " <<G4endl;
	G4cout<<"Voxel container size :"<< ParaContainerHalfSize*2/um<<" um "<<G4endl;
	G4cout<<"Reuse the Hilbert curve "<<fHilbertCurveLayer<<" times to fill each subdomain(voxels)."<<G4endl;
	G4cout<<"HiberterPointDistance = "<<HilbertPointDistance/um<<" um"<<G4endl;
	G4cout<<"*********************************************************************************"<<G4endl;

	//****************************************************************************
	//								Parameterise							  //
	//****************************************************************************
	SetBasicInfo();

	////----- Create parameterisation and set
	param = new TsVoxelParameterisation();
	param->SetVoxelDimensions( fVoxelLength/2, fVoxelLength/2, fVoxelLength/2 ); 
	param->SetNoVoxel( fHilbertCurve3DRepeat, fHilbertCurve3DRepeat, fHilbertCurve3DRepeat ); 
	param->SetContainerDimensions(ParaContainerHalfSize, ParaContainerHalfSize, ParaContainerHalfSize);

	//****************************************************************************
	//							 Build basic geometry						 //
	//****************************************************************************
	// Nucleaus box 
	G4Orb * Nucleaus_solid = new G4Orb("Nucleus", fNucleusRadius);
	fEnvelopeLog		   = CreateLogicalVolume(Nucleaus_solid);
	fEnvelopePhys		  = CreatePhysicalVolume(fEnvelopeLog);
	G4Colour blue (0.0, 0.0, 1.0) ;
	G4VisAttributes* Vis   = new G4VisAttributes( blue );
	Vis->SetVisibility(fShowNucleus);
	fEnvelopeLog->SetVisAttributes(Vis);

	
	// Fiber Envelope
	G4double length = std::sqrt( pow(fiberPosX[2]-fiberPosX[1],2) + pow(fiberPosY[2]-fiberPosY[1],2) + pow(fiberPosZ[2]-fiberPosZ[1],2));
	G4double scaleFactor = HilbertPointDistance/length; // scaleFactor used to get fiber lengths
	length = scaleFactor*length;

	fNumberOfBasePairs = 0;
	// Logic volume of fiber
	fFiberLogic = ConstructFiberLogicalVolume();
//	G4VisAttributes* FiberEnvelopeVis = new G4VisAttributes( blue ); //
//	FiberEnvelopeVis->SetVisibility(true);
//	FiberEnvelopeVis->SetForceWireframe(true);
//	FiberEnvelopeVis->SetForceAuxEdgeVisible(true);
//	fFiberLogic->SetVisAttributes(FiberEnvelopeVis);

	//****************************************************************************
	//					   Build Hilbert loop geometry						//
	//****************************************************************************
	G4int maxSize = fHilbertCurveLayer*TotalPoints*(G4int)pow(fHilbertCurve3DRepeat,3);
	fFiberPhysVolLoop.resize(maxSize);

	//----- Define voxel logical volume
	G4String voxName = "Voxel";
	G4Box		   * voxel_solid = new G4Box( voxName, fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ);
	G4LogicalVolume * voxel_logic = CreateLogicalVolume(voxel_solid); 
	G4Colour  yellow  (1.0, 1.0, 0.0) ;
	G4VisAttributes* voxelVis = new G4VisAttributes( yellow ); //
	voxelVis->SetVisibility(fShowDNAVoxels);
	voxelVis->SetForceWireframe(true);
	voxelVis->SetForceAuxEdgeVisible(true);
	voxel_logic->SetVisAttributes(voxelVis);

	G4int CountFibers = 0;
	G4double xshift = -FiberEnvelopeRadius*(fHilbertCurveLayer-1);
	G4double yshift = -FiberEnvelopeRadius*(fHilbertCurveLayer-1);
	G4double zshift = -FiberEnvelopeRadius*(fHilbertCurveLayer-1);

	for (G4int jj = 0; jj < fHilbertCurveLayer; jj++)
	{		
		G4int max_loopthelayer = fiberPosX.size();
		if (fOnlyBuildOneHistone)
			max_loopthelayer = 2;
		for (G4int loopthelayer = 1; loopthelayer < max_loopthelayer; loopthelayer++) 
		{
			G4double layershift = FiberEnvelopeRadius*2.01;
			G4int i=0;
			if(jj%2==0) i= loopthelayer;
			if(jj%2==1) i= max_loopthelayer-loopthelayer; 
	
			G4double midpoint_x = jj*layershift + xshift + (fiberPosX[i]+fiberPosX[i-1])/2*scaleFactor;
			G4double midpoint_y = jj*layershift + yshift + (fiberPosY[i]+fiberPosY[i-1])/2*scaleFactor;
			G4double midpoint_z = jj*layershift + zshift + (fiberPosZ[i]+fiberPosZ[i-1])/2*scaleFactor;
			G4ThreeVector* midpoint  = new G4ThreeVector(midpoint_x,midpoint_y,midpoint_z);
			G4ThreeVector direction = G4ThreeVector(fiberPosX[i]-fiberPosX[i-1],fiberPosY[i]-fiberPosY[i-1],fiberPosZ[i]-fiberPosZ[i-1]);
			G4RotationMatrix* rotLoop = new G4RotationMatrix();
			if ((direction.x() == 0) && (direction.y() == 0)) rotLoop->rotateZ(0*deg);
			else if ((direction.y() == 0) && (direction.z() == 0)) rotLoop->rotateY(90*deg); //rotate along x
			else if ((direction.x() == 0) && (direction.z() == 0)) rotLoop->rotateX(90*deg);
			else G4cout << "Two consecutive Hilbert paths have the same direction." << G4endl;

			G4String volumeName = "FiberCylinder"+G4UIcommand::ConvertToString(CountFibers);

			fFiberPhysVolLoop[CountFibers] = CreatePhysicalVolume(
														volumeName,	 //name
														CountFibers,	//copy number
														true,		   //use logical volume multiple times many
														fFiberLogic,	 //logical volume
														rotLoop,		// rotation
														midpoint,	   //translation 
														voxel_logic);   // mother logical volume
			
			if(fCheckOverlap)
			{
				G4cout<<"checking CountFibers="<<CountFibers<<G4endl;
				if( fFiberPhysVolLoop[CountFibers]->CheckOverlaps(1000, 0, false)) 
				{
					ofstream outfile;
					outfile.open("overlap.txt",std::ios::out|std::ios::app);
					outfile<<"In "<<" th subdomain "<<jj+1<<" th layer "<<CountFibers<<" th volume detected overlap\n";
					outfile.close(); 
					G4cout << "Detected Overlap, see file overlap.txt for more details." << G4endl;
				}
			}
			CountFibers++;	 
		}
	}

	G4double basePairsInVoxel = fNumberOfBasePairs * CountFibers;
	CreatePhysicalVolume("DNAContent", voxel_logic, fEnvelopePhys, kXAxis, fnVoxels, param);  
	InstantiateChildren(fEnvelopePhys);
	G4cout << "Number of fibers: " << CountFibers <<  " - Number of base pairs in a voxel: " << basePairsInVoxel << " - Number of voxels: " << fnVoxels << G4endl;
	fNumberOfBasePairs = basePairsInVoxel * fnVoxels;
	G4cout << "DNA Construction done. Number of total base pairs: " << fNumberOfBasePairs << G4endl;
	return fEnvelopePhys;
}

void  TsNucleus::SetBasicInfo()
{
	G4double nVoxelX=fHilbertCurve3DRepeat; 
	G4double nVoxelY=fHilbertCurve3DRepeat;
	G4double nVoxelZ=fHilbertCurve3DRepeat; 
	fnVoxels = nVoxelX*nVoxelY*nVoxelZ;
	
	G4double voxelDimX = fVoxelLength; 
	G4double voxelDimY = fVoxelLength; 
	G4double voxelDimZ = fVoxelLength;
	fVoxelHalfDimX = voxelDimX/2;
	fVoxelHalfDimY = voxelDimY/2;
	fVoxelHalfDimZ = voxelDimZ/2;

	fminX = -nVoxelX*voxelDimX/2; fmaxX = nVoxelX*voxelDimX/2;
	fminY = -nVoxelY*voxelDimY/2, fmaxY = nVoxelY*voxelDimY/2;
	fminZ = -nVoxelZ*voxelDimZ/2; fmaxZ = nVoxelZ*voxelDimZ/2;

}


G4LogicalVolume* TsNucleus::ConstructFiberLogicalVolume(){

	//****************************************************************************
	//				   Cylinder containing Histone & DNA (envelope)		   //
	//****************************************************************************

	G4Tubs* gCylinder = new G4Tubs("Chromatin",
								   0,
								   fFiberRadius,
								   fFiberLength/2.,
								   180 * deg,
								   360.0 * deg);

	fFiberLogic = CreateLogicalVolume("ChromatinFiber", gCylinder);

	G4Colour  blue	(0.0, 0.0, 1.0) ;
	G4VisAttributes* EnvelopeVis = new G4VisAttributes( blue ); //
	EnvelopeVis->SetVisibility(fShowChromatinCylinders);
	EnvelopeVis->SetForceAuxEdgeVisible(true);
	fFiberLogic->SetVisAttributes(EnvelopeVis);

	if(fFillCylindersWithDNA || fOnlyBuildOneHistone){
		//****************************************************************************
		//							  Histones									//
		//****************************************************************************
		//Record the position and rotation of each Histone placed
		vector<pair<G4ThreeVector, G4RotationMatrix*>> HistoneDetails;
		BuildHistones(HistoneDetails, fFiberRadius, fFiberLength/2);

		//****************************************************************************
		//								DNA									   //
		//****************************************************************************
		//Builds DNA around and between the Histones
		if(fFillCylindersWithDNA){
			BuildDNA(HistoneDetails);
		}
	}

//	InstantiateChildren(fEnvelopePhys);
	return fFiberLogic;
}


void TsNucleus::BuildHistones(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
							G4double ChromatinRadius,
							G4double ChromatinLength){

	//****************************************************************************
	//							  Histones
	//****************************************************************************

	//place the histones in a left handed solenoid conformation

	G4int nb_HistPerTurn = 6;
	G4double HistoneRadius = 3.3*nm;
	G4double HistoneLength = (5.7 / 2.0)*nm;

	G4String SubComp = "Histone";

	G4Tubs* gHistone = new G4Tubs("Histone",
									  0,
									  HistoneRadius,
									  HistoneLength,
									  180*deg,
									  360*deg);
	G4LogicalVolume* lHistone = CreateLogicalVolume(SubComp, gHistone);

	//if (fScoreOnHistones)
	//	SetLogicalVolumeToBeSensitive(lHistone);

	HistoneRadius+=3.0*nm;  //extra spacing (3nm) for the double helix

	//Generate path for nucleosomes -- solenoid
	G4double x = ((ChromatinRadius-(HistoneRadius+(1.1*nm))) * sin(0.0));
	G4double y = ((ChromatinRadius-(HistoneRadius+(1.1*nm))) * cos(0.0));
	G4double z = (-ChromatinLength + HistoneRadius);

	G4int nNucleosome = (6.0*2.0*(-z)) / (2.0*HistoneRadius); //how many histones fit along chromatin z
	G4double zStep = (2.0*HistoneRadius)/nb_HistPerTurn;
	G4ThreeVector position (x,y,z);

	G4double theta = 0.0;
	G4double thetaStep = (2.0*pi)/(G4double)nb_HistPerTurn;

	if (fOnlyBuildOneHistone) // for visualization and smaller scale simulations.
		nNucleosome = 1;

	G4int built=0;
	for (G4int i=0;i<nNucleosome;i++){
		theta+=thetaStep;
		position[0]=((ChromatinRadius-HistoneRadius) * sin(theta));
		position[1]=((ChromatinRadius-HistoneRadius) * cos(theta));
		position[2]=z;
		z+=zStep;

		if (position[2]+HistoneRadius>=ChromatinLength){
			nNucleosome=built;
			break;
		}

		G4RotationMatrix *HistoneRotation = new G4RotationMatrix();
		HistoneRotation->rotateZ((-120.0+((G4double)(i)*360.0/(G4double)nb_HistPerTurn))*deg);
		HistoneRotation->rotateY(90.0*deg);

		//record the position and rotation of the histone
		pair<G4ThreeVector,G4RotationMatrix*> Details;
		Details.first=position;
		Details.second=HistoneRotation;
		HistoneDetails.push_back(Details);

		if (fAddHistones)
			CreatePhysicalVolume("Histone", i, true, lHistone, HistoneRotation, &position, fFiberLogic);
		built++;
	}

	if (fAddHistones){
		G4VisAttributes* HistoneVis = new G4VisAttributes(G4Colour(0.,0.,1.));
		HistoneVis->SetVisibility(true);
		HistoneVis->SetForceSolid(true);
		lHistone->SetVisAttributes(HistoneVis);
	}
}


//build the DNA
void TsNucleus::BuildDNA(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails)
{
	//****************************************************************************
	//								DNA
	//****************************************************************************

	G4bool BuildHalfCyl=false;
	G4bool BuildQuartCyl=false;
	G4bool BuildSphere=false;

	if (fDNAModel=="HalfCylinder")
		BuildHalfCyl=true;
	else if (fDNAModel=="Sphere")
		BuildSphere=true;
	else if (fDNAModel=="QuarterCylinder")
		BuildQuartCyl=true;

	SetDNAVolumes(BuildHalfCyl, BuildQuartCyl, BuildSphere);

	vector<G4ThreeVector> DNAPath;
	GenerateDNAPath(HistoneDetails, DNAPath);
	SegmentDNAPath(DNAPath);

	if (BuildSphere){
		PlaceDNASphere(DNAPath);
	} else {
		PlaceDNA(DNAPath);
	}
}


//set up DNA volumes
void TsNucleus::SetDNAVolumes(G4bool BuildHalfCyl,
							G4bool BuildQuartCyl,
							G4bool BuildSphere)
{
	// A choice to build the DNA volumes as Half Cylinders, Quarter Cylinders, or Spheres
	// Using Quarter cylinders by default
	// And a denser water for DNA backbone
	
	//sphere DNA
	if (BuildSphere){
		if (fAddBases){
			G4Sphere* gDNA_base1 = new G4Sphere("DNA_base1",
											0*nm,
											0.208*nm,
											0*deg,
											360*deg,
											0*deg,
											180*deg);
			G4Sphere* gDNA_base2 = new G4Sphere("DNA_base2",
											0*nm,
											0.208*nm,
											0*deg,
											360*deg,
											0*deg,
											180*deg);
			
			lBase1 = CreateLogicalVolume("Base1", gDNA_base1);
			lBase2 = CreateLogicalVolume("Base2", gDNA_base2);
		 }
		 if (fAddBackbones){
			 G4Sphere* gDNA_backbone = new G4Sphere("DNA_backbone",
													0.208*nm,
													0.24*nm,
													0*deg,
													360*deg,
													0*deg,
													180*deg);

			 lBack1 = CreateLogicalVolume("Backbone1", gDNA_backbone);
			 lBack2 = CreateLogicalVolume("Backbone2", gDNA_backbone);
		 }

		// ************************** build hydration shell layer **************************
		if (fAddHydrationShell){
			G4Sphere* gWater = new G4Sphere("DNA_WaterLayer",
											   0.24*nm,
											   0.24*nm+fHydrationShellThickness,
											   0*deg,
											   360*deg,
											   0*deg,
											   180*deg);

			lHydrationShell1 = CreateLogicalVolume("HydrationShell1", gWater);
			lHydrationShell2 = CreateLogicalVolume("HydrationShell2", gWater);

		}
	}

	//half cylinder
	else if (BuildHalfCyl){
		G4double CylinderWidth = (0.34-0.05)/2.0 * nm;
		if (fAddBases){
			G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
										0*nm,
										0.5*nm,
											//										((0.34-0.01)/2.0)*nm,
										CylinderWidth,
										180.*deg,
											180.*deg);
			G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
											0*nm,
											0.5*nm,
											CylinderWidth,
											0.*deg,
											180.*deg);

			lBase1 = CreateLogicalVolume("Base1", gDNA_base1);
			lBase2 = CreateLogicalVolume("Base2", gDNA_base2);
		}
			if (fAddBackbones){
				G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_backbone1",
													0.5*nm,
													1.15*nm,
													CylinderWidth,
													180.*deg,
													180.*deg);
				G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_backbone2",
													0.5*nm,
													1.15*nm,
													CylinderWidth,
													0*deg,
													180.*deg);

				lBack1 = CreateLogicalVolume("Backbone1", gDNA_backbone1);
				lBack2 = CreateLogicalVolume("Backbone2", gDNA_backbone2);
			}
		// ************************** build hydration shell layer **************************
		if (fAddHydrationShell){
			G4Tubs* gWater1 = new G4Tubs("WaterLayer1",
											1.15*nm,
											1.15*nm+fHydrationShellThickness,
											CylinderWidth,
											180.*deg,
											180.*deg);
			G4Tubs* gWater2 = new G4Tubs("WaterLayer2",
											1.15*nm,
											1.15*nm+fHydrationShellThickness,
											CylinderWidth,
											0*deg,
											180*deg);

			lHydrationShell1 = CreateLogicalVolume("HydrationShell1", gWater1);
			lHydrationShell2 = CreateLogicalVolume("HydrationShell2", gWater2);
		}
	}
	//quarter cylinder
	else if (BuildQuartCyl){
		G4double CylinderWidth = (0.34-0.05)/2.0 * nm;
		G4double angle = 105 * deg;
		if (fAddBases){
			G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
											0*nm,
											0.5*nm,
											CylinderWidth,
											135.*deg,
											180.*deg);
			G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
											0*nm,
											0.5*nm,
											CylinderWidth,
											-45.*deg,
											180.*deg);

			lBase1 = CreateLogicalVolume("Base1", gDNA_base1);
			lBase2 = CreateLogicalVolume("Base2", gDNA_base2);
		}
		if (fAddBackbones){
			G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_backbone1",
												0.5*nm,
												1.15*nm,
												CylinderWidth,
												180.*deg,
												angle);
			G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_backbone2",
												0.5*nm,
												1.15*nm,
												CylinderWidth,
												0*deg,
												angle);

			lBack1 = CreateLogicalVolume("Backbone1", gDNA_backbone1);
			lBack2 = CreateLogicalVolume("Backbone2", gDNA_backbone2);
		}
		// ************************** build hydration shell layer **************************
		if (fAddHydrationShell){
			G4Tubs* gWater1 = new G4Tubs("WaterLayer1",
											1.15*nm,
											1.15*nm+fHydrationShellThickness,
											CylinderWidth,
											180.*deg,
											angle);
			G4Tubs* gWater2 = new G4Tubs("WaterLayer2",
											1.15*nm,
											1.15*nm+fHydrationShellThickness,
											CylinderWidth,
											0*deg,
											angle);

			lHydrationShell1 = CreateLogicalVolume("HydrationShell1", gWater1);
			lHydrationShell2 = CreateLogicalVolume("HydrationShell2", gWater2);
		}
	}
		
	/*if (fScoreOnBases) {
		SetLogicalVolumeToBeSensitive(lBase1);
		SetLogicalVolumeToBeSensitive(lBase2);
	}
	if (fScoreOnBackbones) {
		SetLogicalVolumeToBeSensitive(lBack1);
		SetLogicalVolumeToBeSensitive(lBack2);
	}
	if (fScoreOnHydrationShell){
		SetLogicalVolumeToBeSensitive(lHydrationShell1);
		SetLogicalVolumeToBeSensitive(lHydrationShell2);
	}*/

	//visualization
	if(fAddBases)	{
		G4VisAttributes * Base1Vis = new G4VisAttributes(G4Colour(0.02, 0.92, 1.0));
		Base1Vis->SetVisibility(true);
		Base1Vis->SetForceSolid(true);
		lBase1 -> SetVisAttributes(Base1Vis);

		G4VisAttributes * Base2Vis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.8));
		Base2Vis->SetVisibility(true);
		Base2Vis->SetForceSolid(true);
		lBase2 -> SetVisAttributes(Base2Vis);	}
	
	if(fAddBackbones)	{
		G4VisAttributes* Back1Vis = new G4VisAttributes(G4Colour(1.0, 0.1, 0.02));
		Back1Vis->SetVisibility(true);
		Back1Vis->SetForceSolid(true);
		lBack1->SetVisAttributes(Back1Vis);

		G4VisAttributes * Back2Vis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.2));
		Back2Vis->SetVisibility(true);
		Back2Vis->SetForceSolid(true);
		lBack2 -> SetVisAttributes(Back2Vis);

	}
	
	if(fAddHydrationShell)	{
		G4VisAttributes * Water1Vis = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
		Water1Vis->SetVisibility(true);
		Water1Vis->SetForceSolid(true);
		lHydrationShell1 -> SetVisAttributes(Water1Vis);

		G4VisAttributes * Water2Vis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
		Water2Vis->SetVisibility(true);
		Water2Vis->SetForceSolid(true);
		lHydrationShell2 -> SetVisAttributes(Water2Vis);
	}
}


//generate a path for DNA around the histones
void TsNucleus::GenerateDNAPath(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
							  vector<G4ThreeVector> &path)
{
	//generate a path of 3vector points that defines the centre of the double helix

	G4double TraceBP=15.0;  //the number of BP to trace back for bezier curve (more = bendier or smoother)
	G4double nCoils = 1.65; //#coils around histone
//	G4double nCoils = 2.38; //#coils around histone
	G4int nHistones=HistoneDetails.size();
	G4double offset = (pi/2.0) - (2.0*pi*nCoils);
	G4double zHist = (5.7 / 2.0)*nm; // Histone Length
	G4double Radius = 3.3*nm + 1.15*nm + 0.16*nm + 0.05*nm; //central radius of double helix around histone = HistRad

	G4int nSteps=200; //go around the histone in 200 steps 

	for (G4int i=0;i<nHistones;i++){
		G4RotationMatrix *rot=HistoneDetails[i].second;

		//histone DNA
		for (G4int j=0;j<nSteps;j++){
			G4double angle = ((G4double)j * 2.0 * pi * nCoils / (G4double)nSteps)+(pi/2.0)+offset;
			//take a pt around the histone
			G4ThreeVector pt( (Radius * sin(angle)),
							 (Radius * cos(angle)),
							 (-zHist + ((G4double)j/(G4double)nSteps)*2.0*(zHist)));

			//rotate the step to the histone
			G4ThreeVector RotPt;
			RotPt[0] = pt.x()*rot->xx() + pt.y()*rot->yx() + pt.z()*rot->zx();
			RotPt[1] = pt.x()*rot->xy() + pt.y()*rot->yy() + pt.z()*rot->zy();
			RotPt[2] = pt.x()*rot->xz() + pt.y()*rot->yz() + pt.z()*rot->zz();

			//translate to the histone position
			RotPt += HistoneDetails[i].first;

			path.push_back(RotPt);
		}


		//link to nextPt -- off histone (start from end of histone i -> start of histone i+1)
		if (i!=nHistones-1){
			G4ThreeVector start = path[path.size()-1]; //last pt on histone

			//calc 1st point on next histone
			G4double angle = (pi/2.0)+offset;
			G4ThreeVector pt( (Radius * sin(angle)),
							 (Radius * cos(angle)),
							 -zHist );
			G4RotationMatrix *nextRot = HistoneDetails[i+1].second;
			G4ThreeVector RotPt;
			RotPt[0] = pt.x()*nextRot->xx() + pt.y()*nextRot->yx() + pt.z()*nextRot->zx();
			RotPt[1] = pt.x()*nextRot->xy() + pt.y()*nextRot->yy() + pt.z()*nextRot->zy();
			RotPt[2] = pt.x()*nextRot->xz() + pt.y()*nextRot->yz() + pt.z()*nextRot->zz();
			G4ThreeVector end = RotPt + HistoneDetails[i+1].first; //1st pt on next histone

			//calc 2nd point on next histone
			angle= (1.0 * 2.0 * pi * nCoils / (G4double)nSteps)+(pi/2.0)+offset;
			G4ThreeVector secpt( (Radius * sin(angle)),
								(Radius * cos(angle)),
								(-zHist + (1.0/(G4double)nSteps)*2.0*(zHist)));
			RotPt[0] = secpt.x()*nextRot->xx() + secpt.y()*nextRot->yx() + secpt.z()*nextRot->zx();
			RotPt[1] = secpt.x()*nextRot->xy() + secpt.y()*nextRot->yy() + secpt.z()*nextRot->zy();
			RotPt[2] = secpt.x()*nextRot->xz() + secpt.y()*nextRot->yz() + secpt.z()*nextRot->zz();
			secpt = RotPt + HistoneDetails[i+1].first; //2nd pt on next histone

			G4ThreeVector StartDir = (start - path[path.size()-2]).unit(); //direction coming off histone
			G4ThreeVector MidPoint1 = ((TraceBP*0.3)*StartDir); //trace it forward x bp
			MidPoint1[0] *= nm;
			MidPoint1[1] *= nm;
			MidPoint1[2] *= nm;
			MidPoint1 += start;

			G4ThreeVector EndDir = (secpt-end).unit();
			G4ThreeVector MidPoint2 = ((TraceBP*0.3)*EndDir); //trace it forward x bp
			MidPoint2[0] *= nm;
			MidPoint2[1] *= nm;
			MidPoint2[2] *= nm;
			MidPoint2 = end - MidPoint2;

			//link the 4 points with Bezier in 100 steps
			Bezier(start, MidPoint1, MidPoint2, end, path, 100);
		}
	}
}


//Bezier curve to smoothly join 2 points (bending through 2 other points)
void TsNucleus::Bezier(G4ThreeVector &start,
					 G4ThreeVector &MidPoint1,
					 G4ThreeVector &MidPoint2,
					 G4ThreeVector &end,
					 vector<G4ThreeVector> &path,
					 G4int nSteps)
{
	G4double j=0.0;
	for (G4int k=0;k<nSteps-1;k++){
		j+=1.0/(G4double)nSteps;
		G4ThreeVector BezPt = pow(1.0-j,3.0)*start + 3.0*pow(1.0-j,2.0)*j*MidPoint1 +
		3.0*(1.0-j)*j*j*MidPoint2 + j*j*j*end;
		path.push_back(BezPt);
	}
}


//split the DNA path into 0.34nm steps
void TsNucleus::SegmentDNAPath(std::vector<G4ThreeVector> &path)
{
	std::vector<G4ThreeVector> newPath;
	double rise=0.34*nm;
	int nPts=path.size();
	int counter=0;
	newPath.push_back(path[0]);
	for (int i=0;i<nPts-1;i++){
		G4ThreeVector vector = (path[i+1]-newPath[counter]);
		double length = vector.mag();
		int nDiv = length/rise;
		G4ThreeVector unit = vector.unit();
		unit *= rise;
		for (int j=0;j<nDiv;j++){
			G4ThreeVector nuwe = newPath[counter] + unit;
			newPath.push_back(nuwe);
			counter++;
		}
	}
}


//use DNA path to place sphere DNA volumes
void TsNucleus::PlaceDNASphere(vector<G4ThreeVector> &newPath)
{
	G4double helixRadius = 1.0*nm;
	G4double rotPair = ((2.0*pi)/10.0);   //10bp per turn
	G4int nBP=newPath.size();
	G4double rBack=helixRadius - 0.24*nm;
	G4double rBase=rBack - 0.24*nm - 0.208*nm;

	for (int bp=0; bp<nBP-1; bp++){
		fNumberOfBasePairs++;
		G4double angle1 = (G4double)bp * rotPair;
		G4double angle2 = angle1+pi;// + (120.0*pi/180.0); //offset for strand2 (major and minor groove)

		//temporary positions
		G4ThreeVector back1temp = G4ThreeVector((rBack*cos(angle1)), (rBack*sin(angle1)), 0.0);
		G4ThreeVector back2temp = G4ThreeVector((rBack*cos(angle2)), (rBack*sin(angle2)), 0.0);
		G4ThreeVector base1temp = G4ThreeVector((rBase*cos(angle1)), (rBase*sin(angle1)), 0.0);
		G4ThreeVector base2temp = G4ThreeVector((rBase*cos(angle2)), (rBase*sin(angle2)), 0.0);

		//Rotation to point to next plane
		G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
		G4ThreeVector norm (0.,0.,-1.); //the normal to the plane (G4 build planes facing -z)
		G4double DotProd = norm.dot(vecNext);
		G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
		G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm

		//set up new 3Vectors for rotated pos
		G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.), hydration1(0.,0.,0.), hydration2(0.,0.,0.);

		//Apply rotation
		G4RotationMatrix *rot = new G4RotationMatrix;
		rot->rotate(AngBetween, cross);

		ApplyRotation(back1, back1temp, rot);
		ApplyRotation(back2, back2temp, rot);
		ApplyRotation(base1, base1temp, rot);
		ApplyRotation(base2, base2temp, rot);
		delete rot;

		//Translate
		base1+=newPath[bp];
		base2+=newPath[bp];
		back1+=newPath[bp];
		back2+=newPath[bp];
		hydration1+=newPath[bp];
		hydration2+=newPath[bp];

		G4int bpID=bp+1;

		if(fAddBases) {
			G4ThreeVector *posBase1 = &base1;
			CreatePhysicalVolume("Base1_", bpID, true, lBase1, rot, posBase1, fFiberLogic);
			G4ThreeVector *posBase2 = &base2;
			CreatePhysicalVolume("Base2_", bpID, true, lBase2, rot, posBase2, fFiberLogic);
		}

		if(fAddBackbones) {
			G4ThreeVector *posBack1 = &back1;
			CreatePhysicalVolume("Backbone1_", bpID, true, lBack1, rot, posBack1, fFiberLogic);  // ???? false?
			G4ThreeVector *posBack2 = &back2;
			CreatePhysicalVolume("Backbone2_", bpID, true, lBack2, rot, posBack2, fFiberLogic);
		}

		if(fAddHydrationShell) {
			G4ThreeVector *posHydration1 = &hydration1;
			G4ThreeVector *posHydration2 = &hydration2;
			CreatePhysicalVolume("HydrationShell1_", bpID, true, lHydrationShell1, rot, posHydration1, fFiberLogic);
			CreatePhysicalVolume("HydrationShell2_", bpID, true, lHydrationShell2, rot, posHydration2, fFiberLogic);
		}
	}
	G4cout<<G4endl<<"Spherical DNA BUILT: "<<(G4double)nBP/1000.0<<" kbp"<<G4endl<<G4endl;
	fFiberDNAContent = (G4double)nBP;
}

//rotate a 3vector by a rot matrix and return new coordinates
void TsNucleus::ApplyRotation(G4ThreeVector &rotated, G4ThreeVector &vector, G4RotationMatrix *rot)
{
	rotated[0] = vector[0]*rot->xx() + vector[1]*rot->yx() + vector[2]*rot->zx();
	rotated[1] = vector[0]*rot->xy() + vector[1]*rot->yy() + vector[2]*rot->zy();
	rotated[2] = vector[0]*rot->xz() + vector[1]*rot->yz() + vector[2]*rot->zz();

}

//Use DNA path to place halfCyl or quartCyl DNA volumes
void TsNucleus::PlaceDNA(vector<G4ThreeVector> &newPath)
{
	G4double rotPair = ((2.0*pi)/10.0);   //10bp per turn
	G4int nBP=newPath.size();

	for (G4int bp=0; bp<nBP-1; bp++){
		fNumberOfBasePairs++;
		//Position of base + back in xy
		//Definitely gives right handed coil (checked) -- left handed in -ve z?
		G4double angle1 = -(G4double)bp * rotPair;

		//Rotation to point to next plane
		G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
		G4ThreeVector norm (0.,0.,1.); //the normal to the plane (G4 build planes facing -z)
		G4double DotProd = vecNext.dot(norm);
		G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
		G4ThreeVector cross = (vecNext.cross(norm)); //vector perp to vecnext and norm

		//set up new 3Vectors for rotated pos
		G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.), hydration1(0.,0.,0.), hydration2(0.,0.,0.);

		//Apply rotation
		G4RotationMatrix *rot1 = new G4RotationMatrix(); 
		rot1->rotate(AngBetween, cross);  // this causes a strange behavior when AngBetween close to pi, to be resolved ... easiest to see by setting angle1 = 0.
		rot1->rotateZ(angle1);

		//Translate
		back1+=newPath[bp];
		back2+=newPath[bp];
		base1+=newPath[bp];
		base2+=newPath[bp];
		hydration1+=newPath[bp];
		hydration2+=newPath[bp];

		int bpID=bp+1;

		if(fAddBases){
			G4ThreeVector *posBase1 = &base1;
			CreatePhysicalVolume("Base1_", bpID, true, lBase1, rot1, posBase1, fFiberLogic);
			G4ThreeVector *posBase2 = &base2;
			CreatePhysicalVolume("Base2_", bpID, true, lBase2, rot1, posBase2, fFiberLogic);
		}

		if(fAddBackbones){
			G4ThreeVector *posBack1 = &back1;
			CreatePhysicalVolume("Backbone1_", bpID, true, lBack1, rot1, posBack1, fFiberLogic);  // ???? false?
			G4ThreeVector *posBack2 = &back2;
			CreatePhysicalVolume("Backbone2_", bpID, true, lBack2, rot1, posBack2, fFiberLogic);
		}

		if(fAddHydrationShell) {
			G4ThreeVector *posHydration1 = &hydration1;
			G4ThreeVector *posHydration2 = &hydration2;
			CreatePhysicalVolume("HydrationShell1_", bpID, true, lHydrationShell1, rot1, posHydration1, fFiberLogic);
			CreatePhysicalVolume("HydrationShell2_", bpID, true, lHydrationShell2, rot1, posHydration2, fFiberLogic);
		}
	}

	G4cout<<G4endl<<"Cylindrical DNA BUILT: "<<(G4double)nBP/1000.0<<" kbp"<<G4endl<<G4endl;
	fFiberDNAContent = (G4double)nBP;
}

