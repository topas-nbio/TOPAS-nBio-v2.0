//
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//

#ifndef TsNucleus_hh
#define TsNucleus_hh

#include "TsVGeometryComponent.hh"
#include "G4PVParameterised.hh"
#include "TsVoxelParameterisation.hh"
#include "G4NistManager.hh"

class TsNucleus : public TsVGeometryComponent
{	
public:
	TsNucleus(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsNucleus();
	
	G4VPhysicalVolume* Construct();
	void SetBasicInfo();
	void UpdateForNewRun(G4bool force);
	
private:
	G4LogicalVolume* ConstructFiberLogicalVolume();
	void BuildHistones(std::vector<std::pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
					   G4double ChromatinRadius,
					   G4double ChromatinLength,
					   G4VPhysicalVolume* vol = NULL);

	void BuildDNA(std::vector<std::pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails);
	void SetDNAVolumes(G4bool BuildHalfCyl,
					   G4bool BuildQuartCyl,
					   G4bool BuildSphere);
	void GenerateDNAPath(std::vector<std::pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
						 std::vector<G4ThreeVector> &path);
	void SegmentDNAPath(std::vector<G4ThreeVector> &path);
	void PlaceDNASphere(std::vector<G4ThreeVector> &path, G4VPhysicalVolume* vol = NULL);
	void PlaceDNA(std::vector<G4ThreeVector> &path, G4VPhysicalVolume* vol = NULL);
	G4double GetFiberDNAContent() {return fFiberDNAContent;}
	void ApplyRotation(G4ThreeVector &Rotated,
					   G4ThreeVector &Position,
					   G4RotationMatrix *Rot);
	void Bezier(G4ThreeVector &start,
				G4ThreeVector &MidPoint1,
				G4ThreeVector & MidPoint2,
				G4ThreeVector &end,
				std::vector<G4ThreeVector> &path,
				G4int nSteps);

	G4LogicalVolume * fFiberLogic;
	std::vector<G4VPhysicalVolume*> fFiberPhysVolLoop;
	std::vector<G4double> fiberPosX;
	std::vector<G4double> fiberPosY;
	std::vector<G4double> fiberPosZ;

	G4LogicalVolume *lBase1, *lBase2, *lBack1, *lBack2;
	G4LogicalVolume *lHydrationShell1, *lHydrationShell2;

	// control parameters
	G4double fFiberRadius;
	G4double fFiberLength;
	G4double fNucleusRadius ;
	G4double fVoxelLength;
	G4int fHilbertCurveLayer;
	G4int fHilbertCurve3DRepeat;   
	G4bool fCheckOverlap;
	G4bool fFillCylindersWithDNA;
	G4bool fAddHistones;
	G4bool fOnlyBuildOneHistone;
	G4bool fAddBases;
	G4bool fAddBackbones;
	G4bool fAddHydrationShell;
	G4double fFiberDNAContent;
	G4double fHydrationShellThickness; 
	G4String fDNAModel;
	G4bool fShowNucleus, fShowDNAVoxels, fShowChromatinCylinders;
	G4double fNumberOfBasePairs;
	G4bool fRotateNucleusForEachRun;

	//ParameterisationInfo
	TsVoxelParameterisation* param;
	G4int fnVoxels;
	G4double fminX, fmaxX, fminY, fmaxY, fminZ, fmaxZ;
	G4double fVoxelHalfDimX,  fVoxelHalfDimY, fVoxelHalfDimZ;

	G4String fWater;
	G4String fBaseMaterialName;
	G4String fBackboneMaterialName;
	G4String fHistoneMaterialName;

	// Selecting Scoring components
	//G4bool fScoreOnBases, fScoreOnBackbones, fScoreOnHydrationShell, fScoreOnHistones;
};

#endif
