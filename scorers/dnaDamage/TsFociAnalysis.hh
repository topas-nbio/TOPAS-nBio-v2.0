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
// Authors: Alejandro Bertolet, Jan Schuemann

#ifndef TsFociAnalysis_hh
#define TsFociAnalysis_hh

#include "TsVNtupleScorer.hh"

class TsFociAnalysis
{
public:
	TsFociAnalysis(TsVGeometryComponent* component);
	~TsFociAnalysis();

	std::vector<G4int> GetNumberOfFoci(std::vector<G4ThreeVector> dsbPositions);
	G4float GetDistance(G4ThreeVector a, G4ThreeVector b);
	G4bool CheckIfAnyIndexIsAvailable(std::vector<G4bool> indexIsAvailable);

	void Produce3DImage(std::vector<G4ThreeVector> dsbPositions);
	void Produce2DImages(std::vector<G4ThreeVector> dsbPositions);

	G4float Gaussian3D(G4float x, G4float y, G4float z, G4float sigma);
	G4double Gaussian2D(G4double x, G4double y, G4double sigma);

	void inline SetFociSizes(std::vector<G4double> v)				{ fFociSizes = v; }
	void inline SetPlanesFor2DFociImages(std::vector<G4String> v)	{ f2DPlanesForFociImage = v; }
	void inline SetPSFShape(G4String v)								{ fMicroscopePSFShape = v; }
	void inline SetPSFWidth(G4double v)								{ fMicroscopePSFWidth = v; }
	void inline Set2DImageResolutions(std::vector<G4double> v)		{ f2DResolutions = v; }
	void inline Set3DImageResolution(G4double v)					{ f3DResolution = v; }
	void inline SetMinX(G4double v)									{ fxmin = v; }
	void inline SetMaxX(G4double v)									{ fxmax = v; }
	void inline SetMinY(G4double v)									{ fymin = v; }
	void inline SetMaxY(G4double v)									{ fymax = v; }
	void inline SetMinZ(G4double v)									{ fzmin = v; }
	void inline SetMaxZ(G4double v)									{ fzmax = v; }

private:
	std::vector<G4double> fFociSizes;
	std::vector<G4String> f2DPlanesForFociImage;
	G4String fMicroscopePSFShape;
	G4double fMicroscopePSFWidth;
	std::vector<G4double> f2DResolutions;
	G4double f3DResolution;
	G4double fxmin, fxmax, fymin, fymax, fzmin, fzmax;

	TsVGeometryComponent* fComponent;
};

#endif
