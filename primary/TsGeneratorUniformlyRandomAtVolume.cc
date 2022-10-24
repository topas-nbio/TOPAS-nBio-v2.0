// Particle Generator for VolumeOrSurface
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

#include "TsGeneratorUniformlyRandomAtVolume.hh"

#include "TsParameterManager.hh"
#include "TsVGenerator.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VSolid.hh"

#include "G4VoxelLimits.hh"
#include "G4BoundingEnvelope.hh"
#include "G4AffineTransform.hh"

#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include <stdlib.h>

TsGeneratorUniformlyRandomAtVolume::TsGeneratorUniformlyRandomAtVolume(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName)
:TsVGenerator(pM, gM, pgM, sourceName)
{
	ResolveParameters();
}


TsGeneratorUniformlyRandomAtVolume::~TsGeneratorUniformlyRandomAtVolume()
{
}


void TsGeneratorUniformlyRandomAtVolume::ResolveParameters() {
	TsVGenerator::ResolveParameters();
	
	fSolid  = fPm->GetStringParameter(GetFullParmName("Component"));
	
	G4String generateOn = "volume";
	if ( fPm->ParameterExists(GetFullParmName("DistributePrimariesAtVolumeRegion")))
		generateOn = fPm->GetStringParameter(GetFullParmName("DistributePrimariesAtVolumeRegion"));
	generateOn.toLower();
	
	if ( generateOn != "volume" && generateOn != "surface" ) {
		G4cerr << "Topas is exiting due to a serious error in source setup." << G4endl;
		G4cerr << "Parameter name: " << GetFullParmName("DistributePrimariesAtVolumeRegion") << G4endl;
		G4cerr << "Define the region to distribute the primaries: volume or surface" << G4endl;
		fPm->AbortSession(1);
	}
	
	fGenerateOnSurface = generateOn == "surface" ? true : false;
	fTheDirection = "isotropic";
	if ( fPm->ParameterExists(GetFullParmName("DirectPrimariesTo")))
		fTheDirection = fPm->GetStringParameter(GetFullParmName("DirectPrimariesTo"));
	fTheDirection.toLower();
	
	if ( fTheDirection != "inside" && ( fTheDirection != "outside" && fTheDirection != "isotropic")) {
		G4cerr << "Topas is exiting due to a serious error in source setup." << G4endl;
		G4cerr << "Parameter name: " << GetFullParmName("DirectPrimariesTo") << G4endl;
		G4cerr << "Specify the direction of the primaries: inside, outside or isotropic" << G4endl;
		fPm->AbortSession(1);
	}
}


void TsGeneratorUniformlyRandomAtVolume::GeneratePrimaries(G4Event* anEvent) {
	if (CurrentSourceHasGeneratedEnough())
		return;
	
	TsPrimaryParticle p;
	
	p.particleDefinition = fParticleDefinition;
	p.isNewHistory = true;
	p.weight = 1;
	
	p.isOpticalPhoton = fIsOpticalPhoton;
	p.isGenericIon = fIsGenericIon;
	p.ionCharge = fIonCharge;
	
	G4VSolid* solid = G4LogicalVolumeStore::GetInstance()->GetVolume(fSolid)->GetSolid();
	
	G4double x, y, z;
	
	G4ThreeVector insideVolume, onSurface;
	
	// A random point inside the volume to determine direction
	G4VoxelLimits voxelLimits;
	G4AffineTransform dontUse;
	G4double xmin, xmax, ymin, ymax, zmin, zmax;
	
	solid->CalculateExtent(kXAxis, voxelLimits, dontUse, xmin, xmax);
	solid->CalculateExtent(kYAxis, voxelLimits, dontUse, ymin, ymax);
	solid->CalculateExtent(kZAxis, voxelLimits, dontUse, zmin, zmax);
	
	while (1) {
		x = G4RandFlat::shoot(xmin, xmax);
		y = G4RandFlat::shoot(ymin, ymax);
		z = G4RandFlat::shoot(zmin, zmax);
		if ( solid->Inside(G4ThreeVector(x,y,z)) == kInside )
			break;
	}
	
	insideVolume = G4ThreeVector(x,y,z);
	onSurface = solid->GetPointOnSurface();
	
	// In case a energy spetrum is defined.
	if ( fUseSpectrum ) {
		if (fSpectrumIsContinuous) {
			G4double aRandom = G4UniformRand() * fSpectrumWeightSums[fSpectrumNBins - 1];
			G4int j = fSpectrumNBins - 2;
			while ((fSpectrumWeightSums[j] > aRandom) && (j > 0)) j--;
			p.kEnergy = fSpectrumEnergies[j];
			
			if (fSpectrumSlopes[j] != 0.) {
				G4double b = fSpectrumWeights[j] / fSpectrumSlopes[j];
				G4double c = 2 * (aRandom - fSpectrumWeightSums[j]) / fSpectrumSlopes[j];
				G4double delta = b * b + c;
				
				G4int sign = 1;
				if (fSpectrumSlopes[j] < 0.) sign = -1;
				
				p.kEnergy += sign * sqrt(delta) - b;
			} else if (fSpectrumWeights[j] > 0.) {
				p.kEnergy += (aRandom - fSpectrumWeightSums[j]) / fSpectrumWeights[j];
			}
		} else {
			G4double aRandom = G4UniformRand();
			G4int j = fSpectrumNBins - 1;
			while ((fSpectrumBinTops[j] >= aRandom) && (j >= 0)) j--;
			p.kEnergy = fSpectrumEnergies[j+1];
		}
	} else {
		p.kEnergy = fEnergy;
	}
	
	G4ThreeVector startDirection;
	
	if ( fTheDirection == "isotropic" ) {
		startDirection = G4RandomDirection();
	} else if ( fTheDirection == "outside") {
		startDirection = (onSurface - insideVolume).unit();
	} else {
		startDirection = (insideVolume - onSurface).unit();
	}
	
	if ( fGenerateOnSurface ) {
		p.posX = onSurface.x();
		p.posY = onSurface.y();
		p.posZ = onSurface.z();
	} else {
		p.posX = insideVolume.x();
		p.posY = insideVolume.y();
		p.posZ = insideVolume.z();
	}
	
	p.dCos1 = startDirection.x();
	p.dCos2 = startDirection.y();
	p.dCos3 = startDirection.z();
	
	TransformPrimaryForComponent(&p);
	GenerateOnePrimary(anEvent, p);
	AddPrimariesToEvent(anEvent);
}
