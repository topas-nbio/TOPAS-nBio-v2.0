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

#ifndef TsGeneratorUniformlyRandomAtVolume_hh
#define TsGeneratorUniformlyRandomAtVolume_hh

#include "TsVGenerator.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class TsGeneratorUniformlyRandomAtVolume : public TsVGenerator
{
public:
	TsGeneratorUniformlyRandomAtVolume(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName);
	~TsGeneratorUniformlyRandomAtVolume();
	
	void ResolveParameters();
	void GeneratePrimaries(G4Event* );
	
private:
	G4bool fGenerateOnSurface;
	G4ThreeVector fPos;
	G4ThreeVector fCenter;
	G4String fSolid;
	G4String fTheDirection;
};
#endif
