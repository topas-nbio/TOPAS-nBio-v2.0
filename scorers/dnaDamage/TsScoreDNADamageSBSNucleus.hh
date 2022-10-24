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

#ifndef TsScoreDNADamageSBSNucleus_hh
#define TsScoreDNADamageSBSNucleus_hh

#include "TsScoreDNADamageSBS.hh"

class TsScoreDNADamageSBSNucleus : public TsScoreDNADamageSBS
{
public:
	TsScoreDNADamageSBSNucleus(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreDNADamageSBSNucleus();

	void UserHookForEndOfRun();
	std::pair<G4int, G4int> CalculateChromosomeAndBasePairID(std::vector<G4int> hids);
	void ReadGeometryAndCopyNumberMaps();

private:
	G4String fGeometryInfoFile;
	G4String fMapOfVoxelCopyNumbersInsideNucleusFile;
	G4String fMapChromosomeVoxelFile;

	G4int fBasePairsInAVoxel;
	G4int fBasePairsInAFiber;
	G4int fnVoxelIs3DRepeated;

	std::vector<G4int> fVoxelIDsInWorld;
	std::vector<G4int> fVoxelIDsInNucleus;

	std::vector<G4int> fChromosomeIDForEachVoxel;
	std::vector<G4int> fVoxelIDs;
};

#endif
