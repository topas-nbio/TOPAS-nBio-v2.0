// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// Authors: Hongyu Zhu, Jan Schuemann, Alejandro Bertolet

#ifndef TsScorerNucleusDNADamage_hh
#define TsScorerNucleusDNADamage_hh

#include "TsVNtupleScorer.hh"
#include "TsHitsRecord.hh"
#include "TsDefineDamage.hh"

class TsScorerNucleusDNADamage : public TsVNtupleScorer
{
public:
	TsScorerNucleusDNADamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScorerNucleusDNADamage();

	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	void AccumulateEvent();
	void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
	void UserHookForEndOfRun();

	void GetGeometryInfo();
	void FindPlaceInChromosome(G4int VoxelID, G4int FiberID, G4int bpIDinFiber,G4int &VoxelNumInSphere, G4int &ChromosomeID, G4int &bpIDinChromosome);
	G4int CalculateVoxelID(G4ThreeVector position);
	G4int Analyze(std::vector<TsHitsRecord*> hits, G4int eventID);

private:
	G4int fNumberOfHistoriesInRun;
	G4double fProbabilityOfOHDamage;
	G4double fDamageThreshold;
	G4bool fUseLinearProbabilityThreshold;
	G4double fLinearProbability_lower_limit;
	G4double fLinearProbability_upper_limit;
	G4int fDSBSeparation;
	G4bool fExcludeShortFragment;
	G4int fLowerFragmentDetectionThreshold;
	G4int fUpperFragmentDetectionThreshold;
	G4bool fOnlyIncludeDSBinSDD;
	G4bool fWriteCSV;
	G4bool fScoreDirectDamages;
	G4bool fScoreIndirectDamages;
	G4bool fScoreOnBases = true;
	G4bool fScoreOnBackbones = true;
	G4bool fScoreOnHydrationShell = true;
	G4double fDosePerExposure;
	G4double fScoringRadius;
	G4bool fHistoneAsScavenger;
	G4String fPrimaryParticle;
	G4bool fMinimalSDDOutput;
	G4int fBasePairDepth;
	std::vector<G4Material*> fStrand1Materials;
	std::vector<G4Material*> fStrand2Materials;

	G4double fTrackAveragedLET;
	G4double fTravelDistance;
	G4double fEdep;
	G4double fDoseInThisExposure;
	G4double fEdepkeV;
	G4double fDoseGy;
	G4double fLETtkeVum;
	G4int fExposureID;
	G4bool fCalZetaBeta;
	G4double fZetaBeta_sq;

	// quantification of damage
    G4int numSSB;
    G4int numSSB_dir;
    G4int numSSB_indir;
    G4int numDSB;
    G4int numDSB_dir;
    G4int numDSB_indir;
    G4int numDSB_hybrid;
    G4int numBaseDam;
    G4int numBaseDam_dir;
    G4int numBaseDam_indir;
    G4int numSSBPlus;
    G4int numDSBPlus;
    G4int numMoreComplex;
    G4int Excluded_numSSB;
    G4int Excluded_numSSB_dir;
    G4int Excluded_numSSB_indir;
    G4int Excluded_numDSB;
    G4int Excluded_numDSB_dir;
    G4int Excluded_numDSB_indir;
    G4int Excluded_numDSB_hybrid;
    G4int Excluded_baseDam;
    G4int Excluded_baseDam_dir;
    G4int Excluded_baseDam_indir;
    G4double yieldSSB;
    G4double yieldDSB;
    G4double yieldBaseDam;
    G4double yieldSSBPlus;
    G4double yieldDSBPlus;
    G4double yieldMoreComplex;

	std::vector<TsHitsRecord*> Hits;

	TsDefineDamage * fDefineDamage;
	std::vector< std::vector<TsHitsRecord*>>HitsOfEvents;
	std::vector<G4double> eventsEdep;
	std::vector<G4double> eventsLength;

	// Geometry information
	G4String fGeometryInfo, fCopyNoTable, fSignedCHVoxel;
	G4double fFiberDNAContent;  //bp
	G4double fVoxelDNAContent;  //bp
	G4double fTotalDNAContent;  //bp
	G4double fNucleusMass;	  //g
	G4int	fNumofVoxel;
	G4int	fVoxel3Drepeat;
	G4double fVoxelSize;
	G4double fVoxelContainerminXYZ;

	std::vector<G4int> fChromosomeDNAContent;
	std::vector<G4double> fChromosomeDNAContentSum;

	std::vector<G4int> fVoxelNumInBox;
	std::vector<G4int> fVoxelNumInSphere;

	std::vector<G4int> fCH_ID;
	std::vector<G4int> fVoxel_ID;
};



#endif
