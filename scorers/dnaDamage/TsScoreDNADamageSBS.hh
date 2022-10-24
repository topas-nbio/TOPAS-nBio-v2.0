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


#ifndef TsScoreDNADamageSBS_hh
#define TsScoreDNADamageSBS_hh

#include "TsVNtupleScorer.hh"

#include "TsDNADamageCalculator.hh"
#include "TsHitInDNA.hh"
#include "TsFociAnalysis.hh"

class TsScoreDNADamageSBS : public TsVNtupleScorer
{
public:
	TsScoreDNADamageSBS(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	virtual ~TsScoreDNADamageSBS();

	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	void AccumulateEvent();
	void AbsorbResultsFromWorkerScorer(TsVScorer*);
	virtual void UserHookForEndOfRun();
	G4int Analyze(std::vector<TsHitInDNA*> hits, G4int eventID);
	void CalculateYields();
	G4double CalculateDoseInGray(G4double edep);

	void inline AddHierarchyLevel(G4String level)	{ fHierarchicalLevels.push_back(level); }
	virtual std::pair<G4int, G4int> CalculateChromosomeAndBasePairID(std::vector<G4int> hids);
	virtual std::pair<G4int, G4int> GetDNAComponentAndStrandID(G4TouchableHistory* touchable);

protected:
	G4int fNumberOfHistoriesInRun;
	// For Material-based filter
	G4int fBasePairDepth;
	std::vector<G4Material*> fStrand1Materials;
	std::vector<G4Material*> fStrand2Materials;

	// For creating a scoring radius
	G4double fScoringRadius;

	// For direct damage
	G4double fDirectDamageThreshold;
	G4bool fUseLinearProbabilityForDirectDamage;
	G4double fLowerLimitLinearProbability;
	G4double fUpperLimitLinearProbability;
	// For quasi-direct damage
	G4double fProbabilityOfChargeTransferFromHydrationShellToBackbone;
	// For indirect damage
	G4bool fAlwaysScavengeSpeciesInDNAComponents;
	G4double fProbabilityOfScavengingInBackbone;
	G4double fProbabilityOfScavengingInBase;
	G4double fProbabilityOfDamageInBackbone;
	G4double fProbabilityOfDamageInBase;
	G4bool fScavengeInHistones;

	// For defining types of damage to be accounted for
	G4int fNumberOfBasePairsForDSB;

	// Foci scoring, creation of foci images
	G4bool fScoreFoci;
	std::vector<G4double> fFociSizes;
	G4bool fGet3DFociImage;
	G4bool fGet2DFociImage;
	std::vector<G4String> f2DPlanesForFociImage;
	G4String fMicroscopePSFShape;
	G4double fMicroscopePSFWidth;
	std::vector<G4double> f2DImageResolutions;
	G4double f3DImageResolution;
	G4double fImXmin, fImXmax, fImYmin, fImYmax, fImZmin, fImZmax;

	G4bool fExcludeShortFragments;
	G4int fLowerThresholdForFragmentDetection;
	G4int fUpperThresholdForFragmentDetection;

	// Options for the output
	G4bool fWriteCSVWithExtensiveDamage;
	G4bool fScoreDirectDamage;
	G4bool fScoreIndirectDamage;
	G4bool fScoreQuasiDirectDamage;
	G4bool fScoreOnBases;
	G4bool fScoreOnBackbones;
	G4bool fBreakDownPerDamageOrigin;
	std::vector<G4int> fChromosomeContents;

	// For SDD specification
	G4double fDosePerExposure;
	G4bool fOnlyIncludeDSBinSDD;
	G4bool fWriteMinimalSDDOutput;
	G4String fPrimaryParticle;

	// Scoring of physical quantities
	G4double fEdep;
	G4double fTrackAveragedLET;
	G4double fDoseInThisExposure;
	G4double fAccumulatedDoseInRun;
	G4int fExposureID;

	// Stop tracking at a given dose
	G4double fStopAtDose;

	// Quantification of damage
	G4int fNumSB;
	G4int fNumSBDirect;
	G4int fNumSBQuasiDirect;
	G4int fNumSBIndirect;
	G4int fNumSSB;
	G4int fNumSSBDirect;
	G4int fNumSSBQuasiDirect;
	G4int fNumSSBIndirect;
	G4int fNumDSB;
	G4int fNumDSBDirect;
	G4int fNumDSBIndirect;
	G4int fNumDSBDirectIndirect;
	G4int fNumDSBDirectQuasiDirect;
	G4int fNumDSBQuasiDirectQuasiDirect;
	G4int fNumDSBIndirectQuasiDirect;
	G4int fNumBaseDamage;
	G4int fNumBaseDamageDirect;
	G4int fNumBaseDamageQuasiDirect;
	G4int fNumBaseDamageIndirect;
	G4int fNumSSBPlus;
	G4int fNumDSBPlus;
	G4int fNumDSBComplex;

	// Quantification for foci (considers up to 5 different foci sizes)
	G4int fNumFoci1, fNumFoci2, fNumFoci3, fNumFoci4, fNumFoci5;
	std::vector<G4ThreeVector> fDSBPositionsInRun;

	// Yields
	G4double fYBaseDam;
	G4double fYSB;
	G4double fYSSB;
	G4double fYDSB;
	G4double fYSSBPlus;
	G4double fYDSBPlus;
	G4double fYDSBComplex;

	// Accumulation of events
	std::vector<G4double> fEventsEdep;

	// Map to check if a track enters in a volume
	std::map<G4double, G4int> fTrackSteps;
	std::map<G4double, G4String> fTrackLastVolume;

	// Hits recorded
	std::vector<TsHitInDNA*> fHits;

	// Collections of hits for each event
	std::vector<std::vector<TsHitInDNA*>> fCollectionsOfHits;

	// Auxiliary class for calculate damage
	TsDNADamageCalculator* fDamageCalculator;
	TsFociAnalysis* fFociAnalyzer;

	// Vector of geometric hierarchical levels
	std::vector<G4String> fHierarchicalLevels;

	// Parameter for SDD header
	G4String fAuthor;
	G4String fSimulationDetails;
	G4String fSourceDetails;
	G4int fSourceType;
	G4double fMeanEnergy;
	G4String fEnergyDist;
	G4String fIrrTarget;
	G4String fCellCycle;
	G4String fDNAStructure;
	G4int fInVitroOrInVivo;
	G4String fProliferationStatus;
	G4String fMicroenvironment;
	G4double fTime;
	G4String fAddInfo;

	// Codes for components ID
	G4int histone = -1;
	G4int base = 0;
	G4int backbone = 1;
	G4int hydrationshell = 2;

	// Codes for types of damage
	G4int nodamage = -1;
	G4int direct = 1;
	G4int indirect = 2;
	G4int multiple = 3;
	G4int quasidirect = 4;
	G4int multiplewithquasidirect = 5;
};



#endif
