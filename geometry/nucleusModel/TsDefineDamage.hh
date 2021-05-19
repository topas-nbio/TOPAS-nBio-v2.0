// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//
// Authors: Hongyu Zhu, Alexander Klapproth, Jan Schuemann, Alejandro Bertolet

#ifndef TsDefineDamage_hh
#define TsDefineDamage_hh

#include "TsVNtupleScorer.hh"
#include "TsHitsRecord.hh"

class TsDefineDamage
{
public:
	TsDefineDamage();
	~TsDefineDamage();

	// -----------
	// New method
	// -----------
	void ComputeStrandBreaks(std::vector<TsHitsRecord*> Hits);
	void QuantifyDamage(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> pairsOfDSB, std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> SSB, std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> BaseDamage);
	std::map<G4int, std::vector<G4int>> IdentifyDamageBlocks();
	std::vector<G4ThreeVector> GetDamageCenterAndBoundaries(std::vector<G4ThreeVector> positions);
	void OutputSDDHeader(G4String author = "", G4String simDetails = "", G4String sourceDetails = "", G4int sourceType = 1, G4String energyDist = "M, 0",
			G4String irrTarget = "", G4String cellCycle = "0", G4String DNAStructure = "0, 1", G4int inVitroOrInVivo = 0, G4String prolStatus = "1", G4String microenv = "20, 0.01", G4double time = 0,
			G4String addInfo = "");

	G4int OutputSDDFile(std::map<G4int, std::vector<G4int>> damageSites, G4int eventID, G4int exposureID);
	void ExcludeShortDNAFragments(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs);
	void UpdateDamageAndPrimaryCount(G4int nLesions, G4int nPrimary);

	// -----------
	// Old methods
	// -----------
	void SeparateHitsOnDifferentDNAStrands(std::vector<TsHitsRecord*> Hits, std::vector<TsHitsRecord*> &HitsBack1, G4double IsStrand1);
	void DefineDSBorSSB(std::vector<TsHitsRecord*> HitsBack1, std::vector<TsHitsRecord*> HitsBack2,
			std::vector<std::pair<TsHitsRecord*, TsHitsRecord*>> &DSB_pairs,
			std::vector<TsHitsRecord*> &SSBonStrand1, std::vector<TsHitsRecord*> &SSBonStrand2);
	void GetSSBonStrand(std::vector<TsHitsRecord*> Hits, std::vector<TsHitsRecord*> & SSBonStrand);
	//void OutputSDDHeader(G4String filename);
	G4bool CauseDirectDamage(G4double edep);
	void DefineDSBType(TsHitsRecord* &hit1, TsHitsRecord* &hit2);
	void ExcludeShortDNAFragments(std::vector <std::pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, std::vector<TsHitsRecord*> SSBonStrand1,std::vector<TsHitsRecord*> SSBonStrand2, G4int ChromosomeID);
	void SeparateDamageSitesOnOneStrand(std::vector <TsHitsRecord*> Hits, std::vector <G4ThreeVector> &XYZPosition, std::vector <double> &ChromosomePosition,
            std::vector <double  > &Cause, std::vector <int > &DamageSiteDSBflag, std::vector< std::vector <std::vector <int>>> &DamageSites);
	void OutputDNAdamageTuple(std::vector<TsHitsRecord*> HitsBack1, std::vector<TsHitsRecord*> HitsBack2, G4String filename);
	void OutputDNAdamageTupleHeader(G4String filename);
	void OutputDNAdamageSummary(G4String filename);
	void OutputSDDFile(bool OnlyIncludeDSBinSDD, std::vector <std::pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs,
            std::vector<TsHitsRecord*> SSBonStrand1,std::vector<TsHitsRecord*> SSBonStrand2,
            G4String filename, G4int EventID, G4int ExposureID, G4int ChromosomeID);
	void WriteDNADamageCSV();

	// Set methods
	inline void SetDamageThreshold(G4double dt)					{ fDamageThreshold = dt; }
	inline void SetUseLinearProbabilityThreshold(G4bool ulpt)	{ fUseLinearProbabilityThreshold = ulpt; }
	inline void SetLinearProbabilityLowerLimit(G4double ll)		{ fLinearProbabilityLowerLimit = ll; }
	inline void SetLinearProbabilityUpperLimit(G4double ul)		{ fLinearProbabilityUpperLimit = ul; }
	inline void SetDSBSeparation(G4int s)						{ fDSBSeparation = s;}
	inline void SetDefineComplexity(G4bool dc)					{ fDefineComplexity = dc; }
	inline void SetComplexitySeparation(G4int s)				{ fComplexitySeparation = s; }
	inline void SetExcludeShortFragment(G4bool esf)				{ fExcludeShortFragment = esf; }
	inline void SetLowerFragmentDetectionThreshold(G4int lfdt)	{ fLowerFragmentDetectionThreshold = lfdt; }
	inline void SetUpperFragmentDetectionThreshold(G4int ufdt)	{ fUpperFragmentDetectionThreshold = ufdt; }
	inline void SetOnlyIncludeDSBinSDD(G4bool oidid)			{ fOnlyIncludeDSBinSDD = oidid; }
	inline void SetMinimalSDD(G4bool msdd)						{ fMinimalSDD = msdd; }
	inline void SetEventID(G4int eID)							{ fEventID = eID; }
	inline void SetEdep(G4double e)								{ fEdep = e; }
	inline void SetLET(G4double LET)							{ fLET = LET; }
	inline void SetNucleusMass(G4double nm)						{ fNucleusMass = nm; }
	inline void SetPrimaryParticle(G4String pp)					{ fPrimaryParticle = pp; }
	inline void SetMeanEnergy(G4double e)						{ fMeanEnergy = e; }
	inline void SetDosePerExposure(G4double d)					{ fDosePerExposure = d; }
	inline void SetScoreDirectDamages(G4bool sdd)				{ fScoreDirectDamages = sdd; }
	inline void SetScoreIndirectDamages(G4bool sid)				{ fScoreIndirectDamages = sid; }
	inline void SetTotalDNAContent(G4double dnac)				{ fTotalDNAContent = dnac; }
	inline void SetZetaBeta_sq(G4double zbsq)					{ fZetaBeta_sq = zbsq; }
	inline void	SetChromosomeDNAContent(std::vector<G4int> DNAc){ fChromosomeDNAContent = DNAc; }
	inline void SetScoreOnBases(G4bool sob)						{ fScoreOnBases = sob; }
	inline void SetScoreOnBackbones(G4bool sob)					{ fScoreOnBackbones = sob; }
	inline void SetScoreOnHydrationShell(G4bool sohs)			{ fScoreOnHydrationShell = sohs; }
	inline void SetWriteCSV(G4bool wcsv)						{ fWriteCSV = wcsv; }
	inline void SetOutputFileName(G4String fn)					{ fOutputFile = fn; }
	inline void SetOutputFileMode(G4String fmode)				{ fOutputMode = fmode; }

	// Get inline methods
	inline G4int GetNumSSB()			{ return numSSB; }
	inline G4int GetNumDirSSB()			{ return numSSB_dir; }
	inline G4int GetNumIndirSSB()		{ return numSSB_indir; }
	inline G4int GetNumDSB()			{ return numDSB; }
	inline G4int GetNumDirDSB()			{ return numDSB_dir; }
	inline G4int GetNumIndirDSB()		{ return numDSB_indir; }
	inline G4int GetNumHybridDSB()		{ return numDSB_hybrid; }
	inline G4int GetNumBaseDam()  		{ return numBaseDam; }
	inline G4int GetNumDirBaseDam()		{ return numBaseDam_dir; }
	inline G4int GetNumIndirBaseDam()	{ return numBaseDam_indir; }
	inline G4int GetNumSSBPlus()		{ return numSSBPlus; }
	inline G4int GetNumDSBPlus()		{ return numDSBPlus; }
	inline G4int GetNumComplex()		{ return numMoreComplex; }
	inline G4int GetExcludedNumSSB()			{ return Excluded_numSSB; }
	inline G4int GetExcludedNumDirSSB()			{ return Excluded_numSSB_dir; }
	inline G4int GetExcludedNumIndirSSB()		{ return Excluded_numSSB_indir; }
	inline G4int GetExcludedNumDSB()			{ return Excluded_numDSB; }
	inline G4int GetExcludedNumDirDSB()			{ return Excluded_numDSB_dir; }
	inline G4int GetExcludedNumIndirDSB()		{ return Excluded_numDSB_indir; }
	inline G4int GetExcludedNumHybridDSB()		{ return Excluded_numDSB_hybrid; }
	inline G4int GetExcludedNumBaseDam()  		{ return Excluded_baseDam; }
	inline G4int GetExcludedNumDirBaseDam()		{ return Excluded_baseDam_dir; }
	inline G4int GetExcludedNumIndirBaseDam()	{ return Excluded_baseDam_indir; }

private:
	G4double	fDamageThreshold;
	G4bool	 	fUseLinearProbabilityThreshold;
	G4double	fLinearProbabilityLowerLimit;
	G4double	fLinearProbabilityUpperLimit;
	G4int		fDSBSeparation;
	G4bool		fDefineComplexity;
	G4int		fComplexitySeparation;
	G4bool		fExcludeShortFragment;
	G4int		fLowerFragmentDetectionThreshold;
	G4int		fUpperFragmentDetectionThreshold;
	G4bool		fOnlyIncludeDSBinSDD;
	G4String	fPrimaryParticle;
	G4double	fMeanEnergy;
	G4double	fDosePerExposure;
	G4bool		fScoreDirectDamages;
	G4bool		fScoreIndirectDamages;
	G4bool		fMinimalSDD;
	G4bool		fWriteCSV;

	G4String 	fOutputFile;
	G4String	fOutputMode;

	G4int fEventID;
	G4double fEdep; // MeV
	G4double fLET; // keV/um
	G4double fNucleusMass; // g
	G4double fTotalDNAContent; // BP
	G4double fZetaBeta_sq;

    G4int numSB;
    G4int numSB_dir;
    G4int numSB_indir;
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
    // complex damages
    G4int numSSB_P;
    G4int numDSB_P;
    G4int numDSB_PP;
    G4int numSSBPlus;
    G4int numDSBPlus;
    G4int numMoreComplex;

    G4int Excluded_numSB;
    G4int Excluded_numSB_dir;
    G4int Excluded_numSB_indir;
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
    // complex damages
    G4int Excluded_numSSB_P;
    G4int Excluded_numDSB_P;
    G4int Excluded_numDSB_PP;

    std::vector<G4int> fChromosomeDNAContent; // In unit of BP

    G4int lastOutputEventID;
    G4int lastOutputExposureID;
    G4int lastDSBonlyOutputEventID;
    G4int lastDSBonlyOutputExposureID;

    G4bool fScoreOnBases, fScoreOnBackbones, fScoreOnHydrationShell;

    std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fVEdepInDNA;
    std::map<G4int, std::map<G4int, std::map<G4int, G4bool>>> fVIndirectDamageInDNA;
    std::map<G4int, std::map<G4int, std::map<G4int, G4bool>>> fVIonizationInHydrationShell;

    std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fDamage;
    std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fDSB;
    std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fSSB;
    std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fBaseDamage;

    std::map<G4int, std::map<G4int, std::map<G4int, G4ThreeVector>>> fDamagePositions;
};



#endif
