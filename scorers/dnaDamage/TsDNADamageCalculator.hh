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

#ifndef TsDNADamageCalculator_hh
#define TsDNADamageCalculator_hh

#include "TsVNtupleScorer.hh"
#include "TsHitInDNA.hh"

class TsDNADamageCalculator
{
public:
	TsDNADamageCalculator();
	~TsDNADamageCalculator();

	void ComputeStrandBreaks(std::vector<TsHitInDNA*> hits);
	void QuantifyDamage(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs);
	std::map<G4int, std::vector<G4int>> GetDamageSites();
	void ExcludeShortDNAFragments(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs);
	G4bool CauseDirectDamage(G4double edep);
	void OutputSDDHeader(G4bool minimalSDD, G4String primaryParticle, G4double energy, G4double dosePerExposure, std::vector<G4int> chromosomeContents, G4bool scoreIndirectDamages,
			G4bool scoreOnBases, G4String author, G4String simDetails, G4String sourceDetails, G4int sourceType, G4String energyDist, G4String irrTarget, G4String cellCycle, G4String DNAStructure,
			G4int inVitroOrInVivo, G4String prolStatus, G4String microenv, G4double time, G4String addInfo);
	G4int OutputSDDFile(std::map<G4int, std::vector<G4int>> damageSites, G4int eventID, G4int exposureID, std::vector<G4int> chromosomeContents);
	std::vector<G4ThreeVector> GetDamageCenterAndBoundaries(std::vector<G4ThreeVector> positions);
	void WriteDNADamageCSV();

	// Set methods
	inline void SetDistanceBasePairsForDSB(G4int v)				{ fNumberOfBasePairForDSB = v; }
	inline void SetDirectDamageThreshold(G4double v)			{ fDirectDamageThreshold = v; }
	inline void SetDirectDamageAsLinearProbability(G4bool v)	{ fUseLinearProbabilityForDirectDamage = v; }
	inline void SetLowerLimitForLinearProbability(G4double v)	{ fLowerLimitLinearProbability = v; }
	inline void SetUpperLimitForLinearProbability(G4double v)	{ fUpperLimitLinearProbability = v; }
	inline void SetExcludeShortFragments(G4bool v)				{ fExcludeShortFragments = v; }
	inline void SetLowerThresholdForFragmentDetection(G4int v)	{ fLowerThresholdForFragmentDetection = v; }
	inline void SetUpperThresholdForFragmentDetection(G4int v)	{ fUpperThresholdForFragmentDetection = v; }
	inline void SetOutputFileName(G4String v)					{ fOutputFileName = v; }
	inline void SetOutputFileMode(G4String v)					{ fOutputMode = v; }
	inline void SetWriteCSVFile(G4bool v)						{ fWriteCSVFile = v; }
	inline void SetMinimalModeForSDD(G4bool v)					{ fMinimalModeForSDD = v; }
	inline void SetReturnOnlyDSBInSDD(G4bool v)					{ fReturnOnlyDSBinSDD = v; }
	inline void SetEventID(G4int v)								{ fEventID = v; }

	// Get methods
	inline G4int GetSB()										{ return fNumSB; }
	inline G4int GetSBDirect()									{ return fNumSBDirect; }
	inline G4int GetSBQuasiDirect()								{ return fNumSBQuasiDirect; }
	inline G4int GetSBIndirect()								{ return fNumSBIndirect; }
	inline G4int GetSSB()										{ return fNumSSB; }
	inline G4int GetSSBDirect()									{ return fNumSSBDirect; }
	inline G4int GetSSBQuasiDirect()							{ return fNumSSBQuasiDirect; }
	inline G4int GetSSBIndirect()								{ return fNumSSBIndirect; }
	inline G4int GetDSB()										{ return fNumDSB; }
	inline G4int GetDSBDirect()									{ return fNumDSBDirect; }
	inline G4int GetDSBIndirect()								{ return fNumDSBIndirect; }
	inline G4int GetDSBDirectIndirect()							{ return fNumDSBDirectIndirect; }
	inline G4int GetDSBDirectQuasiDirect()						{ return fNumDSBDirectQuasiDirect; }
	inline G4int GetDSBQuasiDirectQuasiDirect()					{ return fNumDSBQuasiDirectQuasiDirect; }
	inline G4int GetDSBIndirectQuasiDirect()					{ return fNumDSBIndirectQuasiDirect; }
	inline G4int GetBD()										{ return fNumBaseDamage; }
	inline G4int GetBDDirect()									{ return fNumBaseDamageDirect; }
	inline G4int GetBDQuasiDirect()								{ return fNumBaseDamageQuasiDirect; }
	inline G4int GetBDIndirect()								{ return fNumBaseDamageIndirect; }
	inline G4int GetSSBPlus()									{ return fNumSSBPlus; }
	inline G4int GetDSBPlus()									{ return fNumDSBPlus; }
	inline G4int GetDSBComplex()								{ return fNumDSBComplex; }
	inline std::vector<G4ThreeVector> GetDSB3DPositions()		{ return fDSB3DPositions; }

private:
	// Options to define damage
	G4int fNumberOfBasePairForDSB;
	G4double fDirectDamageThreshold;
	G4bool fUseLinearProbabilityForDirectDamage;
	G4double fLowerLimitLinearProbability;
	G4double fUpperLimitLinearProbability;

	// Options for whether counting types of damage
	G4bool fExcludeShortFragments;
	G4int fLowerThresholdForFragmentDetection;
	G4int fUpperThresholdForFragmentDetection;

	// Options for the output
	G4String fOutputFileName;
	G4String fOutputMode;
	G4bool fWriteCSVFile;
	G4bool fMinimalModeForSDD;
	G4bool fReturnOnlyDSBinSDD;

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

	// Maps of damage. First index is chromosome, second index is base pair, third index is element in the base pair (base, backbone, hydshell)
	std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fMapEdep;
	std::map<G4int, std::map<G4int, std::map<G4int, G4bool>>> fIsThereIndirectDamage;
	std::map<G4int, std::map<G4int, std::map<G4int, G4bool>>> fIsThereQuasiDirectDamage;
	std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fDamageMap;
	// Maps of SB or BD. First index is chromosome, second index is base pair, third index is strand (1 or 2)
	std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fDSBMap;
	std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fSSBMap;
	std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> fBDMap;

	std::map<G4int, std::map<G4int, std::map<G4int, G4ThreeVector>>> fDamagePositions;

	std::vector<G4ThreeVector> fDSB3DPositions;

	// Data from end of the event
	G4double fEventID;

	// Codes for types of damage
	G4int nodamage = -1;
	G4int direct = 1;
	G4int indirect = 2;
	G4int multiple = 3;
	G4int quasidirect = 4;
	G4int multiplewithquasidirect = 5;

	// Codes for components ID
//	G4int base = 0;
//	G4int backbone = 1;
};

#endif
