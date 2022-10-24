#ifndef TsIRTConfiguration_hh
#define TsIRTConfiguration_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <vector>
#include <map>

class TsParameterManager;
class TsIRTUtils;

class TsIRTConfiguration {
public:
	TsIRTConfiguration(G4String, TsParameterManager*);
	
	~TsIRTConfiguration();
	
	void AddMolecule(G4String name, G4int moleculeID, G4double diffusionCoefficient, G4double charge, G4double radius);
	
	void AddMolecule(G4String name, G4double diffusionCoefficient, G4double charge, G4double radius);
	
	void AddMolecule(G4String name);
	
	void AdjustDiffusionCoefficientyForTemperature(G4double);
	
	void AdjustDiffusionCoefficientyForTemperatureArrehniusFit(G4double);
	
	void AdjustReactionRateForTemperature(G4double);
	
	void AdjustReactionRateForPH(G4String);
	
	G4double Arrhenius(G4double A, G4double T, G4double E);
	
	G4double Noyes(G4double kobs, G4double kact, G4double kdiff);
	
	G4double DebyeFactor(G4double T, G4int molA, G4int molB, G4double r);
	
	G4bool MoleculeExists(G4String name);
	
	G4double GetMoleculeRadius(G4int);
	G4int GetMoleculeCharge(G4int);
	
	G4int GetReactionIndex(G4int pdgA, G4int pdgB);
	
	G4double GetOnsagerRadius(G4int molA, G4int molB);
	
	void ResolveReactionParameters(G4int molA, G4int molB, G4double kobs, G4int reactionType);
	
	void ResolveReactionRateCoefficients();
	
	void CalculateContactProbabilities();
	
	void ResolveRemainerReactionParameters();
	
	void InsertReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
						G4double kobs, G4int reactionType);
	
	void InsertReaction(G4int molA, G4int molB, std::vector<G4int> products,
						G4double kobs, G4int reactionType);
	
	void InsertBackgroundReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
								  G4double kobs, G4double concentration, G4bool sampleExponential);
	
	void InsertBackgroundReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
								  G4double scavengingCapacity, G4bool sampleExponential);
	
	void QuitIfMoleculeNotFound(G4String mol);
	
	void Quit(const G4String& name, G4String message);
	
	void SetTimeLimits(G4double, G4double);
	
	G4int GetNumberOfReactions();
	
	std::pair<G4String, G4String> GetReactants(G4int);
	
	std::vector<G4String> GetProducts(G4int);
	
	// PH adjust functions
	
	std::vector<G4double> GetH2SO4ComponentsConcentrationP(G4double);
	
	std::vector<G4double> GetH2SO4ComponentsConcentrationPH(G4double);
	
	G4double IonicRate(G4double, G4double, G4int, G4int);
	
	G4double GetIonicStrength(std::vector<G4double>);
	
	void roots(double*,int, double*, double*);
	
	void deflate(double*,int, double*,double*, double *);
	
	void find_quad(double*, int, double*, double*, double*, int*);
	
	void diff_poly(double *, int, double*);
	
	void recurse(double *, int, double *, int, double*, double *, int *);
	
	void get_quads(double*, int, double*, double*);
	
	std::vector<double> GetRoots(int, std::vector<double>);
	
	// Temperature adjust functions
	
	G4double K_H(G4double);
	
	G4double K_OH(G4double);
	
	G4double K_water(G4double);
	
	G4double K_H2O2(G4double);
	
	G4double lH2Ol(G4double);
	
	G4double k23_inverse(G4double);
	
	G4double k26_inverse(G4double);
	
	G4double k27_inverse(G4double);
	
	G4double k29_inverse(G4double);
	
	G4double k30_inverse(G4double);
	
	G4double k_26(G4double);
	
	G4double k_28(G4double);
	
	G4double k_29(G4double);
	
	G4double k_30(G4double);
	
private:
	TsParameterManager* fPm;
	TsIRTUtils* fUtils;
	G4String fName;
	
	struct TsMolecularReaction {
		G4int    index;
		
		G4int reactorA;
		G4int reactorB;
		std::vector<G4int> products;
		
		G4double kobs;
		G4double kdif;
		G4double kact;
		G4double reactionRadius;
		G4double effectiveReactionRadius;
		G4double effectiveTildeReactionRadius;
		G4double probabilityOfReaction;
		G4double alpha;
		G4int    reactionType;
		G4bool   sampleExponential;
		G4double OnsagerRadius;
		G4double concentration;
		G4double scavengingCapacity;
	};
	
	struct TsMoleculeDefinition {
		G4double diffusionCoefficient;
		G4double charge;
		G4double radius;
	};
	
public:
	struct TsMolecule {
		G4int    id;
		G4double time;
		G4ThreeVector position;
		G4bool reacted;
		G4int trackID;
		G4int spin;
		G4bool isDNA;
		std::vector<G4int> tested;
	};
	
	inline std::map<G4String, G4int> GetMoleculeIDs() {return fMoleculesID;};
	
	inline std::map<G4int, G4String> GetMoleculeNames() { return fMoleculesName;};
	
private:
	std::map<G4int, TsMolecule> fMolecules;
	
	std::map<G4int, TsMoleculeDefinition> fMoleculesDefinition;
	std::map<G4String, G4int> fMoleculesID;
	std::map<G4int, G4String> fMoleculesName;
	
	std::map<G4int, TsMolecularReaction > fReactions;
	std::map<G4int, std::vector<std::pair<G4int,G4int>>> fMoleculeCanReactWith;
	
	std::map<G4String, G4String> fExistingMolecules;
	
	G4double fUpperTime;
	G4double fLowerTime;
	G4int fReactionID;
	G4int fTotalBinaryReaction;
	G4int fLastMoleculeID;
	G4bool fUseSimpleScavengerModel;
	
	G4double fTemperature;
	G4bool fScaleForTemperature;
	G4bool fKick;
	//G4double fH2SO4Concentration;
	
	G4String fpHSolvent;
	G4double fpHSolventConcentration;
	G4double fpHValue;
	
	G4double fFeSO4;
	G4bool   fSampleExponential;
	
	G4double fObservedReactionRate;
	G4double fDiffusionReactionRate;
	G4double fActivationReactionRate;
	G4double fProbability;
	G4double fAlpha;
	G4double fReactionRadius;
	G4double fEffectiveReactionRadius;
	G4double fEffectiveTildeReactionRadius;
	G4double fOnsagerRadius;
	
	G4bool fQualityAssurance;
	
	G4bool fAllTotallyDiffusionControlled;
	
public:
	
	G4double SampleExponentialTime(G4int pdgA, G4int pdgB, G4int indexOfReaction);
	G4double GetIndependentReactionTime(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	G4double SampleIRTTotallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	G4double SampleIRTPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	std::pair<G4int, G4double> SampleIRTFirstOrderAndBackgroundReactions(TsMolecule molA );
	std::vector<std::pair<G4int, G4double>> SampleAllIRTFirstOrderAndBackgroundReactions(TsMolecule molA );
	
	G4int ContactFirstOrderAndBackgroundReactions(TsMolecule molA );
	
	G4double CalculateProbabilityPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double t);
	G4double CalculateProbabilityOfScavenger(TsMolecule molA, G4int indexOfReaction, G4double t);
	
	G4double brents_fun(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double offset);
	G4double brents_fun_scav(TsMolecule molA, G4int indexOfReaction, G4double offset);
	G4double SolveTime(TsMolecule molA, G4int indexOfReaction, G4double offset);
	
	void ResampleReactantsPosition(TsMolecule& molA, TsMolecule& molB, G4int index, G4double time);
	std::vector<G4ThreeVector> GetPositionOfProducts(TsMolecule molA, TsMolecule molB, G4int index);
	G4double GetRCutOff(G4double tCutOff);
	
	G4bool MakeReaction(std::vector<TsMolecule> &initialSpecies,
						std::map<G4int, std::map<G4int, std::map<G4int, std::vector<G4int>>>> &spaceBinned,
						G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin, G4double YMax, G4double ZMin, G4double ZMax,
						std::map<G4int, std::map<G4int, G4int>> &theGvalue,
						std::vector<G4double> timeSteps,
						G4int iM, G4int indexOfReaction, G4double irt, std::vector<G4bool> &used);
	
	G4bool MakeReaction(std::vector<TsMolecule> &initialSpecies,
						std::map<G4int, std::map<G4int, std::map<G4int, std::vector<G4int>>>> &spaceBinned,
						G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin, G4double YMax, G4double ZMin, G4double ZMax,
						std::map<G4int, std::map<G4int, G4int>> &theGvalue,
						std::vector<G4double> timeSteps,
						G4int iM, G4int jM, G4int indexOfReaction, G4double irt, G4double probabilityOfReaction, std::vector<G4bool> &used);
	
	G4bool MakeReaction(std::vector<TsMolecule> &initialSpecies,
						std::map<G4int, std::map<G4int, std::map<G4int, std::vector<G4int>>>> &spaceBinned,
						G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin, G4double YMax, G4double ZMin, G4double ZMax,
						std::map<G4int, std::map<G4int, G4int>> &theGvalue,
						std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
						std::vector<G4double> timeSteps,
						G4int iM, G4int jM, G4int indexOfReaction, G4double irt, G4double probabilityOfReaction, std::vector<G4bool> &used);
	
	G4bool Inside(G4ThreeVector p);
	
	void ScoreGvalue(std::vector<TsMolecule> &initialSpecies,
					 std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
					 std::vector<G4double> timeSteps,
					 G4int iM, G4int jM, G4int indexOfReaction, G4double irt);
	
	inline TsMolecularReaction GetReaction(G4int index) { return fReactions[index];};
	
	void Diffuse(TsMolecule& mol, G4double dt);
	
	
	void TestSampling(G4int indexOfReaction, G4int nHistories);
	
	void PrintMoleculesInformation();
	void PrintReactionsInformation();
	
	inline G4int GetLastMoleculeID() { return fLastMoleculeID; };
	inline G4int GetLastReactionID() { return fReactionID-1; };
	
	G4double Smoluchowski(G4double Beta, TsMolecularReaction Reaction);
	
	G4double Debye(G4double Beta, TsMolecularReaction Reaction, G4double T);
	
	G4double ProbabilityOfReactionT(TsMolecularReaction Reaction, G4double temperatureInKelvin, G4double Xi);
	
	G4double OnsagerRadius(G4double temperatureInKelvin);
	
	G4double IonicRate(G4double, TsMolecularReaction);
	
};
#endif

