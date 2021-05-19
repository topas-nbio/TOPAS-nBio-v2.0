// Scorer for NucleusDNADamage
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

#include "TsScorerNucleusDNADamage.hh"
#include "TsChromosome.hh"
//#include "TsGeometryManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Molecule.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"

#include "G4ParticleTable.hh"

#include <algorithm>
#include <stdint.h>
using namespace std;

TsScorerNucleusDNADamage::TsScorerNucleusDNADamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
									 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	// Initialize
	fTrackAveragedLET 	= 0;
	fTravelDistance		= 0;
	fEdep				= 0;
	fDoseInThisExposure	= 0;
	fExposureID			= 0;
	fCalZetaBeta		= false;
	fZetaBeta_sq		= 0;

	// IS THIS NECESSARY?
	//G4String componentName = fPm->GetStringParameter(GetFullParmName("Component"));
	//TsVGeometryComponent* Component = gM->GetComponent(componentName);
	//
		// Parameters
	fOutFileName = fPm->GetStringParameter(GetFullParmName("OutputFile"));

	fMinimalSDDOutput = false;
	if ( fPm->ParameterExists(GetFullParmName("MinimalSDDOutput")) )
		fMinimalSDDOutput = fPm->GetBooleanParameter(GetFullParmName("MinimalSDDOutput"));

	fNumberOfHistoriesInRun =  fPm->GetIntegerParameter(GetFullParmName("NumberOfHistoriesInRun"));

	fProbabilityOfOHDamage = 0.4;
	if ( fPm->ParameterExists(GetFullParmName("ProbabilityForOHToCauseDamage")) )
		fProbabilityOfOHDamage = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityForOHToCauseDamage"));

	fDamageThreshold = 17.5*eV;
	if ( fPm->ParameterExists(GetFullParmName("DamageThreshold")) )
		fDamageThreshold = fPm->GetDoubleParameter(GetFullParmName("DamageThreshold"),"Energy");

	fUseLinearProbabilityThreshold = false;
	if ( fPm->ParameterExists(GetFullParmName("UseLinearProbabilityThreshold")) )
		fUseLinearProbabilityThreshold = fPm->GetBooleanParameter(GetFullParmName("UseLinearProbabilityThreshold"));

	fLinearProbability_lower_limit= 5*eV;
		if ( fPm->ParameterExists(GetFullParmName("LinearProbability_lower_limit")) )
		fLinearProbability_lower_limit = fPm->GetDoubleParameter(GetFullParmName("LinearProbability_lower_limit"),"Energy");

	fLinearProbability_upper_limit= 37.5*eV;
		if ( fPm->ParameterExists(GetFullParmName("LinearProbability_upper_limit")) )
		fLinearProbability_upper_limit = fPm->GetDoubleParameter(GetFullParmName("LinearProbability_upper_limit"),"Energy");

	fDSBSeparation = 10; //10 BP
	if ( fPm->ParameterExists(GetFullParmName("DSBSeparation")) )
		fDSBSeparation = fPm->GetIntegerParameter(GetFullParmName("DSBSeparation"));

	fExcludeShortFragment = true;
	if ( fPm->ParameterExists(GetFullParmName("ExcludeShortFragment")) )
		fExcludeShortFragment = fPm->GetBooleanParameter(GetFullParmName("ExcludeShortFragment"));

	fLowerFragmentDetectionThreshold = 0; //In unit of basepair
	if ( fPm->ParameterExists(GetFullParmName("LowerFragmentDetectionThreshold")) )
		 fLowerFragmentDetectionThreshold = fPm->GetIntegerParameter(GetFullParmName("LowerFragmentDetectionThreshold"));

	fUpperFragmentDetectionThreshold = 1*3E8; //In unit of basepair, maximum chromosome size
	if ( fPm->ParameterExists(GetFullParmName("UpperFragmentDetectionThreshold")) )
		 fUpperFragmentDetectionThreshold = fPm->GetIntegerParameter(GetFullParmName("UpperFragmentDetectionThreshold"));

	fOnlyIncludeDSBinSDD = true;
	if ( fPm->ParameterExists(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD")) )
		fOnlyIncludeDSBinSDD = fPm->GetBooleanParameter(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD"));

	fWriteCSV = false;
	if ( fPm->ParameterExists(GetFullParmName("WriteCSVOutputWithAllDamageSpecification")))
		fWriteCSV = fPm->GetBooleanParameter(GetFullParmName("WriteCSVOutputWithAllDamageSpecification"));

	fScoreDirectDamages = true;
	if ( fPm->ParameterExists(GetFullParmName("ScoreDirectDamages")) )
		fScoreDirectDamages = fPm->GetBooleanParameter(GetFullParmName("ScoreDirectDamages"));

	fScoreIndirectDamages = true;
	if ( fPm->ParameterExists(GetFullParmName("ScoreIndirectDamages")) )
		fScoreIndirectDamages = fPm->GetBooleanParameter(GetFullParmName("ScoreIndirectDamages"));

	/*fScoreOnBases = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnBases")))
		fScoreOnBases = fPm->GetBooleanParameter(GetFullParmName("ScoreOnBases"));

	fScoreOnBackbones = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnBackbones")))
		fScoreOnBackbones = fPm->GetBooleanParameter(GetFullParmName("ScoreOnBackbones"));

	fScoreOnHydrationShell = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreOnHydrationShell")))
		fScoreOnHydrationShell = fPm->GetBooleanParameter(GetFullParmName("ScoreOnHydrationShell"));*/

	// This is to specify the dose of each exposure for SDD output
	// Please see Schuemann, J., et al. (2019). "A New Standard DNA Damage (SDD) Data Format." Radiat Res 191(1): 76-92.
	// for definition of exposure
	fDosePerExposure = 1 * gray ;//In unit of Gy
	if ( fPm->ParameterExists(GetFullParmName("DosePerExposure") ))
		fDosePerExposure = fPm->GetDoubleParameter(GetFullParmName("DosePerExposure"), "Dose");

	fScoringRadius = 4.65*um ;
	if ( fPm->ParameterExists(GetFullParmName("ScoringRadius")) )
		fScoringRadius = fPm->GetDoubleParameter(GetFullParmName("ScoringRadius"),"Length");

	fHistoneAsScavenger = true;
	if ( fPm->ParameterExists(GetFullParmName("HistoneAsScavenger")) )
		fHistoneAsScavenger = fPm->GetBooleanParameter(GetFullParmName("HistoneAsScavenger"));

	// Parameter to specify the primary particle. If not specified, it's assumed that it is proton
	fPrimaryParticle = "proton";
	if ( fPm->ParameterExists(GetFullParmName("PrimaryParticle")) )
		fPrimaryParticle = fPm->GetStringParameter(GetFullParmName("PrimaryParticle"));

	// Parameter to specify the level of base pairs in the geometry hierarchy
	fBasePairDepth = 0;
	if ( fPm->ParameterExists(GetFullParmName("BasePairPositionAtGeometricHierarchy")) )
		fBasePairDepth = fPm->GetIntegerParameter(GetFullParmName("BasePairPositionAtGeometricHierarchy"));

	// Information about geometry, chromosome and voxel IDs
	fGeometryInfo 	= "";
	if ( fPm->ParameterExists(GetFullParmName("GeometryInfo")))
		fGeometryInfo = fPm->GetStringParameter(GetFullParmName("GeometryInfo"));
	fCopyNoTable 	= "";
	if ( fPm->ParameterExists(GetFullParmName("CopyNoTable")))
		fCopyNoTable = fPm->GetStringParameter(GetFullParmName("CopyNoTable"));
	fSignedCHVoxel	= "";
	if ( fPm->ParameterExists(GetFullParmName("signedCHVoxel")))
		fSignedCHVoxel = fPm->GetStringParameter(GetFullParmName("signedCHVoxel"));

	if (fGeometryInfo.size() == 0)
	{
		fVoxel3Drepeat = 1; // To be updated
		fNucleusMass = 997e3 * 4*CLHEP::pi/3 * pow(fScoringRadius/m, 3);
		fChromosomeDNAContent.push_back(15149);
	}
	G4String* strand1materialNames;
	G4String* strand2materialNames;
	G4int strand1Length, strand2Length;
	if ( fPm->ParameterExists(GetFullParmName("Strand1MaterialNames")) )
	{
		strand1materialNames = fPm->GetStringVector(GetFullParmName("Strand1MaterialNames"));
		strand1Length = fPm->GetVectorLength(GetFullParmName("Strand1MaterialNames"));
	}
	else
		strand1materialNames[0] = "G4_WATER";
	for (G4int i = 0; i < strand1Length; i++)
		fStrand1Materials.push_back(GetMaterial(strand1materialNames[i]));
	if ( fPm->ParameterExists(GetFullParmName("Strand2MaterialNames")) )
	{
		strand2materialNames = fPm->GetStringVector(GetFullParmName("Strand2MaterialNames"));
		strand2Length = fPm->GetVectorLength(GetFullParmName("Strand2MaterialNames"));
	}
	else
		strand2materialNames[0] = "G4_WATER";
	for (G4int i = 0; i < strand2Length; i++)
		fStrand2Materials.push_back(GetMaterial(strand1materialNames[i]));

	// Parameters for the SDD header
	G4String author = "@";
	if ( fPm->ParameterExists(GetFullParmName("AuthorName")) )
		author = fPm->GetStringParameter(GetFullParmName("AuthorName"));
	G4String simulationDetails = "Sim details";
	if ( fPm->ParameterExists(GetFullParmName("SimulationDetails")) )
		simulationDetails = fPm->GetStringParameter(GetFullParmName("SimulationDetails"));
	G4String sourceDetails = "Source details";
	if ( fPm->ParameterExists(GetFullParmName("SourceDetails")) )
		sourceDetails = fPm->GetStringParameter(GetFullParmName("SourceDetails"));
	G4int sourceType = 1;
	if ( fPm->ParameterExists(GetFullParmName("SourceType")) )
		sourceType = fPm->GetIntegerParameter(GetFullParmName("SourceType"));
	G4double meanEnergy = 0.0;
	if ( fPm->ParameterExists(GetFullParmName("MeanEnergy")))
		meanEnergy = fPm->GetDoubleParameter(GetFullParmName("MeanEnergy"), "Energy");
	G4String energyDist = "M, 0";
	if ( fPm->ParameterExists(GetFullParmName("EnergyDistribution")) )
		energyDist = fPm->GetStringParameter(GetFullParmName("EnergyDistribution"));
	G4String irrTarget = "";
	if ( fPm->ParameterExists(GetFullParmName("IrradiationTarget")) )
		irrTarget = fPm->GetStringParameter(GetFullParmName("IrradiationTarget"));
	G4String cellCycle = "0";
	if ( fPm->ParameterExists(GetFullParmName("CellCycleStage")))
		cellCycle = fPm->GetStringParameter(GetFullParmName("CellCycleStage"));
	G4String DNAStructure = "0, 1";
	if ( fPm->ParameterExists(GetFullParmName("DNAStructure")))
		cellCycle = fPm->GetStringParameter(GetFullParmName("DNAStructure"));
	G4int inVitroOrInVivo = 0;
	if ( fPm->ParameterExists(GetFullParmName("InVitroOrInVivo")) )
		inVitroOrInVivo = fPm->GetIntegerParameter(GetFullParmName("InVitroOrInVivo"));
	G4String proliferationStatus = "1";
	if ( fPm->ParameterExists(GetFullParmName("ProliferationStatus")))
		proliferationStatus = fPm->GetStringParameter(GetFullParmName("ProliferationStatus"));
	G4String microenvironment = "20, 0.01";
	if ( fPm->ParameterExists(GetFullParmName("Microenvironment")))
		microenvironment = fPm->GetStringParameter(GetFullParmName("Microenvironment"));
	G4double time = 0;
	if ( fPm->ParameterExists(GetFullParmName("Time")))
		time = fPm->GetDoubleParameter(GetFullParmName("Time"), "Time");
	G4String addInfo = "";
	if ( fPm->ParameterExists(GetFullParmName("AdditionalInfo")))
		addInfo = fPm->GetStringParameter(GetFullParmName("AdditionalInfo"));

	GetGeometryInfo();
	fDefineDamage = new TsDefineDamage();
	fDefineDamage->SetDamageThreshold(fDamageThreshold);
	fDefineDamage->SetUseLinearProbabilityThreshold(fUseLinearProbabilityThreshold);
	fDefineDamage->SetLinearProbabilityLowerLimit(fLinearProbability_lower_limit);
	fDefineDamage->SetLinearProbabilityUpperLimit(fLinearProbability_upper_limit);
	fDefineDamage->SetDSBSeparation(fDSBSeparation);
	fDefineDamage->SetExcludeShortFragment(fExcludeShortFragment);
	fDefineDamage->SetLowerFragmentDetectionThreshold(fLowerFragmentDetectionThreshold);
	fDefineDamage->SetUpperFragmentDetectionThreshold(fUpperFragmentDetectionThreshold);
	fDefineDamage->SetNucleusMass(fNucleusMass);
	fDefineDamage->SetTotalDNAContent(fTotalDNAContent);
	fDefineDamage->SetChromosomeDNAContent(fChromosomeDNAContent);
	fDefineDamage->SetOnlyIncludeDSBinSDD(fOnlyIncludeDSBinSDD);
	fDefineDamage->SetScoreOnBases(fScoreOnBases);
	fDefineDamage->SetScoreOnBackbones(fScoreOnBackbones);
	fDefineDamage->SetScoreOnHydrationShell(fScoreOnHydrationShell);
	fDefineDamage->SetPrimaryParticle(fPrimaryParticle);
	fDefineDamage->SetMeanEnergy(meanEnergy);
	fDefineDamage->SetDosePerExposure(fDosePerExposure);
	fDefineDamage->SetScoreDirectDamages(fScoreDirectDamages);
	fDefineDamage->SetScoreIndirectDamages(fScoreIndirectDamages);
	fDefineDamage->SetMinimalSDD(fMinimalSDDOutput);
	fDefineDamage->SetWriteCSV(fWriteCSV);
	fDefineDamage->SetOutputFileName(fOutFileName);
	fDefineDamage->SetOutputFileMode(fOutFileMode);
	fDefineDamage->OutputSDDHeader(author, simulationDetails, sourceDetails, sourceType, energyDist, irrTarget, cellCycle, DNAStructure,
								   inVitroOrInVivo, proliferationStatus, microenvironment, time, addInfo);

	// ********************************************************************************
	// Print parameters
	// ********************************************************************************
	G4cout<<"*********************************************************************************"<<G4endl;
	G4cout << "fProbabilityOfOHDamage = "<<fProbabilityOfOHDamage<<G4endl;
	if (!fUseLinearProbabilityThreshold)
	G4cout << "fDamageThreshold = "<<fDamageThreshold/eV<<" eV"<<G4endl;
	else
	G4cout << "fDamageThreshold = "<<fLinearProbability_lower_limit/eV<<" - "<<fLinearProbability_upper_limit/eV<<" eV"<<G4endl;
	G4cout << "fDSBSeparation = "<<fDSBSeparation<<" bp"<<G4endl;
	G4cout << "fLowerFragmentDetectionThreshold = "<< fLowerFragmentDetectionThreshold<<G4endl;
	G4cout << "fUpperFragmentDetectionThreshold = "<< fUpperFragmentDetectionThreshold<<G4endl;
	G4cout<<"*********************************************************************************"<<G4endl;

	// Prepare nTuple
	fNtuple->RegisterColumnD(&fEdepkeV, "Energy_imparted_keV", "");
	fNtuple->RegisterColumnD(&fDoseGy, "Dose_per_event_Gy", "");
	fNtuple->RegisterColumnD(&fLETtkeVum, "Track_averaged_LET_keV/um", "");
	fNtuple->RegisterColumnI(&numSSB, "Number_of_SSBs");
	fNtuple->RegisterColumnD(&yieldSSB, "SSB/Gy/Gbp", "");
	fNtuple->RegisterColumnI(&numSSB_dir, "Number_of_direct_SSBs");
	fNtuple->RegisterColumnI(&numSSB_indir, "Number_of_indirect_SSBs");
	fNtuple->RegisterColumnI(&numDSB, "Number_of_DSBs");
	fNtuple->RegisterColumnD(&yieldDSB, "DSB/Gy/Gbp", "");
	fNtuple->RegisterColumnI(&numDSB_dir, "Number_of_direct_DSBs");
	fNtuple->RegisterColumnI(&numDSB_indir, "Number_of_indirect_DSBs");
	fNtuple->RegisterColumnI(&numDSB_hybrid, "Number_of_hybrid_DSBs");
	if (fScoreOnBases)
	{
		fNtuple->RegisterColumnI(&numBaseDam, "Number_of_base_damages");
		fNtuple->RegisterColumnD(&yieldBaseDam, "BDs/Gy/Gbp", "");
		fNtuple->RegisterColumnI(&numBaseDam_dir, "Number_of_direct_base_damages");
		fNtuple->RegisterColumnI(&numBaseDam_indir, "Number_of_indirect_base_damages");
	}
	fNtuple->RegisterColumnI(&numSSBPlus, "Number_of_SSB_plus_SSB");
	fNtuple->RegisterColumnD(&yieldSSBPlus, "SSB_plus/Gy/Gbp", "");
	fNtuple->RegisterColumnI(&numDSBPlus, "Number_of_DSB_plus_DSB");
	fNtuple->RegisterColumnD(&yieldDSBPlus, "DSB_plus/Gy/Gbp", "");
	fNtuple->RegisterColumnI(&numMoreComplex, "Number_of_more_complex_damages");
	fNtuple->RegisterColumnD(&yieldMoreComplex, "More_Complex_Breaks/Gy/Gbp", "");
	if (fExcludeShortFragment)
	{
		fNtuple->RegisterColumnI(&Excluded_numSSB, "Excluded_SSBs");
		fNtuple->RegisterColumnI(&Excluded_numSSB_dir, "Excluded_direct_SSBs");
		fNtuple->RegisterColumnI(&Excluded_numSSB_indir, "Excluded_indirect_SSBs");
		fNtuple->RegisterColumnI(&Excluded_numDSB, "Excluded_DSBs");
		fNtuple->RegisterColumnI(&Excluded_numDSB_dir, "Excluded_direct_DSBs");
		fNtuple->RegisterColumnI(&Excluded_numDSB_indir, "Excluded_indirect_DSBs");
		fNtuple->RegisterColumnI(&Excluded_numDSB_hybrid, "Excluded_hybrid_DSBs");
	    if (fScoreOnBases)
	    {
	    	fNtuple->RegisterColumnI(&Excluded_baseDam, "Excluded_base_damages");
	    	fNtuple->RegisterColumnI(&Excluded_baseDam_dir, "Excluded_direct_base_damages");
	    	fNtuple->RegisterColumnI(&Excluded_baseDam_indir, "Excluded_indirect_base_damages");
	    }
	}
}

TsScorerNucleusDNADamage::~TsScorerNucleusDNADamage() { }

G4bool TsScorerNucleusDNADamage::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive)
	{
		fSkippedWhileInactive++;
		return false;
	}

	// Check if within nucleus
	G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
	G4bool withinNucleus = ( (pow(pos.x(),2) + pow(pos.y(), 2) + pow(pos.z(), 2)) < pow(fScoringRadius, 2) );
	if (!withinNucleus)
		return false;
	// Accumulate energy
	fEdep += aStep->GetTotalEnergyDeposit();

	// Accumulating length in the nucleus for protons and alpha particles
	const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
	G4String particleName = particle->GetParticleName();
	if (fPrimaryParticle == "proton" && (particleName == "proton" || particleName == "hydrogen"))
		fTravelDistance += aStep->GetStepLength();
	if (fPrimaryParticle == "alpha" && (particleName == "alpha" || particleName == "alpha+" || particleName == "helium"))
		fTravelDistance += aStep->GetStepLength();

	// Calculate (Zeff/beta)^2
	if (!fCalZetaBeta)
	{
		G4double particleVelocity = aStep->GetTrack()->GetVelocity();
		G4double lightSpeed = 3e8 * m/s;
		G4double a_beta = particleVelocity/lightSpeed;
		G4double z_eff = 1-exp(-125*a_beta);
		fZetaBeta_sq = pow(z_eff/a_beta, 2);
		G4cout<<"ParticleVelocity= "<< particleVelocity*s/m<<" m/sec;  "
				<<"beta = "<<a_beta<<"  "<<"(Z_eff/β)^2 ="<<fZetaBeta_sq<<G4endl;
		fCalZetaBeta = true;
	}

	// Get volume at the depth of base pairs
	G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	G4String volumeName = touchable->GetVolume(fBasePairDepth)->GetName();

	// Check if hit happens in any of the strand materials
	G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
	G4bool materialMatched = false;
	for (G4int i = 0; i < fStrand1Materials.size(); i++)
	{
		if (material == fStrand1Materials[i])
		{
			materialMatched = true;
			break;
		}
	}
	if (!materialMatched)
	{
		for (G4int i = 0; i < fStrand2Materials.size(); i++)
		{
			if (material == fStrand2Materials[i])
			{
				materialMatched = true;
				break;
			}
		}
	}

	// Filtering by material
	if (materialMatched)
	{
		// IDs for voxel, fiber, base pair in fiber (except for histones)
		G4int bpIDInFiber = touchable->GetVolume(fBasePairDepth)->GetCopyNo();;
		G4int VoxelID = touchable->GetCopyNumber(fBasePairDepth+2);
		G4int FiberID = touchable->GetCopyNumber(fBasePairDepth+1);
		if (strstr(volumeName, "Histone") != NULL)
			bpIDInFiber = -1;

		G4int VoxelNumInSphere = 0;
		G4int ChromosomeID = 0;
		G4int bpIDInChromosome = 0;
		FindPlaceInChromosome(VoxelID, FiberID, bpIDInFiber, VoxelNumInSphere, ChromosomeID, bpIDInChromosome);

		// If VoxelNumInSphere is -1, it means the voxel is outside the spherical shape of nucleus
		if (VoxelNumInSphere == -1 || VoxelID >= pow(fVoxel3Drepeat, 3))
			return false;
		G4double energyDeposit = aStep->GetTotalEnergyDeposit();
		// Set hit record information
		TsHitsRecord* hit = new TsHitsRecord();
		hit->SetEventID(GetEventID() + fNumberOfHistoriesInRun * GetRunID());
		hit->SetPosX(pos.x());
		hit->SetPosY(pos.y());
		hit->SetPosZ(pos.z());
		hit->SetPosition( G4ThreeVector(pos.x(), pos.y(), pos.z()) );
		hit->SetBasePairID(bpIDInChromosome);
		hit->SetChromosomeID(ChromosomeID);
		hit->SetVoxelID(VoxelNumInSphere);
		hit->SetVolumeName(volumeName);
		hit->SetEdep(energyDeposit);
		hit->SetParticleName(particleName);
		hit->SetProcess((G4String)aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());

		// Adding direct damage
		if (aStep->GetTrack()->GetTrackID() >= 0 && energyDeposit > 0 && fScoreDirectDamages)
		{
			Hits.push_back(hit);
			return true;
		}
		// Including indirect damage
		if (aStep->GetTrack()->GetTrackID() < 0 && fScoreIndirectDamages)
		{
			G4String speciesName 	= GetMolecule(aStep->GetTrack())->GetName();
			G4bool isSpeciesToKill 	= (speciesName == "OH^0" || speciesName == "e_aq^-1" || speciesName == "H^0");
			G4bool isHydroxyl		= (speciesName == "OH^0");
			G4bool justEnterVolume	= (aStep->GetPreStepPoint()->GetStepStatus() != fGeomBoundary);
			G4String backName = "Backbone";
			G4String baseName = "Base";
			// Kill all species generated inside DNA volume
			if (aStep->IsFirstStepInVolume())
			{
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				delete hit;
				return false;
			}
			// Make the hydroxyl damage to DNA
			else if (justEnterVolume && isHydroxyl && (strstr(volumeName, backName) != NULL || strstr(volumeName, baseName) != NULL))
			{
				hit->SetEdep(-0.001 * eV);
				hit->SetIsDirectDamage(false);
				G4bool reacted = true;
				if (reacted && G4UniformRand() < fProbabilityOfOHDamage)
				{
					Hits.push_back(hit);
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					return true;
				}
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				delete hit;
				return false;
			}
			// Scavenge {OH*, e_aq- and H*} diffusing into DNA volume
			else if (!justEnterVolume && isSpeciesToKill)
			{
				G4String histoneName = "Histone";
				if (fHistoneAsScavenger && strstr(volumeName, histoneName) != NULL)
				{
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					delete hit;
					return false;
				}
			}
		}
		delete hit;
	}
	return false;
}

void TsScorerNucleusDNADamage::AccumulateEvent()
{
	HitsOfEvents.push_back(Hits);
	Hits.clear();
	G4double edep = fEdep;
	eventsEdep.push_back(edep);
	fEdep = 0;
	G4double length = fTravelDistance;
	eventsLength.push_back(length);
	fTravelDistance = 0;
}

void TsScorerNucleusDNADamage::GetGeometryInfo()
{
	if (fGeometryInfo.size() > 0 && fCopyNoTable.size() > 0 && fSignedCHVoxel.size() > 0)
	{
		// ------------------------------ //
		// Read geometry information file
		// ------------------------------ //
		ifstream infile;
		infile.open(fGeometryInfo);
		G4String input;
		if (!infile.is_open())
		{
			G4cout << "ERROR: Unable to open file " << fGeometryInfo << G4endl;
			exit(0);
		}
		else
			G4cout << "Reading " << fGeometryInfo << G4endl;
			while (infile >> input)
		{
			if (input == "DNAContentperFiber(bp):")			infile >> fFiberDNAContent;	// bp
			if (input == "DNAContentperVoxel(bp):")			infile >> fVoxelDNAContent; // bp
			if (input == "DNAContentinTotal(bp):")			infile >> fTotalDNAContent; // bp
			if (input == "NucleusMass(g):")					infile >> fNucleusMass;		// g
			if (input == "NumberofVoxels(subdomains):")		infile >> fNumofVoxel;
			if (input == "Voxel3DrepeatTimes:")				infile >> fVoxel3Drepeat;
			if (input == "VoxelSize(um):")					infile >> fVoxelSize;		// um
			fVoxelSize = fVoxelSize * um;
			fVoxelContainerminXYZ = -fVoxelSize * fVoxel3Drepeat / 2;
		}
		infile.close();
		// --------------- //
		// Get DNA content
		// --------------- //
		TsChromosome* HumanFibroblastNucleus = new TsChromosome();
		fChromosomeDNAContent = HumanFibroblastNucleus->GetDNAcontentofEachChromosome_BP();
		delete HumanFibroblastNucleus;

		// --------------------- //
		// Read copy number file
		// --------------------- //
		infile.open(fCopyNoTable);

		if (!infile.is_open())
		{
			G4cout << "ERROR: Unable to open file " << fCopyNoTable << G4endl;
			exit(0);
		}
		else
			G4cout << "Reading " << fCopyNoTable << G4endl;

		G4int cpbox = 0, cpsphere = 0;
		G4String line;
		while (getline(infile, line))
		{
			stringstream stream(line.data());
			stream >> cpbox >> cpsphere;
			fVoxelNumInBox.push_back(cpbox);
			fVoxelNumInSphere.push_back(cpsphere);
		}
		infile.close();

		// ------------------------------ //
		// Read chromosomeID and voxel ID
		// ------------------------------ //
		infile.open(fSignedCHVoxel);

		if (!infile.is_open())
		{
			G4cout << "ERROR: Unable to open file " << fSignedCHVoxel << G4endl;
			exit(0);
		}
		else
			G4cout << "Reading " << fSignedCHVoxel << G4endl;

		G4int IDX = 0, IDY = 0, IDZ = 0, CHID = 0, VoxelID = 0;
		while (getline(infile, line))
		{
			stringstream stream(line.data());
			stream >> cpsphere >> IDX >> IDY >> IDZ >> CHID >> VoxelID;
			fCH_ID.push_back(CHID);
			fVoxel_ID.push_back(VoxelID);
		}
		infile.close();
	}
}

void TsScorerNucleusDNADamage::FindPlaceInChromosome(G4int VoxelID, G4int FiberID, G4int bpIDinFiber, G4int &VoxelNumInSphere, G4int &ChromosomeID, G4int &bpIDinChromosome)
{
	if (fVoxelNumInSphere.size() > 0)
	{
		G4int voxelIDInNucleus = fVoxelNumInSphere[VoxelID];
		G4int voxelIDonCh = fVoxel_ID[voxelIDInNucleus];

		VoxelNumInSphere = voxelIDInNucleus;
		ChromosomeID = fCH_ID[voxelIDInNucleus];
		bpIDinChromosome = voxelIDonCh * fVoxelDNAContent + FiberID * fFiberDNAContent + bpIDinFiber;
	}
	else
	{
		VoxelNumInSphere = VoxelID;
		ChromosomeID = 1;
		bpIDinChromosome = bpIDinFiber;
	}
}

G4int TsScorerNucleusDNADamage::CalculateVoxelID(G4ThreeVector position)
{
	G4int calculatedVoxelID = 0, nx = 0, ny = 0, nz = 0;
	G4cout << "Position " << position/um << " fVoxelContainerminXYZ = " << fVoxelContainerminXYZ/um
		   << " fVoxelSize = " << fVoxelSize/um << " fVoxel3Drepeat = " << fVoxel3Drepeat << G4endl;
	nx = floor((position.x() - fVoxelContainerminXYZ) / fVoxelSize);
	ny = floor((position.y() - fVoxelContainerminXYZ) / fVoxelSize);
	nz = floor((position.z() - fVoxelContainerminXYZ) / fVoxelSize);
	calculatedVoxelID = nz * fVoxel3Drepeat * fVoxel3Drepeat + ny * fVoxel3Drepeat + nx;
	return calculatedVoxelID;
}

void TsScorerNucleusDNADamage::UserHookForEndOfRun()
{
	G4cout << "\n\n------------------------------------" << G4endl;
	G4cout << "HitsOfEvents size: " << HitsOfEvents.size() << G4endl;
	G4int distinctLesions = 0;
	for (G4int id = 0; id < HitsOfEvents.size(); id++)
	{
		distinctLesions += Analyze(HitsOfEvents[id], id);
		fNtuple->Fill();
	}
	fDefineDamage->UpdateDamageAndPrimaryCount(distinctLesions, HitsOfEvents.size());
}

G4int TsScorerNucleusDNADamage::Analyze(vector<TsHitsRecord*> hits, G4int eventID)
{
	G4bool useNewMethod = true;

	fDoseInThisExposure = (1.6e-13 * eventsEdep[eventID] / MeV) / (fNucleusMass / 1000);	// Gy
	fTrackAveragedLET	= (eventsEdep[eventID] / keV) / (eventsLength[eventID] / um);
	//G4cout << "Dose in this exposure: " << fDoseInThisExposure << "; Track-averaged LET: " << fTrackAveragedLET << G4endl;
	fDefineDamage->SetEdep(eventsEdep[eventID] / MeV);
	fDefineDamage->SetLET(fTrackAveragedLET);
	fDefineDamage->SetZetaBeta_sq(fZetaBeta_sq);
	fDefineDamage->SetEventID(eventID);

	fEdepkeV = eventsEdep[eventID] / keV;
	fDoseGy = fDoseInThisExposure;
	fLETtkeVum = fTrackAveragedLET * keV/um;

	G4int distinctLesions = 0;

	if (fDoseInThisExposure >= fExposureID * fDosePerExposure / gray)
	{
		fExposureID++;
		G4cout << "Start new exposure (" << fDosePerExposure / gray << " Gy per exposure)!"
				<< " ExposureID = " << fExposureID << " - Event ID = " << eventID << G4endl;
	}

	// NEW METHOD TO COMPUTE DAMAGE
	if (useNewMethod)
	{
		// Separation in chromosomes is now done in fDefineDamage
		fDefineDamage->ComputeStrandBreaks(hits);
		std::map<G4int, std::vector<G4int>> damageInitPos = fDefineDamage->IdentifyDamageBlocks();
		distinctLesions = fDefineDamage->OutputSDDFile(damageInitPos, eventID, fExposureID);

		numSSB 				= fDefineDamage->GetNumSSB();
		numSSB_dir			= fDefineDamage->GetNumDirSSB();
		numSSB_indir		= fDefineDamage->GetNumIndirSSB();
		numDSB				= fDefineDamage->GetNumDSB();
		numDSB_dir			= fDefineDamage->GetNumDirDSB();
		numDSB_indir		= fDefineDamage->GetNumIndirDSB();
		numDSB_hybrid		= fDefineDamage->GetNumHybridDSB();
		if (fScoreOnBases)
		{
			numBaseDam			= fDefineDamage->GetNumBaseDam();
			numBaseDam_dir		= fDefineDamage->GetNumDirBaseDam();
			numBaseDam_indir	= fDefineDamage->GetNumIndirBaseDam();
			yieldBaseDam		= (G4double)numBaseDam / fDoseGy / fTotalDNAContent * 1e9;
		}
		numSSBPlus			= fDefineDamage->GetNumSSBPlus();
		numDSBPlus			= fDefineDamage->GetNumDSBPlus();
		numMoreComplex		= fDefineDamage->GetNumComplex();
		if (fExcludeShortFragment)
		{
		    Excluded_numSSB			= fDefineDamage->GetExcludedNumSSB();
		    Excluded_numSSB_dir		= fDefineDamage->GetExcludedNumDirSSB();
		    Excluded_numSSB_indir	= fDefineDamage->GetExcludedNumIndirSSB();
		    Excluded_numDSB			= fDefineDamage->GetExcludedNumDSB();
		    Excluded_numDSB_dir		= fDefineDamage->GetExcludedNumDirDSB();
		    Excluded_numDSB_indir	= fDefineDamage->GetExcludedNumIndirDSB();
		    Excluded_numDSB_hybrid	= fDefineDamage->GetExcludedNumHybridDSB();
		    if (fScoreOnBases)
		    {
		    	Excluded_baseDam		= fDefineDamage->GetExcludedNumBaseDam();
		    	Excluded_baseDam_dir	= fDefineDamage->GetExcludedNumDirBaseDam();
		    	Excluded_baseDam_indir	= fDefineDamage->GetExcludedNumIndirBaseDam();
		    }
		}
		yieldSSB			= (G4double)numSSB / fDoseGy / fTotalDNAContent * 1e9;
		yieldDSB			= (G4double)numDSB / fDoseGy / fTotalDNAContent * 1e9;
		yieldSSBPlus		= (G4double)numSSBPlus / fDoseGy / fTotalDNAContent * 1e9;
		yieldDSBPlus		= (G4double)numDSBPlus / fDoseGy / fTotalDNAContent * 1e9;
		yieldMoreComplex	= (G4double)numMoreComplex / fDoseGy / fTotalDNAContent * 1e9;
	}
	// PREVIOUS METHOD TO COMPUTE DAMAGE
	else
	{
		// Separate hits in different chromosomes
		vector<G4int> chromosomeCopyNum;
		for (G4int i = 0; i < hits.size(); i++)
			chromosomeCopyNum.push_back(hits[i]->GetChromosomeID());
		sort(chromosomeCopyNum.begin(), chromosomeCopyNum.end());
		chromosomeCopyNum.erase(unique(chromosomeCopyNum.begin(), chromosomeCopyNum.end()), chromosomeCopyNum.end()); // @suppress("Invalid arguments")

		for (G4int i = 0; i < chromosomeCopyNum.size(); i++)
		{
			vector<TsHitsRecord*> hitsOnThisChromosome;
			for (G4int j = 0; j < hits.size(); j++)
			{
				if (hits[j]->GetChromosomeID() == chromosomeCopyNum[i])
					hitsOnThisChromosome.push_back(hits[j]);
			}
			vector<TsHitsRecord*>   SBonStrand1,  SBonStrand2;
			vector<TsHitsRecord*>  SSBonStrand1, SSBonStrand2;
			vector <pair<TsHitsRecord*, TsHitsRecord*>> DSBpairs;
			fDefineDamage->SeparateHitsOnDifferentDNAStrands(hitsOnThisChromosome, SBonStrand1, true);
			fDefineDamage->SeparateHitsOnDifferentDNAStrands(hitsOnThisChromosome, SBonStrand2, false);
			fDefineDamage->DefineDSBorSSB(SBonStrand1,SBonStrand2, DSBpairs, SSBonStrand1, SSBonStrand2);
			fDefineDamage->OutputDNAdamageTuple(SBonStrand1,SBonStrand2, fOutFileName);
			fDefineDamage->OutputSDDFile(fOnlyIncludeDSBinSDD, DSBpairs, SSBonStrand1, SSBonStrand2, fOutFileName, eventID, fExposureID, chromosomeCopyNum[i]);

			if(fExcludeShortFragment)
				fDefineDamage->ExcludeShortDNAFragments(DSBpairs, SSBonStrand1, SSBonStrand2, chromosomeCopyNum[i]);

			hitsOnThisChromosome.clear();
			SBonStrand1.clear();
			SBonStrand2.clear();
			DSBpairs.clear();
			SSBonStrand1.clear();
			SSBonStrand2.clear();
		}
		fDefineDamage->OutputDNAdamageTupleHeader(fOutFileName);
		fDefineDamage->OutputDNAdamageSummary(fOutFileName);
	}
	return distinctLesions;
}

void TsScorerNucleusDNADamage::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
	TsScorerNucleusDNADamage* workerMTScorer = dynamic_cast<TsScorerNucleusDNADamage*>(workerScorer);

	fZetaBeta_sq= workerMTScorer->fZetaBeta_sq;
	for(unsigned int i = 0; i < workerMTScorer->HitsOfEvents.size(); i++)
		HitsOfEvents.push_back(workerMTScorer->HitsOfEvents[i]);
	workerMTScorer->HitsOfEvents.clear();
	for(unsigned int i = 0; i < workerMTScorer->eventsEdep.size(); i++)
		eventsEdep.push_back(workerMTScorer->eventsEdep[i]);
	workerMTScorer->eventsEdep.clear();
	for(unsigned int i = 0; i < workerMTScorer->eventsLength.size(); i++)
		eventsLength.push_back(workerMTScorer->eventsLength[i]);
	workerMTScorer->eventsLength.clear();
}
