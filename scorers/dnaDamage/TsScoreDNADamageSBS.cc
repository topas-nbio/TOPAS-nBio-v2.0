// Scorer for DNADamageStepByStep
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

#include "TsScoreDNADamageSBS.hh"

#include "TsVGeometryComponent.hh"

#include "G4SystemOfUnits.hh"
#include "G4VTouchable.hh"

#include "G4VProcess.hh"
#include "G4Molecule.hh"

TsScoreDNADamageSBS::TsScoreDNADamageSBS(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
							: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	// Initialize physical quantities
	fEdep = 0;
	fTrackAveragedLET = 0;
	fDoseInThisExposure = 0; fAccumulatedDoseInRun = 0;
	fExposureID = 0;

	// Initialize quantification of damage
	fNumSB = 0; fNumSBDirect = 0; fNumSBQuasiDirect = 0; fNumSBIndirect = 0;
	fNumSSB = 0; fNumSSBDirect = 0; fNumSSBQuasiDirect = 0; fNumSSBIndirect = 0;
	fNumDSB = 0; fNumDSBDirect = 0; fNumDSBIndirect = 0; fNumDSBDirectIndirect = 0; fNumDSBDirectQuasiDirect = 0; fNumDSBQuasiDirectQuasiDirect = 0; fNumDSBIndirectQuasiDirect = 0;
	fNumBaseDamage = 0; fNumBaseDamageDirect = 0; fNumBaseDamageQuasiDirect = 0; fNumBaseDamageIndirect = 0;
	fNumSSBPlus = 0; fNumDSBPlus = 0; fNumDSBComplex = 0;
	fYSB = 0; fYSSB = 0; fYDSB = 0; fYSSBPlus = 0; fYDSBPlus = 0; fYDSBComplex = 0;
	fYBaseDam = 0;

	// Add the basic hierarchical level (base pair)
	fHierarchicalLevels.push_back("BasePair");

	//---------------
	// Get parameters
	//---------------
	if (!fPm->ParameterExists(GetFullParmName("NumberOfHistoriesInRun")))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("NumberOfHistoriesInRun") << " has to be specified to manage the events and runs correctly." << G4endl;
		fPm->AbortSession(1);
	}
	fNumberOfHistoriesInRun = fPm->GetIntegerParameter(GetFullParmName("NumberOfHistoriesInRun"));

	// Output filename
	fOutFileName = fPm->GetStringParameter(GetFullParmName("OutputFile"));

	// Parameters for material filter
	fBasePairDepth = 0;
	if (fPm->ParameterExists(GetFullParmName("BasePairPositionAtGeometricHierarchy")))
		fBasePairDepth = fPm->GetIntegerParameter(GetFullParmName("BasePairPositionAtGeometricHierarchy"));
	G4String* strand1Materials = NULL;
	G4String* strand2Materials = NULL;
	G4int strand1Length = 0;
	G4int strand2Length = 0;
	if (fPm->ParameterExists(GetFullParmName("Strand1MaterialNames")))
	{
		strand1Materials = fPm->GetStringVector(GetFullParmName("Strand1MaterialNames"));
		strand1Length = fPm->GetVectorLength(GetFullParmName("Strand1MaterialNames"));
	}
	else
		strand1Materials[0] = "G4_WATER";
	if (fPm->ParameterExists(GetFullParmName("Strand2MaterialNames")))
	{
		strand2Materials = fPm->GetStringVector(GetFullParmName("Strand2MaterialNames"));
		strand2Length = fPm->GetVectorLength(GetFullParmName("Strand2MaterialNames"));
	}
	else
		strand2Materials[0] = "G4_WATER";
	for (G4int i = 0; i < strand1Length; i++)
		fStrand1Materials.push_back(GetMaterial(strand1Materials[i]));
	for (G4int i = 0; i < strand2Length; i++)
		fStrand2Materials.push_back(GetMaterial(strand2Materials[i]));

	fScoringRadius = 0;
	if (fPm->ParameterExists(GetFullParmName("ScoringRadius")))
		fScoringRadius = fPm->GetDoubleParameter(GetFullParmName("ScoringRadius"), "Length");

	// Options for direct damage
	fDirectDamageThreshold = 11.75 * eV; // 17,5 eV for half-cylinder
	if (fPm->ParameterExists(GetFullParmName("DirectDamageThreshold")))
		fDirectDamageThreshold = fPm->GetDoubleParameter(GetFullParmName("DirectDamageThreshold"), "Energy");
	fUseLinearProbabilityForDirectDamage = false;
	if ( fPm->ParameterExists(GetFullParmName("UseLinearProbabilityForDirectDamage")) )
		fUseLinearProbabilityForDirectDamage = fPm->GetBooleanParameter(GetFullParmName("UseLinearProbabilityForDirectDamage"));
	if (!fUseLinearProbabilityForDirectDamage && (fPm->ParameterExists(GetFullParmName("LowerLimitForLinearProbabilityFunction")) || fPm->ParameterExists(GetFullParmName("UpperLimitForLinearProbabilityFunction"))))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("UseLinearProbabilityForDirectDamage") << " is set to False and limits for linear probability function in DNA components were given." << G4endl;
		fPm->AbortSession(1);
	}
	else if (fUseLinearProbabilityForDirectDamage)
	{
		fLowerLimitLinearProbability = 5 * eV;
		if (fPm->ParameterExists(GetFullParmName("LowerLimitForLinearProbabilityFunction")))
			fLowerLimitLinearProbability = fPm->GetDoubleParameter(GetFullParmName("LowerLimitForLinearProbabilityFunction"), "Energy");
		fUpperLimitLinearProbability = 37.5 * eV;
		if (fPm->ParameterExists(GetFullParmName("UpperLimitForLinearProbabilityFunction")))
			fUpperLimitLinearProbability = fPm->GetDoubleParameter(GetFullParmName("UpperLimitForLinearProbabilityFunction"), "Energy");
	}
	// Options for quasi-direct damage
	fProbabilityOfChargeTransferFromHydrationShellToBackbone = 0.33333;
	if (fPm->ParameterExists(GetFullParmName("ProbabilityOfChargeTransferFromHydrationShellToBackbone")))
		fProbabilityOfChargeTransferFromHydrationShellToBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfChargeTransferFromHydrationShellToBackbone"));

	// Options for indirect damage
	fAlwaysScavengeSpeciesInDNAComponents = false;
	if (fPm->ParameterExists(GetFullParmName("AlwaysScavengeSpeciesInDNAComponents")))
		fAlwaysScavengeSpeciesInDNAComponents = fPm->GetBooleanParameter(GetFullParmName("AlwaysScavengeSpeciesInDNAComponents"));

	if (fAlwaysScavengeSpeciesInDNAComponents && (fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBackbone")) || fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBase"))))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("AlwaysScavengeSpeciesInDNAComponents") << " is set to True and probabilities for scavenging in DNA components were given." << G4endl;
		fPm->AbortSession(1);
	}
	if (!fAlwaysScavengeSpeciesInDNAComponents)
	{
		fProbabilityOfScavengingInBackbone = 0.25; // 0.0585 for half-cylinder
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBackbone")))
			fProbabilityOfScavengingInBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfScavengingInBackbone"));
		fProbabilityOfScavengingInBase = 1.0;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBase")))
			fProbabilityOfScavengingInBase = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfScavengingInBase"));
		fProbabilityOfDamageInBackbone = 0.55;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBackbone")))
			fProbabilityOfDamageInBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBackbone"));
		fProbabilityOfDamageInBase = 1.0;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBase")))
			fProbabilityOfDamageInBase = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBase"));
	}
	else
	{
		fProbabilityOfScavengingInBackbone = 1.0;
		fProbabilityOfScavengingInBase = 1.0;
		fProbabilityOfDamageInBackbone = 0.4;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBackbone")))
			fProbabilityOfDamageInBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBackbone"));
		fProbabilityOfDamageInBase = 0.4;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBase")))
			fProbabilityOfDamageInBase = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBase"));
	}

	fScavengeInHistones = true;
	if (fPm->ParameterExists(GetFullParmName("ScavengeInHistones")))
		fScavengeInHistones = fPm->GetBooleanParameter(GetFullParmName("ScavengeInHistones"));

	// Classify damage as SSBs and DSBs
	fNumberOfBasePairsForDSB = 10;
	if (fPm->ParameterExists(GetFullParmName("MaximumBasePairDistanceToConsiderDSB")))
		fNumberOfBasePairsForDSB = fPm->GetIntegerParameter(GetFullParmName("MaximumBasePairDistanceToConsiderDSB"));

	// Foci scoring, creation of foci images
	fScoreFoci = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreNumberOfFoci")))
		fScoreFoci = fPm->GetBooleanParameter(GetFullParmName("ScoreNumberOfFoci"));
	fFociSizes.push_back(500 * nm);
	if (fPm->ParameterExists(GetFullParmName("FociSizes")))
	{
		G4double* sizes = fPm->GetDoubleVector(GetFullParmName("FociSizes"), "Length");
		G4int vectorLength = fPm->GetVectorLength(GetFullParmName("FociSizes"));
		fFociSizes.clear();
		for (G4int i = 0; i < vectorLength; i++)
			fFociSizes.push_back(sizes[i]);
	}
	fGet3DFociImage = false;
	if (fPm->ParameterExists(GetFullParmName("Get3DFociImage")))
		fGet3DFociImage = fPm->GetBooleanParameter(GetFullParmName("Get3DFociImage"));
	fGet2DFociImage = false;
	if (fPm->ParameterExists(GetFullParmName("Get2DFociImages")))
		fGet2DFociImage = fPm->GetBooleanParameter(GetFullParmName("Get2DFociImages"));

	if (fGet3DFociImage || fGet2DFociImage)
	{
		if (fGet2DFociImage)
		{
			G4String* planes;
			G4int vectorLength = 0;
			if (fPm->ParameterExists(GetFullParmName("2DFociImagePlanes")))
			{
				planes = fPm->GetStringVector(GetFullParmName("2DFociImagePlanes"));
				vectorLength = fPm->GetVectorLength(GetFullParmName("2DFociImagePlanes"));
			}
			for (G4int i = 0; i < vectorLength; i++)
				f2DPlanesForFociImage.push_back(planes[i]);
		}
		fMicroscopePSFShape = "Gaussian";
		if (fPm->ParameterExists(GetFullParmName("MicroscopePSFShape")))
			fMicroscopePSFShape = fPm->GetStringParameter(GetFullParmName("MicroscopePSFShape"));
		fMicroscopePSFWidth = 500 * nm;
		if (fPm->ParameterExists(GetFullParmName("MicroscopePSFWidth")))
			fMicroscopePSFWidth = fPm->GetDoubleParameter(GetFullParmName("MicroscopePSFWidth"), "Length");
		f2DImageResolutions.push_back(5.0 * nm);
		if (fPm->ParameterExists(GetFullParmName("Foci2DImageResolutions")))
		{
			G4double* resolutions = fPm->GetDoubleVector(GetFullParmName("Foci2DImageResolutions"), "Length");
			G4int vectorLength = fPm->GetVectorLength(GetFullParmName("Foci2DImageResolutions"));
			f2DImageResolutions.clear();
			for (G4int i = 0; i < vectorLength; i++)
				f2DImageResolutions.push_back(resolutions[i]);
		}
		f3DImageResolution = 500 *nm;
		if (fPm->ParameterExists(GetFullParmName("Foci3DImageResolution")))
			f3DImageResolution = fPm->GetDoubleParameter(GetFullParmName("Foci3DImageResolution"), "Length");
		fImXmin = -5 * um;
		if (fPm->ParameterExists(GetFullParmName("ImageXMinPosition")))
			fImXmin = fPm->GetDoubleParameter(GetFullParmName("ImageXMinPosition"), "Length");
		fImXmax = 5 * um;
		if (fPm->ParameterExists(GetFullParmName("ImageXMaxPosition")))
			fImXmax = fPm->GetDoubleParameter(GetFullParmName("ImageXMaxPosition"), "Length");
		fImYmin = -5 * um;
		if (fPm->ParameterExists(GetFullParmName("ImageYMinPosition")))
			fImYmin = fPm->GetDoubleParameter(GetFullParmName("ImageYMinPosition"), "Length");
		fImYmax = 5 * um;
		if (fPm->ParameterExists(GetFullParmName("ImageYMaxPosition")))
			fImYmax = fPm->GetDoubleParameter(GetFullParmName("ImageYMaxPosition"), "Length");
		fImZmin = -5 * um;
		if (fPm->ParameterExists(GetFullParmName("ImageZMinPosition")))
			fImZmin = fPm->GetDoubleParameter(GetFullParmName("ImageZMinPosition"), "Length");
		fImZmax = 5 * um;
		if (fPm->ParameterExists(GetFullParmName("ImageZMaxPosition")))
			fImZmax = fPm->GetDoubleParameter(GetFullParmName("ImageZMaxPosition"), "Length");
	}

	// Considering fragments
	fExcludeShortFragments = false;
	if (fPm->ParameterExists(GetFullParmName("ExcludeShortFragments")))
		fExcludeShortFragments = fPm->GetBooleanParameter(GetFullParmName("ExcludeShortFragments"));
	if (!fExcludeShortFragments && (fPm->ParameterExists(GetFullParmName("LowerThresholdForFragmentDetection")) || fPm->ParameterExists(GetFullParmName("UpperThresholdForFragmentDetection"))))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("ExcludeShortFragments") << " is set to False and limits for fragment detection were given." << G4endl;
		fPm->AbortSession(1);
	}
	else if (fExcludeShortFragments)
	{
		fLowerThresholdForFragmentDetection = 0.0;
		if (fPm->ParameterExists(GetFullParmName("LowerThresholdForFragmentDetection")))
			fLowerThresholdForFragmentDetection = fPm->GetIntegerParameter(GetFullParmName("LowerThresholdForFragmentDetection"));
		fUpperThresholdForFragmentDetection = 3E8;
		if (fPm->ParameterExists(GetFullParmName("UpperThresholdForFragmentDetection")))
			fUpperThresholdForFragmentDetection = fPm->GetIntegerParameter(GetFullParmName("UpperThresholdForFragmentDetection"));
	}

	// Stop at a given dose
	fStopAtDose = 3E8;
	if (fPm->ParameterExists(GetFullParmName("StopTrackingAtDose")))
		fStopAtDose = fPm->GetDoubleParameter(GetFullParmName("StopTrackingAtDose"), "Dose") / gray;

	// Options for the output
	fWriteCSVWithExtensiveDamage = false;
	if (fPm->ParameterExists(GetFullParmName("WriteCSVOutputWithAllDamageSpecification")))
		fWriteCSVWithExtensiveDamage = fPm->GetBooleanParameter(GetFullParmName("WriteCSVOutputWithAllDamageSpecification"));
	fScoreDirectDamage = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreDirectDamage")))
		fScoreDirectDamage = fPm->GetBooleanParameter(GetFullParmName("ScoreDirectDamage"));
	fScoreIndirectDamage = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreIndirectDamage")))
		fScoreIndirectDamage = fPm->GetBooleanParameter(GetFullParmName("ScoreIndirectDamage"));
	fScoreQuasiDirectDamage = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreQuasiDirectDamage")))
		fScoreQuasiDirectDamage = fPm->GetBooleanParameter(GetFullParmName("ScoreQuasiDirectDamage"));
	fScoreOnBases = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreBaseDamages")))
		fScoreOnBases = fPm->GetBooleanParameter(GetFullParmName("ScoreBaseDamages"));
	fScoreOnBackbones = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreBackboneDamages")))
		fScoreOnBackbones = fPm->GetBooleanParameter(GetFullParmName("ScoreBackboneDamages"));
	fBreakDownPerDamageOrigin = true;
	if (fPm->ParameterExists(GetFullParmName("BreakDownOutputPerDamageOrigin")))
		fBreakDownPerDamageOrigin = fPm->GetBooleanParameter(GetFullParmName("BreakDownOutputPerDamageOrigin"));

	// For SDD specification
	fDosePerExposure = 1 * gray;
	if (fPm->ParameterExists(GetFullParmName("DosePerExposure")))
		fDosePerExposure = fPm->GetDoubleParameter(GetFullParmName("DosePerExposure"), "Dose");
	fOnlyIncludeDSBinSDD = false;
	if (fPm->ParameterExists(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD")))
		fOnlyIncludeDSBinSDD = fPm->GetBooleanParameter(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD"));
	fWriteMinimalSDDOutput = false;
	if (fPm->ParameterExists(GetFullParmName("WriteMinimalSDDOutput")))
		fWriteMinimalSDDOutput = fPm->GetBooleanParameter(GetFullParmName("WriteMinimalSDDOutput"));
	fPrimaryParticle = "proton";
	if (fPm->ParameterExists(GetFullParmName("PrimaryParticle")))
		fPrimaryParticle = fPm->GetStringParameter(GetFullParmName("PrimaryParticle"));

	// Parameters for the SDD header
	fAuthor = "@";
	if ( fPm->ParameterExists(GetFullParmName("AuthorName")) )
		fAuthor = fPm->GetStringParameter(GetFullParmName("AuthorName"));
	fSimulationDetails = "Sim details";
	if ( fPm->ParameterExists(GetFullParmName("SimulationDetails")) )
		fSimulationDetails = fPm->GetStringParameter(GetFullParmName("SimulationDetails"));
	fSourceDetails = "Source details";
	if ( fPm->ParameterExists(GetFullParmName("SourceDetails")) )
		fSourceDetails = fPm->GetStringParameter(GetFullParmName("SourceDetails"));
	fSourceType = 1;
	if ( fPm->ParameterExists(GetFullParmName("SourceType")) )
		fSourceType = fPm->GetIntegerParameter(GetFullParmName("SourceType"));
	fMeanEnergy = 0.0;
	if ( fPm->ParameterExists(GetFullParmName("MeanEnergy")))
		fMeanEnergy = fPm->GetDoubleParameter(GetFullParmName("MeanEnergy"), "Energy");
	fEnergyDist = "M, 0";
	if ( fPm->ParameterExists(GetFullParmName("EnergyDistribution")) )
		fEnergyDist = fPm->GetStringParameter(GetFullParmName("EnergyDistribution"));
	fIrrTarget = "";
	if ( fPm->ParameterExists(GetFullParmName("IrradiationTarget")) )
		fIrrTarget = fPm->GetStringParameter(GetFullParmName("IrradiationTarget"));
	fCellCycle = "0";
	if ( fPm->ParameterExists(GetFullParmName("CellCycleStage")))
		fCellCycle = fPm->GetStringParameter(GetFullParmName("CellCycleStage"));
	fDNAStructure = "0, 1";
	if ( fPm->ParameterExists(GetFullParmName("DNAStructure")))
		fDNAStructure = fPm->GetStringParameter(GetFullParmName("DNAStructure"));
	fInVitroOrInVivo = 0;
	if ( fPm->ParameterExists(GetFullParmName("InVitroOrInVivo")) )
		fInVitroOrInVivo = fPm->GetIntegerParameter(GetFullParmName("InVitroOrInVivo"));
	fProliferationStatus = "1";
	if ( fPm->ParameterExists(GetFullParmName("ProliferationStatus")))
		fProliferationStatus = fPm->GetStringParameter(GetFullParmName("ProliferationStatus"));
	fMicroenvironment = "20, 0.01";
	if ( fPm->ParameterExists(GetFullParmName("Microenvironment")))
		fMicroenvironment = fPm->GetStringParameter(GetFullParmName("Microenvironment"));
	fTime = 0;
	if ( fPm->ParameterExists(GetFullParmName("Time")))
		fTime = fPm->GetDoubleParameter(GetFullParmName("Time"), "Time");
	fAddInfo = "";
	if ( fPm->ParameterExists(GetFullParmName("AdditionalInfo")))
		fAddInfo = fPm->GetStringParameter(GetFullParmName("AdditionalInfo"));

	// =============================
	//       PRINT PARAMETERS
	// =============================
	G4cout << "*********************************************************************************" << G4endl;
	if (fScoreDirectDamage)
	{
		G4cout << "DIRECT DAMAGE" << G4endl;
		G4cout << "--------------" << G4endl;
		if (!fUseLinearProbabilityForDirectDamage)
			G4cout << "Single energy threshold for damage = " << fDirectDamageThreshold/eV << " eV" << G4endl;
		else
			G4cout << "Linearly increasing probability of damage from 0 to 1 from " << fLowerLimitLinearProbability/eV << " eV to " << fUpperLimitLinearProbability/eV << " eV" << G4endl;
	}
	if (fScoreQuasiDirectDamage)
	{
		G4cout << "QUASI-DIRECT DAMAGE" << G4endl;
		G4cout << "-------------------" << G4endl;
		G4cout << "Probability for charge transfer from hydration shell to backbone to produce strand break = " << fProbabilityOfChargeTransferFromHydrationShellToBackbone << G4endl;
	}
	if (fScoreIndirectDamage)
	{
		G4cout << "INDIRECT DAMAGE" << G4endl;
		G4cout << "-------------------" << G4endl;
		if (fAlwaysScavengeSpeciesInDNAComponents)
			G4cout << "Probability of scavenging in bases = 1.0; Probability of scavenging in backbones = 1.0" << G4endl;
		else
			G4cout << "Probability of scavenging in bases = " << fProbabilityOfScavengingInBase << "; Probability of scavenging in backbones = " << fProbabilityOfScavengingInBackbone << G4endl;
		G4cout << "Probability of base damage after scavenging = " << fProbabilityOfDamageInBase << "; Probability of strand break after scavenging = " << fProbabilityOfDamageInBackbone << G4endl;
	}
	G4cout << "*********************************************************************************" << G4endl;

	// Register variables in nTuple
	fNtuple->RegisterColumnD(&fEdep, "Energy_imparted_per_event", "keV");
	fNtuple->RegisterColumnD(&fDoseInThisExposure, "Dose_per_event_Gy", "");
	//fNtuple->RegisterColumnD(&fTrackAveragedLET, "LET_kev/um", "");

	if (fScoreOnBackbones)
	{
		fNtuple->RegisterColumnD(&fYDSB, "DSB/Gy/Gbp", "");
		fNtuple->RegisterColumnD(&fYSSB, "SSB/Gy/Gbp", "");
		fNtuple->RegisterColumnD(&fYSB, "SB/Gy/Gbp", "");
		fNtuple->RegisterColumnD(&fYSSBPlus, "SSB+/Gy/Gbp", "");
		fNtuple->RegisterColumnD(&fYDSBPlus, "DSB+/Gy/Gbp", "");
		fNtuple->RegisterColumnD(&fYDSBComplex, "MoreComplexDamage/Gy/Gbp", "");
	}
	if (fScoreOnBases)
		fNtuple->RegisterColumnD(&fYBaseDam, "BD/Gy/Gbp", "");
	if (fBreakDownPerDamageOrigin)
	{
		if (fScoreOnBackbones)
		{
			fNtuple->RegisterColumnI(&fNumDSB, "DSBs");
			if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumDSBDirect, "DSBs_Direct");
			if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumDSBIndirect, "DSBs_Indirect");
			if (fScoreDirectDamage && fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumDSBDirectIndirect, "DSBs_Hybrid");
			if (fScoreDirectDamage && fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumDSBDirectQuasiDirect, "DSBs_Direct_WithOneQuasiDirect");
			if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumDSBQuasiDirectQuasiDirect, "DSBs_Direct_WithBothQuasiDirect");
			if (fScoreIndirectDamage && fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumDSBIndirectQuasiDirect, "DSBs_Hybrid_WithOneQuasiDirect");
			fNtuple->RegisterColumnI(&fNumSSB, "SSBs");
			if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumSSBDirect, "SSBs_Direct");
			if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumSSBQuasiDirect, "SSBs_QuasiDirect");
			if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumSSBIndirect, "SSBs_Indirect");
			fNtuple->RegisterColumnI(&fNumSB, "SBs");
			if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumSBDirect, "SBs_Direct");
			if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumSBQuasiDirect, "SBs_QuasiDirect");
			if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumSBIndirect, "SBs_Indirect");
			fNtuple->RegisterColumnI(&fNumSSBPlus, "SSB+s");
			fNtuple->RegisterColumnI(&fNumDSBPlus, "DSB+s");
			fNtuple->RegisterColumnI(&fNumDSBComplex, "More complex damages");
		}
		if (fScoreOnBases)
		{
			fNtuple->RegisterColumnI(&fNumBaseDamage, "BDs");
			if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumBaseDamageDirect, "BDs_Direct");
			if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumBaseDamageQuasiDirect, "BDs_QuasiDirect");
			if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumBaseDamageIndirect, "BDs_Indirect");
		}
	}

	if (fScoreFoci)
	{
		fNumFoci1 = 0; fNumFoci2 = 0; fNumFoci3 = 0; fNumFoci4 = 0; fNumFoci5 = 0;
		if (fFociSizes.size() >= 1) fNtuple->RegisterColumnI(&fNumFoci1, "Foci_" + std::to_string(int(fFociSizes[0]*1e6)) + "nm");
		if (fFociSizes.size() >= 2) fNtuple->RegisterColumnI(&fNumFoci2, "Foci_" + std::to_string(int(fFociSizes[1]*1e6)) + "nm");
		if (fFociSizes.size() >= 3) fNtuple->RegisterColumnI(&fNumFoci3, "Foci_" + std::to_string(int(fFociSizes[2]*1e6)) + "nm");
		if (fFociSizes.size() >= 4) fNtuple->RegisterColumnI(&fNumFoci4, "Foci_" + std::to_string(int(fFociSizes[3]*1e6)) + "nm");
		if (fFociSizes.size() >= 5) fNtuple->RegisterColumnI(&fNumFoci5, "Foci_" + std::to_string(int(fFociSizes[4]*1e6)) + "nm");
	}
	// Initialize and setup damage computer
	fDamageCalculator = new TsDNADamageCalculator();
	fDamageCalculator->SetDistanceBasePairsForDSB(fNumberOfBasePairsForDSB);
	fDamageCalculator->SetDirectDamageAsLinearProbability(fUseLinearProbabilityForDirectDamage);
	if (fUseLinearProbabilityForDirectDamage)
	{
		fDamageCalculator->SetLowerLimitForLinearProbability(fLowerLimitLinearProbability);
		fDamageCalculator->SetUpperLimitForLinearProbability(fUpperLimitLinearProbability);
	}
	else
		fDamageCalculator->SetDirectDamageThreshold(fDirectDamageThreshold);
	fDamageCalculator->SetExcludeShortFragments(fExcludeShortFragments);
	if (fExcludeShortFragments)
	{
		fDamageCalculator->SetLowerThresholdForFragmentDetection(fLowerThresholdForFragmentDetection);
		fDamageCalculator->SetUpperThresholdForFragmentDetection(fUpperThresholdForFragmentDetection);
	}
	fDamageCalculator->SetOutputFileName(fOutFileName);
	fDamageCalculator->SetOutputFileMode(fOutFileMode);
	fDamageCalculator->SetWriteCSVFile(fWriteCSVWithExtensiveDamage);
	fDamageCalculator->SetMinimalModeForSDD(fWriteMinimalSDDOutput);
	fDamageCalculator->SetReturnOnlyDSBInSDD(fOnlyIncludeDSBinSDD);
	fDamageCalculator->OutputSDDHeader(fWriteMinimalSDDOutput, fPrimaryParticle, fMeanEnergy, fDosePerExposure, fChromosomeContents, fScoreIndirectDamage, fScoreOnBases,
			fAuthor, fSimulationDetails, fSourceDetails, fSourceType, fEnergyDist, fIrrTarget, fCellCycle, fDNAStructure, fInVitroOrInVivo, fProliferationStatus, fMicroenvironment, fTime, fAddInfo);

	if (fScoreFoci)
	{
		fFociAnalyzer = new TsFociAnalysis(fComponent);
		fFociAnalyzer->SetFociSizes(fFociSizes);
		if (fGet3DFociImage || fGet2DFociImage)
		{
			fFociAnalyzer->SetPSFShape(fMicroscopePSFShape);
			fFociAnalyzer->SetPSFWidth(fMicroscopePSFWidth);
			fFociAnalyzer->Set2DImageResolutions(f2DImageResolutions);
			fFociAnalyzer->Set3DImageResolution(f3DImageResolution);
			fFociAnalyzer->SetMinX(fImXmin);
			fFociAnalyzer->SetMaxX(fImXmax);
			fFociAnalyzer->SetMinY(fImYmin);
			fFociAnalyzer->SetMaxY(fImYmax);
			fFociAnalyzer->SetMinZ(fImZmin);
			fFociAnalyzer->SetMaxZ(fImZmax);
			if (fGet2DFociImage)
				fFociAnalyzer->SetPlanesFor2DFociImages(f2DPlanesForFociImage);
		}
	}
}

TsScoreDNADamageSBS::~TsScoreDNADamageSBS() {}

G4bool TsScoreDNADamageSBS::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive)
	{
		fSkippedWhileInactive++;
		return false;
	}

	// Stops tracking if accumulated dose is higher than specified (only for a new event! A single event is always completed)
	if (fAccumulatedDoseInRun >= fStopAtDose)
	{
		G4cout << "Tracking is stopped because the dose limit (" << fStopAtDose << " Gy) was reached. No more particles will be tracked/scored." << G4endl;
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		return false;
	}

	// Gets position
	G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();

	// Checks with respect to the scoring radius
	if (fScoringRadius > 0)
	{
		G4bool withinScoringRadius = ((pow(pos.x(), 2)+pow(pos.y(), 2)+pow(pos.z(), 2)) < pow(fScoringRadius, 2));
		if (!withinScoringRadius)
			return false;
	}

	// Accumulates energy for this event
	G4double edep = aStep->GetTotalEnergyDeposit();
	fEdep += edep;

	// Increments the number of steps this track has made
	G4int trackID = aStep->GetTrack()->GetTrackID();
	fTrackSteps[trackID] += 1;

	// Gets current volume and last volume for this track. Updates new volume for this track
	G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	G4String volID = touchable->GetVolume(fBasePairDepth)->GetName() + std::to_string(touchable->GetVolume(fBasePairDepth)->GetCopyNo());
	G4bool enteringInNewVolume = false;
	if (volID != fTrackLastVolume[trackID])
		enteringInNewVolume = true;
	fTrackLastVolume[trackID] = volID;

	// Checks materials for determining if step is happening in DNA components
	G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
	G4bool materialMatched = false;
	for (unsigned int i = 0; i < fStrand1Materials.size(); i++)
	{
		if (material == fStrand1Materials[i])
		{
			materialMatched = true;
			break;
		}
	}
	if (!materialMatched)
	{
		for (unsigned int i = 0; i < fStrand2Materials.size(); i++)
		{
			if (material == fStrand2Materials[i])
			{
				materialMatched = true;
				break;
			}
		}
	}
	// Goes on only if DNA materials has been matched
	if (materialMatched)
	{
		std::vector<G4int> hierarchicalIDs;
		// Gets IDs in the different hierarchical levels
		hierarchicalIDs.push_back(touchable->GetVolume(fBasePairDepth)->GetCopyNo());
		for (unsigned int i = 1; i < fHierarchicalLevels.size(); i++)
			hierarchicalIDs.push_back(touchable->GetCopyNumber(fBasePairDepth + i));

		std::pair<G4int, G4int> compAndStrandID = GetDNAComponentAndStrandID(touchable);
		G4int componentID = compAndStrandID.first;
		G4int strandID = compAndStrandID.second;
		// Sets base pair ID to -1 if histone is touched
		if (componentID == -1) hierarchicalIDs[0] = -1;

		// Gets particle and process info
		const G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
		G4String particleName = particle->GetParticleName();
		G4String processName = (G4String)aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

		// Sets hit information
		TsHitInDNA* hit = new TsHitInDNA();
		hit->SetEventID(GetEventID() + fNumberOfHistoriesInRun * GetRunID());
		hit->SetEdep(edep);
		hit->SetParticleName(particleName);
		hit->SetHierarchicalIDs(hierarchicalIDs);
		hit->SetStrandNumber(strandID);
		hit->SetDNAComponentID(componentID);
		hit->SetPosition(pos);
		// Default options
		hit->SetChromosomeID(1);
		hit->SetBasePairID(hierarchicalIDs[0]);

		// Adds direct damage
		if (trackID >= 0 && edep > 0 && fScoreDirectDamage && ((componentID == base && fScoreOnBases) || (componentID == backbone && fScoreOnBackbones)))
		{
			hit->SetDamageType(direct);
			fHits.push_back(hit);
			return true;
		}
		// Adds quasi-direct damage
		if (componentID == hydrationshell && strstr(processName, "Ionisation") != NULL && fScoreQuasiDirectDamage)
		{
			if (G4UniformRand() < fProbabilityOfChargeTransferFromHydrationShellToBackbone)
			{
				hit->SetDNAComponentID(backbone);
				hit->SetDamageType(quasidirect);
			}
			else
			{
				hit->SetDNAComponentID(base);
				hit->SetDamageType(quasidirect);
			}
			fHits.push_back(hit);
			return true;
		}
		// Adds indirect damage
		if (trackID < 0 && fScoreIndirectDamage)
		{
			G4String speciesName = GetMolecule(aStep->GetTrack())->GetName();
			G4bool isSpeciesToKill = (speciesName == "OH^0" || speciesName == "e_aq^-1" || speciesName == "H^0");
			G4bool isHydroxil = (speciesName == "OH^0");
			G4bool isHydElectron =  (speciesName == "e_aq^-1");
			//G4cout << "Species: " << speciesName << " in " << volumeName << " with trackID: " << trackID << " at time " << aStep->GetTrack()->GetLocalTime() << G4endl;
			// Kills all species generated inside DNA volumes except for the hydration shell
			if (fTrackSteps[trackID] == 1 && componentID != hydrationshell)
			{
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				delete hit;
				return false;
			}
			// Makes damage to bases. OH and e_aq induce damage to bases. Only one step per species is considered (otherwise all would end up reacting), so we use the enteringInNewVolume flag
			else if ((isHydroxil || isHydElectron) && componentID == base && enteringInNewVolume && fScoreOnBases)
			{
				hit->SetEdep(-0.001 * eV);
				G4bool scavenged = false;
				if (fAlwaysScavengeSpeciesInDNAComponents)
					scavenged = true;
				else
					if (G4UniformRand() < fProbabilityOfScavengingInBase) scavenged = true;
				if (scavenged)
				{
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					if (G4UniformRand() < fProbabilityOfDamageInBase)
					{
						hit->SetDamageType(indirect);
						fHits.push_back(hit);
						return true;
					}
				}
				delete hit;
				return false;
			}
			// Makes damage to backbones. Only OH induces damage to backbones. Only one step per species is considered (otherwise all would end up reacting), so we use the enteringInNewVolume flag
			else if (isHydroxil && componentID == backbone && enteringInNewVolume && fScoreOnBackbones)
			{
				hit->SetEdep(-0.001 * eV);
				G4bool scavenged = false;
				if (fAlwaysScavengeSpeciesInDNAComponents)
					scavenged = true;
				else
					if (G4UniformRand() < fProbabilityOfScavengingInBackbone) scavenged = true;
				if (scavenged)
				{
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					if (G4UniformRand() < fProbabilityOfDamageInBackbone)
					{
						hit->SetDamageType(indirect);
						fHits.push_back(hit);
						return true;
					}
				}
				delete hit;
				return false;
			}
			// Scavenge species by histones
			else if (isSpeciesToKill && fScavengeInHistones && componentID == histone)
			{
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				delete hit;
				return false;
			}
		}
		delete hit;
	}

	return false;
}

std::pair<G4int,G4int> TsScoreDNADamageSBS::GetDNAComponentAndStrandID(G4TouchableHistory* touchable)
{
	G4String volumeName = touchable->GetVolume(fBasePairDepth)->GetName();
	// Gets strand number and DNA component ID (see header file for component IDs)
	std::pair<G4int, G4int> compAndStandId;
	if (strstr(volumeName, "Base1") != NULL || strstr(volumeName, "BasePair") != NULL) { compAndStandId.first = base; compAndStandId.second = 1; }
	else if (strstr(volumeName, "Base2") != NULL) { compAndStandId.first = base; compAndStandId.second = 2; }
	else if (strstr(volumeName, "Backbone1") != NULL) { compAndStandId.first = backbone; compAndStandId.second = 1; }
	else if (strstr(volumeName, "Backbone2") != NULL) { compAndStandId.first = backbone; compAndStandId.second = 2; }
	else if (strstr(volumeName, "HydrationShell1") != NULL) { compAndStandId.first = hydrationshell; compAndStandId.second = 1; }
	else if (strstr(volumeName, "HydrationShell2") != NULL) { compAndStandId.first = hydrationshell; compAndStandId.second = 2; }
	else if (strstr(volumeName, "Histone") != NULL) { compAndStandId.first = histone; }
	return compAndStandId;
}

void TsScoreDNADamageSBS::AccumulateEvent()
{
	fCollectionsOfHits.push_back(fHits);
	fHits.clear();
	G4double edep = fEdep;
	fEventsEdep.push_back(edep);
	fEdep = 0.;
	fAccumulatedDoseInRun += CalculateDoseInGray(edep);
}

void TsScoreDNADamageSBS::UserHookForEndOfRun()
{
	// Print info
	G4cout << "--------------------" << G4endl;
	G4cout << "Number of events comprised in Run: " << fCollectionsOfHits.size() << G4endl;
	G4int numberOfLesions = 0;
	for (unsigned int i = 0; i < fCollectionsOfHits.size(); i++)
	{
		G4int lesionsThisEvent = Analyze(fCollectionsOfHits[i], i);
		numberOfLesions += lesionsThisEvent;
		// Only fill if there is any damage
		if (lesionsThisEvent > 0) fNtuple->Fill();
	}
	if (fScoreFoci)
	{
		if (fGet3DFociImage) fFociAnalyzer->Produce3DImage(fDSBPositionsInRun);
		if (fGet2DFociImage) fFociAnalyzer->Produce2DImages(fDSBPositionsInRun);
	}
	G4cout << "Number of lesions in Run: " << numberOfLesions << G4endl;
	fCollectionsOfHits.clear();
	fEventsEdep.clear();
}

G4int TsScoreDNADamageSBS::Analyze(std::vector<TsHitInDNA*> hits, G4int eventID)
{
	fEdep = fEventsEdep[eventID];
	fDoseInThisExposure = CalculateDoseInGray(fEventsEdep[eventID]);
	fDamageCalculator->SetEventID(eventID);
	if (fDoseInThisExposure >= fExposureID * fDosePerExposure / gray)
	{
		fExposureID++;
		G4cout << "Start new exposure (" << fDosePerExposure / gray << " Gy per exposure)" << " - Exposure ID: " << fExposureID << " - Event ID: " << eventID << G4endl;
	}

	G4int numberOfLesions = 0;
	fDamageCalculator->ComputeStrandBreaks(hits);
	std::map<G4int, std::vector<G4int>> initialPosDamageSites = fDamageCalculator->GetDamageSites();
	numberOfLesions = fDamageCalculator->OutputSDDFile(initialPosDamageSites, eventID, fExposureID, fChromosomeContents);

	fNumSB = fDamageCalculator->GetSB();
	fNumSBDirect = fDamageCalculator->GetSBDirect();
	fNumSBQuasiDirect = fDamageCalculator->GetSBQuasiDirect();
	fNumSBIndirect = fDamageCalculator->GetSBIndirect();
	fNumSSB = fDamageCalculator->GetSSB();
	fNumSSBDirect = fDamageCalculator->GetSSBDirect();
	fNumSSBQuasiDirect = fDamageCalculator->GetSSBQuasiDirect();
	fNumSSBIndirect = fDamageCalculator->GetSSBIndirect();
	fNumDSB = fDamageCalculator->GetDSB();
	fNumDSBDirect = fDamageCalculator->GetDSBDirect();
	fNumDSBIndirect = fDamageCalculator->GetDSBIndirect();
	fNumDSBDirectIndirect = fDamageCalculator->GetDSBDirectIndirect();
	fNumDSBDirectQuasiDirect = fDamageCalculator->GetDSBDirectQuasiDirect();
	fNumDSBQuasiDirectQuasiDirect = fDamageCalculator->GetDSBQuasiDirectQuasiDirect();
	fNumDSBIndirectQuasiDirect = fDamageCalculator->GetDSBIndirectQuasiDirect();
	fNumBaseDamage = fDamageCalculator->GetBD();
	fNumBaseDamageDirect = fDamageCalculator->GetBDDirect();
	fNumBaseDamageQuasiDirect = fDamageCalculator->GetBDQuasiDirect();
	fNumBaseDamageIndirect = fDamageCalculator->GetBDIndirect();
	fNumSSBPlus = fDamageCalculator->GetSSBPlus();
	fNumDSBPlus = fDamageCalculator->GetDSBPlus();
	fNumDSBComplex = fDamageCalculator->GetDSBComplex();
	if (fChromosomeContents.size() > 0)
		CalculateYields();
	// Get all DSB positions to score foci
	if (fScoreFoci)
	{
		std::vector<G4ThreeVector> dsbPosInEvent = fDamageCalculator->GetDSB3DPositions();
		for (unsigned int i = 0; i < dsbPosInEvent.size(); i++)
			fDSBPositionsInRun.push_back(dsbPosInEvent[i]);
		std::vector<G4int> numFoci = fFociAnalyzer->GetNumberOfFoci(fDSBPositionsInRun);
		if (numFoci.size() >= 1) fNumFoci1 = numFoci[0];
		if (numFoci.size() >= 2) fNumFoci2 = numFoci[1];
		if (numFoci.size() >= 3) fNumFoci3 = numFoci[2];
		if (numFoci.size() >= 4) fNumFoci4 = numFoci[3];
		if (numFoci.size() >= 5) fNumFoci5 = numFoci[4];
	}

	return numberOfLesions;
}

void TsScoreDNADamageSBS::CalculateYields()
{
	G4double totalContentOfDNA = 0;
	for (unsigned int i = 0; i < fChromosomeContents.size(); i++)
		totalContentOfDNA += fChromosomeContents[i];

	fYBaseDam = (G4double)fNumBaseDamage / fDoseInThisExposure / totalContentOfDNA * 1e9;
	fYSB = (G4double)fNumSB / fDoseInThisExposure / totalContentOfDNA * 1e9;
	fYSSB = (G4double)fNumSSB / fDoseInThisExposure / totalContentOfDNA * 1e9;
	fYDSB = (G4double)fNumDSB / fDoseInThisExposure / totalContentOfDNA * 1e9;
	fYSSBPlus = (G4double)fNumSSBPlus / fDoseInThisExposure / totalContentOfDNA * 1e9;
	fYDSBPlus = (G4double)fNumDSBPlus / fDoseInThisExposure / totalContentOfDNA * 1e9;
	fYDSBComplex = (G4double)fNumDSBComplex / fDoseInThisExposure / totalContentOfDNA * 1e9;
}

G4double TsScoreDNADamageSBS::CalculateDoseInGray(G4double edep)
{
	// Calculating dose in this event. First we get volume of the component. Density is assumed to be water's
	G4double componentVolume = fComponent->GetEnvelopeLogicalVolume()->GetSolid()->GetCubicVolume() / pow(m, 3); //m3
	G4double waterDensity = 997; // kg/m3
	G4double componentMass = waterDensity * componentVolume; // kg
	return (1.6e-13 * edep / MeV) / componentMass; // This gives Gy
}

void TsScoreDNADamageSBS::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
	TsScoreDNADamageSBS* workerMTScorer = dynamic_cast<TsScoreDNADamageSBS*>(workerScorer);

	for(unsigned int i=0; i < workerMTScorer->fCollectionsOfHits.size(); i++)
		fCollectionsOfHits.push_back(workerMTScorer->fCollectionsOfHits[i]);
	workerMTScorer->fCollectionsOfHits.clear();

	for(unsigned int i=0; i < workerMTScorer->fEventsEdep.size(); i++)
		fEventsEdep.push_back(workerMTScorer->fEventsEdep[i]);
	workerMTScorer->fEventsEdep.clear();
}

// Default implementation (no more hierarchy levels, everything is chromosome 1)
std::pair<G4int, G4int> TsScoreDNADamageSBS::CalculateChromosomeAndBasePairID(std::vector<G4int> ids)
{
	std::pair<G4int, G4int> chromosomeAndBpID;
	chromosomeAndBpID.first = 1;
	chromosomeAndBpID.second = ids[0];
	return chromosomeAndBpID;
}
