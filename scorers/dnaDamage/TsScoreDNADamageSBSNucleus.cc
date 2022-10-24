// Scorer for DNADamageNucleusStepByStep
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

#include "TsScoreDNADamageSBSNucleus.hh"

TsScoreDNADamageSBSNucleus::TsScoreDNADamageSBSNucleus(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
		G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
		: TsScoreDNADamageSBS(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	// Add hierarchical levels for this geometry
	fHierarchicalLevels.push_back("Fiber");
	fHierarchicalLevels.push_back("Voxel");

	// Push back the total number of base pairs to the vector fChromosomeContents.
	fChromosomeContents = {252823200,  252823200,  248157000,  248157000,  204040200,  204040200,
                          195556200,  195556200,  184951200,  184951200,  174770400,  174770400,
                          162468600,  162468600,  149318400,  149318400,  143379600,  143379600,
                          138289200,  138289200,  137440800,  137440800,  135319800,  135319800,
                          116655000,  116655000,  108595200,  108595200,  102656400,  102656400,
                          90778800, 90778800,  80173800,  80173800,  77628600,  77628600,
                          65326800,  65326800,  63630000,  63630000,  47934600,  47934600,
                          50479800,  50479800,  58963800,  158226600};

    // To have this info updated in the SDD header, rerun:
	fDamageCalculator->OutputSDDHeader(fWriteMinimalSDDOutput, fPrimaryParticle, fMeanEnergy, fDosePerExposure, fChromosomeContents, fScoreIndirectDamage, fScoreOnBases,
		fAuthor, fSimulationDetails, fSourceDetails, fSourceType, fEnergyDist, fIrrTarget, fCellCycle, fDNAStructure, fInVitroOrInVivo, fProliferationStatus, fMicroenvironment, fTime, fAddInfo);

	// Geometry parameters specific for this scorer
	if (fPm->ParameterExists(GetFullParmName("GeometryInfo")))
		fGeometryInfoFile = fPm->GetStringParameter(GetFullParmName("GeometryInfo"));
	else
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("GeometryInfo") << " file has to be specified." << G4endl;
		fPm->AbortSession(1);
	}
	if (fPm->ParameterExists(GetFullParmName("FileVoxelCopyNumbersInsideNucleus")))
		fMapOfVoxelCopyNumbersInsideNucleusFile = fPm->GetStringParameter(GetFullParmName("FileVoxelCopyNumbersInsideNucleus"));
	else
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("FileVoxelCopyNumbersInsideNucleus") << " file has to be specified." << G4endl;
		fPm->AbortSession(1);
	}
	if (fPm->ParameterExists(GetFullParmName("MapChromosomeVoxelsFile")))
		fMapChromosomeVoxelFile = fPm->GetStringParameter(GetFullParmName("MapChromosomeVoxelsFile"));
	else
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("MapChromosomeVoxelsFile") << " file has to be specified." << G4endl;
		fPm->AbortSession(1);
	}
	// Reads files containing geometry and copy number maps
	ReadGeometryAndCopyNumberMaps();
}

TsScoreDNADamageSBSNucleus::~TsScoreDNADamageSBSNucleus() { }

// Calculations of the chromosomeID and basepairID for the specific geometry have to be done in the reimplementation of the UserHookForEndOfRun() method.
// This structure is mandatory:
void TsScoreDNADamageSBSNucleus::UserHookForEndOfRun()
{
	for (unsigned int iEvent = 0; iEvent < fCollectionsOfHits.size(); iEvent++)
	{
		for (unsigned int iHit = 0; iHit < fCollectionsOfHits[iEvent].size(); iHit++)
		{
			std::pair<G4int, G4int> chromosomeBpIDs = CalculateChromosomeAndBasePairID(fCollectionsOfHits[iEvent][iHit]->GetHierarchicalIDs());
			fCollectionsOfHits[iEvent][iHit]->SetChromosomeID(chromosomeBpIDs.first);
			fCollectionsOfHits[iEvent][iHit]->SetBasePairID(chromosomeBpIDs.second);
		}
	}
	TsScoreDNADamageSBS::UserHookForEndOfRun();
}

// Implementation of the actual method to get chromosome and bp ID in chromosome is mandatory
std::pair<G4int, G4int> TsScoreDNADamageSBSNucleus::CalculateChromosomeAndBasePairID(std::vector<G4int> hierarchicalIds)
{
	// Using hierarchical Ids
	G4int basePairIDInFiber = hierarchicalIds[0];
	G4int fiberIdInVoxel = hierarchicalIds[1];
	G4int voxelIdInWorld = hierarchicalIds[2];
	G4int chromosomeID, basePairIDInChromosome;
	if (fVoxelIDsInNucleus.size() > 0)
	{
		G4int voxelIDInNucleus = fVoxelIDsInNucleus[voxelIdInWorld];
		if (voxelIDInNucleus == -1 || voxelIdInWorld >= pow(fnVoxelIs3DRepeated, 3)) // Disregard voxels outside spherical shape of nucleus
		{
			chromosomeID = -1;
			basePairIDInChromosome = -1;
		}
		else
		{
			G4int voxelIDInChromosome = fVoxelIDs[voxelIDInNucleus];
			chromosomeID = fChromosomeIDForEachVoxel[voxelIDInNucleus];
			basePairIDInChromosome = voxelIDInChromosome * fBasePairsInAVoxel + fiberIdInVoxel * fBasePairsInAFiber + basePairIDInFiber;
		}
	}
	else
	{
		chromosomeID = 1;
		basePairIDInChromosome = basePairIDInFiber;
	}
	std::pair<G4int, G4int> ids;
	ids.first = chromosomeID;
	ids.second = basePairIDInChromosome;
	return ids;
}

// Auxiliary methods for Nucleus scorer
void TsScoreDNADamageSBSNucleus::ReadGeometryAndCopyNumberMaps()
{
	// Reads geometry information file. Contains the following components
	// fBasePairsInFiber, fBasePairsInVoxel (in header) 
	G4int nBasePairsInNucleus;
	G4int nVoxels;
	G4double voxelSize;
	//G4double minx, miny, minz;

	std::ifstream infile;
	infile.open(fGeometryInfoFile);
	G4String input;
	if (!infile.is_open())
	{
		G4cerr << "Scoring of DNA Damage in the nucleus. ERROR: Unable to open file: " << fGeometryInfoFile << G4endl;
		fPm->AbortSession(1);
	}
	else G4cout << "Reading geometry information file... " << G4endl;
	while (infile >> input)
	{
		if (input == "DNAContentperFiber(bp):")			infile >> fBasePairsInAFiber;
		if (input == "DNAContentperVoxel(bp):")			infile >> fBasePairsInAVoxel;
		if (input == "DNAContentinTotal(bp):")			infile >> nBasePairsInNucleus;
		if (input == "NumberofVoxels(subdomains)")		infile >> nVoxels;
		if (input == "Voxel3DrepeatTimes:")				infile >> fnVoxelIs3DRepeated;
		if (input == "VoxelSize(um):")					infile >> voxelSize;
		voxelSize *= um;
	//	minx = -voxelSize * fnVoxelIs3DRepeated / 2;
	//	miny = -voxelSize * fnVoxelIs3DRepeated / 2;
	//	minz = -voxelSize * fnVoxelIs3DRepeated / 2;
	}
	infile.close();

	// Reads map file connecting copy numbers. Contains the following elements
	G4int copyNoInWorld, copyNoInNucleus;

	infile.open(fMapOfVoxelCopyNumbersInsideNucleusFile);
	if (!infile.is_open())
	{
		G4cerr << "Scoring of DNA Damage in the nucleus. ERROR: Unable to open file: " << fMapOfVoxelCopyNumbersInsideNucleusFile << G4endl;
		fPm->AbortSession(1);
	}
	else G4cout << "Reading copy number map file... " << G4endl;
	G4String line;
	while (std::getline(infile, line))
	{
		std::stringstream stream(line.data());
		stream >> copyNoInWorld >> copyNoInNucleus;
		fVoxelIDsInWorld.push_back(copyNoInWorld);
		fVoxelIDsInNucleus.push_back(copyNoInNucleus);
	}
	infile.close();

	// Reads assignation of voxels to each chromosome. Contains the following elements
	G4int idx = 0, idy = 0, idz = 0, chId = 0, voxelId = 0;

	infile.open(fMapChromosomeVoxelFile);
	if (!infile.is_open())
	{
		G4cerr << "Scoring of DNA Damage in the nucleus. ERROR: Unable to open file: " << fMapChromosomeVoxelFile << G4endl;
		fPm->AbortSession(1);
	}
	else G4cout << "Reading chromosome-voxel map file... " << G4endl;
	while (std::getline(infile, line))
	{
		std::stringstream stream(line.data());
		stream >> copyNoInNucleus >> idx >> idy >> idz >> chId >> voxelId;
		fChromosomeIDForEachVoxel.push_back(chId);
		fVoxelIDs.push_back(voxelId);
	}
	infile.close();
}
