// Extra Class for TsNucleus
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
// Authors: Hongyu Zhu, Alexander Klapproth, Jan Schuemann, Alejandro Bertolet

#include "TsDefineDamage.hh"
#include "TsHitsRecord.hh"
#include "TsChromosome.hh"

#include "Randomize.hh"

using namespace std;

TsDefineDamage::TsDefineDamage()
{
	TsChromosome* HumanFibroblastNucleus = new TsChromosome();
	fChromosomeDNAContent = HumanFibroblastNucleus->GetDNAcontentofEachChromosome_BP(); // In unit of BP
	fTotalDNAContent 	  = HumanFibroblastNucleus->GetTotalDNAContent_MBP()*1e6; // In unit of BP
	delete HumanFibroblastNucleus;

	fScoreDirectDamages = true;
	fScoreIndirectDamages = true;
	fScoreOnBases = true;
	fScoreOnBackbones = true;
	fScoreOnHydrationShell = true;

	fEdep 			= 0;
	fLET  			= 0;
	fZetaBeta_sq 	= 0;
	fEventID		= 0;
	fNucleusMass	= 0;
	fPrimaryParticle = "proton";
	fDosePerExposure = 0;

	numSB			= 0;
	numSB_dir		= 0;
	numSB_indir		= 0;
	numSSB			= 0;
	numSSB_dir		= 0;
	numSSB_indir	= 0;
	numDSB			= 0;
	numDSB_dir		= 0;
	numDSB_indir	= 0;
	numDSB_hybrid	= 0;
	numSSB_P		= 0;
	numDSB_P		= 0;
	numDSB_PP		= 0;
	numBaseDam = 0;
	numBaseDam_dir = 0;
	numBaseDam_indir = 0;

	Excluded_numSB			= 0;
	Excluded_numSB_dir 		= 0;
	Excluded_numSB_indir	= 0;
	Excluded_numSSB			= 0;
	Excluded_numSSB_dir		= 0;
	Excluded_numSSB_indir	= 0;
	Excluded_numDSB			= 0;
	Excluded_numDSB_dir		= 0;
	Excluded_numDSB_indir	= 0;
	Excluded_numDSB_hybrid	= 0;
	Excluded_numSSB_P		= 0;
	Excluded_numDSB_P		= 0;
	Excluded_numDSB_PP		= 0;

	lastOutputEventID			= -1;
	lastOutputExposureID		= -1;
	lastDSBonlyOutputEventID	= -1;
	lastDSBonlyOutputExposureID	= -1;

}

TsDefineDamage::~TsDefineDamage() { }

// ******************************************************************************
//         NEW METHODS TO COMPUTE DAMAGE (Alejandro Bertolet, 03-29-21)
// ******************************************************************************
void TsDefineDamage::ComputeStrandBreaks(std::vector<TsHitsRecord*> Hits)
{
	// Clear all structures from previous computations
	fVEdepInDNA.clear();
	fVIndirectDamageInDNA.clear();
	fVIonizationInHydrationShell.clear();
	fDamagePositions.clear();
	fDamage.clear();
	fDSB.clear();
	fSSB.clear();
	fBaseDamage.clear();

	// Get info from these hits
	for (G4int i = 0; i < Hits.size(); i++)
	{
		G4int basePairID = Hits[i]->GetBasePairID();
		G4double edep = Hits[i]->GetEdep();
		G4int chromosomeID = Hits[i]->GetChromosomeID();
		G4int iElement;
		if (strstr(Hits[i]->GetVolumeName(), "Base1") != NULL) iElement = 0;
		else if (strstr(Hits[i]->GetVolumeName(), "Base2") != NULL) iElement = 1;
		else if (strstr(Hits[i]->GetVolumeName(), "Backbone1") != NULL) iElement = 2;
		else if (strstr(Hits[i]->GetVolumeName(), "Backbone2") != NULL) iElement = 3;
		else if (strstr(Hits[i]->GetVolumeName(), "HydrationShell1") != NULL) iElement = 4;
		else if (strstr(Hits[i]->GetVolumeName(), "HydrationShell2") != NULL) iElement = 5;
		else iElement = -1;
		if (iElement >= 0)
		{
			if (Hits[i]->GetIsDirectDamage()) fVEdepInDNA[chromosomeID][basePairID][iElement] += edep;
			else fVIndirectDamageInDNA[chromosomeID][basePairID][iElement] = true;
			if (strstr(Hits[i]->GetProcess(), "Ionisation"))
			{
				fVIonizationInHydrationShell[chromosomeID][basePairID][iElement] = true;
				fDamagePositions[chromosomeID][basePairID][iElement-2] = Hits[i]->GetPosition();
			}
			fDamagePositions[chromosomeID][basePairID][iElement] = Hits[i]->GetPosition();
		}
	}
	// Classify damages into structures
	for (auto& chromosome : fVEdepInDNA)
	{
		G4int iChr = chromosome.first;
		for (auto& pairbpElem : fVEdepInDNA[iChr])
		{
			G4int ibp = pairbpElem.first;
			for (auto& pairElemEnergy : fVEdepInDNA[iChr][ibp])
			{
				G4int iElem = pairElemEnergy.first;
				if (iElem >= 0 && (((iElem == 0 || iElem == 1) && fScoreOnBases) || ((iElem == 2 || iElem == 3) && fScoreOnBackbones)))
				{
					if (CauseDirectDamage(fVEdepInDNA[iChr][ibp][iElem]) && fScoreDirectDamages)
						fDamage[iChr][ibp][iElem] = 1;
					else
						fDamage[iChr][ibp][iElem] = -1;
				}
			}
		}
	}
	for (auto& chromosome : fVIndirectDamageInDNA)
	{
		G4int iChr = chromosome.first;
		for (auto& pairbpElem : fVIndirectDamageInDNA[iChr])
		{
			G4int ibp = pairbpElem.first;
			for (auto& pairElemEnergy : fVIndirectDamageInDNA[iChr][ibp])
			{
				G4int iElem = pairElemEnergy.first;
				if (iElem >= 0 && (((iElem == 0 || iElem == 1) && fScoreOnBases) || ((iElem == 2 || iElem == 3) && fScoreOnBackbones)) && fScoreIndirectDamages)
				{
					if (fDamage[iChr][ibp][iElem] == 1)
						fDamage[iChr][ibp][iElem] = 3; // Multiple damage if previous direct
					else
						fDamage[iChr][ibp][iElem] = 2; // Indirect if not
				}
			}
		}
	}
	for (auto& chromosome : fVIonizationInHydrationShell)
	{
		G4int iChr = chromosome.first;
		for (auto& pairbpElem : fVIonizationInHydrationShell[iChr])
		{
			G4int ibp = pairbpElem.first;
			for (auto& pairElemEnergy : fVIonizationInHydrationShell[iChr][ibp])
			{
				G4int iElem = pairElemEnergy.first;
				if (iElem >= 0 && (((iElem == 4 || iElem == 5) && fScoreOnHydrationShell)) && fScoreDirectDamages)
				{
					// Score damage in backbone due to ionization of the corresponding hydration shell
					if (fDamage[iChr][ibp][iElem - 2] < 2)
						fDamage[iChr][ibp][iElem - 2] = 1; // Direct if no previous indirect
					else
						fDamage[iChr][ibp][iElem - 2] = 3; // Multiple damage if previous indirect
				}
			}
		}
	}
	if (fWriteCSV) WriteDNADamageCSV();

	// Look for DSB first
	std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs;
	for (auto& chromosome : fDamage)
	{
		G4int iChr = chromosome.first;
		for (auto& pairbpElem : fDamage[iChr])
		{
			G4int ibp = pairbpElem.first;
			if (fDamage[iChr][ibp][2] > 0 && fDSB[iChr][ibp][1] == 0) // Check if there is damage in backbone 1 (fDamage index 2) and if there is not already a DSB identified in strand 1 (fDSB index 1)
			{
				G4bool dsbFound = false;
				for (G4int inc2 = 0; inc2 < fDSBSeparation; inc2++)
				{
					if (fDamage[iChr][ibp+inc2][3] > 0 && fDSB[iChr][ibp+inc2][2] == 0) // Go from ibp to ibp+10 to look for damage in backbone 2 (fDamage index 3) and not already a DSB identified in strand 2 (fDSB index 2)
					{
						G4int typeDamage2 = 2;
						if (fDamage[iChr][ibp+inc2][3] == 1) typeDamage2 = 1; // Check if damage in backbone 2 is direct (1) or indirect (2)
						// Once identified damage in strand 2 (and ensuring damage in strand 1), look for closest damage in strand 1 (just in case there are multiple damages)
						G4int closestS1pos = ibp;
						G4int typeDamage1 = 2;
						for (G4int inc1 = 0; inc1 < fDSBSeparation; inc1++)
						{
							if (fDamage[iChr][ibp+inc2+inc1][2] > 0 && fDSB[iChr][ibp+inc2+inc1][1] == 0) // Go from ibp+inc2 increasing 1 bp....
							{
								closestS1pos = ibp + inc2 + inc1;
								if (fDamage[iChr][ibp+inc2+inc1][2] == 1) typeDamage1 = 1;
								break;
							}
							if (ibp+inc2-inc1 >= 0 && inc1 > 0 && fDamage[iChr][ibp+inc2-inc1][2] > 0 && fDSB[iChr][ibp+inc2-inc1][1] == 0) // and decreasing 1 bp alternatively, breaking when closest damage is found
							{
								closestS1pos = ibp + inc2 - inc1;
								if (fDamage[iChr][ibp+inc2-inc1][2] == 1) typeDamage1 = 1;
								break;
							}
						}
						fDSB[iChr][closestS1pos][1] += typeDamage1;
						fDSB[iChr][ibp + inc2][2] += typeDamage2;
						std::pair<G4int, G4int> pos(closestS1pos, ibp+inc2);
						DSBPairs[iChr][pos] = typeDamage1 + typeDamage2; // 2->Direct, 3->Hybrid, 4->Indirect
						//G4cout << "DSB found at " << closestS1pos << ", " << ibp+inc2 << " of type " << DSBPairs[pos] << G4endl;
						dsbFound = true;
					}
					if (ibp-inc2 >= 0 && inc2 > 0 && fDamage[iChr][ibp-inc2][3] > 0 && fDSB[iChr][ibp-inc2][2] == 0) // Go from ibp to ibp-10 to look for damage in backbone 2 (alternatively, if dsb is found then break)
					{
						G4int typeDamage2 = 2;
						if (fDamage[iChr][ibp-inc2][3] == 1) typeDamage2 = 1; // Check if damage in backbone 2 is direct (1) or indirect (2)
						// Look for closest damage in strand 1
						G4int closestS1pos = ibp;
						G4int typeDamage1 = 2;
						for (G4int inc1 = 0; inc1 < fDSBSeparation; inc1++)
						{
							if (fDamage[iChr][ibp-inc2+inc1][2] > 0 && fDSB[iChr][ibp-inc2+inc1][1] == 0)
							{
								closestS1pos = ibp - inc2 + inc1;
								if (fDamage[iChr][ibp-inc2+inc1][2] == 1) typeDamage1 = 1;
								break;
							}
							if (ibp-inc2-inc1 >= 0 && inc1 > 0 && fDamage[iChr][ibp-inc2-inc1][2] > 0 && fDSB[iChr][ibp-inc2-inc1][1] == 0)
							{
								closestS1pos = ibp - inc2 - inc1;
								if (fDamage[iChr][ibp-inc2-inc1][2] == 1) typeDamage1 = 1;
								break;
							}
						}
						fDSB[iChr][closestS1pos][1] += typeDamage1;
						fDSB[iChr][ibp - inc2][2] += typeDamage2;
						std::pair<G4int, G4int> pos(closestS1pos, ibp-inc2);
						DSBPairs[iChr][pos] = typeDamage1 + typeDamage2; // 2->Direct, 3->Hybrid, 4->Indirect
						//G4cout << "DSB found at " << closestS1pos << ", " << ibp-inc2 << " of type " << DSBPairs[pos] << G4endl;
						dsbFound = true;
					}
					if (dsbFound) break;
				}
			}
		}
	}
	// Look for SSBs excluding damage already classified as DSB
	for (auto& chromosome : fDamage)
	{
		G4int iChr = chromosome.first;
		for (auto& pairbpElem : fDamage[iChr])
		{
			G4int ibp = pairbpElem.first;
			if ((fDamage[iChr][ibp][2] > 0 && !(fDSB[iChr][ibp][1] > 0)) || (fDamage[iChr][ibp][3] > 0 && !(fDSB[iChr][ibp][2] > 0)))
			{
				G4bool s1 = false, s2 = false; G4int strand;
				if (fDamage[iChr][ibp][2] > 0) strand = 1;
				if (fDamage[iChr][ibp][3] > 0) strand = 2;
				fSSB[iChr][ibp][strand] = fDamage[iChr][ibp][strand+1];
			}
		}
	}
	// Base damages
	if (fScoreOnBases)
	{
		for (auto& chromosome : fDamage)
		{
			G4int iChr = chromosome.first;
			for (auto& pairbpElem : fDamage[iChr])
			{
				G4int ibp = pairbpElem.first;
				if (fDamage[iChr][ibp][0] > 0 || fDamage[iChr][ibp][1] > 0)
				{
					G4int strand;
					if (fDamage[iChr][ibp][0] > 0) strand = 1;
					if (fDamage[iChr][ibp][1] > 0) strand = 2;
					fBaseDamage[iChr][ibp][strand] = fDamage[iChr][ibp][strand-1];
				}
			}
		}
	}
	if (fExcludeShortFragment) ExcludeShortDNAFragments(DSBPairs);
	QuantifyDamage(DSBPairs, fSSB, fBaseDamage);
}

void TsDefineDamage::QuantifyDamage(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs, std::map<G4int, std::map<G4int, std::map<G4int, G4int>>> SSB, std::map<G4int,
									std::map<G4int, std::map<G4int, G4int>>> BaseDamage)
{
	// Reinitialize all quantities from previous quantifications
	numDSB = 0; numDSB_dir = 0; numDSB_hybrid = 0; numDSB_indir = 0;
	numSSB = 0; numSSB_dir = 0; numSSB_indir = 0;
	numSB = 0; numSB_dir = 0; numSB_indir = 0;
	numBaseDam = 0; numBaseDam_dir = 0; numBaseDam_indir = 0;
	for (auto& chromosome : DSBPairs)
	{
		G4int iChr = chromosome.first;
		for (auto& pairPosType : DSBPairs[iChr])
		{
			numDSB++;
			if (pairPosType.second == 2)
				numDSB_dir++;
			else if (pairPosType.second == 3)
				numDSB_hybrid++;
			else if (pairPosType.second == 4)
				numDSB_indir++;
		}
	}
	for (auto& chromosome : SSB)
	{
		G4int iChr = chromosome.first;
		for (auto& bp_strand_type : SSB[iChr])
		{
			G4int bpPos = bp_strand_type.first;
			for (auto& strand_type : SSB[iChr][bpPos])
			{
				numSSB++;
				if (strand_type.second == 1)
					numSSB_dir++;
				else if (strand_type.second == 2)
					numSSB_indir++;
			}
		}
	}
	for (auto& chromosome : BaseDamage)
	{
		G4int iChr = chromosome.first;
		for (auto& bp_strand_type : BaseDamage[iChr])
		{
			G4int bpPos = bp_strand_type.first;
			for (auto& strand_type : BaseDamage[iChr][bpPos])
			{
				numBaseDam++;
				if (strand_type.second == 1)
					numBaseDam_dir++;
				else if (strand_type.second == 2)
					numBaseDam_indir++;
			}
		}
	}
	numSB = numSSB + numDSB * 2;
	numSB_dir = numSSB_dir + numDSB_dir * 2 + numDSB_hybrid;
	numSB_indir = numSSB_indir + numDSB_indir * 2 + numDSB_hybrid;
	//G4cout << "numBaseDam: " << numBaseDam << " - numBaseDam_dir: " << numBaseDam_dir << " - numBaseDam_indir: " << numBaseDam_indir << G4endl;
	//G4cout << "numSB: " << numSB << " - numSB_dir: " << numSB_dir << " - numSB_indir: " << numSB_indir << G4endl;
	//G4cout << "numSSB: " << numSSB << " - numSSB_dir: " << numSSB_dir << " - numSSB_indir: " << numSSB_indir << G4endl;
	//G4cout << "numDSB: " << numDSB << " - numDSB_dir: " << numDSB_dir << " - numDSB_indir: " << numDSB_indir <<  " - numDSB_hybrid: " << numDSB_hybrid << G4endl;
}

std::map<G4int, std::vector<G4int>> TsDefineDamage::IdentifyDamageBlocks()
{
	numSSBPlus = 0;
	numDSBPlus = 0;
	numMoreComplex = 0;
	// Look for 10-bp (or specified by fDSBSeparation) length blocks. Returns a list of the initial base pair number for each block
	// Assign 0.001 to each base damage, 0.1 to each SSB, 10 to each DSB (5 per break)
	std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> damageAssignation;

	for (auto& chromosome : fBaseDamage)
	{
		G4int iChr = chromosome.first;
		for (auto& bp_strand_type : fBaseDamage[iChr])
		{
			G4int bpPos = bp_strand_type.first;
			for (auto& strand_type : fBaseDamage[iChr][bpPos])
			{
				G4int strand = strand_type.first;
				if (fBaseDamage[iChr][bpPos][strand] > 0)
					damageAssignation[iChr][bpPos][strand] += 0.001;
			}
		}
	}
	for (auto& chromosome : fSSB)
	{
		G4int iChr = chromosome.first;
		for (auto& bp_strand_type : fSSB[iChr])
		{
			G4int bpPos = bp_strand_type.first;
			for (auto& strand_type : fSSB[iChr][bpPos])
			{
				G4int strand = strand_type.first;
				if (fSSB[iChr][bpPos][strand] > 0)
					damageAssignation[iChr][bpPos][strand] += 0.1;
			}
		}
	}
	for (auto& chromosome : fDSB)
	{
		G4int iChr = chromosome.first;
		for (auto& bp_strand_type : fDSB[iChr])
		{
			G4int bpPos = bp_strand_type.first;
			for (auto& strand_type : fDSB[iChr][bpPos])
			{
				G4int strand = strand_type.first;
				if (fDSB[iChr][bpPos][strand] > 0)
					damageAssignation[iChr][bpPos][strand] += 5;
			}
		}
	}
	// Generate 10-length blocks of damage
	std::map<G4int, std::vector<G4int>> chromBpInitiatingDamageBlocks;
	for (auto& chromosome : damageAssignation)
	{
		G4int iChr = chromosome.first;
		std::vector<G4double> linearAccumDamages;
		std::vector<G4int> initialBps;
		std::vector<G4int> bpInitiatingDamageBlocks;
		G4double maxDamage = 10 * fDSBSeparation; // Maximum possible damage
		G4double minDamage = 1e-8;
		if (fOnlyIncludeDSBinSDD) minDamage = 9.999999; // Only include damages with at least one DSB (assigned 10)
		while (maxDamage > minDamage)
		{
			for (auto& bp_strand_type : damageAssignation[iChr])
			{
				G4int bpPos = bp_strand_type.first;
				G4double sumDamage = 0;
				for (G4int inc = 0; inc < fDSBSeparation; inc++)
				{
					if (damageAssignation[iChr].find(bpPos + inc) != damageAssignation[iChr].end())
						sumDamage += damageAssignation[iChr][bpPos + inc][1] + damageAssignation[iChr][bpPos + inc][2];
				}
				linearAccumDamages.push_back(sumDamage);
				initialBps.push_back(bpPos);
			}
			G4double m = 0; G4int index = 0;
			for (G4int i = 0; i < linearAccumDamages.size(); i++)
			{
				if (m < linearAccumDamages[i])
				{
					m = linearAccumDamages[i];
					index = initialBps[i];
				}
			}
			maxDamage = m;
			if (maxDamage >= 0.2 && maxDamage < 10) numSSBPlus++;
			if (maxDamage >= 20 && maxDamage < 30) numDSBPlus++;
			if (maxDamage >= 30) numMoreComplex++;
			if (maxDamage < minDamage) { break; }
			bpInitiatingDamageBlocks.push_back(index);
			// Exclude these damages for next steps
			for (G4int inc = 0; inc < fDSBSeparation; inc++)
			{
				if (damageAssignation[iChr].find(index + inc) != damageAssignation[iChr].end())
				{
					damageAssignation[iChr][index + inc][1] = 0;
					damageAssignation[iChr][index + inc][2] = 0;
				}
			}
			linearAccumDamages.clear();
			initialBps.clear();
		}
		chromBpInitiatingDamageBlocks[iChr] = bpInitiatingDamageBlocks;
	}
	return chromBpInitiatingDamageBlocks;
}

void TsDefineDamage::OutputSDDHeader(G4String author, G4String simDetails, G4String sourceDetails, G4int sourceType, G4String energyDist, G4String irrTarget,
		G4String cellCycle, G4String DNAStructure, G4int inVitroOrInVivo, G4String prolStatus, G4String microenv, G4double time, G4String addInfo)
{
	// Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al.
	// A New Standard DNA Damage (SDD) Data Format. Radiat. Res. 191, 76-92 (2019)
	G4String dataEntries = "1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0";
	if (fMinimalSDD) dataEntries = "1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0";

	G4String SDDVersion = "v1.0";
	G4String incParticle;
	if (fPrimaryParticle == "proton") incParticle = "2212";
	else if (fPrimaryParticle == "alpha") incParticle = "100002004";
	G4String meanEnergy = (G4String)to_string(fMeanEnergy) + " MeV";

	G4String scoringIndirect = "0";
	if (fScoreIndirectDamages) scoringIndirect = "1";
	G4String scoringOnBases = "-1";
	if (fScoreOnBases) scoringOnBases = (G4String)std::to_string(fDSBSeparation);

	std::ofstream outFile;
	outFile.open(fOutputFile + "_sdd.txt");
	outFile << "SDD Version, SDD" << SDDVersion << ";\n";
	outFile << "Software, TOPAS-nBio;\n";
	outFile << "Author, " << author << ";\n";
	outFile << "Simulation Details, " << simDetails << ";\n";
	outFile << "Source, " << sourceDetails << ";\n";
	outFile << "Source type, " << to_string(sourceType) << ";\n";
	outFile << "Incident particles, " << incParticle << ";\n";
	outFile << "Mean particle energy, " << meanEnergy << ";\n";
	outFile << "Energy distribution, " << energyDist << ";\n";
	outFile << "Particle fraction, 1.0;\n";
	outFile << "Dose or fluence, 1, " << fDosePerExposure << ";\n";
	outFile << "Irradiation target, " << irrTarget << ";\n";
	outFile << "Volumes, 0,5,5,5,0,0,0,1,4.65,4.65,4.65,0,0,0;\n";
	outFile << "Chromosome sizes, " << fChromosomeDNAContent.size() << ",";
	for (G4int i = 0; i < fChromosomeDNAContent.size(); i++)
	{
		outFile << (G4double)fChromosomeDNAContent[i] / 1e6;
		if (i < fChromosomeDNAContent.size() - 1) outFile << ",";
		else outFile << ";\n";
	}
	outFile << "DNA Density, " << fTotalDNAContent / 1e6 / 421.15 << ";\n";
	outFile << "Cell Cycle Phase, " << cellCycle << ";\n";
	outFile << "DNA Structure, " << DNAStructure << ";\n";
	outFile << "In vitro / in vivo, " << to_string(inVitroOrInVivo) << ";\n";
	outFile << "Proliferation status, " << prolStatus << ";\n";
	outFile << "Microenvironment, " << microenv << ";\n";
	outFile << "Damage definition, " << scoringIndirect << ", 0, " << fDSBSeparation << ", " << scoringOnBases << ", " << fDamageThreshold/eV;
	if (fUseLinearProbabilityThreshold) outFile << ", " << fLinearProbabilityUpperLimit;
	outFile << ";\n";
	outFile << "Time, " << to_string(time) << ";\n";
	outFile << "Damage and primary count, ;\n";
	outFile << "Data entries, " << dataEntries << ";\n";
	outFile << "Additional information, " << addInfo << ";\n";
	outFile << "***EndOfHeader***;\n\n";
	outFile.close();
}

void TsDefineDamage::UpdateDamageAndPrimaryCount(G4int nLesions, G4int nPrimary)
{
	std::ifstream outFile;
	outFile.open(fOutputFile + "_sdd.txt");
	G4String line;
	vector<G4String> lines;

	while (getline(outFile, line))
	{
		if (line.find("Damage and primary count") != std::string::npos)
			lines.push_back("Damage and primary count, " + to_string(nLesions) + ", " + to_string(nPrimary) + ";");
		else
			lines.push_back(line);
	}
	outFile.close();
	std::ofstream file;
	file.open(fOutputFile + "_sdd.txt", std::ios::out | std::ios::trunc);
	for (G4int i = 0; i < lines.size(); i++)
		file << lines[i] << "\n";
	file.close();
}

G4int TsDefineDamage::OutputSDDFile(std::map<G4int, std::vector<G4int>> damageSites, G4int eventID, G4int exposureID)
{
	// Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al.
	// A New Standard DNA Damage (SDD) Data Format. Radiat. Res. 191, 76-92 (2019)

	G4int lastExposureID = -1000;
	G4int lastEventID = -1000;

	G4int distinctLesions = 0;

	std::ofstream outFile;
	outFile.open(fOutputFile + "_sdd.txt", std::ios::out|std::ios::app);
	for (auto& chromosome : damageSites)
	{
		G4int iChr = chromosome.first;
		for (G4int i = 0; i < damageSites[iChr].size(); i++)
		{
			G4int initialBp = damageSites[iChr][i];
			// Determining number and type of damages in block
			G4int numDirDamage = 0;
			G4int numIndirDamage = 0;
			G4int numBaseDamagesinBlock = 0;
			G4int numSBinBlock = 0;
			G4int DSBinBlock = 0;
			for (G4int j = 0; j < fDSBSeparation; j++)
			{
				if (fSSB[iChr][initialBp + j][1] == 1) { numDirDamage++; numSBinBlock++; }
				if (fSSB[iChr][initialBp + j][2] == 1) { numDirDamage++; numSBinBlock++; }
				if (fSSB[iChr][initialBp + j][1] >= 2) { numIndirDamage++; numSBinBlock++; }
				if (fSSB[iChr][initialBp + j][2] >= 2) { numIndirDamage++; numSBinBlock++; }
				if (fDSB[iChr][initialBp + j][1] == 1) { numDirDamage++; numSBinBlock++; DSBinBlock = 1; }
				if (fDSB[iChr][initialBp + j][2] == 1) { numDirDamage++; numSBinBlock++; DSBinBlock = 1; }
				if (fDSB[iChr][initialBp + j][1] >= 2) { numIndirDamage++; numSBinBlock++; DSBinBlock = 1; }
				if (fDSB[iChr][initialBp + j][2] >= 2) { numIndirDamage++; numSBinBlock++; DSBinBlock = 1; }
				if (fScoreOnBases)
				{
					if (fBaseDamage[iChr][initialBp + j][1] == 1) { numDirDamage++; numBaseDamagesinBlock++; }
					if (fBaseDamage[iChr][initialBp + j][2] == 1) { numDirDamage++; numBaseDamagesinBlock++; }
					if (fBaseDamage[iChr][initialBp + j][1] >= 2) { numIndirDamage++; numBaseDamagesinBlock++; }
					if (fBaseDamage[iChr][initialBp + j][2] >= 2) { numIndirDamage++; numBaseDamagesinBlock++; }
				}
			}

			distinctLesions++;

			// Field 1
			// -------
			// Determining exposure status
			G4int newExposureFlag = 0;
			if (lastEventID != eventID)
			{
				lastEventID = eventID;
				newExposureFlag = 1;
			}
			if (lastExposureID != exposureID)
			{
				lastExposureID = exposureID;
				newExposureFlag = 2;
			}
			outFile << newExposureFlag << ", " << eventID << "; ";

			// Field 2
			// --------
			// Add damage positions (backbones)
			std::vector<G4ThreeVector> damagePositions;
			for (G4int j = 0; j < fDSBSeparation; j++)
			{
				if (fDamage[iChr][initialBp + j][0] > 0) damagePositions.push_back(fDamagePositions[iChr][initialBp + j][0]);
				if (fDamage[iChr][initialBp + j][1] > 0) damagePositions.push_back(fDamagePositions[iChr][initialBp + j][1]);
				if (fDamage[iChr][initialBp + j][2] > 0) damagePositions.push_back(fDamagePositions[iChr][initialBp + j][2]);
				if (fDamage[iChr][initialBp + j][3] > 0) damagePositions.push_back(fDamagePositions[iChr][initialBp + j][3]);
			}
			std::vector<G4ThreeVector> damageCenterMaxMin = GetDamageCenterAndBoundaries(damagePositions);
			G4ThreeVector center = damageCenterMaxMin[0];
			outFile << center.x() / um << ", " << center.y() / um << ", " << center.z() / um;
			if (!fMinimalSDD)
			{
				 G4ThreeVector max = damageCenterMaxMin[1];
				 G4ThreeVector min = damageCenterMaxMin[2];
				 outFile << " / " << max.x() / um << ", " << max.y() / um << ", " << max.z() / um;
				 outFile << " / " << min.x() / um << ", " << min.y() / um << ", " << min.z() / um;
			}
			outFile << "; ";

			if (!fMinimalSDD)
			{
				// Field 3
				// -------
				// Nucleus considered as heterochromatin
				outFile << "1, ";
				// Chromatid number set to 1 (unduplicated), chromosome arm set to 0 (short)
				outFile << iChr - 1 << ", 1, 0; ";

				// Field 4
				// -------
				G4int chromosomeLength = fChromosomeDNAContent[iChr - 1];
				G4double damageChromPos = (G4double) initialBp / chromosomeLength;
				outFile.precision(12);
				outFile << fixed << damageChromPos << "; ";
				outFile.precision(6);

				// Field 5
				// -------
				G4int typeDamage = 0;
				if (numDirDamage == 0 && numIndirDamage > 0) typeDamage = 1;
				if (numDirDamage > 0 && numIndirDamage > 0) typeDamage = 2;
				outFile << typeDamage << ", " << numDirDamage << ", " << numIndirDamage << "; ";
			}


			// Field 6
			// -------
			outFile << numBaseDamagesinBlock << ", " << numSBinBlock << ", " << DSBinBlock << "; ";
			if (!fMinimalSDD)
			{
				G4String damageSpec;
				// Field 7
				// -------
				for (G4int j = 0; j < fDSBSeparation; j++)
				{
					// Starting with bases of strand 1
					if (fDamage[iChr][initialBp + j][0] != 0)
					{
						G4int typeDamage = fDamage[iChr][initialBp + j][0];
						if (typeDamage == -1) typeDamage = 0;
						damageSpec += "1, " + (G4String)to_string(j+1) + ", " + (G4String)to_string(typeDamage) + " / ";
					}
				}
				for (G4int j = 0; j < fDSBSeparation; j++)
				{
					// Backbones of strand 1
					if (fDamage[iChr][initialBp + j][2] != 0)
					{
						G4int typeDamage = fDamage[iChr][initialBp + j][2];
						if (typeDamage == -1) typeDamage = 0;
						damageSpec += "2, " + (G4String)to_string(j+1) + ", " + (G4String)to_string(typeDamage) + " / ";
					}
				}
				for (G4int j = 0; j < fDSBSeparation; j++)
				{
					// Backbones of strand 2
					if (fDamage[iChr][initialBp + j][3] != 0)
					{
						G4int typeDamage = fDamage[iChr][initialBp + j][3];
						if (typeDamage == -1) typeDamage = 0;
						damageSpec += "3, " + (G4String)to_string(j+1) + ", " + (G4String)to_string(typeDamage) + " / ";
					}
				}
				for (G4int j = 0; j < fDSBSeparation; j++)
				{
					// Bases of strand 2
					if (fDamage[iChr][initialBp + j][1] != 0)
					{
						G4int typeDamage = fDamage[iChr][initialBp + j][1];
						if (typeDamage == -1) typeDamage = 0;
						damageSpec += "4, " + (G4String)to_string(j+1) + ", " + (G4String)to_string(typeDamage) + " / ";
					}
				}
				damageSpec = damageSpec.substr(0, damageSpec.size()-2);
				outFile << damageSpec << ";";
			}
			outFile << "\n";
		}
	}

	outFile.close();
	return distinctLesions;
}

std::vector<G4ThreeVector> TsDefineDamage::GetDamageCenterAndBoundaries(std::vector<G4ThreeVector> positions)
{
	G4double xcenter, ycenter, zcenter;
	G4double xmin = positions[0].x();
	G4double ymin = positions[0].y();
	G4double zmin = positions[0].z();
	G4double xmax = positions[0].x();
	G4double ymax = positions[0].y();
	G4double zmax = positions[0].z();
	for (auto& p : positions)
	{
		if (p.x() > xmax) xmax = p.x();
		if (p.x() < xmin) xmin = p.x();
		if (p.y() > ymax) ymax = p.y();
		if (p.y() < ymin) ymin = p.y();
		if (p.z() > zmax) zmax = p.z();
		if (p.z() < zmin) zmin = p.z();
	}
	xcenter = (xmax + xmin) / 2.0;
	ycenter = (ymax + ymin) / 2.0;
	zcenter = (zmax + zmin) / 2.0;
	std::vector<G4ThreeVector> centerMaxMin;
	centerMaxMin.push_back(G4ThreeVector(xcenter, ycenter, zcenter));
	centerMaxMin.push_back(G4ThreeVector(xmax, ymax, zmax));
	centerMaxMin.push_back(G4ThreeVector(xmin, ymin, zmin));
	return centerMaxMin;
}

G4bool TsDefineDamage::CauseDirectDamage(G4double edep)
{
	if (!fUseLinearProbabilityThreshold)
	{
		if (edep > fDamageThreshold)
			return true;
	}
	else
	{
		if (edep >= fLinearProbabilityUpperLimit)
			return true;
		else if (edep <= fLinearProbabilityLowerLimit)
			return false;
		else
		{
			G4double samplePro = G4UniformRand();
			G4double acceptPro = (edep - fLinearProbabilityLowerLimit) / (fLinearProbabilityUpperLimit - fLinearProbabilityLowerLimit);
			if (samplePro < acceptPro)
				return true;
		}
	}
	return false;
}

void TsDefineDamage::WriteDNADamageCSV()
{
	G4String damFileName = fOutputFile + "_full.csv";
	fstream fout;
	// Creates a new file, appending results
	if (fOutputMode == "overwrite")
		fout.open(damFileName, ios::out | ios::trunc);
	else
		fout.open(damFileName, ios::out | ios::app);
	// Insert headers
	fout << "New event\n";
	fout << "Chromosome,Bp,Backbone1,Base1,Base2,Backbone2\n";

	for (auto& chromosome : fDamage)
	{
		G4int iChr = chromosome.first;
		for (auto& pairbpElem : fDamage[iChr])
		{
			G4String base1 = "0", base2 = "0", bb1 = "0", bb2 = "0";
			G4int ibp = pairbpElem.first;
			if (fDamage[iChr][ibp][0] == 1) base1 = "D";
			if (fDamage[iChr][ibp][0] == 2) base1 = "I";
			if (fDamage[iChr][ibp][0] == 3) base1 = "M";
			if (fDamage[iChr][ibp][0] == -1) base1 = "*";
			if (fDamage[iChr][ibp][1] == 1) base2 = "D";
			if (fDamage[iChr][ibp][1] == 2) base2 = "I";
			if (fDamage[iChr][ibp][1] == 3) base2 = "M";
			if (fDamage[iChr][ibp][1] == -1) base2 = "*";
			if (fDamage[iChr][ibp][2] == 1) bb1 = "D";
			if (fDamage[iChr][ibp][2] == 2) bb1 = "I";
			if (fDamage[iChr][ibp][2] == 3) bb1 = "M";
			if (fDamage[iChr][ibp][2] == -1) bb1 = "*";
			if (fDamage[iChr][ibp][3] == 1) bb2 = "D";
			if (fDamage[iChr][ibp][3] == 2) bb2 = "I";
			if (fDamage[iChr][ibp][3] == 3) bb2 = "M";
			if (fDamage[iChr][ibp][3] == -1) bb2 = "*";
			fout << iChr << ", " << ibp << "," << bb1 << "," << base1 << "," << base2 << "," << bb2 << "\n";
		}
	}
	fout.close();
}

void TsDefineDamage::ExcludeShortDNAFragments(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs)
{
	// Split chromosomes according to a DSB position
	TsChromosome* processChromosome = new TsChromosome();
	std::vector<std::vector<G4int>> splitedChromosome = processChromosome->SplitChromosome(DSBPairs);
	delete processChromosome;

	for (auto& chromosome : fDamage)
	{
		G4int iChr = chromosome.first;
		// Get maximum bp possible to avoid segmentation fault when not using actual chromosomes
		G4int maxbp = 0;
		for (auto& bp_strand : fDamage[iChr])
		{
			if (bp_strand.first > maxbp) maxbp = bp_strand.first;
		}
		G4int inviableFragments = 0;
		std::vector<G4int> allFreeEndsOfThisChromosome = splitedChromosome[iChr - 1];
		for (G4int i = 1; i < allFreeEndsOfThisChromosome.size(); i++)
		{
			G4int fragmentSize = min(maxbp, allFreeEndsOfThisChromosome[i]) - allFreeEndsOfThisChromosome[i-1];
			G4bool isInviableFragment = (fragmentSize < fLowerFragmentDetectionThreshold || fragmentSize > fUpperFragmentDetectionThreshold);
			if (isInviableFragment)
			{
				inviableFragments++;
				// Exclude damage to bases, SSB and DSB on the short fragments
				for (G4int ibp = allFreeEndsOfThisChromosome[i-1]; ibp < allFreeEndsOfThisChromosome[i]; ibp++)
				{
					if (ibp <= maxbp)
					{
						if (fBaseDamage[iChr][ibp][1] == 1) { fBaseDamage[iChr][ibp][1] = 0; fDamage[iChr][ibp][0] = 0; Excluded_baseDam++; Excluded_baseDam_dir++;}
						if (fBaseDamage[iChr][ibp][1] >= 2) { fBaseDamage[iChr][ibp][1] = 0; fDamage[iChr][ibp][0] = 0; Excluded_baseDam++; Excluded_baseDam_indir++;}
						if (fBaseDamage[iChr][ibp][2] == 1) { fBaseDamage[iChr][ibp][2] = 0; fDamage[iChr][ibp][1] = 0; Excluded_baseDam++; Excluded_baseDam_dir++;}
						if (fBaseDamage[iChr][ibp][2] >= 2) { fBaseDamage[iChr][ibp][2] = 0; fDamage[iChr][ibp][1] = 0; Excluded_baseDam++; Excluded_baseDam_indir++;}
						if (fSSB[iChr][ibp][1] == 1) { fSSB[iChr][ibp][1] = 0; fDamage[iChr][ibp][2] = 0; Excluded_numSSB++; Excluded_numSSB_dir++;}
						if (fSSB[iChr][ibp][1] >= 2) { fSSB[iChr][ibp][1] = 0; fDamage[iChr][ibp][2] = 0; Excluded_numSSB++; Excluded_numSSB_indir++;}
						if (fSSB[iChr][ibp][2] == 1) { fSSB[iChr][ibp][2] = 0; fDamage[iChr][ibp][3] = 0; Excluded_numSSB++; Excluded_numSSB_dir++;}
						if (fSSB[iChr][ibp][2] >= 2) { fSSB[iChr][ibp][2] = 0; fDamage[iChr][ibp][3] = 0; Excluded_numSSB++; Excluded_numSSB_indir++;}
						// We use the break in the first break so we only check for damage in that to count DSBs. However, damage in second breaks are also disregarded. Also avoid first base pair
						if (ibp > allFreeEndsOfThisChromosome[i-1])
						{
							if (fDSB[iChr][ibp][1] == 1) { fDSB[iChr][ibp][1] = 0; fDamage[iChr][ibp][2] = 0; Excluded_numDSB++; Excluded_numDSB_dir++;}
							if (fDSB[iChr][ibp][1] == 2) { fDSB[iChr][ibp][1] = 0; fDamage[iChr][ibp][2] = 0; Excluded_numDSB++; Excluded_numDSB_indir++;}
							if (fDSB[iChr][ibp][1] >= 3) { fDSB[iChr][ibp][1] = 0; fDamage[iChr][ibp][2] = 0; Excluded_numDSB++; Excluded_numDSB_hybrid++;}
						}
						fDSB[iChr][ibp][2] = 0; fDamage[iChr][ibp][3] = 0;
						fDSB[iChr][ibp][2] = 2; fDamage[iChr][ibp][3] = 0;
						fDSB[iChr][ibp][2] = 0; fDamage[iChr][ibp][3] = 0;
					}
				}
			}
		}
		allFreeEndsOfThisChromosome.clear();
	}
	splitedChromosome.clear();
}

// ****************************************************************************** //
//                        OLD METHODS TO COMPUTE DAMAGE                           //
// ****************************************************************************** //
void TsDefineDamage::SeparateHitsOnDifferentDNAStrands(std::vector<TsHitsRecord*> Hits, std::vector<TsHitsRecord*> &HitsBack, G4double IsStrand1)
{
	// Merge different hits on the same backbone (and neighboring hydration shell)
	// Accumulate energy deposition and only record one effective hit for each damaged backbone

	G4cout << "Strand " << IsStrand1 << " " << Hits.size() << G4endl;
	if (Hits.size() <= 0)
		return;

	G4String baseName = "Base1";
	G4String backName = "Backbone1";
	G4String hydrationName = "HydrationShell1";
	if (!IsStrand1)
	{
		baseName = "Base2";
		backName = "Backbone2";
		hydrationName = "HydrationShell2";
	}

	std::vector<TsHitsRecord*> HitsOnThisStrand;
	std::vector<G4int> DamageBasepairID;
	for (G4int index = 0; index < Hits.size(); index++)
	{
		G4cout << " huh " << Hits[index]->GetVolumeName() << G4endl;
		G4bool foundBase = false, foundBackbone = false, foundHydrationShell = false;
		if (fScoreOnBases)
			foundBase = (strstr(Hits[index]->GetVolumeName(), baseName) != NULL);
		if (fScoreOnBackbones)
			foundBackbone = (strstr(Hits[index]->GetVolumeName(), backName) != NULL);
		if (fScoreOnHydrationShell)
			foundHydrationShell = (strstr(Hits[index]->GetVolumeName(), hydrationName) != NULL);
		if (foundBase || foundBackbone || foundHydrationShell)
		{
			HitsOnThisStrand.push_back(Hits[index]);
			G4cout << " DSB " << Hits[index]->GetBasePairID() << G4endl;
			DamageBasepairID.push_back(Hits[index]->GetBasePairID());
		}
	}
	sort(DamageBasepairID.begin(), DamageBasepairID.end());
	//DamageBasepairID.erase(unique( DamageBasepairID.begin(), DamageBasepairID.end()), DamageBasepairID.end());
	G4double totalEdepInOneCopy = 0;
	G4int index = 0;
	for (G4int iter = 0; iter < DamageBasepairID.size(); iter++)
	{
		for (G4int i = 0; i < HitsOnThisStrand.size(); i++)
		{
			if (HitsOnThisStrand[i]->GetBasePairID() == DamageBasepairID[iter])
			{
				totalEdepInOneCopy += HitsOnThisStrand[i]->GetEdep();
				index = i;
			}
		}
		if (CauseDirectDamage(totalEdepInOneCopy) || totalEdepInOneCopy < 0)
		{
			TsHitsRecord* effectiveDamage = HitsOnThisStrand[index];
			effectiveDamage->SetEdep(totalEdepInOneCopy);
			HitsBack.push_back(effectiveDamage);
		}
		totalEdepInOneCopy = 0;
	}
	HitsOnThisStrand.clear();
	DamageBasepairID.clear();
}

void AddCountOfDifferentTypeSSB(std::vector<TsHitsRecord*> SSBonStrand, G4int &NumTotal, G4int &NumDir, G4int &NumIndir)
{
	NumTotal = SSBonStrand.size();
	for (G4int i = 0; i < SSBonStrand.size(); i++)
	{
		if (SSBonStrand[i]->GetIsDirectDamage())
			NumDir++;
		else
			NumIndir++;
	}
}

void AddCountOfDifferentTypeDSB(std::vector<std::pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, G4int &NumTotal, G4int &NumDir, G4int &NumIndir, G4int &NumHybrid)
{
	NumTotal += DSB_pairs.size();
	for (G4int i = 0; i < DSB_pairs.size(); i++)
	{
		G4bool isDir1 = DSB_pairs[i].first->GetIsDirectDamage();
		G4bool isDir2 = DSB_pairs[i].second->GetIsDirectDamage();
		if (isDir1 && isDir2)
			NumDir++;
		else if (!isDir1 && !isDir2)
			NumIndir++;
		else
			NumHybrid++;
	}
}

void TsDefineDamage::DefineDSBorSSB(std::vector<TsHitsRecord*> HitsBack1, std::vector<TsHitsRecord*> HitsBack2,
									std::vector<std::pair<TsHitsRecord*, TsHitsRecord*>> &DSB_pairs,
									std::vector<TsHitsRecord*> &SSBonStrand1, std::vector<TsHitsRecord*> &SSBonStrand2)
{
	G4cout << " TEST " << HitsBack1.size() << " " << HitsBack2.size() << G4endl;
	if (HitsBack1.size() <= 0 && HitsBack2.size() <= 0)
		return;
	else
	{
		G4int cp1 = 0, cp2 = 0;
		// Find DSB
		for (G4int iter = 0; iter < HitsBack1.size(); iter++)
		{
			cp1 = HitsBack1[iter]->GetBasePairID();
			for (G4int loop = 0; loop < HitsBack2.size(); loop++)
			{
				cp2 = HitsBack2[loop]->GetBasePairID();
				if (abs(cp1 - cp2) <= fDSBSeparation)
				{
					G4bool DSBflag1 = HitsBack1[iter]->GetFlagDSB();
					G4bool DSBflag2 = HitsBack2[loop]->GetFlagDSB();
					if (!DSBflag1 && !DSBflag2)
					{
						HitsBack1[iter]->SetFlagDSB(true);
						HitsBack2[iter]->SetFlagDSB(true);
						DefineDSBType(HitsBack1[iter], HitsBack2[loop]);
						std::pair<TsHitsRecord*, TsHitsRecord*> OnePairDSB;
						OnePairDSB.first = HitsBack1[iter];
						OnePairDSB.second = HitsBack2[loop];
						DSB_pairs.push_back(OnePairDSB);
					}
				}
			}
		}
		// Find SSB
		GetSSBonStrand(HitsBack1, SSBonStrand1);
		GetSSBonStrand(HitsBack2, SSBonStrand2);

		// Calculate damage number
		AddCountOfDifferentTypeDSB(DSB_pairs, numDSB, numDSB_dir, numDSB_indir, numDSB_hybrid);
		AddCountOfDifferentTypeSSB(SSBonStrand1, numSSB, numSSB_dir, numSSB_indir);
		AddCountOfDifferentTypeSSB(SSBonStrand2, numSSB, numSSB_dir, numSSB_indir);
		numSB = numSSB + numDSB * 2;
		numSB_dir = numSSB_dir + numDSB_dir * 2 + numDSB_hybrid;
		numSB_indir = numSSB_indir + numDSB_indir * 2 + numDSB_hybrid;
	}
	G4cout << "NumDSB = " << numDSB << ". DSB_pairs = " << DSB_pairs.size()
			<< ". NumSSB = " << numSSB << ". NumSSB_dir = " << numSSB_dir << ". NumSSB_indir = " << numSSB_indir
			<< ". SSBonStrand1.size() = " << SSBonStrand1.size()
			<< ". SSBonStrand2.size() = " << SSBonStrand2.size() << G4endl;
}

void TsDefineDamage::GetSSBonStrand(std::vector<TsHitsRecord*> Hits, std::vector<TsHitsRecord*> &SSBonStrand)
{
	for (G4int iter = 0; iter < Hits.size(); iter++)
	{
		G4bool isDSB = Hits[iter]->GetFlagDSB();
		if (!isDSB)
		{
			Hits[iter]->SetFlagSSB(true);
			SSBonStrand.push_back(Hits[iter]);
		}
	}
}

void TsDefineDamage::DefineDSBType(TsHitsRecord* &hit1, TsHitsRecord* &hit2)
{
	G4bool isDir1 = hit1->GetIsDirectDamage();
	G4bool isDir2 = hit2->GetIsDirectDamage();
	if (isDir1 && isDir2)
	{
		hit1->SetDSBType(0);
		hit2->SetDSBType(0);
	}
	else if (!isDir1 && !isDir2)
	{
		hit1->SetDSBType(1);
		hit2->SetDSBType(1);
	}
	else
	{
		hit1->SetDSBType(2);
		hit2->SetDSBType(2);
	}
}

void TsDefineDamage::ExcludeShortDNAFragments(std::vector <std::pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs, std::vector<TsHitsRecord*> SSBonStrand1, std::vector<TsHitsRecord*> SSBonStrand2, G4int ChromosomeID)
{
	// Split chromosome according to DSB pairs positions
	TsChromosome* processChromosome = new TsChromosome();
	std::vector<std::vector<G4int>> splitedChromosome = processChromosome->SplitChromosome(DSB_pairs);
	delete processChromosome;

	G4int indivisibleFragmentNum = 0;
	std::vector<TsHitsRecord*> SSBDamages;
	SSBDamages.insert(SSBDamages.end(), SSBonStrand1.begin(), SSBonStrand1.end()); // @suppress("Invalid arguments")
	SSBDamages.insert(SSBDamages.end(), SSBonStrand2.begin(), SSBonStrand2.end()); // @suppress("Invalid arguments")
	std::vector<G4int> AllFreeEndsOfThisChromosome = splitedChromosome[ChromosomeID - 1];

	for (G4int i = 1; i < AllFreeEndsOfThisChromosome.size(); i++)
	{
		G4int fragmentSize = AllFreeEndsOfThisChromosome[i] - AllFreeEndsOfThisChromosome[i - 1];
		G4bool isIndivisibleFragment = (fragmentSize < fLowerFragmentDetectionThreshold || fragmentSize > fUpperFragmentDetectionThreshold);
		if (isIndivisibleFragment)
			indivisibleFragmentNum++;

		// Exclude SSB on the short fragments
		std::vector<TsHitsRecord*>::iterator it;
		for (it = SSBDamages.begin(); it != SSBDamages.end();)
		{
			G4int basePairID = (*it)->GetBasePairID();
			G4bool SSBonThisFragment = (basePairID >= AllFreeEndsOfThisChromosome[i - 1] && basePairID <= AllFreeEndsOfThisChromosome[i]);
			if (isIndivisibleFragment && SSBonThisFragment)
				it = SSBDamages.erase(it); // @suppress("Invalid arguments")
			else
				it++;
		}
		// Exclude DSB on the short fragments
		std::vector<std::pair<TsHitsRecord*, TsHitsRecord*>>::iterator iter;
		for (iter = DSB_pairs.begin(); iter != DSB_pairs.end();)
		{
			G4int basePairID = (*iter).first->GetBasePairID();
			G4bool DSBonThisFragment = (basePairID == AllFreeEndsOfThisChromosome[i]);
			if (isIndivisibleFragment && DSBonThisFragment)
				iter = DSB_pairs.erase(iter); // @suppress("Invalid arguments")
			else
				iter++;
		}
	}
	AddCountOfDifferentTypeDSB(DSB_pairs, Excluded_numDSB, Excluded_numDSB_dir, Excluded_numDSB_indir, Excluded_numDSB_hybrid);
	AddCountOfDifferentTypeSSB(SSBDamages, Excluded_numSSB, Excluded_numSSB_dir, Excluded_numSSB_indir);
	Excluded_numSB = Excluded_numSSB + 2 * Excluded_numDSB;
	Excluded_numSB_dir = Excluded_numSSB_dir + 2 * Excluded_numDSB_dir + Excluded_numDSB_hybrid;
	Excluded_numSB_indir = Excluded_numSSB_indir + 2 * Excluded_numDSB_indir + Excluded_numDSB_hybrid;

	splitedChromosome.clear();
	AllFreeEndsOfThisChromosome.clear();
	SSBDamages.clear();
}

void TsDefineDamage::OutputDNAdamageTuple(std::vector<TsHitsRecord*> HitsBack1, std::vector<TsHitsRecord*> HitsBack2, G4String filename)
{
    std::vector<TsHitsRecord*> Hit;
	std::ofstream outfile;
    outfile.open(filename+"_damage.txt",std::ios::out|std::ios::app);


	for(int backind=1; backind<=2; backind++)
	{
		if (backind==1)
			Hit=HitsBack1;
		else
			Hit=HitsBack2;

		for(unsigned i=0; i<Hit.size(); i++)
		{
			if(Hit[i]->GetFlagSSB()==true || Hit[i]->GetFlagDSB()==true)
            {
                outfile << std::setiosflags(std::ios::fixed)<<std::setprecision(0)<<std::setw(3) <<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetEventID() <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(0)<<std::setw(3) <<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetChromosomeID()<<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(0)<<std::setw(12)<<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetBasePairID() <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(6) <<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetVolumeName() <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(12)<<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetParticleName() <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(18)<<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetParticleEnergy()/eV  <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(12)<<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetEdep()/eV <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(2) <<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagSSB() <<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(2) <<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagDSB()<<"    "
                        << std::setiosflags(std::ios::fixed)<<std::setprecision(4)<<std::setw(2) <<std::setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetDSBType()<<"    ";
                // if(fDefineComplexity)
                // outfile << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(2) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagSSB_P() <<"    "
                //         << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(2) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagDSB_P() <<"    "
                //         << setiosflags(ios::fixed)<<setprecision(4)<<std::setw(2) <<setfill(' ')<<std::setiosflags(std::ios::right)<< Hit[i]->GetFlagDSB_PP();

                outfile << "\n";
            }
		}
	}
	outfile.close();
    Hit.clear();
}

void TsDefineDamage::OutputDNAdamageTupleHeader(G4String filename)
{
	std::ofstream outfile;
    outfile.open(filename+"_damage.header",std::ios::out);

    outfile << "Columns of data are as follows:\n"
            << "1: Event ID" <<"\n"
            << "2: Chromosom ID " <<"\n"
            << "3: bp ID on chromosome "<<"\n"
            << "4: Volume name" <<"\n"
            << "5: Particel name" <<"\n"
            << "6: Particel Energy [eV]" <<"\n"
            << "7: Energy deposit [eV]" <<"\n"
            << "8: FlagSSB    (0-Not SSB, 1-Is SSB)" <<"\n"
            << "9: FlagDSB    (0-Not DSB, 1-Is DSB)" <<"\n"
            << "10: DSBType   (-1-Not DSB, 0-dir DSB, 1-indir DSB, 2-hybrid DSB)" <<"\n";
    // if(fDefineComplexity)
    // outfile << "11: FlagSSB+  (0-Not SSB+, 1-Is SSB+)" <<"\n"
    //         << "12: FlagDSB+  (0-Not DSB+, 1-Is DSB+)" <<"\n"
    //         << "13: FlagDSB++ (0-Not DSB++, 1-Is DSB++)" <<"\n\n";
	outfile.close();

}

void TsDefineDamage::OutputDNAdamageSummary(G4String filename)
{
    G4double NucleusDose = (fEdep*1.6E-13)/(fNucleusMass/1000); // MeV/g —> Gy
	std::ofstream outfile;
    outfile.open(filename+"_damage.header",std::ios::out);

    outfile <<"Recoreded Energy deposition(MeV): "<<fEdep<<"\n"
            <<"(Zeff/β)^2 :                      "<<fZetaBeta_sq<<"\n"
            <<"LET(keV/um):                      "<<fLET<<"\n"
            <<"Dose in nucleusmass(Gy):          "<<NucleusDose<<"\n"
            <<"DNA content in nucleus (Gbp):     "<<fTotalDNAContent/1E9<<"\n"
            <<"Simulated events:                 "<<fEventID+1<<"\n\n";
    outfile <<"************************************************"<<"\n"
            <<"************** Yied of DNA damage **************" <<"\n"
            <<"Total SB:           "<<  numSB       <<"\n"
            <<"Total SB_direct:    "<<  numSB_dir   <<"\n"
            <<"Total SB_indirect:  "<<  numSB_indir <<"\n"
            <<"Total SSB:          "<<  numSSB      <<"\n"
            <<"Total SSB_direct:   "<<  numSSB_dir  <<"\n"
            <<"Total SSB_indirect: "<<  numSSB_indir<<"\n"
            <<"Total DSB:          "<<  numDSB    <<"\n"
            <<"Total DSB_direct:   "<<  numDSB_dir    <<"\n"
            <<"Total DSB_indirect: "<<  numDSB_indir  <<"\n"
            <<"Total DSB_hybrid:   "<<  numDSB_hybrid <<"\n"
            <<"Total SSB+:  "<<  numSSB_P  <<"\n"
            <<"Total DSB+:  "<<  numDSB_P  <<"\n"
            <<"Total DSB++: "<<  numDSB_PP <<"\n\n"
            <<"SSB/Gy/Gbp : "<<   numSSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"DSB/Gy/Gbp : "<<   numDSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"************************************************" <<"\n\n";
    if(fExcludeShortFragment)
    outfile <<"************************************************"<<"\n"
            <<"**  Yied of DNA damage (ExcludeShortFragment) **" <<"\n"
            <<"Total SB:           "<<  Excluded_numSB       <<"\n"
            <<"Total SB_direct:    "<<  Excluded_numSB_dir   <<"\n"
            <<"Total SB_indirect:  "<<  Excluded_numSB_indir <<"\n"
            <<"Total SSB:          "<<  Excluded_numSSB  <<"\n"
            <<"Total SSB_direct:   "<<  Excluded_numSSB_dir  <<"\n"
            <<"Total SSB_indirect: "<<  Excluded_numSSB_indir  <<"\n"
            <<"Total DSB:          "<<  Excluded_numDSB    <<"\n"
            <<"Total DSB_direct:   "<<  Excluded_numDSB_dir    <<"\n"
            <<"Total DSB_indirect: "<<  Excluded_numDSB_indir  <<"\n"
            <<"Total DSB_hybrid:   "<<  Excluded_numDSB_hybrid <<"\n"
            <<"Total SSB+:  "<<   Excluded_numSSB_P  <<"\n"
            <<"Total DSB+:  "<<   Excluded_numDSB_P  <<"\n"
            <<"Total DSB++: "<<   Excluded_numDSB_PP <<"\n"
            <<"SSB/Gy/Gbp : "<<   Excluded_numSSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"DSB/Gy/Gbp : "<<   Excluded_numDSB/NucleusDose/(fTotalDNAContent/1E9)<<"\n"
            <<"************************************************" <<"\n";
	outfile.close();
}

int GetStranID(G4String name)
{
    if (name == "Backbone1" || name == "HydrationShell1")
        return 1;
	//else if (name == "back2"|| name == "water2")
    return 4;
}

void TsDefineDamage::SeparateDamageSitesOnOneStrand(std::vector <TsHitsRecord*> Hits, std::vector <G4ThreeVector> &XYZPosition, std::vector <double  > &ChromosomePosition,
                                              std::vector <double  > &Cause, std::vector <int > &DamageSiteDSBflag, std::vector< std::vector <std::vector <int > > >  &DamageSites)
{
    if (Hits.size()==0)
        return;

    int StartIDofDamageSite = Hits[0]->GetBasePairID();
    int EndIDofDamageSite   = StartIDofDamageSite+fDSBSeparation-1;
    int DamageType     = 0;
    int localBasePairID = 0;
    int BasePairID     = 0;
    bool DirectDamageFlag   = false;
    bool IndirectDamageFlag = false;
    G4ThreeVector aXYZPos ;

    unsigned iter=0;
    while( iter<Hits.size()  )
    {

        std::vector <std::vector <int > > aDamageSite;
        for(unsigned i=iter; i<Hits.size(); i++)
        {
            BasePairID = Hits[i]->GetBasePairID();
            if( BasePairID >= StartIDofDamageSite && BasePairID<= EndIDofDamageSite )
            {
                if(Hits[i]->GetIsDirectDamage() ) {DamageType = 1; DirectDamageFlag =true;}
                if(!Hits[i]->GetIsDirectDamage()) {DamageType = 2; IndirectDamageFlag =true;}
                localBasePairID =BasePairID-StartIDofDamageSite +1;
                std::vector <int > aDamage = {GetStranID(Hits[i]->GetVolumeName()), localBasePairID,  DamageType }; // @suppress("Invalid arguments")
                aDamageSite.push_back(aDamage);
                aXYZPos    = Hits[i]->GetPosition();
                iter++;
            }
            else
                break;
        }

        int ChromoID = Hits[0]->GetChromosomeID();
        int ChromoLength = fChromosomeDNAContent[ChromoID-1];
        double ChromoPos = (double)BasePairID/ChromoLength;

        bool aCause =0;
             if (DirectDamageFlag == true && IndirectDamageFlag == false)  aCause =0;
        else if (DirectDamageFlag == false && IndirectDamageFlag == true)  aCause =1;
        else if (DirectDamageFlag == true  && IndirectDamageFlag == true)  aCause =2;


        XYZPosition.push_back(aXYZPos);            // Field 2
        ChromosomePosition.push_back(ChromoPos);   // Field 4
        Cause.push_back(aCause);                   // Field 5
        DamageSiteDSBflag.push_back(0);            // Field 6, DSB flag
        DamageSites.push_back(aDamageSite);        // Field 7

        if (iter == Hits.size())
            break;

        StartIDofDamageSite = Hits[iter]->GetBasePairID();
        EndIDofDamageSite   = StartIDofDamageSite + fDSBSeparation-1;
        G4cout<<"iter="<<iter<<", DamageSites.size()="<<DamageSites.size()<<G4endl;
    }
}

void TsDefineDamage::OutputSDDFile(bool OnlyIncludeDSBinSDD, std::vector <std::pair<TsHitsRecord*, TsHitsRecord*>> DSB_pairs,
                            std::vector<TsHitsRecord*> SSBonStrand1,std::vector<TsHitsRecord*> SSBonStrand2,
                            G4String filename, G4int EventID, G4int ExposureID, G4int ChromosomeID)
{
    // Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al. A New Standard DNA Damage (SDD) Data Format[J]. Radiation research, 2018.
    // Data was generated with 7 SDD filed.
    // Filed 1: Calssification (2 Int) : Is damage associated with new primary particle or new exposure? And event ID.
    // Filed 2: X,Y,Z (3x3 Floats): Spatial X, Y, Z coordinates and extent (unit: um)
    // Filed 3: Chromosome IDs (4 Int) : ID of chromosome/chromatid where damage occurred and on which arm (long/short) or specification of non-nuclear DNA type.
    // Filed 4: Chromosome position (Float): Location of damage within chromosome
    // Filed 5: Cause (3 Int): Cause of damage - direct or indirect and number
    // Filed 6: Damage type (3 int): Types of damage at site (Base damage, SSB, DSB)
    // Filed 7: Full break spec (.../3 int/...) Full description of strand break structure

    std::vector <G4ThreeVector> XYZPosition;             // For field 2
    std::vector <double  > ChromosomePosition;           // For filed 4
    std::vector <double  > Cause;                        // For filed 5
    std::vector <int > DamageSiteDSBflag;                // For field 6
    std::vector< std::vector <std::vector <int > > > DamageSites;  // For field 7

    std::vector<TsHitsRecord*> SSBDamages;
    SSBDamages.insert(SSBDamages.end(),SSBonStrand1.begin(),SSBonStrand1.end()); // @suppress("Invalid arguments")
    SSBDamages.insert(SSBDamages.end(),SSBonStrand2.begin(),SSBonStrand2.end()); // @suppress("Invalid arguments")

    for( unsigned i=0; i<DSB_pairs.size(); i++)
    {
        std::vector <std::vector <int > > aDamageSite;
        std::pair<TsHitsRecord*, TsHitsRecord*> onePairDSB = DSB_pairs[i];
        int DamageSiteStartID = std::min(onePairDSB.first->GetBasePairID(), onePairDSB.second->GetBasePairID());
        int DamageSiteEndID   = DamageSiteStartID + fDSBSeparation-1;
        int BasePairID = 0;
        int DamageType = 0;
        bool DirectDamageFlag   = false;
        bool IndirectDamageFlag = false;

        // ------ Get all SSB/DSB damages on this site
        std::vector<TsHitsRecord*> HitsOnThisSite ;
        HitsOnThisSite.push_back(onePairDSB.first);
        HitsOnThisSite.push_back(onePairDSB.second);

        std::vector<TsHitsRecord*>::iterator it;
        for(it=SSBDamages.begin();it!=SSBDamages.end();)
        {
            BasePairID = (*it)->GetBasePairID();
            if(BasePairID>= DamageSiteStartID && BasePairID<=DamageSiteEndID )
            {
                HitsOnThisSite.push_back(*it);
                it=SSBDamages.erase(it); // @suppress("Invalid arguments")
            }
            else
                ++it;
        }

        // ------ Get damage spectrum in this site
        for (unsigned ihits=0; ihits<HitsOnThisSite.size(); ihits++)
        {
            if ( HitsOnThisSite[ihits]->GetIsDirectDamage() )
                {DamageType = 1; DirectDamageFlag =true;}
            else
                {DamageType = 2; IndirectDamageFlag =true;}

            int localBasePairID = HitsOnThisSite[ihits]->GetBasePairID() - DamageSiteStartID +1;
            std::vector<G4int> OneHit = {GetStranID(HitsOnThisSite[ihits]->GetVolumeName()), localBasePairID,  DamageType }; // @suppress("Invalid arguments")

            aDamageSite.push_back(OneHit);
        }

        // ------ sort damage spectrum, make sure damages are recorded from stran1 to stand 2
        sort(aDamageSite.begin(), aDamageSite.end());

        int ChromoLength      = fChromosomeDNAContent[ChromosomeID-1];
        double ChromoPos      = (double)DamageSiteStartID/ChromoLength;
        G4ThreeVector aXYZPos = HitsOnThisSite[0]->GetPosition();

        int aCause =0;
             if (DirectDamageFlag == true  && IndirectDamageFlag == false)  aCause=0;
        else if (DirectDamageFlag == false && IndirectDamageFlag == true )  aCause=1;
        else if (DirectDamageFlag == true  && IndirectDamageFlag == true )  aCause=2;


        XYZPosition.push_back(aXYZPos);            // Field 2
        ChromosomePosition.push_back(ChromoPos);   // Field 4
        Cause.push_back(aCause);                   // Field 5
        DamageSiteDSBflag.push_back(1);            // Field 6, DSB flag
        DamageSites.push_back(aDamageSite);        // Field 7
    }

    if(!OnlyIncludeDSBinSDD )
        SeparateDamageSitesOnOneStrand(SSBDamages, XYZPosition, ChromosomePosition, Cause, DamageSiteDSBflag, DamageSites);

    // ------ output
    std::ofstream outfile;
    outfile.open(filename+"_sdd.txt",std::ios::out|std::ios::app);

    for(unsigned i=0; i<DamageSites.size(); i++)
    {
        std::vector <std::vector <int > >  aDamageSite = DamageSites[i];

        G4int newExposoureFlag =0;
        if(EventID != lastOutputEventID)
        {
            newExposoureFlag = 1;
            lastOutputEventID = EventID;
        }
        if (ExposureID !=  lastOutputExposureID)
        {
            newExposoureFlag = 2;
            lastOutputExposureID = ExposureID;
        }

        outfile << newExposoureFlag<<","<< EventID<<"; ";                                                 // Filed 1
        outfile << XYZPosition[i].x()/um<<","<<XYZPosition[i].y()/um<<","<<XYZPosition[i].z()/um<<"; ";   // Filed 2
        outfile << "1," << ChromosomeID-1<<","  << "1,"<< "0" <<"; ";                                     // Filed 3, strat from 0
        outfile << ChromosomePosition[i]<<"; ";                                                           // Filed 4
        outfile << Cause[i]<<"; ";                                                                        // Filed 5
        outfile <<"0,"<<aDamageSite.size()<<","<<DamageSiteDSBflag[i]<<"; ";                              // Filed 6, num of base damage, num of strand break, DSB flag

        for(unsigned j=0; j<aDamageSite.size(); j++)
            outfile <<aDamageSite[j][0]<<","<<aDamageSite[j][1]<<","<<aDamageSite[j][2]<<"/";             //  Filed 7: Full break spec

        outfile << "; \n";
    }
    outfile.close();

    XYZPosition.clear();
    ChromosomePosition.clear();
    Cause.clear();
    DamageSiteDSBflag.clear();
    DamageSites.clear();
}
