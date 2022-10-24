// Extra Class for TsScoreDNADamageSBS
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

#include "TsDNADamageCalculator.hh"

#include "Randomize.hh"

TsDNADamageCalculator::TsDNADamageCalculator()
{

	fNumberOfBasePairForDSB = 10;
	fDirectDamageThreshold = 1E6;
	fUseLinearProbabilityForDirectDamage = false;
	fLowerLimitLinearProbability = 1E6;
	fUpperLimitLinearProbability = 1E6;

	fExcludeShortFragments = false;
	fLowerThresholdForFragmentDetection = 0;
	fUpperThresholdForFragmentDetection = 3E8;

	fNumSB = 0;
	fNumSBDirect = 0;
	fNumSBQuasiDirect = 0;
	fNumSBIndirect = 0;
	fNumSSB = 0;
	fNumSSBDirect = 0;
	fNumSSBQuasiDirect = 0;
	fNumSSBIndirect = 0;
	fNumDSB = 0;
	fNumDSBDirect = 0;
	fNumDSBIndirect = 0;
	fNumDSBDirectIndirect = 0;
	fNumDSBDirectQuasiDirect = 0;
	fNumDSBQuasiDirectQuasiDirect = 0;
	fNumDSBIndirectQuasiDirect = 0;
	fNumBaseDamage = 0;
	fNumBaseDamageDirect = 0;
	fNumBaseDamageQuasiDirect = 0;
	fNumBaseDamageIndirect = 0;
	fNumSSBPlus = 0; fNumDSBPlus = 0; fNumDSBComplex = 0;

	fEventID = 0;
	fOutputFileName = "";
	fOutputMode = "";
	fWriteCSVFile = false;

	fMinimalModeForSDD = false;
	fReturnOnlyDSBinSDD = false;
}

TsDNADamageCalculator::~TsDNADamageCalculator() {}

void TsDNADamageCalculator::ComputeStrandBreaks(std::vector<TsHitInDNA*> hits)
{
	// Clears maps from previous calculations
	fMapEdep.clear();
	fIsThereIndirectDamage.clear();
	fIsThereQuasiDirectDamage.clear();
	fDamageMap.clear();
	fDSBMap.clear();
	fSSBMap.clear();
	fBDMap.clear();
	fDamagePositions.clear();
	fDSB3DPositions.clear();

	// Gets types of damage from the hits vector
	for (unsigned int i = 0; i < hits.size(); i++)
	{
		G4int basePairID = hits[i]->GetBasePairID();
		G4double edep = hits[i]->GetEdep();
		G4int chromosomeID = hits[i]->GetChromosomeID();
		if (chromosomeID >= 0)
		{
			G4int strand = hits[i]->GetStrandNumber();
			G4int dnacomp = hits[i]->GetDNAComponentID();
			G4int iComp = dnacomp * 2 + (strand - 1); // This formula assigns 0 to base1, 1 to base2, 2 to backbone1, 3 to backbone2
			if (hits[i]->GetDamageType() == direct) fMapEdep[chromosomeID][basePairID][iComp] += edep;
			else if (hits[i]->GetDamageType() == indirect) fIsThereIndirectDamage[chromosomeID][basePairID][iComp] = true;
			else if (hits[i]->GetDamageType() == quasidirect) fIsThereQuasiDirectDamage[chromosomeID][basePairID][iComp] = true;
			fDamagePositions[chromosomeID][basePairID][iComp] = hits[i]->GetPosition();
		}
	}

	// Classifies damages. First loop through direct damage
	for (auto& chrMap : fMapEdep)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fMapEdep[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& energyMap : fMapEdep[iChr][iBp])
			{
				G4int iComp = energyMap.first;
				if (CauseDirectDamage(fMapEdep[iChr][iBp][iComp]))
					fDamageMap[iChr][iBp][iComp] = direct;
				else
					fDamageMap[iChr][iBp][iComp] = nodamage;
			}
		}
	}
	// Classifies damages. Second loop through indirect damage
	for (auto& chrMap : fIsThereIndirectDamage)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fIsThereIndirectDamage[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& indirectDamageMap : fIsThereIndirectDamage[iChr][iBp])
			{
				G4int iComp = indirectDamageMap.first;
				if (fDamageMap[iChr][iBp][iComp] == direct)
					fDamageMap[iChr][iBp][iComp] = multiple;
				else
					fDamageMap[iChr][iBp][iComp] = indirect;
			}
		}
	}
	// Classifies damages. Third loop through quasi-direct damage
	for (auto& chrMap : fIsThereQuasiDirectDamage)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fIsThereQuasiDirectDamage[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& quasidirectDamageMap : fIsThereQuasiDirectDamage[iChr][iBp])
			{
				G4int iComp = quasidirectDamageMap.first;
				if (fDamageMap[iChr][iBp][iComp] < direct)
					fDamageMap[iChr][iBp][iComp] = quasidirect;
				else if (fDamageMap[iChr][iBp][iComp] == indirect)
					fDamageMap[iChr][iBp][iComp] = multiplewithquasidirect;
			}
		}
	}
	if (fWriteCSVFile) WriteDNADamageCSV();

	// Loops through the damage to obtain DSBs
	std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs;
	for (auto& chrMap : fDamageMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fDamageMap[iChr])
		{
			G4int iBp = bpMap.first;
			// Checks if there is damage in backbone 1 (iComp=2) and there is not already a DSB identified in strand 1 (fDSBMap index 1)
			if (fDamageMap[iChr][iBp][2] > 0 && fDSBMap[iChr][iBp][1] == 0)
			{
				G4bool dsbFound = false;
				// Checks nearby base pairs in strand 2
				for (G4int i2 = 0; i2 < fNumberOfBasePairForDSB; i2++)
				{
					// Goes from iBp to iBp+10 to find damage in backbone 2 (iComp=3) and there is not already a DSB in strand 2 (fDSBMap index 2)
					G4int typeDamageInStrand2 = fDamageMap[iChr][iBp+i2][3];
					if (typeDamageInStrand2 > 0 && fDSBMap[iChr][iBp+i2][2] == 0)
					{
						G4int adjustedTypeDamageInStrand2 = indirect; // To convert multiple in direct or quasi-direct. Indirect if no direct or multiple effects are found
						// Checks if damage in backbone 2 is direct or quasi-direct and keeps it in that case
						if (typeDamageInStrand2 == direct || typeDamageInStrand2 == quasidirect) adjustedTypeDamageInStrand2 = typeDamageInStrand2;
						// If it is multiple damage, considers 'direct' or 'quasidirect' (since it would happen without chemistry)
						if (typeDamageInStrand2 == multiple) adjustedTypeDamageInStrand2 = direct;
						if (typeDamageInStrand2 == multiplewithquasidirect) adjustedTypeDamageInStrand2 = quasidirect;

						// Once damage is identified in strand 2 (and damage ensured in strand 1), looks for the closest damage in strand 1 (just in case there are multiple close damages)
						G4int closestPosInStrand1 = iBp;
						G4int adjustedTypeDamageInStrand1 = indirect;
						for (G4int i1 = 0; i1 < fNumberOfBasePairForDSB; i1++)
						{
							// Goes from ibp+i2 increasing one base pair...
							G4int typeDamageInStrand1 = fDamageMap[iChr][iBp+i2+i1][2];
							if (typeDamageInStrand1 > 0 && fDSBMap[iChr][iBp+i2+i1][1] == 0)
							{
								closestPosInStrand1 = iBp + i2 + i1;
								if (typeDamageInStrand1 == direct || typeDamageInStrand1 == quasidirect) adjustedTypeDamageInStrand1 = typeDamageInStrand1;
								if (typeDamageInStrand1 == multiple) adjustedTypeDamageInStrand1 = direct;
								if (typeDamageInStrand1 == multiplewithquasidirect) adjustedTypeDamageInStrand1 = quasidirect;
								break;
							}
							// ... and decreasing one base pair alternatively (avoiding index < 0)
							typeDamageInStrand1 = fDamageMap[iChr][iBp+i2-i1][2];
							if (iBp+i2-i1 >= 0 && i1 > 0 && typeDamageInStrand1 > 0 && fDSBMap[iChr][iBp+i2-i1][1] == 0)
							{
								closestPosInStrand1 = iBp + i2 - i1;
								if (typeDamageInStrand1 == direct || typeDamageInStrand1 == quasidirect) adjustedTypeDamageInStrand1 = typeDamageInStrand1;
								if (typeDamageInStrand1 == multiple) adjustedTypeDamageInStrand1 = direct;
								if (typeDamageInStrand1 == multiplewithquasidirect) adjustedTypeDamageInStrand1 = quasidirect;
								break;
							}
						}
						fDSBMap[iChr][closestPosInStrand1][1] = adjustedTypeDamageInStrand1;
						fDSBMap[iChr][iBp+i2][2] = adjustedTypeDamageInStrand2;
						std::pair<G4int, G4int> pos(closestPosInStrand1, iBp+i2);
						// Combines damages in both strands: 2->Direct, 3->Hybrid, 4->Indirect, 5->Direct with 1 quasi, 6->Hybrid with 1 quasi, 8->Direct with 2 quasi
						DSBPairs[iChr][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2;
						G4ThreeVector posSB1 = fDamagePositions[iChr][closestPosInStrand1][2];
						G4ThreeVector posSB2 = fDamagePositions[iChr][iBp+i2][3];
						G4ThreeVector center = (posSB1 + posSB2) / 2;
						fDSB3DPositions.push_back(center);
						dsbFound = true;
					}
					typeDamageInStrand2 = fDamageMap[iChr][iBp-i2][3];
					// Goes from iBp to iBp-10 to find damage in backbone 2 (alternatively, if dsb is found then breaks) (avoiding index < 0)
					if (iBp - i2 >= 0 && i2 > 0 && typeDamageInStrand2 > 0 && fDSBMap[iChr][iBp-i2][2] == 0)
					{
						G4int adjustedTypeDamageInStrand2 = indirect; // To convert multiple in direct or quasi-direct. Indirect if no direct or multiple effects are found
						// Checks if damage in backbone 2 is direct or quasi-direct and keeps it in that case
						if (typeDamageInStrand2 == direct || typeDamageInStrand2 == quasidirect) adjustedTypeDamageInStrand2 = typeDamageInStrand2;
						// If it is multiple damage, considers 'direct' or 'quasidirect' (since it would happen without chemistry)
						if (typeDamageInStrand2 == multiple) adjustedTypeDamageInStrand2 = direct;
						if (typeDamageInStrand2 == multiplewithquasidirect) adjustedTypeDamageInStrand2 = quasidirect;

						// Once damage is identified in strand 2 (and damage ensured in strand 1), looks for the closest damage in strand 1 (just in case there are multiple close damages)
						G4int closestPosInStrand1 = iBp;
						G4int adjustedTypeDamageInStrand1 = indirect;
						for (G4int i1 = 0; i1 < fNumberOfBasePairForDSB; i1++)
						{
							// Goes from ibp-i2 increasing one base pair...
							G4int typeDamageInStrand1 = fDamageMap[iChr][iBp-i2+i1][2];
							if (typeDamageInStrand1 > 0 && fDSBMap[iChr][iBp-i2+i1][1] == 0)
							{
								closestPosInStrand1 = iBp - i2 + i1;
								if (typeDamageInStrand1 == direct || typeDamageInStrand1 == quasidirect) adjustedTypeDamageInStrand1 = typeDamageInStrand1;
								if (typeDamageInStrand1 == multiple) adjustedTypeDamageInStrand1 = direct;
								if (typeDamageInStrand1 == multiplewithquasidirect) adjustedTypeDamageInStrand1 = quasidirect;
								break;
							}
							// ... and decreasing one base pair alternatively (avoiding index < 0)
							typeDamageInStrand1 = fDamageMap[iChr][iBp-i2-i1][2];
							if (iBp-i2-i1 >= 0 && i1 > 0 && typeDamageInStrand1 > 0 && fDSBMap[iChr][iBp-i2-i1][1] == 0)
							{
								closestPosInStrand1 = iBp - i2 - i1;
								if (typeDamageInStrand1 == direct || typeDamageInStrand1 == quasidirect) adjustedTypeDamageInStrand1 = typeDamageInStrand1;
								if (typeDamageInStrand1 == multiple) adjustedTypeDamageInStrand1 = direct;
								if (typeDamageInStrand1 == multiplewithquasidirect) adjustedTypeDamageInStrand1 = quasidirect;
								break;
							}
						}
						fDSBMap[iChr][closestPosInStrand1][1] = adjustedTypeDamageInStrand1;
						fDSBMap[iChr][iBp-i2][2] = adjustedTypeDamageInStrand2;
						std::pair<G4int, G4int> pos(closestPosInStrand1, iBp-i2);
						// Combines damages in both strands: 2->Direct, 3->Hybrid, 4->Indirect, 5->Direct with 1 quasi, 6->Hybrid with 1 quasi, 8->Direct with 2 quasi
						DSBPairs[iChr][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2;
						G4ThreeVector posSB1 = fDamagePositions[iChr][closestPosInStrand1][2];
						G4ThreeVector posSB2 = fDamagePositions[iChr][iBp-i2][3];
						G4ThreeVector center = (posSB1 + posSB2) / 2;
						fDSB3DPositions.push_back(center);
						dsbFound = true;
					}
					if (dsbFound) break;
				}
			}
		}
	}

	// Loops through damage map to exclude damage already classified as DSB, to populate SSB and base damage maps
	for (auto& chrMap : fDamageMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fDamageMap[iChr])
		{
			G4int iBp = bpMap.first;
			// Single strand breaks
			if ((fDamageMap[iChr][iBp][2] > 0 && !(fDSBMap[iChr][iBp][1] > 0)) || (fDamageMap[iChr][iBp][3] > 0 && !(fDSBMap[iChr][iBp][2] > 0)))
			{
				if (fDamageMap[iChr][iBp][2] > 0) fSSBMap[iChr][iBp][1] = fDamageMap[iChr][iBp][2];
				if (fDamageMap[iChr][iBp][3] > 0) fSSBMap[iChr][iBp][2] = fDamageMap[iChr][iBp][3];
			}
			// Base damages
			if (fDamageMap[iChr][iBp][0] > 0) fBDMap[iChr][iBp][1] = fDamageMap[iChr][iBp][0];
			if (fDamageMap[iChr][iBp][1] > 0) fBDMap[iChr][iBp][2] = fDamageMap[iChr][iBp][1];
		}
	}
	if (fExcludeShortFragments) ExcludeShortDNAFragments(DSBPairs);
	QuantifyDamage(DSBPairs);
}

void TsDNADamageCalculator::QuantifyDamage(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs)
{
	// Reinitializes all quantities from previous quantifications
	fNumSB = 0; fNumSBDirect = 0; fNumSBQuasiDirect = 0; fNumSBIndirect = 0;
	fNumSSB = 0; fNumSSBDirect = 0; fNumSSBQuasiDirect = 0; fNumSSBIndirect = 0;
	fNumDSB = 0; fNumDSBDirect = 0; fNumDSBIndirect = 0; fNumDSBDirectIndirect = 0; fNumDSBDirectQuasiDirect = 0; fNumDSBQuasiDirectQuasiDirect = 0; fNumDSBIndirectQuasiDirect = 0;
	fNumBaseDamage = 0; fNumBaseDamageDirect = 0; fNumBaseDamageQuasiDirect = 0; fNumBaseDamageIndirect = 0;
	fNumSSBPlus = 0; fNumDSBPlus = 0; fNumDSBComplex = 0;

	// Loops through DSB pairs
	for (auto& chrMap : DSBPairs)
	{
		G4int iChr = chrMap.first;
		for (auto& dsbMap : DSBPairs[iChr])
		{
			fNumDSB++; fNumSB++; fNumSB++;
			G4int typeDamage = dsbMap.second;
			if (typeDamage == direct + direct) 					{ fNumDSBDirect++; fNumSBDirect++; fNumSBDirect++;}
			else if (typeDamage == direct + indirect) 			{ fNumDSBDirectIndirect++; fNumSBDirect++; fNumSBIndirect++; }
			else if (typeDamage == indirect + indirect)			{ fNumDSBIndirect++; fNumSBIndirect++; fNumSBIndirect++; }
			else if (typeDamage == direct + quasidirect)		{ fNumDSBDirectQuasiDirect++; fNumDSBDirect++; fNumSBDirect++; fNumSBDirect++; fNumSBQuasiDirect++; } // Quasi-direct also counts as direct
			else if (typeDamage == indirect + quasidirect)		{ fNumDSBIndirectQuasiDirect++; fNumDSBDirectIndirect++; fNumSBDirect++; fNumSBQuasiDirect++; fNumSBIndirect++; }
			else if (typeDamage == quasidirect + quasidirect)	{ fNumDSBQuasiDirectQuasiDirect++; fNumDSBDirect++; fNumSBDirect++; fNumSBDirect++; fNumSBQuasiDirect++; fNumSBQuasiDirect++;}
		}
	}

	// Loops through SSB
	for (auto& chrMap : fSSBMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fSSBMap[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& ssbMap : fSSBMap[iChr][iBp])
			{
				fNumSSB++; fNumSB++;
				G4int typeDamage = ssbMap.second;
				if (typeDamage == direct || typeDamage == multiple)		{ fNumSSBDirect++; fNumSBDirect++; }
				else if (typeDamage == indirect)						{ fNumSSBIndirect++; fNumSBIndirect++; }
				else if (typeDamage == quasidirect || typeDamage == multiplewithquasidirect)
				{
					fNumSSBDirect++; fNumSSBQuasiDirect++;
					fNumSBDirect++; fNumSBQuasiDirect++;
				}
			}
		}
	}

	// Loops through base damage
	for (auto& chrMap : fBDMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fBDMap[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& bdMap : fBDMap[iChr][iBp])
			{
				fNumBaseDamage++;
				G4int typeDamage = bdMap.second;
				if (typeDamage == direct || typeDamage == multiple)		fNumBaseDamageDirect++;
				else if (typeDamage == indirect)						fNumBaseDamageIndirect++;
				else if (typeDamage == quasidirect || typeDamage == multiplewithquasidirect)	{ fNumBaseDamageDirect++; fNumBaseDamageQuasiDirect++; }
			}
		}
	}
}

std::map<G4int, std::vector<G4int>> TsDNADamageCalculator::GetDamageSites()
{
	// Uses a scoring system to sort damage sites in descendent complexity order
	// Base damage is assigned 0.001 value, SSB is assigned 0.1 value and DSB is assigned 10 value in total (5 per break)
	std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> damageScore;
	// Loops through base damages
	for (auto& chrMap : fBDMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fBDMap[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& bdmap : fBDMap[iChr][iBp])
			{
				G4int strand = bdmap.first;
				if (fBDMap[iChr][iBp][strand] > 0)	damageScore[iChr][iBp][strand] += 0.001;
			}
		}
	}
	// Loops through SSB
	for (auto& chrMap : fSSBMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fSSBMap[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& ssbmap : fSSBMap[iChr][iBp])
			{
				G4int strand = ssbmap.first;
				if (fSSBMap[iChr][iBp][strand] > 0)	damageScore[iChr][iBp][strand] += 0.1;
			}
		}
	}
	// Loops through DSB
	for (auto& chrMap : fDSBMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fDSBMap[iChr])
		{
			G4int iBp = bpMap.first;
			for (auto& dsbmap : fDSBMap[iChr][iBp])
			{
				G4int strand = dsbmap.first;
				if (fDSBMap[iChr][iBp][strand] > 0)	damageScore[iChr][iBp][strand] += 5.0;
			}
		}
	}

	// Separates damage in 10-bp length (or the number specified by fNumberOfBasePairForDSB) damage sites
	std::map<G4int, std::vector<G4int>> chromosomeAndInitialBpIdOfDamageSites;
	for (auto& chrMap : damageScore)
	{
		G4int iChr = chrMap.first;
		std::vector<G4double> linearAccumulatedDamage;
		std::vector<G4int> initialBpIds;
		std::vector<G4int> startingBpIdForDamageSites;
		G4double maxDamage = 10 * fNumberOfBasePairForDSB;
		G4double minDamage = 1e-8;
		if (fReturnOnlyDSBinSDD) minDamage = 9.99999999; // Only includes damage sites with at least one DSB
		while (maxDamage > minDamage)
		{
			for (auto& bpMap : damageScore[iChr])
			{
				G4int iBp = bpMap.first;
				G4double sumDamage = 0;
				for (G4int i = 0; i < fNumberOfBasePairForDSB; i++)
				{
					if (damageScore[iChr].find(iBp+i) != damageScore[iChr].end())
						sumDamage += damageScore[iChr][iBp+i][1] + damageScore[iChr][iBp+i][2];
				}
				linearAccumulatedDamage.push_back(sumDamage);
				initialBpIds.push_back(iBp);
			}
			// Gets maximum damage and its bp index
			G4double m = 0;
			G4int index = 0;
			for (unsigned int i = 0; i < linearAccumulatedDamage.size(); i++)
			{
				if (m < linearAccumulatedDamage[i])
				{
					m = linearAccumulatedDamage[i];
					index = initialBpIds[i];
				}
			}
			maxDamage = m;
			if (maxDamage >= 0.2 && maxDamage < 10) fNumSSBPlus++;
			if (maxDamage >= 20 && maxDamage < 30) fNumDSBPlus++;
			if (maxDamage >= 30) fNumDSBComplex++;
			if (maxDamage < minDamage) break;
			// Populates a new vector sorted by maximum damage
			startingBpIdForDamageSites.push_back(index);
			// Excludes this damage site for next steps
			for (G4int i = 0; i < fNumberOfBasePairForDSB; i++)
			{
				if (damageScore[iChr].find(index + i) != damageScore[iChr].end())
				{
					damageScore[iChr][index + i][1] = 0;
					damageScore[iChr][index + i][2] = 0;
				}
			}
			linearAccumulatedDamage.clear();
			initialBpIds.clear();
		}
		chromosomeAndInitialBpIdOfDamageSites[iChr] = startingBpIdForDamageSites;
	}
	return chromosomeAndInitialBpIdOfDamageSites;
}

void TsDNADamageCalculator::ExcludeShortDNAFragments(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs)
{
	// Get initial positions for DSB to define fragments
	std::map<G4int, std::vector<G4int>> initialBasePairsInEachChromosome;
	for (auto& chrDSBPair : DSBPairs)
	{
		G4int iChr = chrDSBPair.first;
		for (auto& pairDamage : DSBPairs[iChr])
		{
			std::pair<G4int, G4int> pairPos = pairDamage.first;
			initialBasePairsInEachChromosome[iChr].push_back(std::min(pairPos.first, pairPos.second));
		}
	}
	// Exclude damages in short fragments
	for (auto& chrMap : fDamageMap)
	{
		G4int iChr = chrMap.first;
		// Get maximum bp possible to avoid segmentation fault when not using actual chromosomes
		G4int maxBp = 0;
		for (auto& bpMap : fDamageMap[iChr]) { if (bpMap.first > maxBp) maxBp = bpMap.first; }

		G4int inviableFragments = 0;
		std::vector<G4int> allFreeEndsOfThisChromosome = initialBasePairsInEachChromosome[iChr];
		for (unsigned int i = 1; i < allFreeEndsOfThisChromosome.size(); i++)
		{
			G4int fragmentSize = std::min(maxBp, allFreeEndsOfThisChromosome[i] - allFreeEndsOfThisChromosome[i-1]);
			G4bool isInviableFragment = (fragmentSize < fLowerThresholdForFragmentDetection || fragmentSize > fUpperThresholdForFragmentDetection);
			if (isInviableFragment)
			{
				inviableFragments++;
				// Exclude damage to bases, SSB and DSB on the short fragments
				for (G4int ibp = allFreeEndsOfThisChromosome[i-1]; ibp < allFreeEndsOfThisChromosome[i]; ibp++)
				{
					if (ibp <= maxBp)
					{
						if (fBDMap[iChr][ibp][1] == 1) { fBDMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][0] = 0; }
						if (fBDMap[iChr][ibp][1] >= 2) { fBDMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][0] = 0; }
						if (fBDMap[iChr][ibp][2] == 1) { fBDMap[iChr][ibp][2] = 0; fDamageMap[iChr][ibp][1] = 0; }
						if (fBDMap[iChr][ibp][2] >= 2) { fBDMap[iChr][ibp][2] = 0; fDamageMap[iChr][ibp][1] = 0; }
						if (fSSBMap[iChr][ibp][1] == 1) { fSSBMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][2] = 0; }
						if (fSSBMap[iChr][ibp][1] >= 2) { fSSBMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][2] = 0; }
						if (fSSBMap[iChr][ibp][2] == 1) { fSSBMap[iChr][ibp][2] = 0; fDamageMap[iChr][ibp][3] = 0; }
						if (fSSBMap[iChr][ibp][2] >= 2) { fSSBMap[iChr][ibp][2] = 0; fDamageMap[iChr][ibp][3] = 0; }
						// We use the break in the first break so we only check for damage in that to count DSBs. However, damage in second breaks are also disregarded. Also avoid first base pair
						if (ibp > allFreeEndsOfThisChromosome[i-1])
						{
							if (fDSBMap[iChr][ibp][1] == 1) { fDSBMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][2] = 0; }
							if (fDSBMap[iChr][ibp][1] == 2) { fDSBMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][2] = 0; }
							if (fDSBMap[iChr][ibp][1] >= 3) { fDSBMap[iChr][ibp][1] = 0; fDamageMap[iChr][ibp][2] = 0; }
						}
					}
				}
			}
		}

	}
}

G4bool TsDNADamageCalculator::CauseDirectDamage(G4double edep)
{
	if (!fUseLinearProbabilityForDirectDamage)
	{
		if (edep > fDirectDamageThreshold)
			return true;
	}
	else
	{
		if (edep >= fUpperLimitLinearProbability)
			return true;
		else if (edep <= fLowerLimitLinearProbability)
			return false;
		else
		{
			G4double accept = (edep - fLowerLimitLinearProbability) / (fUpperLimitLinearProbability - fLowerLimitLinearProbability);
			if (G4UniformRand() < accept)
				return true;
		}
	}
	return false;
}

void TsDNADamageCalculator::OutputSDDHeader(G4bool minimalSDD, G4String primaryParticle, G4double energy, G4double dosePerExposure, std::vector<G4int> chromosomeContents, G4bool scoreIndirectDamages,
		G4bool scoreOnBases, G4String author, G4String simDetails, G4String sourceDetails, G4int sourceType, G4String energyDist, G4String irrTarget, G4String cellCycle, G4String DNAStructure,
		G4int inVitroOrInVivo, G4String prolStatus, G4String microenv, G4double time, G4String addInfo)
{
	// Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al.
	// A New Standard DNA Damage (SDD) Data Format. Radiat. Res. 191, 76-92 (2019)
	G4String dataEntries = "1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0";
	if (minimalSDD) dataEntries = "1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0";

	G4String SDDVersion = "v1.0";
	G4String incParticle;
	if (primaryParticle == "proton") incParticle = "2212";
	else if (primaryParticle == "alpha") incParticle = "100002004";
	G4String meanEnergy = (G4String)std::to_string(energy) + " MeV";

	G4String scoringIndirect = "0";
	if (scoreIndirectDamages) scoringIndirect = "1";
	G4String scoringOnBases = "-1";
	if (scoreOnBases) scoringOnBases = (G4String)std::to_string(fNumberOfBasePairForDSB);

	std::ofstream outFile;
	outFile.open(fOutputFileName + "_sdd.txt");
	outFile << "SDD Version, SDD" << SDDVersion << ";\n";
	outFile << "Software, TOPAS-nBio;\n";
	outFile << "Author, " << author << ";\n";
	outFile << "Simulation Details, " << simDetails << ";\n";
	outFile << "Source, " << sourceDetails << ";\n";
	outFile << "Source type, " << std::to_string(sourceType) << ";\n";
	outFile << "Incident particles, " << incParticle << ";\n";
	outFile << "Mean particle energy, " << meanEnergy << ";\n";
	outFile << "Energy distribution, " << energyDist << ";\n";
	outFile << "Particle fraction, 1.0;\n";
	outFile << "Dose or fluence, 1, " << dosePerExposure << ";\n";
	outFile << "Dose rate, 0.0" << ";\n";
	outFile << "Irradiation target, " << irrTarget << ";\n";
	outFile << "Volumes, 0,5,5,5,0,0,0,1,4.65,4.65,4.65,0,0,0;\n";
	outFile << "Chromosome sizes, " << chromosomeContents.size() << ",";
	G4double totalDNAContent = 0;
	for (unsigned int i = 0; i < chromosomeContents.size(); i++)
	{
		outFile << (G4double)chromosomeContents[i] / 1e6;
		totalDNAContent += chromosomeContents[i];
		if (i < chromosomeContents.size() - 1) outFile << ",";
		else outFile << ";\n";
	}
	outFile << "DNA Density, " << totalDNAContent / 1e6 / 421.15 << ";\n";
	outFile << "Cell Cycle Phase, " << cellCycle << ";\n";
	outFile << "DNA Structure, " << DNAStructure << ";\n";
	outFile << "In vitro / in vivo, " << std::to_string(inVitroOrInVivo) << ";\n";
	outFile << "Proliferation status, " << prolStatus << ";\n";
	outFile << "Microenvironment, " << microenv << ";\n";
	outFile << "Damage definition, " << scoringIndirect << ", 0, " << fNumberOfBasePairForDSB << ", " << scoringOnBases << ", " << fDirectDamageThreshold/eV;
	if (fUseLinearProbabilityForDirectDamage) outFile << ", " << fUpperLimitLinearProbability;
	outFile << ";\n";
	outFile << "Time, " << std::to_string(time) << ";\n";
	outFile << "Damage and primary count, ;\n";
	outFile << "Data entries, " << dataEntries << ";\n";
	outFile << "Additional information, " << addInfo << ";\n";
	outFile << "***EndOfHeader***;\n";
	outFile.close();
}

G4int TsDNADamageCalculator::OutputSDDFile(std::map<G4int, std::vector<G4int>> damageSites, G4int eventID, G4int exposureID, std::vector<G4int> chromosomeContents)
{
	// Ref. Schuemann J, McNamara A L, Warmenhoven J W, et al. A New Standard DNA Damage (SDD) Data Format. Radiat. Res. 191, 76-92 (2019)
	G4int lastExposureID = -1000; G4int lastEventID = -1000;

	G4int numSites = 0;
	std::ofstream outFile;
	outFile.open(fOutputFileName + "_sdd.txt", std::ios::out | std::ios::app);
	for (auto& chrMap : damageSites)
	{
		G4int iChr = chrMap.first;
		std::vector<G4int> ibpsTakenForThisChromosome;
		for (unsigned int i = 0; i < damageSites[iChr].size(); i++)
		{
			G4int initialBpId = damageSites[iChr][i];
			// Determines number and type of damages in block
			G4int dir = 0; G4int indir = 0;
			G4int bd = 0; G4int sb = 0; G4double dsb = 0;
			for (G4int j = 0; j < fNumberOfBasePairForDSB; j++)
			{
				 // Avoiding using same damage in two sites
				if (std::find(ibpsTakenForThisChromosome.begin(), ibpsTakenForThisChromosome.end(), initialBpId + j) == ibpsTakenForThisChromosome.end())
				{
						if ((fSSBMap[iChr][initialBpId + j][1] == direct || fSSBMap[iChr][initialBpId + j][1] == quasidirect))	{ dir++; sb++; }
						if ((fSSBMap[iChr][initialBpId + j][2] == direct || fSSBMap[iChr][initialBpId + j][2] == quasidirect))	{ dir++; sb++; }
						if (fSSBMap[iChr][initialBpId + j][1] == indirect)														{ indir++; sb++; }
						if (fSSBMap[iChr][initialBpId + j][2] == indirect)														{ indir++; sb++; }
						if ((fDSBMap[iChr][initialBpId + j][1] == direct || fDSBMap[iChr][initialBpId + j][1] == quasidirect))	{ dir++; sb++; dsb++; }
						if ((fDSBMap[iChr][initialBpId + j][2] == direct || fDSBMap[iChr][initialBpId + j][2] == quasidirect))	{ dir++; sb++; }
						if (fDSBMap[iChr][initialBpId + j][1] == indirect)														{ indir++; sb++; dsb++; }
						if (fDSBMap[iChr][initialBpId + j][2] == indirect)														{ indir++; sb++; }
						if ((fBDMap[iChr][initialBpId + j][1] == direct || fBDMap[iChr][initialBpId + j][1] == quasidirect))	{ dir++; bd++; }
						if ((fBDMap[iChr][initialBpId + j][2] == direct || fBDMap[iChr][initialBpId + j][2] == quasidirect))	{ dir++; bd++; }
						if (fBDMap[iChr][initialBpId + j][1] == indirect)														{ indir++; bd++; }
						if (fBDMap[iChr][initialBpId + j][2] == indirect)														{ indir++; bd++; }
				}
				ibpsTakenForThisChromosome.push_back(initialBpId+j);
			}
			numSites++;
			// Field 1: Determines exposure status
			G4int newExposureFlag = 0;
			if (lastEventID != eventID)			{ lastEventID = eventID; newExposureFlag = 1; }
			if (lastExposureID != exposureID)	{ lastExposureID = exposureID; newExposureFlag = 2; }
			outFile << newExposureFlag << ", " << eventID << "; ";

			// Field 2: Adds damage positions (for backbone)
			std::vector<G4ThreeVector> damagePositions;
			for (G4int j = 0; j < fNumberOfBasePairForDSB; j++)
			{
				if (fDamageMap[iChr][initialBpId + j][0] > 0)	damagePositions.push_back(fDamagePositions[iChr][initialBpId + j][0]);
				if (fDamageMap[iChr][initialBpId + j][1] > 0)	damagePositions.push_back(fDamagePositions[iChr][initialBpId + j][1]);
				if (fDamageMap[iChr][initialBpId + j][2] > 0)	damagePositions.push_back(fDamagePositions[iChr][initialBpId + j][2]);
				if (fDamageMap[iChr][initialBpId + j][3] > 0)	damagePositions.push_back(fDamagePositions[iChr][initialBpId + j][3]);
			}
			G4ThreeVector center;
			if (damagePositions.size() > 0)
			{
				std::vector<G4ThreeVector> damageCenterMaxMin = GetDamageCenterAndBoundaries(damagePositions);
				center = damageCenterMaxMin[0];
				outFile << center.x() / um << ", " << center.y() / um << ", " << center.z() / um;
				if (!fMinimalModeForSDD)
				{
					G4ThreeVector max = damageCenterMaxMin[1];
					G4ThreeVector min = damageCenterMaxMin[2];
					outFile << " / " << max.x() / um << ", " << max.y() / um << ", " << max.z() / um;
					outFile << " / " << min.x() / um << ", " << min.y() / um << ", " << min.z() / um;
				}
			}

			outFile << "; ";
			if (!fMinimalModeForSDD)
			{
				// Field 3: Nucleus considered as heterochromatin. Chromatied number set to 1 (unduplicated), chromosme arm set to 0 (short)
				outFile << "1, " << iChr << ", 1, 0; ";

				// Field 4: Chromosome length
				G4int chromosomeLength = chromosomeContents[iChr-1];
				G4double damageChromPos = (G4double) initialBpId / chromosomeLength;
				outFile.precision(12);
				outFile << std::ios::fixed << damageChromPos << "; ";
				outFile.precision(6);

				// Field 5. Type damage
				G4int typeDamage = 0;
				if (dir == 0 && indir > 0) typeDamage = 1;
				if (dir > 0 && indir > 0) typeDamage = 2;
				outFile << typeDamage << ", " << dir << ", " << indir << "; ";
			}
			// Field 6. Number damages in site
			outFile << bd << ", " << sb << ", " << dsb << "; ";
			if (!fMinimalModeForSDD)
			{
				G4String damageSpec;
				// Field 7. Damage specification
				for (G4int j = 0; j < fNumberOfBasePairForDSB; j++)
				{
					// Starting with bases of strand 1
					if (fDamageMap[iChr][initialBpId + j][0] != 0)
					{
						G4int typeDamage = fDamageMap[iChr][initialBpId + j][0];
						if (typeDamage == nodamage) typeDamage = 0;
						if (typeDamage == 4) typeDamage = 1;
						if (typeDamage == 5) typeDamage = 3;
						damageSpec += "1, " + (G4String)std::to_string(j+1) + ", " + (G4String)std::to_string(typeDamage) + " / ";
					}
				}
				for (G4int j = 0; j < fNumberOfBasePairForDSB; j++)
				{
					// Backbones of strand 1
					if (fDamageMap[iChr][initialBpId + j][2] != 0)
					{
						G4int typeDamage = fDamageMap[iChr][initialBpId + j][2];
						if (typeDamage == nodamage) typeDamage = 0;
						if (typeDamage == 4) typeDamage = 1;
						if (typeDamage == 5) typeDamage = 3;
						damageSpec += "2, " + (G4String)std::to_string(j+1) + ", " + (G4String)std::to_string(typeDamage) + " / ";
					}
				}
				for (G4int j = 0; j < fNumberOfBasePairForDSB; j++)
				{
					// Backbones of strand 2
					if (fDamageMap[iChr][initialBpId + j][3] != 0)
					{
						G4int typeDamage = fDamageMap[iChr][initialBpId + j][3];
						if (typeDamage == nodamage) typeDamage = 0;
						if (typeDamage == 4) typeDamage = 1;
						if (typeDamage == 5) typeDamage = 3;
						damageSpec += "3, " + (G4String)std::to_string(j+1) + ", " + (G4String)std::to_string(typeDamage) + " / ";
					}
				}
				for (G4int j = 0; j < fNumberOfBasePairForDSB; j++)
				{
					// Bases of strand 2
					if (fDamageMap[iChr][initialBpId + j][1] != 0)
					{
						G4int typeDamage = fDamageMap[iChr][initialBpId + j][1];
						if (typeDamage == nodamage) typeDamage = 0;
						if (typeDamage == 4) typeDamage = 1;
						if (typeDamage == 5) typeDamage = 3;
						damageSpec += "4, " + (G4String)std::to_string(j+1) + ", " + (G4String)std::to_string(typeDamage) + " / ";
					}
				}
				damageSpec = damageSpec.substr(0, damageSpec.size()-2);
				outFile << damageSpec << ";";
			}
			outFile << "\n";
		}
	}
	outFile.close();
	return numSites;
}

std::vector<G4ThreeVector> TsDNADamageCalculator::GetDamageCenterAndBoundaries(std::vector<G4ThreeVector> positions)
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

void TsDNADamageCalculator::WriteDNADamageCSV()
{
	G4String damageFileName = fOutputFileName + "_full.csv";
	std::fstream out;
	// Creates a new file, appending results depending on whether overwrite is the mode
	if (fOutputMode == "overwrite" && fEventID == 0)
		out.open(damageFileName, std::ios::out | std::ios::trunc);
	else
		out.open(damageFileName, std::ios::out | std::ios::app);

	// Inserts headers
	if (fEventID == 0)
		out << "Event,Chromosome,Bp,Backbone1,Base1,Base2,Backbone2\n";
	for (auto& chrMap : fDamageMap)
	{
		G4int iChr = chrMap.first;
		for (auto& bpMap : fDamageMap[iChr])
		{
			G4int iBp = bpMap.first;
			G4String base1 = "0", base2 = "0", back1 = "0", back2 = "0";
			if (fDamageMap[iChr][iBp][0] == direct || fDamageMap[iChr][iBp][0] == quasidirect)  base1 = "D";
			if (fDamageMap[iChr][iBp][0] == indirect) base1 = "I";
			if (fDamageMap[iChr][iBp][0] == multiple || fDamageMap[iChr][iBp][0] == multiplewithquasidirect) base1 = "M";
			if (fDamageMap[iChr][iBp][0] == nodamage) base1 = "*";
			if (fDamageMap[iChr][iBp][1] == direct || fDamageMap[iChr][iBp][1] == quasidirect) base2 = "D";
			if (fDamageMap[iChr][iBp][1] == indirect) base2 = "I";
			if (fDamageMap[iChr][iBp][1] == multiple || fDamageMap[iChr][iBp][1] == multiplewithquasidirect) base2 = "M";
			if (fDamageMap[iChr][iBp][1] == nodamage) base2 = "*";
			if (fDamageMap[iChr][iBp][2] == direct || fDamageMap[iChr][iBp][2] == quasidirect) back1 = "D";
			if (fDamageMap[iChr][iBp][2] == indirect) back1 = "I";
			if (fDamageMap[iChr][iBp][2] == multiple || fDamageMap[iChr][iBp][2] == multiplewithquasidirect) back1 = "M";
			if (fDamageMap[iChr][iBp][2] == nodamage) back1 = "*";
			if (fDamageMap[iChr][iBp][3] == direct || fDamageMap[iChr][iBp][3] == quasidirect) back2 = "D";
			if (fDamageMap[iChr][iBp][3] == indirect) back2 = "I";
			if (fDamageMap[iChr][iBp][3] == multiple || fDamageMap[iChr][iBp][3] == multiplewithquasidirect) back2 = "M";
			if (fDamageMap[iChr][iBp][3] == nodamage) back2 = "*";
			out << fEventID << "," << iChr << "," << iBp << "," << back1 << "," << base1 << "," << base2 << "," << back2 << "\n";
		}
	}
	out.close();
}

