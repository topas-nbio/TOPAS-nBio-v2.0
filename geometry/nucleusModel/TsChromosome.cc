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
// Authors: Hongyu Zhu, Jan Schuemann

#include "TsChromosome.hh"

 TsChromosome:: TsChromosome()
{

    fDNAcontentofEachChromosome_BP = {252823200,  252823200,  248157000,  248157000,  204040200,  204040200,  
                                    195556200,  195556200,  184951200,  184951200,  174770400,  174770400,  
                                    162468600,  162468600,  149318400,  149318400,  143379600,  143379600,  
                                    138289200,  138289200,  137440800,  137440800,  135319800,  135319800,  
                                    116655000,  116655000,  108595200,  108595200,  102656400,  102656400,  
                                    90778800, 90778800,  80173800,  80173800,  77628600,  77628600,  
                                    65326800,  65326800,  63630000,  63630000,  47934600,  47934600, 
                                    50479800,  50479800,  58963800,  158226600};
    
    fTotalDNAContent =0;
    for(unsigned i=0; i<fDNAcontentofEachChromosome_BP.size(); i++ )
    {
        fTotalDNAContent += (G4double) fDNAcontentofEachChromosome_BP[i]/1E6;
        std::vector<G4int> oneMat = {0, fDNAcontentofEachChromosome_BP[i]};
        fChromosomeMatrix.push_back(oneMat);
    }
}

    
TsChromosome::~ TsChromosome()
{}


std::vector<std::vector<G4int> > TsChromosome::SplitChromosome(std::vector <std::pair <TsHitInDNA*, TsHitInDNA*>> DSB_pairs)
{
    std::vector<std::vector<G4int> > SplitChromo = GetChromosomeMatrix();
    for(unsigned i =0; i<SplitChromo.size(); i++)
    {
        G4int ChromosomeID = i+1;
        for (unsigned iter =0; iter<DSB_pairs.size(); iter++)
        {
            G4int DamagedChromosomeID = DSB_pairs[iter].first->GetChromosomeID();
            G4int DamagedBasepairID   = DSB_pairs[iter].first->GetBasePairID();
            if( ChromosomeID == DamagedChromosomeID )
                SplitChromo[i].push_back( DamagedBasepairID);
        } 
        sort(SplitChromo[i].begin(), SplitChromo[i].end());
    }

    return SplitChromo;
}

std::vector<std::vector<G4int>> TsChromosome::SplitChromosome(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs)
{
	std::vector<std::vector<G4int> > SplitChromo = GetChromosomeMatrix();
	for (unsigned int i = 0; i < SplitChromo.size(); i++)
	{
		for (auto& chrom : DSBPairs)
		{
			unsigned int iChr = chrom.first;
			for (auto& pairOfBPDSB : DSBPairs[iChr])
			{
				G4int initialDamagedBp = pairOfBPDSB.first.first;
				if (i+1 == iChr)
					SplitChromo[i].push_back(initialDamagedBp);
			}
		}
		sort(SplitChromo[i].begin(), SplitChromo[i].end());
	}
	return SplitChromo;
}

G4int TsChromosome::CountDNAFrangmentsWithSize( std::vector<std::vector<G4int> > DNAfragments, G4int LowerFragmentSizeThreshold, G4int UpperFragmentSizeThreshold)
{
    // LowerFragmentSizeThreshold, UpperFragmentSizeThreshold  in the unit of BP
    G4int NumberOfInterestedFragments = 0;
    for(unsigned i =0; i<DNAfragments.size(); i++)
    {
        G4double FrangmentSize =0;
        std::vector<G4int> AllFreeEndofChromosome = DNAfragments[i];
        for(unsigned j =1; j<AllFreeEndofChromosome.size(); j++)
        {
            FrangmentSize = AllFreeEndofChromosome[j] - AllFreeEndofChromosome[j-1];
            if(FrangmentSize>=LowerFragmentSizeThreshold && FrangmentSize<=UpperFragmentSizeThreshold )
                NumberOfInterestedFragments++; 
        }
    }

    return NumberOfInterestedFragments;
}

