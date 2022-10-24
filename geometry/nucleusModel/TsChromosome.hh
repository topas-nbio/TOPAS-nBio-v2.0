// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//
// Authors: Hongyu Zhu

#ifndef  TsChromosome_hh
#define  TsChromosome_hh

#include "G4UnitsTable.hh"
#include "TsHitInDNA.hh"

class  TsChromosome 
{
public:
     TsChromosome();
    ~ TsChromosome();
    
    G4double GetTotalDNAContent_MBP() {return fTotalDNAContent;}
    std::vector<G4int> GetDNAcontentofEachChromosome_BP()  {return fDNAcontentofEachChromosome_BP;}
    std::vector<std::vector<G4int> > GetChromosomeMatrix() { return fChromosomeMatrix;}
    std::vector<std::vector<G4int> > SplitChromosome(std::vector <std::pair < TsHitInDNA*,  TsHitInDNA*>> DSB_pairs);
    std::vector<std::vector<G4int> > SplitChromosome(std::map<G4int, std::map<std::pair<G4int, G4int>, G4int>> DSBPairs);
    G4int CountDNAFrangmentsWithSize(std::vector<std::vector<G4int> > DNAfragments, G4int LowerThreshold, G4int UpperThreshold);
    
protected:
    
private:
    G4double fTotalDNAContent;
    std::vector<G4int> fDNAcontentofEachChromosome_BP; // In unit of BP
    std::vector<std::vector<G4int> > fChromosomeMatrix;     // ChromosomeID, 0, Chromosome DNA content(Mbp)
};

#endif
