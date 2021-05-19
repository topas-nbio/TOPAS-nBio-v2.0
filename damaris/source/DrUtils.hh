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
// Created by John-William Warmenhoven.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//

#pragma once

#include <G4Types.hh>
#include <G4String.hh>
#include <vector>

class DrUtils {
public:
    DrUtils();
    ~DrUtils();

    void ExtractBiologicallyRelevantParameters();

    void PrintToScreenBioParam();
    void PrintBiologicallyRelevantParameters();
    void PrintMSD(std::ofstream &file);
    void PrintFRAP(std::ofstream &file);
    void PrintRecruitment(std::ofstream &file);
    void PrintRawMolecules(std::ofstream &file);
    void PrintResiduals(std::ofstream &file);
    void PrintPerRepBiologicallyRelevantParameters();
    void PrintMRPercent();
    //@@@@ Bins and then prints out to file a vector of doubles, arguments are:
    //@@@@ Name of file, list to print, units to convert to (1 if no conversion)
    //@@@@ number of bins, start value of bins, end value of bins (if unsure of
    //@@@@ the end of range then start = 0, end = 0 will automatically calculate)
    //@@@@ last bool is to set append or not.
    void PrintBinList(G4String,std::vector<G4double>,G4double,G4int,G4double,G4double,G4bool);

    //used in #DEBUG function in DrEndRunForBiology
    std::vector<std::pair<G4double,G4double>> BinDoubleList(std::vector<G4double>,G4int,G4double,G4double);

private:
    unsigned fVerbose;
    void NormaliseMoleculeStore();

    G4double StdErrMeanMultiColumn(G4double, std::vector<std::vector<G4double> >,G4int);
    G4double StdErrMeanList(G4double, std::vector<G4double>);
    G4double StdErrMeanList(G4double, std::vector<G4int>);

};
