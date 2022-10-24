// Extra Class for use by DaMaRiS
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
#include "DrDefinitions.hh"
#include "TsParameterManager.hh"
#include "DrUtils.hh"
#include "DrDSBEnd_Generic.hh"
#include <G4Scheduler.hh>
#include <G4SystemOfUnits.hh>

using namespace std;

DrUtils::DrUtils() {}

DrUtils::~DrUtils() {}

void DrUtils::ExtractBiologicallyRelevantParameters() {

    NormaliseMoleculeStore();

    auto currentData = DrDefinitions::Instance()->fMoleculesRecord[1];

    if(DrDefinitions::Instance()->fMoleculesRecord[0].size() == 0){
        DrDefinitions::Instance()->fMoleculesRecord[0] = currentData;
    }
    else{
        for(auto timePoint: DrDefinitions::Instance()->fMoleculesRecord[1]){
            G4double time = timePoint.first;
            auto oldData = DrDefinitions::Instance()->fMoleculesRecord[0];
            for(auto molecule: timePoint.second){
                const G4MoleculeDefinition* molDef = molecule.first;

                //@@@@ NTS = Values are normalised to the maximum value of the series
                //@@@@ NTI = Values are normalised to the initial number of break(s/ ends).
                G4double oldAverage_NTS = oldData[time][molDef][0];
                G4double oldVariance_NTS = oldData[time][molDef][1];
                G4double oldNumberOfRunsToHere = oldData[time][molDef][2];
                G4double oldAverage_NTI = oldData[time][molDef][3];
                G4double oldVariance_NTI = oldData[time][molDef][4];

                G4double newValue_NTS = molecule.second[0];
                G4double newValue_NTI = molecule.second[3];

                G4double repeatCount = oldNumberOfRunsToHere + 1;

                G4double deltaOldAverage_NTS = newValue_NTS - oldAverage_NTS;
                G4double newAverage_NTS = oldAverage_NTS + (deltaOldAverage_NTS / repeatCount);
                G4double deltaNewAverage_NTS = newValue_NTS - newAverage_NTS;
                G4double newVariance_NTS = oldVariance_NTS + (deltaOldAverage_NTS*deltaNewAverage_NTS);

                G4double deltaOldAverage_NTI = newValue_NTI - oldAverage_NTI;
                G4double newAverage_NTI = oldAverage_NTI + (deltaOldAverage_NTI / repeatCount);
                G4double deltaNewAverage_NTI = newValue_NTI - newAverage_NTI;
                G4double newVariance_NTI = oldVariance_NTI + (deltaOldAverage_NTI*deltaNewAverage_NTI);

                DrDefinitions::Instance()->fMoleculesRecord[0][time][molDef][0] = newAverage_NTS;
                DrDefinitions::Instance()->fMoleculesRecord[0][time][molDef][1] = newVariance_NTS;
                DrDefinitions::Instance()->fMoleculesRecord[0][time][molDef][2] = repeatCount;
                DrDefinitions::Instance()->fMoleculesRecord[0][time][molDef][3] = newAverage_NTI;
                DrDefinitions::Instance()->fMoleculesRecord[0][time][molDef][4] = newVariance_NTI;
            }
        }
    }
    DrDefinitions::Instance()->fMoleculesRecord[1].clear();
}

void DrUtils::NormaliseMoleculeStore(){
    for (const auto& molecule: DrDefinitions::Instance()->GetNameMap()){
        G4double max{0};
        G4double initial{0};
        if(molecule.first.substr(0,6) == "DSBEnd" ||
           molecule.first.substr(0,12) == "DSB_Fixed_HR"){
            initial = DrDefinitions::Instance()->fInitialBreakNumber*2;
        }
        else if(molecule.first.substr(0,6) == "DSBSyn" ||
                molecule.first.substr(0,6) == "DSB_Fi"){
            initial = DrDefinitions::Instance()->fInitialBreakNumber;
        }

        const G4MoleculeDefinition* molDef = molecule.second;

        for (const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[1]){
            auto molCountMap = timePoint.second;
            if(molCountMap[molDef][0] > max) {
                max = molCountMap[molDef][0];
            }
        }
        for (const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[1]){
            if(max > 0){
                DrDefinitions::Instance()->fMoleculesRecord[1][timePoint.first][molDef][0] /= max;
            }
            else if(max == 0){
                DrDefinitions::Instance()->fMoleculesRecord[1][timePoint.first][molDef][0] = 0;
            }
            else if(max < 0){
                G4cout << "ERROR: max value of molecule " << molDef->GetName()
                        << " was calculated to be negative whilst trying to normalise."
                        << G4endl;
                DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
            }
            DrDefinitions::Instance()->fMoleculesRecord[1][timePoint.first][molDef][3] /= initial;
        }
    }
}

void DrUtils::PrintBiologicallyRelevantParameters() {

    //----------------------------------------------------------------------------
    //@@@@ The following writes out to a file the binned values for number of the
    //@@@@ specified states present in the simulation. Residual DSBs left is
    //@@@@ written out with the percentage incorrect and unrepaired DSBs in the
    //@@@@ title.
    G4int alternativeRunID = DrDefinitions::Instance()->GetAlternativeRunID();
    G4String outFileName = "PerRunResults" + to_string(alternativeRunID) + ".out";
    ofstream file(outFileName);
    file << "#Time |"
          << " Normalised to Self | SEM_NTS |"
          << " Normalised to Initial DSBs | SEM_NTI"
         << G4endl;

    //----------------------------------------------------------------------------
    //@@@@ Write accumulation data first so that it is always the same index
    PrintRecruitment(file);

    //----------------------------------------------------------------------------
    //@@@@ Write repair data next so that it is always the same index
    PrintResiduals(file);

    //----------------------------------------------------------------------------
    //@@@@ Write out the raw molecule numbers next
    PrintRawMolecules(file);

    //----------------------------------------------------------------------------
    //@@@@ Write out the FRAP data as PKcs - Bleached
    PrintFRAP(file);

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ The following writes out to file the separation at initiation of the
    //@@@@ two double strand break ends involved in an incorrectly repaired break
    //@@@@ site.
    G4double radius = DrDefinitions::Instance()->GetBoundingCellOrNucleusRadius();
    vector<pair<G4double, G4double>> list;
    file << "Original_Seperation_of_Misrepair" << G4endl;
    list = BinDoubleList(DrDefinitions::Instance()->fMisrepairSeparationsStore,
                         2500, 0, radius);
    file << "# Separation (nm) | Normalised Frequency" << G4endl;
    for (auto data : list) {
        file << data.first / nm << " " << data.second << G4endl;
    }
    file << G4endl << G4endl;
    //----------------------------------------------------------------------------
    file << "Time_of_Misrepair_Interactions" << G4endl;
    file << "# Time of Interaction (s) | Normalised Frequency" << G4endl;
    list = BinDoubleList(DrDefinitions::Instance()->fMisrepairTimeStore,
                         (G4int)(floor(G4Scheduler::Instance()->GetEndTime() / s)), 0,
                         G4Scheduler::Instance()->GetEndTime());
    for (auto data : list) {
        file << data.first / s << " " << data.second << G4endl;
    }
    file << G4endl << G4endl;

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ The following prints out the MSD of all DSB objects in the simulation
    file << "MSD_of_DSB" << G4endl;
    file << "# Time (s) | MSD" << G4endl;
    PrintMSD(file);
    file << G4endl << G4endl;

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ The following writes out the displacement of the particle at the end of
    //@@@@ run from it's original location at t=0
    file << "Displacement_at_End" << G4endl;
    file << "# Displacement (nm) | Normalised Frequency" << G4endl;
    list = BinDoubleList(DrDefinitions::Instance()->fDisplacementTrackingStore,
                         2500, 0, radius);
    for (auto data : list) {
        file << data.first / nm << " " << data.second << G4endl;
    }
    file << G4endl << G4endl;
    //----------------------------------------------------------------------------
    file << "Displacement_of_ResidualEnds" << G4endl;
    list = BinDoubleList(DrDefinitions::Instance()->fResidualEndDisplacementStore,
                         2500, 0, radius);
    for (auto data : list) {
        file << data.first / nm << " " << data.second << G4endl;
    }
    file << G4endl << G4endl;
    //----------------------------------------------------------------------------
    file << "Displacement_of_ResidualSyn" << G4endl;
    list = BinDoubleList(DrDefinitions::Instance()->fResidualSynDisplacementStore,
                         2500, 0, radius);
    for (auto data : list) {
        file << data.first / nm << " " << data.second << G4endl;
    }
    file << G4endl << G4endl;

    file.close();
    //----------------------------------------------------------------------------
}

void DrUtils::PrintPerRepBiologicallyRelevantParameters() {
    //----------------------------------------------------------------------------
    //@@@@ The following calculates the percentage of wrong and unrepaired double
    //@@@@ strand breaks at the end of each chemistry repeat.
    DrDefinitions* definitions = DrDefinitions::Instance();
    G4int RunID = definitions->fCurrentBiologyRepeatNumber;

    std::ofstream file;
    G4int altRunID = DrDefinitions::Instance()->GetAlternativeRunID();
    file.open("PerRepResults" + std::to_string(altRunID) + ".out",
              std::ios_base::app);
    file << definitions->fInitialBreakNumber << " "
            << definitions->fCheckBreakStore[RunID].fMisrepairCounter << " "
            << definitions->fCheckBreakStore[RunID].fNumberUnrepaired << " "
            << definitions->fCheckBreakStore[RunID].fInterChromosomeAberration << " "
            << definitions->fCheckBreakStore[RunID].fInterChromatidAberration << " "
            << definitions->fCheckBreakStore[RunID].fInterArmAberration << " "
            << definitions->fCheckBreakStore[RunID].fIntraArmAberration << " "
            << G4endl;
    file.close();

    //----------------------------------------------------------------------------
}

void DrUtils::PrintMRPercent() {
    auto initDSB = (G4double)DrDefinitions::Instance()->fInitialBreakNumber;
    G4double averagePercentMR{0.0};
    for (auto read : DrDefinitions::Instance()->fMisrepairNumberStore) {
        averagePercentMR += G4double(read) / initDSB;
    }
    averagePercentMR /= DrDefinitions::Instance()->GetBiologyRepeatNumber();

    G4double err = StdErrMeanList(averagePercentMR,
                                  DrDefinitions::Instance()->fMisrepairNumberStore);
    std::ofstream file;
    file.open("Misrepair.out", std::ios_base::app);
    file << DrDefinitions::Instance()->GetDSBTimeDelay() << " "
         << averagePercentMR << " " << err << G4endl;
    file.close();
}

void DrUtils::PrintBinList(G4String fileName, vector<G4double> inList,
                           G4double units, G4int numBins, G4double start,
                           G4double end, G4bool append) {
    G4int altRunID = DrDefinitions::Instance()->GetAlternativeRunID();
    ofstream file;
    if (append) {
        file.open(fileName, std::ios_base::app);
    } else {
        file.open(fileName + to_string(altRunID) + ".out");
    }
    vector<pair<G4double, G4double>> list =
            BinDoubleList(inList, numBins, start, end);
    for (auto read : list) {
        file << read.first / units << " " << read.second << G4endl;
    }
    file.close();
}

void DrUtils::PrintToScreenBioParam() {
    DrDefinitions* definitions = DrDefinitions::Instance();
    G4int RunID = definitions->fCurrentBiologyRepeatNumber;
    G4cout.precision(3);
    G4cout << "Initial number of DSBs: " << definitions->fInitialBreakNumber << G4endl
           << "Misrepaired: " << definitions->fCheckBreakStore[RunID].fMisrepairCounter << G4endl
           << "Chromosome Aberrations: "
           << definitions->fCheckBreakStore[RunID].fInterChromosomeAberration << " / "
           << definitions->fCheckBreakStore[RunID].fInterChromatidAberration << " / "
           << definitions->fCheckBreakStore[RunID].fInterArmAberration << " / "
            << definitions->fCheckBreakStore[RunID].fIntraArmAberration
           << G4endl << "Still unrepaired: " << definitions->fCheckBreakStore[RunID].fNumberUnrepaired
           << G4endl << "Molecules left in simulation at end: "
           << definitions->fCheckBreakStore[RunID].fNumberBreakMoleculeLeft << G4endl;
}

G4double DrUtils::StdErrMeanMultiColumn(G4double average,
                                        vector<vector<G4double>> list,
                                        G4int column) {
    G4int numberOfRuns = DrDefinitions::Instance()->GetBiologyRepeatNumber();
    G4double sum{0.0};
    for (auto data : list) {
        sum += pow(data[column] - average, 2);
    }
    return sqrt(sum / (numberOfRuns - 1)) / sqrt(numberOfRuns);
}

G4double DrUtils::StdErrMeanList(G4double average, vector<G4double> list) {
    G4int numberOfRuns = DrDefinitions::Instance()->GetBiologyRepeatNumber();
    G4double sum{0.0};
    for (auto data : list) {
        sum += pow(data - average, 2);
    }
    return sqrt(sum / (numberOfRuns - 1)) / sqrt(numberOfRuns);
}

G4double DrUtils::StdErrMeanList(G4double average, vector<G4int> list) {
    G4int numberOfRuns = DrDefinitions::Instance()->GetBiologyRepeatNumber();
    G4double sum{0.0};
    for (auto data : list) {
        sum += pow(G4double(data) - average, 2);
    }
    return sqrt(sum / (numberOfRuns - 1)) / sqrt(numberOfRuns);
}

vector<pair<G4double, G4double>> DrUtils::BinDoubleList(vector<G4double> list,
                                                        G4int numBins,
                                                        G4double start,
                                                        G4double end) {

    G4double max{DBL_MIN};
    G4double sizeBins{0.0};
    if (start == 0 && end == 0) {
        for (auto read : list) {
            if (read > max)
                max = read;
        }
        sizeBins = (max * 1.01) / (G4double)numBins;
    } else {
        sizeBins = (end - start) / (G4double)numBins;
    }
    map<G4int, G4double> binValues;
    for (G4int i{0}; i < numBins + 1; i++) {
        binValues[i] = sizeBins * i;
    }
    map<G4int, G4double> binnedData;
    for (auto toBin : list) {
        for (auto bin : binValues) {
            if (toBin <= bin.second) {
                binnedData[bin.first]++;
                break;
            }
        }
    }
    max = DBL_MIN;
    for (auto read : binnedData) {
        if (read.second >= max)
            max = read.second;
    }
    for (auto toNorm : binnedData) {
        binnedData[toNorm.first] /= max;
    }
    vector<pair<G4double, G4double>> tempVec;
    for (auto read : binnedData) {
        tempVec.push_back(make_pair(binValues[read.first], read.second));
    }
    return tempVec;
}

void DrUtils::PrintMSD(std::ofstream& file) {

    map<G4double, vector<G4double>> displacementAtTimeStore =
            DrDefinitions::Instance()->fMSDTrackingStore;

    for (const auto& timePoint : displacementAtTimeStore) {

        auto listSize = G4double(timePoint.second.size());
        //@@@@ SAFETY
        if (listSize <= 0) {
            G4cerr << "(DrUtils::PrintMSD) ERROR: no displacements recorded for this"
                   << " time step!" << G4endl << "Time step: " << timePoint.first
                   << G4endl;
            DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
        }
        G4double sqrd_displacement{0.};
        G4double sum_sqrd_displacement{0.};
        G4double MSD{0.};
        G4double residuals{0.};
        G4double stdev{0.};
        G4double SEM{0.};

        for (auto displacement : timePoint.second) {
            //@@@@ SAFETY
            if (displacement == 0.)
                sqrd_displacement = 0.;
            else
                sqrd_displacement = displacement * displacement;
            sum_sqrd_displacement += sqrd_displacement;
        }

        MSD = sum_sqrd_displacement / listSize;

        for (auto displacement : timePoint.second) {
            //@@@@ SAFETY
            if (displacement == 0.)
                sqrd_displacement = 0.;
            else
                sqrd_displacement = displacement * displacement;
            residuals += (sqrd_displacement - MSD) * (sqrd_displacement - MSD);
        }
        //@@@@ SAFETY
        if (listSize == 1)
            stdev = MSD;
        else
            stdev = sqrt(residuals / (listSize - 1));
        SEM = stdev / sqrt(listSize);

        file << timePoint.first << " " << MSD << " " << SEM << G4endl;
    }
}

void DrUtils::PrintFRAP(std::ofstream &file) {

    DrDefinitions* definitions = DrDefinitions::Instance();
    G4int repeats = DrDefinitions::Instance()->GetBiologyRepeatNumber();
    std::map<G4double,std::pair<G4double,G4double>> storeBleached;
    std::map<G4double,std::pair<G4double,G4double>> storePKcs;
    std::map<G4double,std::pair<G4double,G4double>> storeFRAP;

    for(const auto& read: definitions->fBleachedStoreRun){
        G4double average{0.0};
        for(const auto& read2: read.second){
            average += read2;
        }
        if(!read.second.empty()) average /= read.second.size();
        else average = 0;
        G4double SEM = StdErrMeanList(average,read.second);
        storeBleached[read.first] = std::make_pair(average,SEM);
    }

    for(const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[0]){
        G4double time = timePoint.first;
        auto data = timePoint.second;

        G4double average{0.};
        G4double SEM{0.};
        for(const auto& name: DrDefinitions::Instance()->GetHasPKcs()){
            auto nameMap = DrDefinitions::Instance()->GetNameMap();
            if(nameMap.find(name) != nameMap.end()){
                auto mol = DrDSBEnd_Generic::Definition(name);
                G4double temp_average = data[mol][3];
                G4double temp_residual = data[mol][4];
                G4double temp_SEM = sqrt(temp_residual/(repeats - 1)) / sqrt(repeats);
                average += temp_average;
                SEM += temp_SEM*temp_SEM;
            }
        }
        SEM = sqrt(SEM);
        storePKcs[time] = make_pair(average,SEM);
    }

    file << "Time || Bleached_RelInit SEM_RelInit" << G4endl;
    for(const auto& read: storeBleached){
        G4double time = read.first;
        G4double cutTime = time - 30;
        if(cutTime > 0.0){
            G4double combinedAverage = storePKcs[time].first - storeBleached[time].first;
            if(combinedAverage < 0.0) combinedAverage = 0.0;
            G4double combinedSEM = sqrt(storePKcs[time].second*storePKcs[time].second
                                        + storeBleached[time].second*storeBleached[time].second);
            storeFRAP[cutTime] = std::make_pair(combinedAverage,combinedSEM);
        }
    }

    G4double max{0.0};
    for(const auto& read: storeFRAP){
        if (read.second.first > max) max = read.second.first;
    }

    for(const auto& read: storeFRAP){
        file << read.first << " " << read.second.first/max << " " << read.second.second/max << G4endl;
    }
    file << G4endl << G4endl;

}

void DrUtils::PrintRecruitment(std::ofstream &file) {

    G4int repeats = DrDefinitions::Instance()->GetBiologyRepeatNumber();

    //Storage Structures
    std::map<G4double,std::pair<G4double,G4double> > KuData;
    std::map<G4double,std::pair<G4double,G4double> > PKcsData;
    std::map<G4double,std::pair<G4double,G4double> > XRCC4Data;

    //-------------------------------------------------------------------------------------
    // Ku70/80
    file << "Ku_Recruitment" << G4endl;
    for(const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[0]){
        G4double time = timePoint.first;
        auto data = timePoint.second;
        G4double average{0.};
        G4double SEM{0.};
        for(const auto& name: DrDefinitions::Instance()->GetHasKu()){
            auto nameMap = DrDefinitions::Instance()->GetNameMap();
            if(nameMap.find(name) != nameMap.end()){
                auto mol = DrDSBEnd_Generic::Definition(name);
                G4double temp_average = data[mol][3];
                G4double temp_residual = data[mol][4];
                G4double temp_SEM = sqrt(temp_residual/(repeats - 1)) / sqrt(repeats);
                average += temp_average;
                SEM += temp_SEM*temp_SEM;
            }
        }
        SEM = sqrt(SEM);
        KuData[time] = make_pair(average,SEM);
        if(time >= 60) break;
    }
    // normalise
    G4double max{0.0};
    for(const auto& read: KuData){
        if(read.second.first > max) max = read.second.first;
    }
    // write
    for(const auto& read: KuData){
        file << read.first/60.0 << " " << read.second.first/max << " " << read.second.second/max << G4endl;
    }
    file << G4endl << G4endl;
    //-------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------
    // DNA-PKcs
    file << "PKcs_Recruitment" << G4endl;
    for(const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[0]){
        G4double time = timePoint.first;
        auto data = timePoint.second;
        G4double average{0.};
        G4double SEM{0.};
        for(const auto& name: DrDefinitions::Instance()->GetHasPKcs()){
            auto nameMap = DrDefinitions::Instance()->GetNameMap();
            if(nameMap.find(name) != nameMap.end()){
                auto mol = DrDSBEnd_Generic::Definition(name);
                G4double temp_average = data[mol][3];
                G4double temp_residual = data[mol][4];
                G4double temp_SEM = sqrt(temp_residual/(repeats - 1)) / sqrt(repeats);
                average += temp_average;
                SEM += temp_SEM*temp_SEM;
            }
        }
        SEM = sqrt(SEM);
        PKcsData[time] = make_pair(average,SEM);
        if(time >= 60) break;
    }
    // normalise
    max = 0.0;
    for(const auto& read: PKcsData){
        if(read.second.first > max) max = read.second.first;
    }
    // write
    for(const auto& read: PKcsData){
        file << read.first/60.0 << " " << read.second.first/max << " " << read.second.second/max << G4endl;
    }
    file << G4endl << G4endl;
    //-------------------------------------------------------------------------------------

    //-------------------------------------------------------------------------------------
    // XRCC4
    file << "XRCC4_Recruitment" << G4endl;
    for(const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[0]){
        G4double time = timePoint.first;
        auto data = timePoint.second;
        G4double average{0.};
        G4double SEM{0.};
        for(const auto& name: DrDefinitions::Instance()->GetHasXRCC4()){
            auto nameMap = DrDefinitions::Instance()->GetNameMap();
            if(nameMap.find(name) != nameMap.end()){
                auto mol = DrDSBEnd_Generic::Definition(name);
                G4double temp_average = data[mol][3];
                G4double temp_residual = data[mol][4];
                G4double temp_SEM = sqrt(temp_residual/(repeats - 1)) / sqrt(repeats);
                average += temp_average;
                SEM += temp_SEM*temp_SEM;
            }
        }
        SEM = sqrt(SEM);
        XRCC4Data[time] = make_pair(average,SEM);
        if(time >= 300) break;
    }
    // normalise
    max = 0.0;
    for(const auto& read: XRCC4Data){
        if(read.second.first > max) max = read.second.first;
    }
    // write
    for(const auto& read: XRCC4Data){
        file << read.first/60.0 << " " << read.second.first/max << " " << read.second.second/max << G4endl;
    }
    file << G4endl << G4endl;
    //-------------------------------------------------------------------------------------
}

void DrUtils::PrintRawMolecules(std::ofstream &file){

    G4int repeats = DrDefinitions::Instance()->GetBiologyRepeatNumber();

    for(const auto& molList: DrDefinitions::Instance()->GetNameMap()){
        const G4MoleculeDefinition* molDef = molList.second;
        if(molList.first.substr(0,9) == "DSB_Fixed"){
            file << "DSB_Residual" << G4endl;
        }
        else{
            file << molList.first << G4endl;
        }

        for(const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[0]){
            G4double time = timePoint.first;
            auto data = timePoint.second;

            G4double average_NTS = data[molDef][0];
            G4double residuals_NTS = data[molDef][1];
            G4double numberOfRepeats = data[molDef][2];
            G4double average_NTI = data[molDef][3];
            G4double residuals_NTI = data[molDef][4];

            G4double SEM_NTS = sqrt(residuals_NTS / (numberOfRepeats - 1)) / sqrt(repeats);
            G4double SEM_NTI = sqrt(residuals_NTI / (numberOfRepeats - 1)) / sqrt(repeats);

            if(molList.first.substr(0,9) == "DSB_Fixed"){
                file << time / 60 << " "
                     << 1 - average_NTS << " "
                     << SEM_NTS << " "
                     << 1 - average_NTI << " "
                     << SEM_NTI << G4endl;
            }
            else{
                file << time / 60 << " "
                     << average_NTS << " "
                     << SEM_NTS << " "
                     << average_NTI << " "
                     << SEM_NTI << G4endl;
            }
        }
        file << G4endl << G4endl;
    }
}

void DrUtils::PrintResiduals(std::ofstream &file) {

    G4int repeats = DrDefinitions::Instance()->GetBiologyRepeatNumber();

    G4double numberOfRepeats{0.};
    G4double average_NTI{0.};
    G4double residuals_NTI{0.};
    G4double SEM_NTI{0.};

    std::map<G4double, std::pair<G4double, G4double> > tempStore;

    file << "DSB_Residual" << G4endl;

    for(const auto& molList: DrDefinitions::Instance()->GetNameMap()){
        const G4MoleculeDefinition* molDef = molList.second;
        if(molList.first.substr(0,9) == "DSB_Fixed"){
            for(const auto& timePoint: DrDefinitions::Instance()->fMoleculesRecord[0]){
                G4double time = timePoint.first; // time is in seconds
                auto data = timePoint.second;

                numberOfRepeats = data[molDef][2];
                average_NTI = data[molDef][3];
                residuals_NTI = data[molDef][4];

                SEM_NTI = sqrt(residuals_NTI / (numberOfRepeats - 1)) / sqrt(repeats);

                if(tempStore.find(time) != tempStore.end()){
                    tempStore[time].first += 1 - average_NTI;
                    G4double pSEM = tempStore[time].second;
                    tempStore[time].second = sqrt(pSEM*pSEM+SEM_NTI*SEM_NTI);
                }
                else{
                    tempStore[time].first = 1- average_NTI;
                    tempStore[time].second = SEM_NTI;
                }
            }
        }
    }

    G4double relVal = tempStore[24*60*60].first;
    G4double absVal{0.};
    G4double initVal{0.};

    if(relVal != 0.){
        for(auto read: DrDefinitions::Instance()->fCheckBreakStore){
            absVal += read.second.fNumberUnrepaired;
        }
        absVal /= repeats;
        initVal = absVal/relVal;
    }
    else{
        initVal = 1.;
    }

    for(const auto& read: tempStore){
        file << read.first/60 << " "
             << read.second.first * initVal << " "
             << read.second.second * initVal << " "
             << read.second.first << " "
             << read.second.second << G4endl;
    }
    file << G4endl << G4endl;
}
