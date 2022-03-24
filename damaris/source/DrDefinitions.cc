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
#include "DrDSBEnd_Generic.hh"
#include <G4MoleculeTable.hh>
#include "DrProteinKinetics_Generic.hh"
#include <regex>
#include <G4SystemOfUnits.hh>
#include <TsParameterManager.hh>

using namespace std;

DrDefinitions* DrDefinitions::definitionsInstance(nullptr);

DrDefinitions* DrDefinitions::Instance(){
    if (!definitionsInstance) definitionsInstance = new DrDefinitions;
    return definitionsInstance;
}

DrDefinitions::DrDefinitions():
        fVerbose(1), fPm(nullptr), fCurrentExplicitBinNuber(-1), fExplicitBinning("Null"), fDSBPattern({0}),
        fName("Default"), fObserveDurationForMSD(300*s),fObserveStepSizeForMSD(1*s), fTurnOffTime(true),
        fDSBOriginNumber(-1), fDSBOffset(0.0), fDSBColumnNumber(-1), fDSBColumnRadius(0), fDSBSeparation(-1.0),
        fDSBTimeDelay(0.0), fBiologyEndTime(-1.0), fIsSubDiffusion(true), fDiffusionCoefficientForJump(2.4e11*nm*nm/s),
        fDiffusionCoefficientForTrapped(0.0), fMinWaitingTime(1e9*ps), fIsReadFromFile(false), fSelectFromExposure(-1),
        fAlternativeRunID(1), fBiologyRepeatNumber(30), fBoundingRadius(-1.0){

    // Define which molecule names has which proteins attached.
    fHasPKcs = {"DSBEnd_PKcs",
                "DSBEnd_PKcs_XRCC4",
                "DSBSynaptic",
                "DSBSynaptic_Stable",
                "DSBSynaptic_PKcs",
                "DSBSynaptic_PKcs_XRCC4"};

    fHasKu = {"DSBEnd_Ku",
              "DSBEnd_PKcs",
              "DSBEnd_PKcs_XRCC4",
              "DSBEnd_XRCC4",
              "DSBSynaptic",
              "DSBSynaptic_Stable",
              "DSBSynaptic_PKcs",
              "DSBSynaptic_XRCC4",
              "DSBSynaptic_PKcs_XRCC4"};

    fHasXRCC4 = {"DSBEnd_XRCC4",
                 "DSBEnd_PKcs_XRCC4",
                 "DSBSynaptic_Stable",
                 "DSBSynaptic_XRCC4",
                 "DSBSynaptic_PKcs_XRCC4"};
}




DrDefinitions::~DrDefinitions()=default;

void DrDefinitions::Initialise(TsParameterManager* pM){
    fPm = pM;
}

void DrDefinitions::InitialiseParameters(){
    DrDefinitions::SetParameters();
}

void DrDefinitions::SetParameters(){

    //---------------------------------------------
    //  Diffusion Parameters
    //
    G4int selectDiffusion = fPm->GetIntegerParameter(GetFullParmName("DiffusionMode"));
    if(selectDiffusion == 0) fIsSubDiffusion = false;
    else if (selectDiffusion == 1) fIsSubDiffusion = true;
    else{
        G4cerr <<"WARNING: Improper selection of diffusion mode; "
               <<"assuming subdiffusive motion"
               <<G4endl;
        fIsSubDiffusion = true;
    }

    if(fIsSubDiffusion){
        fDiffusionCoefficientForJump = fPm->GetUnitlessParameter(GetFullParmName("DiffusionCoefficientForJump"))*nm*nm/s;
        fDiffusionCoefficientForTrapped = fPm->GetUnitlessParameter(GetFullParmName("DiffusionCoefficientForTrapped"))*nm*nm/s;
        fMinWaitingTime = fPm->GetDoubleParameter(GetFullParmName("MinWaitingTime"),"Time");
    }
    if(fPm->ParameterExists(GetFullParmName("ObserveDurationForMSD"))){
        fObserveDurationForMSD = fPm->GetDoubleParameter(GetFullParmName("ObserveDurationForMSD"),"Time");
    }
    if(fPm->ParameterExists(GetFullParmName("ObserveStepSizeForMSD"))){
        fObserveStepSizeForMSD = fPm->GetDoubleParameter(GetFullParmName("ObserveStepSizeForMSD"),"Time");
    }


    //---------------------------------------------
    // Simulation Parameters
    //
    fAlternativeRunID = fPm->GetIntegerParameter(GetFullParmName("AlternativeRunID"));
    fBiologyRepeatNumber = fPm->GetIntegerParameter(GetFullParmName("BiologyRepeatNumber"));

    if(fBiologyEndTime <= 0){
        fBiologyEndTime = fPm->GetDoubleParameter(GetFullParmName("DaMaRiSStageTimeEnd"), "Time");
        //@@@@ extends time by 0.1% to ensure last time point is
        //@@@@ captured by DrReportSystem
        fBiologyEndTime *= 1.001;
    }
    if(fPm->ParameterExists(GetFullParmName("TurnOffTime"))) {
        fTurnOffTime = fPm->GetBooleanParameter(GetFullParmName("TurnOffTime"));
    }
    if(fPm->ParameterExists(GetFullParmName("ExplicitBinningFileName"))){
        fExplicitBinning = fPm->GetStringParameter(GetFullParmName("ExplicitBinningFileName"));
    }

    //---------------------------------------------
    // Geometry Parameters
    //
    fBoundingRadius = fPm->GetDoubleParameter(GetFullParmName("BoundingCellOrNucleusRadius"),"Length");

    //---------------------------------------------
    // Damage Parameters
    //
    if(fPm->ParameterExists(GetFullParmName("STDFormatDamageFileName"))){
        fDamageFileName = fPm->GetStringParameter(GetFullParmName("STDFormatDamageFileName"));
        fIsReadFromFile = true;
        if(fPm->ParameterExists(GetFullParmName("SelectFromExposureNumber"))){
            fSelectFromExposure = fPm->GetIntegerParameter(GetFullParmName("SelectFromExposureNumber"));
        }
        else{
            fSelectFromExposure = -1;
        }

    }
    else if(fPm->ParameterExists(GetFullParmName("DSBOriginNumber"))) {
        fDSBOriginNumber = fPm->GetIntegerParameter(GetFullParmName("DSBOriginNumber"));
    }
    else if(fPm->ParameterExists(GetFullParmName("DSBOffset"))) {
        fDSBOffset = fPm->GetDoubleParameter(GetFullParmName("DSBOffset"), "Length");
    }
    else if(fPm->ParameterExists(GetFullParmName("DSBColumnNumber"))) {
        fDSBColumnNumber = fPm->GetIntegerParameter(GetFullParmName("DSBColumnNumber"));
        fDSBColumnRadius = fPm->GetDoubleParameter(GetFullParmName("DSBColumnRadius"), "Length");
    }
    else if(fPm->ParameterExists(GetFullParmName("DSBSeparation"))){
        fDSBSeparation = fPm->GetDoubleParameter(GetFullParmName("DSBSeparation"),"Length");
        if(fPm->ParameterExists(GetFullParmName("DSBTimeDelay"))){
            fDSBTimeDelay = fPm->GetDoubleParameter(GetFullParmName("DSBTimeDelay"),"Time");
        }
        else{
            fDSBTimeDelay = 0;
        }
    }
}

void DrDefinitions::InitialiseBinning(){
    unsigned repeatIt;
    G4double timeInterval;
    G4double previousTimeExtent = 0.0*s;
    G4double toThisTime;

    if (fExplicitBinning != "Null"){
        G4double previousBin = 0.;
        std::ifstream binFile(DrDefinitions::Instance()->GetExplicitBinning());
        vector<G4double> binValues;
        if(!binFile) {
            G4cout << "ERROR: could not load in " << DrDefinitions::Instance()->GetExplicitBinning() << G4endl;
            fPm->AbortSession(1);
        }
        else {
            string value;
            string line;

            while (binFile.good()) {
                getline(binFile, value, ',');
                double data;
                data = (stod(string(value, 0, value.length())))*s;
                if(data!=0.) fExplicitBins.push_back(data);
            }
        }
        for (auto bin: fExplicitBins){
            fExplicitBinDiff.push_back((bin-previousBin));
            previousBin = bin;
        }
        previousTimeExtent = fExplicitBins.back();
    }
    if (previousTimeExtent < 60.0*s) {
        //@@@@ From 0 second to 1 minute in 1 second steps
        timeInterval = 1.0*s;
        toThisTime = 60.0*s;
        repeatIt = (unsigned)floor((toThisTime-previousTimeExtent)/timeInterval);
        for(unsigned i=0; i < repeatIt; i++){
            fExplicitBinDiff.push_back(timeInterval);
        }
        previousTimeExtent = toThisTime;
    }

    if (previousTimeExtent < 600.0*s) {
        //@@@@ From 1 minute to 10 minutes in 10 second steps
        timeInterval = 10.0*s;
        toThisTime = 600.0*s;
        repeatIt = (unsigned)floor((toThisTime-previousTimeExtent)/timeInterval);
        for(unsigned i=0; i < repeatIt; i++){
            fExplicitBinDiff.push_back(timeInterval);
        }
        previousTimeExtent = toThisTime;
    }

    if (previousTimeExtent < 3600.0*s) {
        //@@@@ From 10 minutes to 1 hour in 1 minute steps
        timeInterval = 60.0*s;
        toThisTime = 3600.0*s;
        repeatIt = (unsigned)floor((toThisTime-previousTimeExtent)/timeInterval);
        for(unsigned i=0; i < repeatIt; i++){
            fExplicitBinDiff.push_back(timeInterval);
        }
        previousTimeExtent = toThisTime;
    }

    if (previousTimeExtent < fBiologyEndTime) {
        //@@@@ From 1 hour till the end in 10 minute steps
        timeInterval = 600.0*s;
        toThisTime = fBiologyEndTime;
        repeatIt = (unsigned)floor((toThisTime-previousTimeExtent)/timeInterval);
        for(unsigned i=0; i < repeatIt+1; i++){
            fExplicitBinDiff.push_back(timeInterval);
        }
    }

    fExplicitBinDiff.push_back(DBL_MAX);
}

void DrDefinitions::InitialisePathway(){

        G4int numberOfConnections = fPm->GetIntegerParameter(GetFullParmName("numberOfConnections"));

        for (G4int i = 1 ; i <= numberOfConnections ; i++){

            G4String part = "ProteinKinetic"+std::to_string(i)+"/";

            G4MoleculeDefinition* fromMoleculeDef{nullptr};
            std::vector<G4MoleculeDefinition*> toMolDefList;


            G4String fromName = fPm->GetStringParameter(GetFullParmName(part+"changeFrom"));
            if(fNameMap.find(fromName) == fNameMap.end()){
                fNameMap[fromName] = DrDSBEnd_Generic::Definition(fromName);
                G4MoleculeTable::Instance()->CreateConfiguration(fromName, DrDSBEnd_Generic::Definition(fromName));
            }
            fromMoleculeDef = DrDSBEnd_Generic::Definition(fromName);

            G4int sizeOfToNameList = fPm->GetVectorLength(GetFullParmName(part+"changeTo"));
            G4String* toNameList = fPm->GetStringVector(GetFullParmName(part+"changeTo"));
            for (G4int j = 0 ; j < sizeOfToNameList ; j++ ){
                if(fNameMap.find(toNameList[j]) == fNameMap.end()){
                    fNameMap[toNameList[j]] = DrDSBEnd_Generic::Definition(toNameList[j]);
                    G4MoleculeTable::Instance()->CreateConfiguration(toNameList[j], DrDSBEnd_Generic::Definition(toNameList[j]));
                }
                toMolDefList.emplace_back(DrDSBEnd_Generic::Definition(toNameList[j]));
            }

            G4double changeTime = fPm->GetDoubleParameter(GetFullParmName(part+"changeTime"), "Time");
            G4bool requiresCleaning = fPm->GetBooleanParameter(GetFullParmName(part+"requireClean"));
            G4bool cleanBackbone = fPm->GetBooleanParameter(GetFullParmName(part+"cleanBackbone"));
            G4bool cleanBase = fPm->GetBooleanParameter(GetFullParmName(part+"cleanBase"));
            G4int cleaningWhat{0};
            if(cleanBackbone) cleaningWhat = 1;
            else if(cleanBase) cleaningWhat = 2;

            DrProteinKinetics_Generic* proc = new DrProteinKinetics_Generic("ProteinKinetic"+std::to_string(i),
                                                                            fromMoleculeDef,
                                                                            toMolDefList,
                                                                            changeTime,
                                                                            requiresCleaning,
                                                                            cleaningWhat);

            fProcList.emplace_back(fromMoleculeDef,proc);

            if(fVerbose > 1){
                G4cout << fromName << " " << "--> ";
                for (G4int j = 0 ; j < sizeOfToNameList ; j++ ){
                    G4cout << toNameList[j] << " ";
                }
                G4cout << " with time constant of " << changeTime/s <<" s" << G4endl;
            }
        }
    if(fVerbose >1){
        G4cout << "Press any button to continue" << G4endl;
        G4cin.get();
    }
}

vector<pair<G4double,vector<G4String> > > DrDefinitions::InitialiseReactions(){

    vector<pair<G4double,vector<G4String> > > tempReactionStore;

    G4int numberOfReactions = fPm->GetIntegerParameter(GetFullParmName("numberOfReactions"));

    for (G4int i = 1 ; i <= numberOfReactions ; i++){
        G4String part = "Reaction"+std::to_string(i)+"/";

        std::vector<G4MoleculeDefinition*> fromMolDefList;
        std::vector<G4MoleculeDefinition*> toMolDefList;

        G4int sizeOfFromNameList = fPm->GetVectorLength(GetFullParmName(part+"reactants"));
        G4String* fromNameList = fPm->GetStringVector(GetFullParmName(part+"reactants"));
        for (G4int j = 0 ; j < sizeOfFromNameList ; j++ ){
            if(fNameMap.find(fromNameList[j]) == fNameMap.end()){
                fNameMap[fromNameList[j]] = DrDSBEnd_Generic::Definition(fromNameList[j]);
                G4MoleculeTable::Instance()->CreateConfiguration(fromNameList[j], DrDSBEnd_Generic::Definition(fromNameList[j]));
            }
            fromMolDefList.emplace_back(DrDSBEnd_Generic::Definition(fromNameList[j]));
        }

        G4int sizeOfToNameList = fPm->GetVectorLength(GetFullParmName(part+"products"));
        G4String* toNameList = fPm->GetStringVector(GetFullParmName(part+"products"));
        for (G4int j = 0 ; j < sizeOfToNameList ; j++ ){
            if(fNameMap.find(toNameList[j]) == fNameMap.end()){
                fNameMap[toNameList[j]] = DrDSBEnd_Generic::Definition(toNameList[j]);
                G4MoleculeTable::Instance()->CreateConfiguration(toNameList[j], DrDSBEnd_Generic::Definition(toNameList[j]));
            }
            toMolDefList.emplace_back(DrDSBEnd_Generic::Definition(toNameList[j]));
        }

        G4double reactionRange = fPm->GetDoubleParameter(GetFullParmName(part+"reactionRange"),"Length");

        vector<G4String> tempReactionString(3);
        tempReactionString[0]=fromNameList[0];
        tempReactionString[1]=fromNameList[1];
        tempReactionString[2]=toNameList[0];

        tempReactionStore.emplace_back(make_pair(reactionRange,tempReactionString));
    }

    return tempReactionStore;
}

G4String DrDefinitions::GetFullParmName(const G4String& parameterName){
    // Set the biology name
//    if (fName != "Default") return fName;
    if ( fPm->ParameterExists("Ch/ChemistryName") ) fName = fPm->GetStringParameter("Ch/ChemistryName");
    return "Ch/"+fName+"/"+parameterName;
}
