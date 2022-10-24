// EndRun2 for TOPAS
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
#include "DrEndRunForBiology.hh"
#include "DrDefinitions.hh"
#include "DrDSBGun.hh"
#include "DrBiologyTimeStepAction.hh"
#include "DrCheckBreaks.hh"
#include "TsDamagePhaseSpaceStore.hh"
#include "DrUtils.cc"
#include "DrPrecompiler.hh"
#include <G4Scheduler.hh>
#include <G4DNAChemistryManager.hh>
#include <G4SystemOfUnits.hh>
#include <TsParameterManager.hh>
#include <G4UnitsTable.hh>
#include <G4MoleculeFinder.hh>      //Used to clean memory between runs

DrEndRunForBiology::DrEndRunForBiology(TsParameterManager* fPm) {

    if(fPm->ParameterExists("Ch/DaMaRiS/Bool_DaMaRiS")){
        if(fPm->GetBooleanParameter("Ch/DaMaRiS/Bool_DaMaRiS")){
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            //@@   Initialising Definitions
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            DrDefinitions::Instance()->SetVerbosity(0);
            DrDefinitions::Instance()->InitialiseBinning();

            G4PhysicsModelCatalog::Register("DrBreakMolecule");
            G4int numberOfRepeats = DrDefinitions::Instance()->GetBiologyRepeatNumber();
            G4Scheduler::Instance()->SetEndTime(DrDefinitions::Instance()->GetBiologyEndTime());

            G4cout << G4endl
                   << "-----------------------------------------------------------"
                   << G4endl << "Physics stage ends" << G4endl
                   << "DrBreakMolecule registered as physics process "
                   << G4PhysicsModelCatalog::GetIndex("DrBreakMolecule") << G4endl
                   << "Number of biology repeats selected is " << numberOfRepeats
                   << G4endl << "Biology run time selected is "
                   << G4Scheduler::Instance()->GetEndTime() / 1e9 << " s" << G4endl
                   << G4endl;

            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            //@@@@ Overwrite the PerRepResults file and insert header
            G4int altRunID = DrDefinitions::Instance()->GetAlternativeRunID();
            std::ofstream file("PerRepResults" + std::to_string(altRunID) + ".out");
            file << "Initial DSBs | No. MisRep | No. UnRep | "
                 << "No. InterCosomeAb | No. InterCatidAb | "
                 << "No. InterArmAb | No. ItraArmAb" << G4endl;
            file.close();
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            for (G4int i{0}; i < numberOfRepeats; i++) {
                //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                //@@   Placement Of Molecules Before Activating Chemistry   @@
                //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                G4cout << "Placing DNA break molecules..." << G4endl;

                DrDefinitions::Instance()->fCurrentBiologyRepeatNumber = i + 1;
                auto *placeDSBMolecules = new DrDSBGun();
                auto *definitions = DrDefinitions::Instance();
                G4int goAhead{0};

                //@@@@ Decide which method to use to place DSBs
                G4bool origin = (definitions->GetDSBOriginNumber() >= 0);
                G4bool offset = (definitions->GetDSBOffset() != 0.0);
                G4bool column = (definitions->GetDSBColumnNumber() >= 0);
                G4bool separation = (definitions->GetDSBSeparation() >= 0.0);
                G4bool damageStore = (!TsDamagePhaseSpaceStore::Instance()->GetDamageEventStore().empty());
                G4bool fromFile = definitions->GetIsReadFromFile();

                G4int checkPlacements = origin + offset + column + separation + fromFile;

                if(checkPlacements > 1){
                    G4cerr << "ERROR: Only one DSB placement method may be used" << G4endl
                            << "Origin = " << origin << G4endl
                            << "Offset = " << offset << G4endl
                            << "Column = " << column << G4endl
                            << "Separation = " << separation << G4endl
                            << "From File = " << fromFile << G4endl << G4endl;
                    fPm->AbortSession(1);
                }
                else if (checkPlacements == 0 && !damageStore){
                    G4cerr << "ERROR: No DSB placement method specified therefore expecting" << G4endl
                            << "DamagePhaseSpaceStore to be populated by another simulation." << G4endl
                            << "However, DamagePhaseSPaceStore was empty." << G4endl << G4endl;
                    fPm->AbortSession(1);
                }

                if (definitions->GetDSBOriginNumber() >= 0) {
                    G4int numberOf = definitions->GetDSBOriginNumber();
                    goAhead = placeDSBMolecules->BuildDSBEnd_Origin(numberOf,"DSBEnd");
                }
                else if(definitions->GetDSBOffset() != 0.0){
                    G4double offset_nm = definitions->GetDSBOffset();
                    goAhead = placeDSBMolecules->BuildDSBEnd_Offset("DSBEnd", offset_nm);
                }
                else if(definitions->GetDSBColumnNumber() >= 0){
                    G4int numberOf = definitions->GetDSBColumnNumber();
                    G4double radiusOf = definitions->GetDSBColumnRadius();
                    goAhead = placeDSBMolecules->BuildDSB_Column(numberOf,radiusOf);
                }
                else if (definitions->GetDSBSeparation() >= 0.0) {
                    //@@@@ to place two double strand breaks a set distance apart (first
                    //@@@@ argument). The complexity can be set with the next two arguments,
                    //@@@@ backbones first then bases, but please note this complexity will
                    //@@@@ be copied for all four break ends placed. The last argument
                    //@@@@ allows for placing eith two double strand breaks (0) or instead
                    //@@@@ place two double strand break ends (1).
                    G4double DSBSep = definitions->GetDSBSeparation();
                    G4double DSBdelay = definitions->GetDSBTimeDelay();
                    goAhead = placeDSBMolecules->BuildDSB_SepSpaceAndTime(DSBSep * nm, 0, 0, 0, DSBdelay);
                }
                else if (false) {
                    //@@@@ to place a custom distribution of double strand breaks as defined
                    //@@@@ by the vector dsbPattern. Each element represents a group of DSBs
                    //@@@@ to be placed. The groups will be separated by 100 nm
                    //@@@@ symmetrically about the origin on the x axis. Make sure you are
                    //@@@@ not trying to place more groups than can fit in your cell nucleus
                    std::vector<G4int> dsbPattern = {1};
                    goAhead = placeDSBMolecules->BuildDSB_Pattern(0, 0, dsbPattern);
                }
                else if (!TsDamagePhaseSpaceStore::Instance()->GetDamageEventStore().empty()){
                    //@@@@ Here we have previously read in a file and populated the store
                    //@@@@ of damages. To save time we should use this store rather than
                    //@@@@ reading in the file again.
                    if (definitions->GetIsReadFromFile()) {
                        std::vector<std::vector<std::vector<DrDamageEvent*> > > store;
                        store = TsDamagePhaseSpaceStore::Instance()->GetDamageEventStore();
                        G4int selection = placeDSBMolecules->SelectRandomFromSTDInput(store);
                        goAhead = placeDSBMolecules->PlaceFromSTDInput(store[selection]);
                    }
                    else {
                        //@@@@ The store of damages has been populated by another simulation
                        //@@@@ and we should use this as the source of damages rather than
                        //@@@@ require that simulation to write out to file.
                        std::vector<std::vector<std::vector<DrDamageEvent*> > > store;
                        store = TsDamagePhaseSpaceStore::Instance()->GetDamageEventStore();
                        goAhead = placeDSBMolecules->PlaceFromSTDInput(store[0]);
                        //@@@@ To clean up the phase space before next run
                        TsDamagePhaseSpaceStore::Instance()->Destroy();
                    }
                }
                else {
                    //@@@@ to place using new STD format file. Can take either a single new
                    //@@@@ event (ie. one run of an experiment) or can take multiple new
                    //@@@@ events and will select one at random to simulate. Can also select
                    //@@@@ a specific wingle new event from the multiple new events by
                    //@@@@ setting SelectFromExposureNumber.
                    goAhead = placeDSBMolecules->PlaceFromSTDFile(
                            definitions->GetDamageFile());
                }

                //@@@@ Place the clock molecule into the world to force time steps
                //@@@@ of a certain size.
                placeDSBMolecules->PlaceClockMolecule();

                if (goAhead != 0) {

                    placeDSBMolecules->DefineAllTracks();

                    //@@@@ Avoid error from large amount of zero steps which could be
                    //@@@@ caused by stationary (waiting) particles on the cell boundary
                    G4Scheduler::Instance()->SetMaxZeroTimeAllowed(10000);

                    G4cout << "-----------------------------------------------------------"
                           << G4endl << "Biology stage begins..." << G4endl << G4endl;
                    G4Scheduler::Instance()->SetVerbose(1);

                    //@@@@ Switches to biology user time step action in order to handle
                    //@@@@ the reactions between DSB ends.
                    auto *TSA = new DrBiologyTimeStepAction(fPm,"DaMaRiS");
                    G4Scheduler::Instance()->SetUserAction(TSA);
                    G4Scheduler::Instance()->SetVerbose(0);

                    G4cout << "*** Biology Simulation Starts ***" << G4endl;
                    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    G4DNAChemistryManager::Instance()->Run(); // starts chemistry
                    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    G4cout << "*** Biology Simulation Finishes At: "
                           << G4BestUnit(G4Scheduler::Instance()->GetEndTime(), "Time")
                           << " ***" << G4endl;

                    //@@@@ checking if the correct partner ends have been joined
                    //@@@@ will output a file as well as write to screen a summary
                    DrUtils utilities = DrUtils();

                    utilities.PrintPerRepBiologicallyRelevantParameters();
                    utilities.PrintToScreenBioParam();
                    //@@@@ Writes to file each change of molecule type and it's time.
                    //@@@@ Used to track evolution of system in time
                    utilities.ExtractBiologicallyRelevantParameters();

                    //@@@@ cleans memory between runs?
                    G4MoleculeFinder::Instance()->Clear();
                    definitions->ResetCurrentExplicitBinNuber();

                    G4cout << "Biology Repeat " << i + 1 << " Ended" << G4endl
                           << "-----------------------------------------------------------"
                           << G4endl << G4endl;
                }
                delete placeDSBMolecules;
            }
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            DrUtils *utilities = new DrUtils();
            //@@@@ Writes to file the kinetics
            utilities->PrintBiologicallyRelevantParameters();

            if(DrDefinitions::Instance()->GetDSBSeparation() >= 0.0){
                utilities->PrintMRPercent();
            }

        #ifdef DEBUG_DAMARIS

            std::ofstream debugFile("debugProcesses.out", std::ios_base::app);

            for(auto procMap: DrDefinitions::Instance()->debugProcMap){
                G4String name = procMap.first;
                std::vector<std::pair<G4double, G4double>> outList;
                std::vector<G4double> inList = DrDefinitions::Instance()->debugProcMap[name];
                outList = utilities->BinDoubleList(inList,100,0,0);
                debugFile << name << G4endl;
                for (auto read: outList){
                    debugFile << read.first << " " << read.second << G4endl;
                }
                debugFile << G4endl << G4endl;
            }

            debugFile << G4endl << G4endl << "DEBUG_SpaceSteps" << G4endl;
            utilities->PrintBinList("debugProcesses.out",
                                    DrDefinitions::Instance()->spaceStepStore, 1 * nm,
                                    100000, 0, 10000 * nm, true);

            debugFile << G4endl << G4endl << "DEBUG_ActualReactionRange" << G4endl;
            utilities->PrintBinList("debugProcesses.out",
                                    DrDefinitions::Instance()->actualReactionRangeStore, 1 * nm,
                                    1000, 0, 1000 * nm, true);

            debugFile << G4endl << G4endl << "DEBUG_SuggestedReactionRange" << G4endl;
            utilities->PrintBinList("debugProcesses.out",
                                    DrDefinitions::Instance()->suggestedReactionRangeStore,
                                    1*nm, 1000, 0, 1000 * nm, true);
            debugFile << G4endl << G4endl << "DEBUG_DiffusionTime" << G4endl;
            utilities->PrintBinList("debugProcesses.out",
                                    DrDefinitions::Instance()->diffTimeStore,
                                    1*s, 1000, 0, 1000 * s, true);

        #endif /*DEBUG_DAMARIS*/

            delete utilities;

            if (!TsDamagePhaseSpaceStore::Instance()->GetDamageEventStore().empty() &&
                DrDefinitions::Instance()->GetIsReadFromFile()) {
                TsDamagePhaseSpaceStore::Instance()->Destroy();
            }
        }
    }
}
