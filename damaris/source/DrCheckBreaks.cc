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
#include "DrCheckBreaks.hh"
#include "DrDefinitions.hh"
#include "DrBreakMolecule.hh"
#include <G4Scheduler.hh>
#include <G4Molecule.hh>


DrCheckBreaks::DrCheckBreaks() {}
DrCheckBreaks::~DrCheckBreaks() {}

void DrCheckBreaks::CheckRepairFidelity() {
    G4cout << G4endl << "Checking Fixed Double Stand Breaks" << G4endl;

    DrDefinitions* definitions = DrDefinitions::Instance();
    G4int misrepairCounter{0};
    G4int correctRepairCounter{0};
    G4int interChromosomeAberration{0};
    G4int interChromatidAberration{0};
    G4int interArmAberration{0};
    G4int intraArmAberration{0};

    {
        G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
        G4TrackManyList::iterator it = mainList->begin();
        G4TrackManyList::iterator end = mainList->end();

        for (; it != end; ++it) {
            G4Track *track = *it;

            G4Molecule *molecule = GetMolecule(track);
            const G4MoleculeDefinition *moleculeDefinition =
                    molecule->GetDefinition();

            if (moleculeDefinition->GetName().substr(0, 3) == "DSB") {

                DrBreakMolecule *breakMol = (DrBreakMolecule *) (track->GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));

                if (!breakMol) {
                    G4cerr << "WARNING: Simulation had no double strand break ends" << G4endl;
                    break;
                } else {
                    //!! Molecule exists
                    if (breakMol->sBreakEndB->fOriginalBreakMoleculeID != -1) {
                        //!! Has a synaptic partner
                        if (breakMol->sBreakEndA->fOriginalBreakMoleculeID !=
                            breakMol->sBreakEndB->fCorrectPartnerBreakMoleculeID) {
                            //!! Synaptic partner is incorrect end

                            misrepairCounter++;
                            //----------------------------------------------------------
                            G4double mag = (breakMol->sBreakEndA->fOriginalPosition -
                                            breakMol->sBreakEndB->fOriginalPosition)
                                    .mag();
                            definitions->fMisrepairSeparationsStore.push_back(mag);
                            //----------------------------------------------------------
                            //[DNA Structure];[Chromosome #];[Chromatid #];[Chromosome Arm]
                            std::vector<G4int> breakIDA = breakMol->sBreakEndA->fChromasomeID;
                            G4double breakPosA = breakMol->sBreakEndA->fChromasomePosition;

                            std::vector<G4int> breakIDB = breakMol->sBreakEndB->fChromasomeID;
                            G4double breakPosB = breakMol->sBreakEndB->fChromasomePosition;

                            if (breakIDA[1] != breakIDB[1]) {
                                interChromosomeAberration++;
                            } else if (breakIDA[2] != breakIDB[2]) {
                                interChromatidAberration++;
                            } else if (breakIDA[3] != breakIDB[3]) {
                                interArmAberration++;
                            } else if (breakPosA != breakPosB) {
                                intraArmAberration++;
                            }
                            //----------------------------------------------------------
                        } else {
                            correctRepairCounter++;
                        }
                    }
                }
            }
        }
    }

    G4int initB = definitions->fInitialBreakNumber;
    G4int absoluteUnrepaired = initB - misrepairCounter - correctRepairCounter;
    G4int RunID = definitions->fCurrentBiologyRepeatNumber;
    definitions->fCheckBreakStore[RunID].fMisrepairCounter = misrepairCounter;
    definitions->fCheckBreakStore[RunID].fNumberUnrepaired = absoluteUnrepaired;
    definitions->fCheckBreakStore[RunID].fInterChromosomeAberration = interChromosomeAberration;
    definitions->fCheckBreakStore[RunID].fInterChromatidAberration = interChromatidAberration;
    definitions->fCheckBreakStore[RunID].fInterArmAberration = interArmAberration;
    definitions->fCheckBreakStore[RunID].fIntraArmAberration = intraArmAberration;
}



void DrCheckBreaks::CheckNumberMoleculesLeft() {

    G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
    G4TrackManyList::iterator it = mainList->begin();
    G4TrackManyList::iterator end = mainList->end();

    G4int breakMoleculeCounter{0};
    for (; it != end; ++it) {
        breakMoleculeCounter++;
    }

    G4int RunID = DrDefinitions::Instance()->fCurrentBiologyRepeatNumber;
    DrDefinitions::Instance()->fCheckBreakStore[RunID].fNumberBreakMoleculeLeft = breakMoleculeCounter;
}

void DrCheckBreaks::StoreMisrepair() {
  G4int RunID = DrDefinitions::Instance()->fCurrentBiologyRepeatNumber;
    DrDefinitions::Instance()->fMisrepairNumberStore.push_back(
          DrDefinitions::Instance()->fCheckBreakStore[RunID].fMisrepairCounter);
}

