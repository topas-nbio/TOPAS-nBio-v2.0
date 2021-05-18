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
#include "DrBreakTable.hh"

DrCheckBreaks::DrCheckBreaks() {}
DrCheckBreaks::~DrCheckBreaks() {}

void DrCheckBreaks::CheckRepairFidelity() {
  G4cout << G4endl << "Checking Fixed Double Stand Breaks" << G4endl;

  DrBreakTable* bTable = DrBreakTable::Instance();
  G4int misrepairCounter{0};
  G4int correctRepairCounter{0};
  G4int interChromosomeAberration{0};
  G4int interChromatidAberration{0};
  G4int interArmAberration{0};
  G4int intraArmAberration{0};

  DrBreakMoleculeIterator it = bTable->GetBreakIterator();
  it.reset();
  while (it()) {
    DrBreakMolecule *breakMol = it.value();

    if (!breakMol) {
      G4cerr << "WARNING: Simulation had no double strand break ends" << G4endl;
      break;
    }
    else {
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
          bTable->fMisrepairSeparationsStore.push_back(mag);
          //----------------------------------------------------------
          //[DNA Structure];[Chromosome #];[Chromatid #];[Chromosome Arm]
          std::vector<G4int> breakIDA =  breakMol->sBreakEndA->fChromasomeID;
          G4double breakPosA = breakMol->sBreakEndA->fChromasomePosition;

          std::vector<G4int> breakIDB =  breakMol->sBreakEndB->fChromasomeID;
          G4double breakPosB = breakMol->sBreakEndB->fChromasomePosition;

          if(breakIDA[1] != breakIDB[1]){
            interChromosomeAberration++;
          }
          else if(breakIDA[2] != breakIDB[2]){
            interChromatidAberration++;
          }
          else if(breakIDA[3] != breakIDB[3]){
            interArmAberration++;
          }
          else if(breakPosA != breakPosB){
            intraArmAberration++;
          }
          //----------------------------------------------------------
        }
        else{
          correctRepairCounter++;
        }
      }
    }
  }

  G4int initB = bTable->fInitialBreakNumber;
  G4int absoluteUnrepaired = initB - misrepairCounter - correctRepairCounter;
  G4int RunID = bTable->fCurrentBiologyRepeatNumber;
  bTable->fCheckBreakStore[RunID].fMisrepairCounter = misrepairCounter;
  bTable->fCheckBreakStore[RunID].fNumberUnrepaired = absoluteUnrepaired;
  bTable->fCheckBreakStore[RunID].fInterChromosomeAberration = interChromosomeAberration;
  bTable->fCheckBreakStore[RunID].fInterChromatidAberration = interChromatidAberration;
  bTable->fCheckBreakStore[RunID].fInterArmAberration = interArmAberration;
  bTable->fCheckBreakStore[RunID].fIntraArmAberration = intraArmAberration;
}

void DrCheckBreaks::CheckNumberMoleculesLeft() {
  G4int breakMoleculeCounter{0};
  DrBreakMoleculeIterator it = DrBreakTable::Instance()->GetBreakIterator();
  it.reset();
  while (it()) {
    breakMoleculeCounter++;
  }
  G4int RunID = DrBreakTable::Instance()->fCurrentBiologyRepeatNumber;
  DrBreakTable::Instance()->fCheckBreakStore[RunID].fNumberBreakMoleculeLeft = breakMoleculeCounter;
}

void DrCheckBreaks::StoreMisrepair() {
  G4int RunID = DrBreakTable::Instance()->fCurrentBiologyRepeatNumber;
  DrBreakTable::Instance()->fMisrepairNumberStore.push_back(
          DrBreakTable::Instance()->fCheckBreakStore[RunID].fMisrepairCounter);
}

