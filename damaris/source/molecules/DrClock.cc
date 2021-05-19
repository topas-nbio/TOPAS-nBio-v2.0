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
#include "DrClock.hh"
#include <G4ParticleTable.hh>

DrClock *DrClock::theInstance = nullptr;

DrClock *DrClock::Definition() {
  if (theInstance != nullptr) return theInstance;
  const G4String name = "Clock";
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr) anInstance = new G4MoleculeDefinition(name, 0., 0., 0, -1, 0., 0, 0);
  theInstance = reinterpret_cast<DrClock *>(anInstance);
  return theInstance;
}
