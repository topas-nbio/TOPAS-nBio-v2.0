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

#include <G4DNASmoluchowskiReactionModel.hh>

/**
  * G4DNASmoluchowskiReactionModel should be used
  * for very fast reactions (high reaction rate) : the reactions between
  * reactants occuring at encounter.
  * When the time step is constrained this model
  * uses brownian bridge : "Absorbing (Smoluchowski) boundary condition"
  * Reference : On the simulation of diffusion processes close to boundaries,
  * N. J. B. Green, Molecular Physics, 65: 6, 1399 — 1408(1988)
  */

class DrDNASmoluchowskiReactionModel : public G4DNASmoluchowskiReactionModel{
public :
    virtual G4bool FindReaction(const G4Track&, const G4Track&, const G4double, G4double&, const G4bool);
};