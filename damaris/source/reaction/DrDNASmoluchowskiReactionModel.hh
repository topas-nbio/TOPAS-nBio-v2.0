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

#include "G4VDNAReactionModel.hh"
#include <vector>

class G4DNAMolecularReactionData;

/**
  * G4DNASmoluchowskiReactionModel should be used
  * for very fast reactions (high reaction rate) : the reactions between
  * reactants occuring at encounter.
  * When the time step is constrained this model
  * uses brownian bridge : "Absorbing (Smoluchowski) boundary condition"
  * Reference : On the simulation of diffusion processes close to boundaries,
  * N. J. B. Green, Molecular Physics, 65: 6, 1399 — 1408(1988)
  */

class DrDNASmoluchowskiReactionModel : public G4VDNAReactionModel{
public :
    DrDNASmoluchowskiReactionModel();
    virtual ~DrDNASmoluchowskiReactionModel();

    DrDNASmoluchowskiReactionModel(const DrDNASmoluchowskiReactionModel&) = delete;
    DrDNASmoluchowskiReactionModel& operator=(const DrDNASmoluchowskiReactionModel&) = delete;

    virtual void Initialise(const G4MolecularConfiguration*, const G4Track&) ;
    virtual void InitialiseToPrint(const G4MolecularConfiguration*) ;
    virtual G4double GetReactionRadius(const G4MolecularConfiguration*,
                                       const G4MolecularConfiguration*);

    virtual G4double GetReactionRadius(const G4int);

    virtual G4bool FindReaction(const G4Track&,
                                const G4Track&,
                                G4double /*reactionRadius*/,
                                G4double& /*separationDistance*/,
                                G4bool /*alongStepInteraction*/) ;

private :
    const std::vector<const G4DNAMolecularReactionData*>* fpReactionData ;
};