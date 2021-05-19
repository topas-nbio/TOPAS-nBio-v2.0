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
#include "DrDNAMolecularReaction.hh"
#include <G4DNAMolecularReactionTable.hh>
#include <G4MoleculeFinder.hh>
#include <G4ITReactionChange.hh>
#include <G4Scheduler.hh>

std::unique_ptr<G4ITReactionChange> DrDNAMolecularReaction::MakeReaction(const G4Track& trackA, const G4Track& trackB){

    std::unique_ptr<G4ITReactionChange> pChanges(new G4ITReactionChange());
    pChanges->Initialize(trackA, trackB);

    const G4MolecularConfiguration* pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
    const G4MolecularConfiguration* pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();

    const G4DNAMolecularReactionData* pReactionData = fMolReactionTable->GetReactionData(pMoleculeA, pMoleculeB);

    const G4int nbProducts = pReactionData->GetNbProducts();

    if (nbProducts){

        const G4double DiffusionCoefficientMoleculeA = pMoleculeA->GetDiffusionCoefficient();
        const G4double DiffusionCoefficientMoleculeB = pMoleculeB->GetDiffusionCoefficient();
        const G4double sqrtDA = sqrt(DiffusionCoefficientMoleculeA);
        const G4double sqrtDB = sqrt(DiffusionCoefficientMoleculeB);
        const G4double numerator = sqrt(DiffusionCoefficientMoleculeA) + sqrt(DiffusionCoefficientMoleculeB);
        const G4ThreeVector reactionSite = sqrtDA / numerator * trackA.GetPosition() + sqrtDB / numerator * trackB.GetPosition();

        for (G4int j = 0; j < nbProducts; j++){
            //@@@@-----------------------------------------------------------
            //@@@@ setting time with G4Scheduler GlobalTime rather than with
            //@@@@ track.GetGlobalTime() as is the case in similar code by M.K.. This is
            //@@@@ done because at larger time scales the track global time diverges
            //@@@@ from the scheduler time. Have been unable to determine why. Interest-
            //@@@@ ingly it seems to incrementally grow in deviation size by the same
            //@@@@ .something ps amount across all simulations.
            G4double time = G4Scheduler::Instance()->GetGlobalTime();
            //@@@@ Should be:
            //G4double time = trackA.GetGlobalTime();
            //@@@@ But this throws exception ITStepManager015
            //@@@@-----------------------------------------------------------

            auto product = new G4Molecule(pReactionData->GetProduct(j));
            G4Track* productTrack = product->BuildTrack(time, reactionSite);

            productTrack->SetTrackStatus(fAlive);

            pChanges->AddSecondary(productTrack);
            G4MoleculeFinder::Instance()->Push(productTrack);
        }
    }
    pChanges->KillParents(true);
    return pChanges;
}
