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
#include "DrDNASmoluchowskiReactionModel.hh"
#include <G4UnitsTable.hh>
#include <G4Molecule.hh>

//@@@@
#include "DrBreakTable.hh"
#include "DrDefinitions.hh"

G4bool DrDNASmoluchowskiReactionModel::FindReaction(const G4Track& trackA, const G4Track& trackB, const G4double reactionRange, G4double& separationDistance, const G4bool isAlongStepReaction){
    G4double postStepSeparation = 0;
    bool doBreak = false;
    G4double reactionRangeSq = reactionRange * reactionRange;
    int k = 0;

    for (; k < 3; k++){
        postStepSeparation += std::pow(trackA.GetPosition()[k] - trackB.GetPosition()[k], 2);

        if (postStepSeparation > reactionRangeSq){
            doBreak = true;
            break;
        }
    }

    //@@@@ If molecules are inside reaction radius at end of step
    //@@@@ then the reaction happens.
    if (!doBreak){
        separationDistance = std::sqrt(postStepSeparation);
        return true;
    }
        //@@@@ Else calculate the probability of meeting each other in
        //@@@@ that step
    else if (isAlongStepReaction){
        //Along step check and the loop has broken
        // Continue loop
        //@@@@ Finish calculating the separation at end of step
        for (; k < 3; k++){
            postStepSeparation += std::pow(trackA.GetPosition()[k] - trackB.GetPosition()[k], 2);
        }

        // Use Green approach : the Brownian bridge
        separationDistance = (postStepSeparation = std::sqrt(postStepSeparation));

        G4Molecule* moleculeA = GetMolecule(trackA);
        G4Molecule* moleculeB = GetMolecule(trackB);

        G4double combinedDiffusionCoefficient = moleculeA->GetDiffusionCoefficient() + moleculeB->GetDiffusionCoefficient();

        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        //this could be made more efficient but I'm doing it this way to keep the codes separate
        G4String MoleculeAName = GetMolecule(trackA)->GetDefinition()->GetName();
        G4String MoleculeBName = GetMolecule(trackB)->GetDefinition()->GetName();

        if(MoleculeAName.substr(0,3) == "DSB" && MoleculeBName.substr(0,3) == "DSB"){

            G4bool isWaitingMoleculeA = DrBreakTable::Instance()->GetBreakMolecule(trackA, "SmolReact")->fIsWaiting;
            G4double diffusionCoefficientMoleculeA = (isWaitingMoleculeA) ? DrDefinitions::Instance()->GetTrapDiff() : DrDefinitions::Instance()->GetJumpDiff();

            G4bool isWaitingMoleculeB = DrBreakTable::Instance()->GetBreakMolecule(trackB, "SmolReact")->fIsWaiting;
            G4double diffusionCoefficientMoleculeB = (isWaitingMoleculeB) ? DrDefinitions::Instance()->GetTrapDiff() : DrDefinitions::Instance()->GetJumpDiff();

            combinedDiffusionCoefficient = diffusionCoefficientMoleculeA + diffusionCoefficientMoleculeB;
        }
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        G4ThreeVector preStepPositionA = trackA.GetStep()->GetPreStepPoint()->GetPosition();
        G4ThreeVector preStepPositionB = trackB.GetStep()->GetPreStepPoint()->GetPosition();

        if (preStepPositionA == trackA.GetPosition()
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            && GetMolecule(trackA)->GetDefinition()->GetName().substr(0,3) != "DSB"
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                )
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription << "The molecule : " << moleculeA->GetName();
            exceptionDescription << " with track ID :" << trackA.GetTrackID();
            exceptionDescription << " did not move since the previous step." << G4endl;
            exceptionDescription << "Current position : "
                                 << G4BestUnit(trackA.GetPosition(), "Length")
                                 << G4endl;
            exceptionDescription << "Previous position : "
                                 << G4BestUnit(preStepPositionA, "Length") << G4endl;
            G4Exception("DrDNASmoluchowskiReactionModel::FindReaction",
                        "DrDNASmoluchowskiReactionModel", FatalErrorInArgument,
                        exceptionDescription);
        }

        G4double preStepSeparation = (preStepPositionA - preStepPositionB).mag();

        //===================================
        // Brownian bridge
        G4double probabilityOfEncounter = std::exp(-(preStepSeparation - reactionRange) * (postStepSeparation - reactionRange) / (combinedDiffusionCoefficient * (trackB.GetStep()->GetDeltaTime())));
        G4double selectedProbabilityOfEncounter = G4UniformRand();

        if (selectedProbabilityOfEncounter <= probabilityOfEncounter) return true;
        //===================================
    }
    return false;
}
