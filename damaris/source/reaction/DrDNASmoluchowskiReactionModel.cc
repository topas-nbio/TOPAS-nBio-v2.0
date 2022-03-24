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
#include "DrBreakMolecule.hh"
#include "DrDefinitions.hh"
#include "TsParameterManager.hh"
#include <G4UnitsTable.hh>
#include <G4Molecule.hh>
#include <G4Track.hh>
#include <Randomize.hh>
#include <G4SystemOfUnits.hh>
#include <G4DNAMolecularReactionTable.hh>

DrDNASmoluchowskiReactionModel::DrDNASmoluchowskiReactionModel()
        : G4VDNAReactionModel()
        , fpReactionData(nullptr)
{}

DrDNASmoluchowskiReactionModel::~DrDNASmoluchowskiReactionModel() = default;

G4double DrDNASmoluchowskiReactionModel::GetReactionRadius(const G4int __i)
{
    G4double __output = (*fpReactionData)[__i]->GetEffectiveReactionRadius();
    return __output;
}

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

            DrBreakMolecule* breakMoleculeA = (DrBreakMolecule*)(trackA.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));
            DrBreakMolecule* breakMoleculeB = (DrBreakMolecule*)(trackB.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));

            G4bool isWaitingMoleculeA = breakMoleculeA->fIsWaiting;
            G4double jumpDiffusionCoefficient = DrDefinitions::Instance()->GetJumpDiff();
            G4double trappedDiffusionCoefficient = DrDefinitions::Instance()->GetTrapDiff();
            G4double diffusionCoefficientMoleculeA = (isWaitingMoleculeA) ? trappedDiffusionCoefficient : jumpDiffusionCoefficient;

            G4bool isWaitingMoleculeB = breakMoleculeB->fIsWaiting;
            G4double diffusionCoefficientMoleculeB = (isWaitingMoleculeB) ? trappedDiffusionCoefficient : jumpDiffusionCoefficient;

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
          G4cerr
          << "DrDNASmoluchowskiReactionModel::FindReaction"
          << "The molecule : " << moleculeA->GetName() << G4endl
          << " with track ID :" << trackA.GetTrackID() << G4endl
          << " did not move since the previous step." << G4endl
          << "Current position : " << G4endl
          << G4BestUnit(trackA.GetPosition(), "Length") << G4endl
          << G4endl << G4endl
          << "Previous position : " << G4endl
          << G4BestUnit(preStepPositionA, "Length") << G4endl;
            DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
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

void DrDNASmoluchowskiReactionModel::Initialise(const G4MolecularConfiguration* pMolecule,
                                                const G4Track&)
{
    fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

void DrDNASmoluchowskiReactionModel::InitialiseToPrint(const G4MolecularConfiguration* pMolecule)
{
    fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

G4double DrDNASmoluchowskiReactionModel::GetReactionRadius(const G4MolecularConfiguration* pMol1,
                                                           const G4MolecularConfiguration* pMol2)
{
    G4double __output = fpReactionTable->GetReactionData(pMol1, pMol2)->GetEffectiveReactionRadius();
    return __output;
}
