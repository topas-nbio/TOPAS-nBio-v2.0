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
#include <G4VDNAReactionModel.hh>

DrDNAMolecularReaction::DrDNAMolecularReaction()
        : G4VITReactionProcess()
        , fMolReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
        , fpReactionModel(nullptr)
{
}

DrDNAMolecularReaction::DrDNAMolecularReaction(G4VDNAReactionModel* pReactionModel)
        : DrDNAMolecularReaction()
{
    fpReactionModel = pReactionModel;
}

G4bool DrDNAMolecularReaction::TestReactibility(const G4Track &trackA,
                                                const G4Track &trackB,
                                                double currentStepTime,
                                                bool userStepTimeLimit) /*const*/
{
    const auto pMoleculeA = GetMolecule(trackA)->GetMolecularConfiguration();
    const auto pMoleculeB = GetMolecule(trackB)->GetMolecularConfiguration();

    const G4double reactionRadius = fpReactionModel->GetReactionRadius(pMoleculeA, pMoleculeB);

    G4double separationDistance = -1.;

    if (currentStepTime == 0.)
    {
        userStepTimeLimit = false;
    }

    G4bool output = fpReactionModel->FindReaction(trackA, trackB, reactionRadius,
                                                  separationDistance, userStepTimeLimit);
    return output;
}

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

void DrDNAMolecularReaction::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fpReactionModel = pReactionModel;
}

std::vector<std::unique_ptr<G4ITReactionChange>> DrDNAMolecularReaction::FindReaction(
        G4ITReactionSet* pReactionSet,
        const double currentStepTime,
        const double /*fGlobalTime*/,
        const bool reachedUserStepTimeLimit)
{
    std::vector<std::unique_ptr<G4ITReactionChange>> fReactionInfo;
    fReactionInfo.clear();

    if (pReactionSet == nullptr)
    {
        return fReactionInfo;
    }

    G4ITReactionPerTrackMap& reactionPerTrackMap = pReactionSet->GetReactionMap();
    for (auto tracks_i = reactionPerTrackMap.begin();
         tracks_i != reactionPerTrackMap.end();
         tracks_i = reactionPerTrackMap.begin())
    {
        G4Track* pTrackA = tracks_i->first;
        if (pTrackA->GetTrackStatus() == fStopAndKill)
        {
            continue;
        }

        G4ITReactionPerTrackPtr reactionPerTrack = tracks_i->second;
        G4ITReactionList& reactionList = reactionPerTrack->GetReactionList();

        assert(reactionList.begin() != reactionList.end());

        for (auto it = reactionList.begin(); it != reactionList.end(); it = reactionList.begin())
        {
            G4ITReactionPtr reaction(*it);
            G4Track* pTrackB = reaction->GetReactant(pTrackA);
            if (pTrackB->GetTrackStatus() == fStopAndKill)
            {
                continue;
            }

            if (pTrackB == pTrackA)
            {
                G4ExceptionDescription exceptionDescription;
                exceptionDescription
                        << "The IT reaction process sent back a reaction between trackA and trackB. ";
                exceptionDescription << "The problem is trackA == trackB";
                G4Exception("G4ITModelProcessor::FindReaction",
                            "ITModelProcessor005",
                            FatalErrorInArgument,
                            exceptionDescription);
            }

            pReactionSet->SelectThisReaction(reaction);

            if (TestReactibility(*pTrackA, *pTrackB, currentStepTime, reachedUserStepTimeLimit))
            {
                auto pReactionChange = MakeReaction(*pTrackA, *pTrackB);

                if (pReactionChange)
                {
                    fReactionInfo.push_back(std::move(pReactionChange));
                    break;
                }
            }
        }
    }

    pReactionSet->CleanAllReaction();
    return fReactionInfo;
}