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
#include "DrDNAMoleculeEncounterStepper.hh"
#include <G4VDNAReactionModel.hh>
#include <G4DNAMolecularReactionTable.hh>
#include <G4UnitsTable.hh>
#include <G4MoleculeFinder.hh>
#include <G4MolecularConfiguration.hh>
#include "DrDefinitions.hh"
#include "DrBreakTable.hh"
#include "DrPrecompiler.hh"

using namespace std;
using namespace CLHEP;


#ifdef DEBUG_MEM
#include <G4MemStat.hh>
using namespace G4MemStat;
#endif

DrDNAMoleculeEncounterStepper::Utils::Utils(const G4Track& tA, const G4MolecularConfiguration* pMoleculeB) : fpTrackA(tA), fpMoleculeB(pMoleculeB){
  fpMoleculeA = GetMolecule(tA);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  G4String substringOfName = fpMoleculeA->GetDefinition()->GetName().substr(0,3);

  if(substringOfName=="DSB" && DrDefinitions::Instance()->GetIsSubDiffusion()){
    fDiffusionCoefficientMoleculeA= DrDefinitions::Instance()->GetJumpDiff();
    fDiffusionCoefficientMoleculeB= DrDefinitions::Instance()->GetJumpDiff();
  }
  else{
    fDiffusionCoefficientMoleculeA= fpMoleculeA->GetDiffusionCoefficient();
    fDiffusionCoefficientMoleculeB= fpMoleculeB->GetDiffusionCoefficient();
  }

  G4double diffusionCoefficientSum = fDiffusionCoefficientMoleculeA+ fDiffusionCoefficientMoleculeB;
  G4double diffusionCoefficientProduct = fDiffusionCoefficientMoleculeA * fDiffusionCoefficientMoleculeB;
  fConstant = 8 * ( diffusionCoefficientSum + 2 * sqrt(diffusionCoefficientProduct));
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

DrDNAMoleculeEncounterStepper::DrDNAMoleculeEncounterStepper()
        :G4DNAMoleculeEncounterStepper()
        , fHasAlreadyReachedNullTime(false)
        , fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
        , fReactionModel(nullptr)
        , fVerbose(0)
{
}

void DrDNAMoleculeEncounterStepper::InitializeForNewTrack()
{
    if (fReactants)
    {
        fReactants.reset();
    }
    fSampledMinTimeStep = DBL_MAX;
    fHasAlreadyReachedNullTime = false;
}

void DrDNAMoleculeEncounterStepper::CheckAndRecordResults(const Utils& utils,
#ifdef G4VERBOSE
                                                          const G4double R,
#endif
                                                          G4KDTreeResultHandle& results)
{
    if (results == 0)
    {
#ifdef G4VERBOSE
        if (fVerbose > 1)
        {
            G4cout << "No molecule " << utils.fpMoleculeB->GetName()
                   << " found to react with " << utils.fpMoleculeA->GetName()
                   << G4endl;
        }
#endif
        return;
    }

    for (results->Rewind(); !results->End(); results->Next())
    {
        G4IT* reactiveB = results->GetItem<G4IT>();

        if (reactiveB == 0)
        {
            continue;
        }

        G4Track *trackB = reactiveB->GetTrack();

        if (trackB == 0)
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                    << "The reactant B found using the MoleculeFinder does not have a valid "
                       "track attached to it. If this is done on purpose, please do "
                       "not record this molecule in the MoleculeFinder."
                    << G4endl;
            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults",
                        "MoleculeEncounterStepper001", FatalErrorInArgument,
                        exceptionDescription);
            continue;
        }

        if (trackB->GetTrackStatus() != fAlive)
        {
            continue;
        }

        if (trackB == &utils.fpTrackA)
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                    << "A track is reacting with itself (which is impossible) ie fpTrackA == trackB"
                    << G4endl;
            exceptionDescription << "Molecule A (and B) is of type : "
                                 << utils.fpMoleculeA->GetName() << " with trackID : "
                                 << utils.fpTrackA.GetTrackID() << G4endl;

            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults",
                        "MoleculeEncounterStepper003", FatalErrorInArgument,
                        exceptionDescription);

        }

        if (fabs(trackB->GetGlobalTime() - utils.fpTrackA.GetGlobalTime())
            > utils.fpTrackA.GetGlobalTime() * (1 - 1 / 100))
        {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                    << "The interacting tracks are not synchronized in time" << G4endl;
            exceptionDescription
                    << "trackB->GetGlobalTime() != fpTrackA.GetGlobalTime()" << G4endl;

            exceptionDescription << "fpTrackA : trackID : " << utils.fpTrackA.GetTrackID()
                                 << "\t Name :" << utils.fpMoleculeA->GetName()
                                 << "\t fpTrackA->GetGlobalTime() = "
                                 << G4BestUnit(utils.fpTrackA.GetGlobalTime(), "Time") << G4endl;

            exceptionDescription << "trackB : trackID : " << trackB->GetTrackID()
                                 << "\t Name :" << utils.fpMoleculeB->GetName()
                                 << "\t trackB->GetGlobalTime() = "
                                 << G4BestUnit(trackB->GetGlobalTime(), "Time") << G4endl;

            G4Exception("G4DNAMoleculeEncounterStepper::RetrieveResults",
                        "MoleculeEncounterStepper004", FatalErrorInArgument,
                        exceptionDescription);
        }

#ifdef G4VERBOSE
        if (fVerbose > 1)
        {
            G4double r2 = results->GetDistanceSqr();
            G4cout << "\t ************************************************** " << G4endl;
            G4cout << "\t Reaction between "
                   << utils.fpMoleculeA->GetName() << " (" << utils.fpTrackA.GetTrackID() << ") "
                   << " & " << utils.fpMoleculeB->GetName() << " (" << trackB->GetTrackID() << "), "
                   << "Interaction Range = "
                   << G4BestUnit(R, "Length") << G4endl;
            G4cout << "\t Real distance between reactants  = "
                   << G4BestUnit((utils.fpTrackA.GetPosition() - trackB->GetPosition()).mag(), "Length") << G4endl;
            G4cout << "\t Distance between reactants calculated by nearest neighbor algorithm = "
                   << G4BestUnit(sqrt(r2), "Length") << G4endl;
        }
#endif

        fReactants->push_back(trackB);
    }
}

G4double
DrDNAMoleculeEncounterStepper::CalculateStep(const G4Track& trackA,
                                             const G4double& userMinTimeStep)
{
    auto pMoleculeA = GetMolecule(trackA);
    InitializeForNewTrack();
    fUserMinTimeStep = userMinTimeStep;

#ifdef G4VERBOSE
    if (fVerbose)
    {
        G4cout
                << "_______________________________________________________________________"
                << G4endl;
        G4cout << "G4DNAMoleculeEncounterStepper::CalculateStep" << G4endl;
        G4cout << "Check done for molecule : " << pMoleculeA->GetName()
               << " (" << trackA.GetTrackID() << ") "
               << G4endl;
    }
#endif

    //__________________________________________________________________
    // Retrieve general informations for making reactions
    auto pMolConfA = pMoleculeA->GetMolecularConfiguration();

    const auto pReactantList = fMolecularReactionTable->CanReactWith(pMolConfA);

    if (!pReactantList)
    {
#ifdef G4VERBOSE
        if (fVerbose > 1)
        {
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
            G4cout << "!!! WARNING" << G4endl;
            G4cout << "G4MoleculeEncounterStepper::CalculateStep will return infinity "
                      "for the reaction because the molecule "
                   << pMoleculeA->GetName()
                   << " does not have any reactants given in the reaction table."
                   << G4endl;
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
        }
#endif
        return DBL_MAX;
    }

    G4int nbReactives = pReactantList->size();

    if (nbReactives == 0)
    {
#ifdef G4VERBOSE
        if (fVerbose)
        {
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
            G4cout << "!!! WARNING" << G4endl;
            G4cout << "G4MoleculeEncounterStepper::CalculateStep will return infinity "
                      "for the reaction because the molecule "
                   << pMoleculeA->GetName()
                   << " does not have any reactants given in the reaction table."
                   << "This message can also result from a wrong implementation of the reaction table."
                   << G4endl;
            G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
        }
#endif
        return DBL_MAX;
    }

    fReactants.reset(new vector<G4Track*>());
    fReactionModel->Initialise(pMolConfA, trackA);

    //__________________________________________________________________
    // Start looping on possible reactants
    for (G4int i = 0; i < nbReactives; i++)
    {
        auto pMoleculeB = (*pReactantList)[i];

        //______________________________________________________________
        // Retrieve reaction range
        const G4double R = fReactionModel->GetReactionRadius(i);

        //______________________________________________________________
        // Use KdTree algorithm to find closest reactants
        G4KDTreeResultHandle resultsNearest(
                G4MoleculeFinder::Instance()->FindNearest(pMoleculeA,
                                                          pMoleculeB->GetMoleculeID()));

        if (resultsNearest == 0) continue;

        G4double r2 = resultsNearest->GetDistanceSqr();
        Utils utils(trackA, pMoleculeB);

        if (r2 <= R * R) // ==> Record in range
        {
            // Entering in this condition may due to the fact that molecules are very close
            // to each other
            // Therefore, if we only take the nearby reactant into account, it might have already
            // reacted. Instead, we will take all possible reactants that satisfy the condition r<R

            if (fHasAlreadyReachedNullTime == false)
            {
                fReactants->clear();
                fHasAlreadyReachedNullTime = true;
            }

            fSampledMinTimeStep = 0.;
            G4KDTreeResultHandle resultsInRange(
                    G4MoleculeFinder::Instance()->FindNearestInRange(pMoleculeA,
                                                                     pMoleculeB->GetMoleculeID(),
                                                                     R));
            CheckAndRecordResults(utils,
#ifdef G4VERBOSE
                                  R,
#endif
                                  resultsInRange);
        }
        else
        {
            G4double r = sqrt(r2);
            G4double tempMinET = pow(r - R, 2) / utils.fConstant;
            // constant = 16 * (fDA + fDB + 2*sqrt(fDA*fDB))

            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if(pMoleculeA->GetDefinition()->GetName().substr(0,3) == "DSB"){
                //@@@@ Here the minimum time required for the particles to diffuse
                //@@@@ within interaction range is set depending on sub-diffusive
                //@@@@ parameters
                if(DrDefinitions::Instance()->GetIsSubDiffusion()){
                    G4IT* reactiveB = resultsNearest->GetItem<G4IT>();
                    G4Track *trackB = reactiveB->GetTrack();

                    G4bool isWaitingMoleculeA = DrBreakTable::Instance()->GetBreakMolecule(trackA)->fIsWaiting;
                    G4bool isWaitingMoleculeB = DrBreakTable::Instance()->GetBreakMolecule(*trackB)->fIsWaiting;

                    G4double JumpDCMoleculeA = DrDefinitions::Instance()->GetJumpDiff();
                    G4double JumpDCMoleculeB = DrDefinitions::Instance()->GetJumpDiff();
                    G4double TrappedDCMoleculeA = DrDefinitions::Instance()->GetTrapDiff();
                    G4double TrappedDCMoleculeB = DrDefinitions::Instance()->GetTrapDiff();

                    if(isWaitingMoleculeA && isWaitingMoleculeB){
                        //@@@@ This means that both particles, although within the
                        //@@@@ maximum range that could lead to a reaction were
                        //@@@@ they both freely moving, are in fact trapped and thus
                        //@@@@ the estimated time for them to react is infinity
                        G4double calc = pow(r - R, 2) /(8 * (TrappedDCMoleculeA + TrappedDCMoleculeB + 2 * sqrt(TrappedDCMoleculeA * TrappedDCMoleculeB)));
                        (calc == 0)? tempMinET = DBL_MAX : tempMinET = calc;


                    }
                    else if( isWaitingMoleculeA || !isWaitingMoleculeB ){
                        //@@@@ This means that one particle is trapped and the other
                        //@@@@ is not so the estimated time for them to react is
                        //@@@@ the time for one to move to the other
                        tempMinET = pow(r - R, 2) /(8 * (JumpDCMoleculeB + TrappedDCMoleculeA + 2 * sqrt(JumpDCMoleculeB * TrappedDCMoleculeA)));
                    }
                    else if( !isWaitingMoleculeA || isWaitingMoleculeB ){
                        //@@@@ This means that one particle is trapped and the other
                        //@@@@ is not so the estimated time for them to react is
                        //@@@@ the time for one to move to the other
                        tempMinET = pow(r - R, 2) /(8 * (JumpDCMoleculeA + TrappedDCMoleculeB + 2 * sqrt(JumpDCMoleculeA * TrappedDCMoleculeB)));
                    }
                }
            }
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            if (tempMinET <= fSampledMinTimeStep)
            {
                if (fUserMinTimeStep < DBL_MAX/*IsInf(fUserMinTimeStep) == false*/
                    && tempMinET <= fUserMinTimeStep) // ==> Record in range
                {
                    if (fSampledMinTimeStep > fUserMinTimeStep)
                    {
                        fReactants->clear();
                    }

                    fSampledMinTimeStep = fUserMinTimeStep;

                    G4double range = R + sqrt(fUserMinTimeStep*utils.fConstant);

                    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                    if(pMoleculeA->GetDefinition()->GetName().substr(0,3) == "DSB" && DrDefinitions::Instance()->GetIsSubDiffusion()){

                        G4IT* reactiveB = resultsNearest->GetItem<G4IT>();
                        G4Track *trackB = reactiveB->GetTrack();

                        G4bool isWaitingMoleculeA = DrBreakTable::Instance()->GetBreakMolecule(trackA, "MolEnc")->fIsWaiting;
                        G4bool isWaitingMoleculeB = DrBreakTable::Instance()->GetBreakMolecule(*trackB, "MolEnc")->fIsWaiting;

                        G4double JumpDCMoleculeA = DrDefinitions::Instance()->GetJumpDiff();
                        G4double JumpDCMoleculeB = DrDefinitions::Instance()->GetJumpDiff();
                        G4double TrappedDCMoleculeA = DrDefinitions::Instance()->GetTrapDiff();
                        G4double TrappedDCMoleculeB = DrDefinitions::Instance()->GetTrapDiff();

                        if(isWaitingMoleculeA && isWaitingMoleculeB){
                            //@@@@ This means that both particles, although within the
                            //@@@@ maximum range that could lead to a reaction were
                            //@@@@ they both freely moving, are in fact trapped and thus
                            //@@@@ the range should be 0
                            range = R + sqrt(fUserMinTimeStep * (8 * (TrappedDCMoleculeA + TrappedDCMoleculeB + 2 * sqrt(TrappedDCMoleculeA * TrappedDCMoleculeB))));
                        }
                        else if( isWaitingMoleculeA || !isWaitingMoleculeB ){
                            //@@@@ This means that one particle is trapped and the other
                            //@@@@ is not so the estimated time for them to react is
                            //@@@@ the time for one to move to the other
                            range = R + sqrt(fUserMinTimeStep * (8 * (JumpDCMoleculeB + TrappedDCMoleculeA + 2 * sqrt(JumpDCMoleculeB * TrappedDCMoleculeA))));
                        }
                        else if( !isWaitingMoleculeA || isWaitingMoleculeB ){
                            //@@@@ This means that one particle is trapped and the other
                            //@@@@ is not so the estimated time for them to react is
                            //@@@@ the time for one to move to the other
                            range = R + sqrt(fUserMinTimeStep * (8 * (JumpDCMoleculeA + TrappedDCMoleculeB + 2 * sqrt(JumpDCMoleculeA * TrappedDCMoleculeB))));
                        }

#ifdef DEBUG_DAMARIS
                        DrBreakTable::Instance()->suggestedReactionRangeStore.push_back(range);
#endif /*DEBUG_DAMARIS*/

                    }
                    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                    G4KDTreeResultHandle resultsInRange(
                            G4MoleculeFinder::Instance()->
                                    FindNearestInRange(pMoleculeA,
                                                       pMoleculeB->GetMoleculeID(),
                                                       range));

                    CheckAndRecordResults(utils,
#ifdef G4VERBOSE
                                          range,
#endif
                                          resultsInRange);
                }
                else // ==> Record nearest
                {
                    if (tempMinET < fSampledMinTimeStep)
                        // to avoid cases where fSampledMinTimeStep == tempMinET
                    {
                        fSampledMinTimeStep = tempMinET;
                        fReactants->clear();
                    }

                    CheckAndRecordResults(utils,
#ifdef G4VERBOSE
                                          R,
#endif
                                          resultsNearest);
                }
            }
        }
    }

#ifdef G4VERBOSE
    if (fVerbose)
    {
        G4cout << "G4MoleculeEncounterStepper::CalculateStep will finally return :"
               << G4BestUnit(fSampledMinTimeStep, "Time") << G4endl;

        if (fVerbose > 1)
        {
            G4cout << "Selected reactants for trackA: " << pMoleculeA->GetName()
                   << " (" << trackA.GetTrackID() << ") are: ";

            vector<G4Track*>::iterator it;
            for (it = fReactants->begin(); it != fReactants->end(); it++)
            {
                G4Track* trackB = *it;
                G4cout << GetMolecule(trackB)->GetName() << " ("
                       << trackB->GetTrackID() << ") \t ";
            }
            G4cout << G4endl;
        }
    }
#endif
    return fSampledMinTimeStep;
}