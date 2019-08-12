// Extra Class for TsEmDNAChemistry
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
// $Id: TsDNASmoluchowskiReactionModel.cc 95948 2016-03-03 10:40:33Z gcosmo $

#include "TsDNASmoluchowskiReactionModel.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4UnitsTable.hh"
#include "G4Molecule.hh"
//#include "G4Scheduler.hh"
#include "G4Exp.hh"

TsDNASmoluchowskiReactionModel::TsDNASmoluchowskiReactionModel() :
    G4VDNAReactionModel()
{
  fReactionData = 0;
}

TsDNASmoluchowskiReactionModel::TsDNASmoluchowskiReactionModel(const TsDNASmoluchowskiReactionModel& __right) :
    G4VDNAReactionModel(__right)
{
  fReactionData = 0;
}

TsDNASmoluchowskiReactionModel& TsDNASmoluchowskiReactionModel::operator=(const TsDNASmoluchowskiReactionModel& right)
{
  if (this == &right) return *this;
  fReactionData = 0;
  return *this;
}

TsDNASmoluchowskiReactionModel::~TsDNASmoluchowskiReactionModel()
{
  fReactionData = 0;
}

void TsDNASmoluchowskiReactionModel::Initialise(G4MolecularConfiguration* __molecule,
                                                const G4Track&)
{
  fReactionData = fReactionTable->GetReactionData(__molecule);
}

void
TsDNASmoluchowskiReactionModel::
InitialiseToPrint(G4MolecularConfiguration* __molecule)
{
  fReactionData = fReactionTable->GetReactionData(__molecule);
}

G4double
TsDNASmoluchowskiReactionModel::GetReactionRadius(G4MolecularConfiguration* __mol1,
                                                  G4MolecularConfiguration* __mol2)
{
  G4double __output = fReactionTable->GetReactionData(__mol1, __mol2)
      ->GetEffectiveReactionRadius();
  return __output;
}

G4double TsDNASmoluchowskiReactionModel::GetReactionRadius(const G4int __i)
{
  G4double __output = (*fReactionData)[__i]->GetEffectiveReactionRadius();
  return __output;
}

G4bool TsDNASmoluchowskiReactionModel::FindReaction(const G4Track& __trackA,
                                                    const G4Track& __trackB,
                                                    const G4double __R,
                                                    G4double& __r,
                                                    const G4bool __alongStepReaction)
{
  G4double postStepSeparation = 0;
  bool do_break = false;
  G4double R2 = __R * __R;
  int k = 0;

  for (; k < 3; k++)
  {
    postStepSeparation += std::pow(
        __trackA.GetPosition()[k] - __trackB.GetPosition()[k], 2);

    if (postStepSeparation > R2)
    {
      do_break = true;
      break;
    }
  }

  if (do_break == false)
  {
    // The loop was not break
    // => __r^2 < __R^2
    __r = std::sqrt(postStepSeparation);
    return true;
  }
  else if (__alongStepReaction == true)
  {
    //Along step cheack and
    // the loop has break
    
    // Continue loop
    for (; k < 3; k++)
    {
      postStepSeparation += std::pow(
          __trackA.GetPosition()[k] - __trackB.GetPosition()[k], 2);
    }
    // Use Green approach : the Brownian bridge
    __r = (postStepSeparation = std::sqrt(postStepSeparation));

    G4Molecule* __moleculeA = GetMolecule(__trackA);
    G4Molecule* __moleculeB = GetMolecule(__trackB);

    G4double __D = __moleculeA->GetDiffusionCoefficient()
        + __moleculeB->GetDiffusionCoefficient();

    G4ThreeVector __preStepPositionA = __trackA.GetStep()->GetPreStepPoint()
        ->GetPosition();
    G4ThreeVector __preStepPositionB = __trackB.GetStep()->GetPreStepPoint()
        ->GetPosition();
    
    /*if (__preStepPositionA == __trackA.GetPosition())
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "The molecule : " << __moleculeA->GetName();
      exceptionDescription << " with track ID :" << __trackA.GetTrackID();
      exceptionDescription << " did not move since the previous step." << G4endl;
      exceptionDescription << "Current position : "
                           << G4BestUnit(__trackA.GetPosition(), "Length")
                           << G4endl;
      exceptionDescription << "Previous position : "
                           << G4BestUnit(__preStepPositionA, "Length") << G4endl;
      exceptionDescription << "Times: " << G4BestUnit(__trackA.GetStep()->GetPreStepPoint()->GetGlobalTime(), "Time") << " "
                                        << G4BestUnit(__trackA.GetGlobalTime(), "Time") <<G4endl;

      G4Exception("TsDNASmoluchowskiReactionModel::FindReaction",
                  "TsDNASmoluchowskiReactionModel", FatalErrorInArgument,
                  exceptionDescription);
    }*/

    G4double __preStepSeparation =
        (__preStepPositionA - __preStepPositionB).mag();

    //===================================
    // Brownian bridge


//    if(G4Scheduler::Instance()->GetTimeStep() != __trackB.GetStep()->GetDeltaTime())
//    {
//      G4cout << G4Scheduler::Instance()->GetTimeStep() << G4endl;
//      G4cout << __trackB.GetStep()->GetDeltaTime() << G4endl;
//      assert(G4Scheduler::Instance()->GetTimeStep() == __trackB.GetStep()->GetDeltaTime());
//    }

    G4double __probabiltyOfEncounter = G4Exp(
        -(__preStepSeparation - __R) * (postStepSeparation - __R) / (__D
            * (__trackB.GetStep()->GetDeltaTime())));
    G4double __selectedPOE = G4UniformRand();

    if (__selectedPOE <= __probabiltyOfEncounter) return true;

    //===================================
  }

  return false;
}
