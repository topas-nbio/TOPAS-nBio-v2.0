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
// $Id: TsDNAMolecularStepByStepModel.cc 94218 2015-11-09 08:24:48Z gcosmo $
//

#include "TsDNAMolecularStepByStepModel.hh"
#include "TsDNASmoluchowskiReactionModel.hh"

#include "globals.hh"
#include "G4DNAMolecularReaction.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ExceptionSeverity.hh"
#include "G4Molecule.hh"
#include "G4ReferenceCast.hh"
#include "G4VDNAReactionModel.hh"

TsDNAMolecularStepByStepModel::TsDNAMolecularStepByStepModel(const G4String& name) :
    G4VITStepModel(name),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
  fpTimeStepper = new G4DNAMoleculeEncounterStepper();
  fpReactionProcess = new G4DNAMolecularReaction();

  fType1 = G4Molecule::ITType();
  fType2 = G4Molecule::ITType();
  fReactionModel = 0;
}

TsDNAMolecularStepByStepModel::~TsDNAMolecularStepByStepModel()
{
  if(fReactionModel) delete fReactionModel;
}

TsDNAMolecularStepByStepModel& TsDNAMolecularStepByStepModel::operator=(const TsDNAMolecularStepByStepModel& right)
{
  G4ExceptionDescription exceptionDescription("Use copy constructor rather than assignement operator.");
  G4Exception("TsDNAMolecularStepByStepModel::operator=(const TsDNAMolecularStepByStepModel&)",
              "TsDNAMolecularStepByStepModel001",
              FatalErrorInArgument,
              exceptionDescription);

  if(&right == this) return *this;
  return *this; // avoid warnings
}

TsDNAMolecularStepByStepModel::TsDNAMolecularStepByStepModel(const TsDNAMolecularStepByStepModel& right) :
    G4VITStepModel(right),
    fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
  fpReactionTable = right.fpReactionTable;
  if(right.fReactionModel)
  {
    fReactionModel = right.fReactionModel->Clone();
    ((G4DNAMolecularReaction*) fpReactionProcess)->SetReactionModel(fReactionModel);
    ((G4DNAMoleculeEncounterStepper*) fpTimeStepper)->SetReactionModel(fReactionModel);
  }
  else fReactionModel = 0;
}

void TsDNAMolecularStepByStepModel::Initialize()
{
  if(fpReactionTable == 0)
  {
    SetReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
  }

  if(fReactionModel == 0)
  {
    fReactionModel = new TsDNASmoluchowskiReactionModel();
  }

  fReactionModel->SetReactionTable((const G4DNAMolecularReactionTable*) fpReactionTable);

  ((G4DNAMolecularReaction*) fpReactionProcess)->SetReactionModel(fReactionModel);
  ((G4DNAMoleculeEncounterStepper*) fpTimeStepper)->SetReactionModel(fReactionModel);

  G4VITStepModel::Initialize();
}

void TsDNAMolecularStepByStepModel::PrintInfo()
{
#ifdef G4VERBOSE
  G4cout << "DNAMolecularStepByStepModel will be used" << G4endl;
#endif
}
