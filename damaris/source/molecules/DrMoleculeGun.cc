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
#include "DrMoleculeGun.hh"
#include "DrDSBMoleculeManager.hh"
#include <G4Molecule.hh>
#include <G4MoleculeGunMessenger.hh>
#include <G4MoleculeTable.hh>
#include <Randomize.hh>

//------------------------------------------------------------------------------

template <>
void TDrMoleculeShoot<G4Track>::ShootAtRandomPosition(DrMoleculeGun *gun) {
  G4ThreeVector positionInLocalCoordinate;

  for (int i = 0; i < fNumber; ++i) {
    RandomPosInBox(*fBoxSize, positionInLocalCoordinate);
    gun->BuildAndPushTrack(fMoleculeName, fPosition + positionInLocalCoordinate,
                           fTime);
  }
}

//------------------------------------------------------------------------------

template <>
void TDrMoleculeShoot<G4Track>::ShootAtFixedPosition(DrMoleculeGun *gun) {
  for (int i = 0; i < fNumber; ++i) {
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if (fpBreakMolecule != nullptr){
      gun->BuildAndPushTrack(fpBreakMolecule,fMoleculeName);
    }
    else if (fMoleculeName.substr(0, 3) == "DSB") {
      gun->BuildAndPushTrack(fChromosomeID, fInsults, fMoleculeName, fPosition,
                             fTime);
    } else {
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      gun->BuildAndPushTrack(fMoleculeName, fPosition, fTime);
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }
}

//------------------------------------------------------------------------------

template <> void TDrMoleculeShoot<G4Track>::Shoot(DrMoleculeGun *gun) {
  if (fBoxSize)
    ShootAtRandomPosition(gun);
  else
    ShootAtFixedPosition(gun);
}

//------------------------------------------------------------------------------

// template<>
// void TG4MoleculeShoot<G4ContinuousMedium>::Shoot(G4MoleculeGun*)
//{
//    G4MolecularConfiguration* conf = G4MoleculeTable::Instance()
//        ->GetConfiguration(fMoleculeName);
//    G4MIWorkspace::GetWorldWorkspace()->GetSpeciesInCM().Add(conf,
//                                                             fNumber);
//}

//------------------------------------------------------------------------------

DrMoleculeGun::DrMoleculeGun() {
  fpMessenger = new G4MoleculeGunMessenger(this);
}

//------------------------------------------------------------------------------

DrMoleculeGun::~DrMoleculeGun() {
  if (fpMessenger)
    delete fpMessenger;
}

//------------------------------------------------------------------------------

void DrMoleculeGun::DefineTracks() {
  for (size_t i = 0; i < fShoots.size(); i++) {
    fShoots[i]->Shoot(this);
  }
}

//------------------------------------------------------------------------------

void DrMoleculeGun::AddMolecule(const G4String &name,
                                const G4ThreeVector &position, double time) {
  G4shared_ptr<DrMoleculeShoot> shoot(new TDrMoleculeShoot<G4Track>());
  shoot->fMoleculeName = name;
  shoot->fPosition = position;
  shoot->fTime = time;
  fShoots.push_back(shoot);
}
//------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void DrMoleculeGun::AddMolecule(std::pair<G4int, G4int> insults,
                                const G4String &name,
                                const G4ThreeVector &position, double time) {
  G4shared_ptr<DrMoleculeShoot> shoot(new TDrMoleculeShoot<G4Track>());
  shoot->fMoleculeName = name;
  shoot->fPosition = position;
  shoot->fTime = time;
  shoot->fInsults = insults;
  fShoots.push_back(shoot);
}

void DrMoleculeGun::AddMolecule(G4int chromosomeID,
                                std::pair<G4int, G4int> insults,
                                const G4String &name,
                                const G4ThreeVector &position, double time) {
  G4shared_ptr<DrMoleculeShoot> shoot(new TDrMoleculeShoot<G4Track>());
  shoot->fMoleculeName = name;
  shoot->fPosition = position;
  shoot->fTime = time;
  shoot->fInsults = insults;
  shoot->fChromosomeID = chromosomeID;
  fShoots.push_back(shoot);
}

void DrMoleculeGun::AddMolecule(DrBreakMolecule* breakMolecule,
                                const G4String &name) {
  G4shared_ptr<DrMoleculeShoot> shoot(new TDrMoleculeShoot<G4Track>());
  shoot->fMoleculeName = name;
  shoot->fPosition = breakMolecule->sBreakEndA->fOriginalPosition;
  shoot->fTime = breakMolecule->sBreakEndA->fLesionTime[0];
  shoot->fpBreakMolecule = breakMolecule;
  fShoots.push_back(shoot);
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//------------------------------------------------------------------------------

void DrMoleculeGun::AddNMolecules(size_t n, const G4String &moleculeName,
                                  const G4ThreeVector &position, double time) {
  G4shared_ptr<DrMoleculeShoot> shoot(new TDrMoleculeShoot<G4Track>());
  shoot->fNumber = n;
  shoot->fMoleculeName = moleculeName;
  shoot->fPosition = position;
  shoot->fTime = time;
  fShoots.push_back(shoot);
}

//------------------------------------------------------------------------------

void DrMoleculeGun::AddMoleculesRandomPositionInBox(
    size_t n, const G4String &moleculeName, const G4ThreeVector &boxCenter,
    const G4ThreeVector &boxSize, double time) {
  G4shared_ptr<DrMoleculeShoot> shoot(new TDrMoleculeShoot<G4Track>());
  shoot->fNumber = n;
  shoot->fMoleculeName = moleculeName;
  shoot->fPosition = boxCenter;
  shoot->fBoxSize = new G4ThreeVector(boxSize);
  shoot->fTime = time;
  fShoots.push_back(shoot);
}

//------------------------------------------------------------------------------

void DrMoleculeGun::BuildAndPushTrack(const G4String &name,
                                      const G4ThreeVector &position,
                                      double time) {
  G4MolecularConfiguration *conf =
      G4MoleculeTable::Instance()->GetConfiguration(name);
  assert(conf != 0);
  G4Molecule *molecule = new G4Molecule(conf);

  PushTrack(molecule->BuildTrack(time, position));
}

//------------------------------------------------------------------------------
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void DrMoleculeGun::BuildAndPushTrack(G4int /*chromosomeID*/,
                                      std::pair<G4int, G4int> insults,
                                      const G4String &name,
                                      const G4ThreeVector &position,
                                      double time) {
    G4MolecularConfiguration *conf = G4MoleculeTable::Instance()->GetConfiguration(name);
    assert(conf != 0);
    G4Molecule *molecule = new G4Molecule(conf);
    G4Track *track = molecule->BuildTrack(time, position);
    PushTrack(track);

    auto DSBMoleculeManager = new DrDSBMoleculeManager();

    DrBreakMolecule *breakMolecule = DSBMoleculeManager->LinkNewBreakMolecule(track);
    DSBMoleculeManager->NewBreakID(breakMolecule);

    breakMolecule->sBreakEndA->fDamageTypes[1] = insults.first;
    breakMolecule->sBreakEndA->fDamageTypes[0] = insults.second;
    breakMolecule->sBreakEndA->fOriginalPosition = position;
    G4int breakID = breakMolecule->sBreakEndA->fOriginalBreakMoleculeID;
    if (breakID % 2 != 0) breakMolecule->sBreakEndA->fCorrectPartnerBreakMoleculeID = breakID - 1; // odd
    else breakMolecule->sBreakEndA->fCorrectPartnerBreakMoleculeID = breakID + 1; // even

    delete DSBMoleculeManager;
}

void DrMoleculeGun::BuildAndPushTrack(DrBreakMolecule* breakMol,
                                      const G4String &name) {
    G4MolecularConfiguration *conf = G4MoleculeTable::Instance()->GetConfiguration(name);
    assert(conf != 0);
    auto *molecule = new G4Molecule(conf);
    G4Track *track = molecule->BuildTrack(breakMol->sBreakEndA->fLesionTime[0],
            breakMol->sBreakEndA->fOriginalPosition);
    PushTrack(track);
    auto DSBMoleculeManager = new DrDSBMoleculeManager();
    DSBMoleculeManager->LinkThisBreakMolecule(track,breakMol);
    delete DSBMoleculeManager;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//------------------------------------------------------------------------------

void DrMoleculeGun::GetNameAndNumber(DrMoleculeGun::NameNumber &output) {
  for (size_t i = 0; i < fShoots.size(); ++i) {
    output[fShoots[i]->fMoleculeName] += fShoots[i]->fNumber;
  }
}

//------------------------------------------------------------------------------

void DrMoleculeShoot::RandomPosInBox(const G4ThreeVector &boxSize,
                                     G4ThreeVector &output) {
  output[0] = boxSize.x() * G4UniformRand() - boxSize.x() / 2;
  output[1] = boxSize.y() * G4UniformRand() - boxSize.y() / 2;
  output[2] = boxSize.z() * G4UniformRand() - boxSize.z() / 2;
}

//------------------------------------------------------------------------------

DrMoleculeShoot::DrMoleculeShoot()
    : G4enable_shared_from_this<DrMoleculeShoot>() {
  fMoleculeName = "";
  fTime = 0;
  fNumber = 1;
  fBoxSize = 0;
  //@@@@
  fInsults = {0, 0};
  //@@@@
}

//------------------------------------------------------------------------------

DrMoleculeShoot::~DrMoleculeShoot() {
  if (fBoxSize)
    delete fBoxSize;
}

//------------------------------------------------------------------------------

void DrMoleculeGun::AddMoleculeShoot(G4shared_ptr<DrMoleculeShoot> shoot) {
  fShoots.push_back(shoot);
}

void DrMoleculeGun::AddMoleculeInCMRepresentation(size_t n,
                                                  const G4String &moleculeName,
                                                  double time) {
  G4shared_ptr<DrMoleculeShoot> shoot(
      new TDrMoleculeShoot<G4ContinuousMedium>());
  shoot->fNumber = n;
  shoot->fMoleculeName = moleculeName;
  shoot->fTime = time;
  fShoots.push_back(shoot);
}
