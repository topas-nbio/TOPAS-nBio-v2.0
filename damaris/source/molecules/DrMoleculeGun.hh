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

//@@@@
#include "DrBreakMolecule.hh"
#include <G4MoleculeGun.hh>

class G4Track;
class DrMoleculeGun;
typedef int G4ContinuousMedium;

//------------------------------------------------------------------------------

/*
 * Define a specific species shoot
 * Multiple shoots maybe be defined
 */

class DrMoleculeShoot : public G4enable_shared_from_this<DrMoleculeShoot> {
public:
  DrMoleculeShoot();
  virtual ~DrMoleculeShoot();
  virtual void Shoot(DrMoleculeGun *) = 0;

  template <typename TYPE> G4shared_ptr<DrMoleculeShoot> ChangeType();

  G4String fMoleculeName;
  G4ThreeVector fPosition;
  G4double fTime;
  G4int fNumber;
  G4ThreeVector *fBoxSize;
  //@@@@
  std::pair<G4int, G4int> fInsults{-1, -1};
  G4int fChromosomeID{0};
  DrBreakMolecule *fpBreakMolecule{nullptr};
  //@@@@
  static void RandomPosInBox(const G4ThreeVector &boxSize,
                             G4ThreeVector &output);
};

//------------------------------------------------------------------------------

/*
 * Define a shoot type =
 *   track (used by the Smoluchowski code)
 *   or continuous medium (used by the gillespie code)
 */

template <typename TYPE> class TDrMoleculeShoot : public DrMoleculeShoot {
public:
  TDrMoleculeShoot() : DrMoleculeShoot() { ; }
  virtual ~TDrMoleculeShoot() { ; }
  void Shoot(DrMoleculeGun *) {}

protected:
  void ShootAtRandomPosition(DrMoleculeGun *) {}
  void ShootAtFixedPosition(DrMoleculeGun *) {}
};

template <typename TYPE>
G4shared_ptr<DrMoleculeShoot> DrMoleculeShoot::ChangeType() {
  G4shared_ptr<DrMoleculeShoot> output(new TDrMoleculeShoot<TYPE>);
  output->fMoleculeName = fMoleculeName;
  output->fPosition = fPosition;
  output->fTime = fTime;
  output->fNumber = fNumber;
  output->fBoxSize = fBoxSize;
  //@@@@
  output->fInsults = fInsults;
  output->fChromosomeID = fChromosomeID;
  output->fpBreakMolecule = fpBreakMolecule;
  //@@@@
  return output;
}

//------------------------------------------------------------------------------

class DrMoleculeGun : public G4MoleculeGun {
public:
  DrMoleculeGun();
  virtual ~DrMoleculeGun();

  virtual void DefineTracks();

  /*
   * Create a single molecule
   * @param moleculeName name of the molecule such as recorded in molecule table
   * @param position position where the molecule should pop up
   * @param time time at which the molecule should pop up
   */
  void AddMolecule(const G4String &moleculeName, const G4ThreeVector &position,
                   double time = 0);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  void AddMolecule(std::pair<G4int, G4int>, const G4String &moleculeName,
                   const G4ThreeVector &position, double time = 0);

  void AddMolecule(G4int, std::pair<G4int, G4int>, const G4String &moleculeName,
                   const G4ThreeVector &position, double time = 0);
  void AddMolecule(DrBreakMolecule *, const G4String &moleculeName);

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  /*
   * Create N molecules at a single point
   * @param n number of molecules to create
   * @param moleculeName name of the molecules as recorded in molecule table
   * @param position position where the molecules should pop up
   * @param time time at which the molecules should pop up
   */
  void AddNMolecules(size_t n, const G4String &moleculeName,
                     const G4ThreeVector &position, double time = 0);

  /*
   * Create N molecules in a box
   * @param n number of molecules to create
   * @param moleculeName name of the molecules as recorded in molecule table
   * @param boxCenter center of the box
   * @param boxExtension size of the box
   * @param time time at which the molecules should pop up
   */
  void AddMoleculesRandomPositionInBox(size_t n, const G4String &moleculeName,
                                       const G4ThreeVector &boxCenter,
                                       const G4ThreeVector &boxExtension,
                                       double time = 0);

  /*
   * Create N molecules as component of the continuous medium in a box
   * @param n number of molecules to create
   * @param moleculeName name of the molecules as recorded in molecule table
   * @param boxCenter center of the box
   * @param boxExtension size of the box
   * @param time time at which the molecules should pop up
   */
  //  void AddMoleculeInCMRepresentation(size_t n,
  //                                     const G4String& moleculeName,
  //                                     const G4ThreeVector& boxCenter,
  //                                     const G4ThreeVector& boxExtension,
  //                                     double time = 0);

  void AddMoleculeInCMRepresentation(size_t n, const G4String &moleculeName,
                                     double time = 0);

  const std::vector<G4shared_ptr<DrMoleculeShoot>> &GetMoleculeShoot() {
    return fShoots;
  }

  typedef std::map<G4String, int> NameNumber;
  void GetNameAndNumber(NameNumber &);

  void AddMoleculeShoot(G4shared_ptr<DrMoleculeShoot>);

protected:
  void BuildAndPushTrack(const G4String &name, const G4ThreeVector &position,
                         double time = 0);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  void BuildAndPushTrack(G4int, std::pair<G4int, G4int>, const G4String &name,
                         const G4ThreeVector &position, double time = 0);
  void BuildAndPushTrack(DrBreakMolecule*, const G4String &name);
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  G4MoleculeGunMessenger *fpMessenger;

  std::vector<G4shared_ptr<DrMoleculeShoot>> fShoots;
  friend class DrMoleculeShoot;
  template <class T> friend class TDrMoleculeShoot;
};