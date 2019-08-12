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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// --------------------------------------------------------------
// Authors: E. Delage
// november 2013
// --------------------------------------------------------------
//
// $Id$
//
/// \file PDBmolecule.hh
/// \brief Definition of the PDBmolecule class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MOLECULE_H
#define MOLECULE_H

#include "PDBresidue.hh"

#include <iostream>

using namespace std;

class Residue;

//! Molecule Class
/*!
 * This Class define Molecule model ... 
 */
class Molecule
{
public:
  //! First constructor
  Molecule();
  //! Second constructor
  Molecule(string resName,int mNum);
  //! Destructor
  ~Molecule() {};

  //! information about molecule (not implemented)
  //void PrintInfo();
  //! Get the next molecule
  Molecule *GetNext();
  //! Get the first Residue
  Residue *GetFirst();
  //! Get number Molecule
  int GetID();
  //! Set the next Molecule
  void SetNext(Molecule *);
  //! Set the first Residue
  void SetFirst(Residue *);

  string fMolName;   //!< Molecule name
  int fMolNum;       //!< Molecule number

  double fMinGlobZ;   //Cylinder length => min Z
  double fMaxGlobZ;
  double fMinGlobX;   //Radius => min X
  double fMaxGlobX;
  double fMinGlobY;   //=> min Y
  double fMaxGlobY;

  int fCenterX;      //!< "X center" of this Molecule (for rotation...)
  int fCenterY;      //!< "Y center" of this Molecule (for rotation...)
  int fCenterZ;//!< "Z center" of this Molecule (for rotation...)
  int fDistCenterMax;//!< dist from center to most away most of the molecule
  int fNbResidue;        //!< Number of residue inside the molecule

private:
  Molecule *fpNext;//!< Header of the next Molecule (usage before vector)
  Residue *fpFirst;//!< Header of the first Residue (usage before vector)
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
