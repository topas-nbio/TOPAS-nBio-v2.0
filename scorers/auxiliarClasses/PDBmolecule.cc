// Extra Class for PDBlib
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
// $Id$
// 
/// \file PDBmolecule.cc
/// \brief Implementation file for PDBmolecule class

#include "PDBmolecule.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule::Molecule():
fMolName(""),fMolNum(0),fMinGlobZ(0),fMaxGlobZ(0),
fMinGlobX(0),fMaxGlobX(0),fMinGlobY(0),fMaxGlobY(0),
fCenterX(0),fCenterY(0),fCenterZ(0),fDistCenterMax(0),fNbResidue(0),
fpNext(0),fpFirst(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule::Molecule(string mN,int mNum)
{
  fMolName=mN;  //Molecule name
  fMolNum=mNum; //Molecule number
  fMinGlobZ=0;
  fMaxGlobZ=0;
  fMinGlobX=0;
  fMaxGlobX=0;
  fMinGlobY=0;
  fMaxGlobY=0;
  fCenterX=0;
  fCenterY=0;
  fCenterZ=0;
  fDistCenterMax=0;
  fNbResidue=0;
  fpNext=0;
  fpFirst=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Molecule *Molecule::GetNext()
{
  return fpNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue *Molecule::GetFirst()
{
  return fpFirst;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Molecule::GetID()
{
  return fMolNum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Molecule::SetNext(Molecule *moleculeNext)
{
  fpNext=moleculeNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Molecule::SetFirst(Residue *resFirst)
{
  fpFirst=resFirst;
}

