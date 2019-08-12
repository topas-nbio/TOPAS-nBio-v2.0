// Extra Class for PDBmolecule
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
/// \file PDBresidue.cc
/// \brief Implementation of the PDBresidue class

#include "PDBresidue.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue::Residue():fResName(""),fResSeq(0),fVisible(false),fSelected(false),
fCenterX(0),fCenterY(0),fCenterZ(0),fNbAtom(0),fpNext(0),fpFirst(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue::Residue(string rN,int rS)
{
  fResName=rN;// Residue name
  fResSeq=rS; // Residue sequence number
  fVisible=false;
  fSelected=false;
  fCenterX=0;
  fCenterY=0;
  fCenterZ=0;
  fNbAtom=0;
  fpNext=0;
  fpFirst=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Residue *Residue::GetNext()
{
  return fpNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Atom *Residue::GetFirst()
{
  return fpFirst;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Residue::GetID()
{
  return fResSeq;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Residue::SetNext(Residue *residueNext)
{
  fpNext=residueNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Residue::SetFirst(Atom *atomFirst)
{
  fpFirst=atomFirst;
}
