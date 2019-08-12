// Extra Class for PDBresidue
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
/// \file PDBatom.cc
/// \brief Implementation of the PDBatom class

#include "PDBatom.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Atom::Atom(int s,string n,string rN,int numInRes,int rS,
      double xInit,double yInit,double zInit,
      double radius,
      double o, double tF, string e)
{
  fSerial=s;
  fName=n;//!< Atom name
  fResName=rN;//!< Residue name
  fNumInRes=numInRes;
  fResSeq=rS;//!< Residue sequence number
  fX=xInit;
  fY=yInit;
  fZ=zInit;
  fVdwRadius=radius;
  fOccupancy=o;//!< occupancy
  fTempFactor=tF;
  fElement=e;
  fpNext=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Atom *Atom::GetNext()
{
  return fpNext;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Atom::GetX()
{
  return fX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Atom::GetY()
{
  return fY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Atom::GetZ()
{
  return fZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Atom::GetID()
{
  return fSerial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

string Atom::GetName()
{
  return fName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

string Atom::GetElementName()
{
  return fElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Atom::GetVanDerWaalsRadius()
{
  return fVdwRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Atom::SetNext(Atom *AtomNext)
{
  fpNext=AtomNext;
}
