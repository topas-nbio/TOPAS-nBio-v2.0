// Extra Class for ClusteringAlgo ClusterSBPoints
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
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Authors: Henri Payno and Yann Perrot
//
// $Id$
//
/// \file SBPoint.cc
/// \brief Implementation of the SBPoint class

#include "SBPoint.hh"

#include <Randomize.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

SBPoint::SBPoint(unsigned int pId, G4ThreeVector pPos, G4double pEdep, G4int pStrand ):
fId(pId),
fPosition(pPos),
fEdep(pEdep),
fpCluster(0),
fStrand(pStrand)
{
  // pick randomly one strand
  if ( fStrand < 0 )
  	fStrand = (G4UniformRand() < 0.5) ? 0 : 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SBPoint::~SBPoint()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator != (const SBPoint& pPt) const
{
  return pPt.fId != fId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator == (const SBPoint& pPt) const
{
  return pPt.fId == fId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator < (const SBPoint& pPt) const
{
  return pPt.fId < fId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator > (const SBPoint& pPt) const
{
  return pPt.fId > fId;
}
