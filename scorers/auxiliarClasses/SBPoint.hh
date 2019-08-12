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
/// \file SBPoint.hh
/// \brief Definition of the SBPoint class

#ifndef SB_POINT_HH
#define SB_POINT_HH

#include <assert.h>
#include <G4ThreeVector.hh>

class ClusterSBPoints;
/// \brief defines a point of energy deposition which defines a damage to the DNA.
class SBPoint
{
public:
  /// \brief constructor
  SBPoint(unsigned int, G4ThreeVector pPos, G4double pEdep, G4int pStrand );
  /// \brief destructor
  ~SBPoint();

  // Get methods
  G4int GetID() const
  {
    return fId;
  }
  G4ThreeVector GetPosition() const {
    return fPosition;
  }
  G4double GetEdep() const {
    return fEdep;
  }
  ClusterSBPoints* GetCluster() const {
    return fpCluster;
  }
  G4int GetTouchedStrand() const {
    return fStrand;
  }

  // Set methods
  void SetCluster(ClusterSBPoints* pCluster)
  {
    assert(pCluster); fpCluster = pCluster;
  }

  void SetPosition(G4ThreeVector pPos)
  {
    fPosition=pPos;
  }

  bool HasCluster() const {
    return fpCluster != 0;
  }
  void CleanCluster() {
    fpCluster = 0;
  }

  bool operator != (const SBPoint& ) const;
  bool operator == (const SBPoint& ) const;
  bool operator < (const SBPoint& ) const;
  bool operator > (const SBPoint& ) const;

private:

  unsigned int fId;             //ID
  G4ThreeVector fPosition;      //Position
  G4double fEdep;               //Edp
  ClusterSBPoints* fpCluster;   // Associated clustered points
  G4int fStrand;                // Strand

};

#endif // SB_POINT_HH
