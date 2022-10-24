// Extra Class for TsNucleus
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
// Author: Hongyu Zhu
// Created date: 12/07/2018
// Last edit   : 07/08/2019 by Alexander Klapproth

#include "TsVoxelParameterisation.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "TsParameterManager.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"


TsVoxelParameterisation::TsVoxelParameterisation(TsParameterManager* pM)
    : G4VPVParameterisation(), fPm(pM)
{}

TsVoxelParameterisation::~TsVoxelParameterisation()
{}

void TsVoxelParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4ThreeVector trans = GetTranslation( copyNo );
  physVol->SetTranslation( trans );
}


G4ThreeVector TsVoxelParameterisation::GetTranslation(const G4int copyNo ) const
{
  CheckCopyNo( copyNo );

  size_t nx;
  size_t ny;
  size_t nz;

  ComputeVoxelIndices( copyNo, nx, ny, nz );

  G4ThreeVector trans( (2*nx+1)*fVoxelHalfX - fContainerWallX,
                       (2*ny+1)*fVoxelHalfY - fContainerWallY,
                       (2*nz+1)*fVoxelHalfZ - fContainerWallZ);

  return trans;
}

void TsVoxelParameterisation::ComputeVoxelIndices(const G4int copyNo, size_t& nx,
                            size_t& ny, size_t& nz ) const
{
  CheckCopyNo( copyNo );
  nx = size_t(copyNo%fNoVoxelX);
  ny = size_t( (copyNo/fNoVoxelX)%fNoVoxelY );
  nz = size_t(copyNo/fNoVoxelXY);
}

void TsVoxelParameterisation::SetNoVoxel( size_t nx, size_t ny, size_t nz )
{
  fNoVoxelX = nx; 
  fNoVoxelY = ny; 
  fNoVoxelZ = nz; 
  fNoVoxelXY = nx*ny; 
  fNoVoxel = nx*ny*nz;
}

void TsVoxelParameterisation::SetVoxelDimensions( G4double halfx, G4double halfy, G4double halfz )
{
  fVoxelHalfX = halfx; 
  fVoxelHalfY = halfy; 
  fVoxelHalfZ = halfz; 
}

void TsVoxelParameterisation::SetContainerDimensions( G4double halfx, G4double halfy, G4double halfz )
{
    fContainerWallX = halfx; 
    fContainerWallY = halfy; 
    fContainerWallZ = halfz; 
}

void TsVoxelParameterisation::CheckCopyNo( const G4int copyNo ) const
{ 
  if( copyNo < 0 || copyNo >= G4int(fNoVoxel) )
  {
    G4cerr << "Copy number is negative or too big!" << G4endl;
    G4cerr << "        Copy number: " << copyNo << G4endl;
    G4cerr << "        Total number of voxels: " << fNoVoxel << G4endl;
    fPm->AbortSession(1);
  }
}

