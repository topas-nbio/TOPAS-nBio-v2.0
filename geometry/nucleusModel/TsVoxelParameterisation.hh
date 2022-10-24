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

#ifndef TsVOXELPARAMETERISATION_HH
#define TsVOXELPARAMETERISATION_HH

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"
#include "TsVGeometryComponent.hh"
#include "TsParameterManager.hh"

#include "TsVoxelParameterisation.hh"

class TsVoxelParameterisation : public G4VPVParameterisation
{
public:
    TsVoxelParameterisation(TsParameterManager*);
    virtual ~TsVoxelParameterisation();
    virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

    
    void ComputeVoxelIndices(const G4int copyNo, size_t& nx,size_t& ny, size_t& nz ) const;// Convert the copyNo to voxel numbers in x, y and z.
    void SetVoxelDimensions( G4double halfx, G4double halfy, G4double halfz );
    void SetContainerDimensions( G4double halfx, G4double halfy, G4double halfz );
    void SetNoVoxel( size_t nx, size_t ny, size_t nz );
    void CheckCopyNo( const G4int copyNo ) const;
    G4ThreeVector GetTranslation(const G4int copyNo ) const;


protected:

    G4double fVoxelHalfX,fVoxelHalfY,fVoxelHalfZ;
    G4double fContainerWallX, fContainerWallY, fContainerWallZ;
    size_t fNoVoxelX,fNoVoxelY,fNoVoxelZ;
    size_t fNoVoxelXY;
    size_t fNoVoxel;

private:
    TsParameterManager* fPm;

};

#endif // TsVOXELPARAMETERISATION_HH
