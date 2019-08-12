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
/// \file PDBbarycenter.hh
/// \brief Definition of the Barycenter class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BARY_H
#define BARY_H

//! Molecule Class
/*!
 * This Class define Molecule model ... 
 */
class Barycenter
{
public:
  //! First constructor
  Barycenter();
  //! Second constructor
  Barycenter(int bNum,double x,double y, double z,//Nucleotide bar. coordinates
      double Bx,double By, double Bz, //Base bar. coordinates
      double Sx,double Sy, double Sz, //Sugar bar. coordinates
      double Px,double Py, double Pz);//Phosphate bar. coordinates
  //! Destructor
  ~Barycenter() {};

  //! Get the next Barycenter
  Barycenter *GetNext();
  //! Get the first
  //Residue *GetFirst();
  //! Get number Barycenter
  int GetID();
  //! Set the next Barycenter
  void SetNext(Barycenter *);
  //! Set the distance between atom i and nucleotide barycenter
  void SetDistance(int i, double);
  //! Get the distance between atom i and nucleotide barycenter
  double GetDistance(int i);
  //! Set the distance between the farther atom and nucleotide barycenter
  void SetRadius(double );
  //! Get the distance between the farther atom and nucleotide barycenter
  double GetRadius();

  int fBaryNum;//!< Barycenter number
  double fDistanceTab[33];//!< distance table [0..32] (11 hydrogens!)
  double fRadius;

  double fCenterX;            //!< "X coordinate" of this nucelotide Barycenter
  double fCenterY;            //!< "Y coordinate" of this nucelotide Barycenter
  double fCenterZ;            //!< "Z coordinate" of this nucelotide Barycenter

  double fCenterBaseX;        //!< "X coordinate" of this Base Barycenter
  double fCenterBaseY;        //!< "Y coordinate" of this Base Barycenter
  double fCenterBaseZ;        //!< "Z coordinate" of this Base Barycenter

  double fCenterSugarX;       //!< "X coordinate" of this Sugar Barycenter
  double fCenterSugarY;       //!< "Y coordinate" of this Sugar Barycenter
  double fCenterSugarZ;       //!< "Z coordinate" of this Sugar Barycenter

  double fCenterPhosphateX;   //!< "X coordinate" of this Phosphate Barycenter
  double fCenterPhosphateY;   //!< "Y coordinate" of this Phosphate Barycenter
  double fCenterPhosphateZ;   //!< "Z coordinate" of this Phosphate Barycenter

private:
  Barycenter *fpNext;    //!< Header of the next Molecule (usage before vector)
};
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

