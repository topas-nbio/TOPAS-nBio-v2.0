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

#include <map>
#include <G4Types.hh>
#include <G4String.hh>

template<typename MOLECULE>
class DrBreakIterator
{
protected:
  typedef std::map<G4int, MOLECULE*> MAP;
  MAP* fMap;
  G4bool fDefined;
  typename MAP::iterator fIt;

public:
  DrBreakIterator(MAP& _map): fMap(&_map){
    fDefined = false;
  }

  virtual ~DrBreakIterator(){}

  DrBreakIterator(const DrBreakIterator& right){
    fMap = right.fMap;
    fDefined = right.fDefined;
    fIt = right.fIt;
  }

  DrBreakIterator& operator=(const DrBreakIterator& right){
    if (this == &right) return *this;
    fMap = right.fMap;
    fDefined = right.fDefined;
    fIt = right.fIt;
    return *this;
  }

  G4bool operator++(int){
    if (!fDefined) return false;
    fIt++;
    return fIt != fMap->end() ? true : false;
  }

  G4bool operator++(){
    if (!fDefined) return false;
    fIt++;
    return fIt != fMap->end() ? true : false;
  }

  void reset(){
    fDefined = false;
  }

  G4bool operator()(){
    if (fDefined == false){
      fDefined = true;
      fIt = fMap->begin();
      return true;
    }
    else{
      fIt++;
    }
    if (fIt == fMap->end()) return false;
    return true;
  }

  const G4String& Name(){
    return fIt->first;
  }

  MOLECULE* value(){
    return fIt->second;
  }
};
