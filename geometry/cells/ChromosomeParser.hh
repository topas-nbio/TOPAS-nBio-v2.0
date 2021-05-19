/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#ifndef ChromosomeParser_h
#define ChromosomeParser_h

#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <fstream>
#include <vector>

struct ChromObj{
    G4ThreeVector Pos;
    G4int ChromID=-1;
    G4int ChromCopy=-1;
    G4int VolumeCopyNum=-1;
    G4double Radius=-1.0;
    G4double GenomeLength=-1.0;
    G4double FractionOfChrom=-1.0;
    G4String name;
    G4double MinFractAlong=-1.0;
    G4double MaxFractAlong=-1.0;
    G4int MinBPAlong=-1;
    G4int MaxBPAlong=-1;
    G4double FractionHeterochromatic=-1.0;
};

class ChromosomeParser{
public:
    ChromosomeParser(){}
    ~ChromosomeParser(){}
    
    void ReadBeads(G4String filename);
    void RecentreBeads();
    void ResizeNucleus();
    
    std::map<G4int,std::vector<ChromObj>> GetBeads(){return Beads;}
    G4ThreeVector GetDimensions(){return Dimensions;}
    G4RotationMatrix* GetRotation(){return rot;}
    std::map<G4int,G4int> GetChromosomeSize(){return ChromSize;}
    
    
private:
    
    
    std::map<G4int,std::vector<ChromObj>> Beads;
    G4ThreeVector Dimensions;
    G4RotationMatrix*rot;
    std::map<G4int,G4int> ChromSize;
};

#endif /* ChromosomeObject_h */
