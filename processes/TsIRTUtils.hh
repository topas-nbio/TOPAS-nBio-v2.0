#ifndef TsIRTUtils_hh
#define TsIRTUtils_hh

#include "globals.hh"
#include <vector>

class TsIRTUtils {

public:
    TsIRTUtils();
   ~TsIRTUtils();

    G4double NormQuantile(G4double x);
    G4double erfcx_y100(G4double x);
    G4double erfcx(G4double x);
    G4double erfc(G4double x);
    G4double erfcInv(G4double x);
    G4double erfcWxy(G4double c, G4double x, G4double y);
    
    G4int FindBin(G4double x, std::vector<G4double> vx);
    G4int FindBin(G4int n, G4double xmin, G4double xmax, G4double);

    std::vector<G4double> CreateTimeSteps(G4double t0, G4double tf, G4int tbins, G4bool isLogSpace);
 
    G4double SampleTypeII(G4double alpha, G4double sigma, G4double r0, G4double D);
    G4double Lambda(G4double x, G4double beta, G4double alphatilde);
    
};

#endif

