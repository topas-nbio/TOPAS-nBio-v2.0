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

#ifndef TsDNAOneStepThermalizationModel_
#define TsDNAOneStepThermalizationModel_

#include "G4VEmModel.hh"

class G4ITNavigator;
class G4Navigator;

namespace DNA{
    namespace Penetration{
        //-----------------------
        /*
         * Article: Jintana Meesungnoen, Jean-Paul Jay-Gerin,
         *          Abdelali Filali-Mouhim, and Samlee Mankhetkorn (2002)
         *          Low-Energy Electron Penetration Range in Liquid Water.
         *          Radiation Research: November 2002, Vol. 158, No. 5, pp.657-660.
         */
        struct Meesungnoen2002A{
            static void GetPenetration(G4double energy,
                                       G4ThreeVector& displacement);
            static double GetRmean(double energy);
            //-----
            // Polynomial fit of Meesungnoen, 2002
            static const double gCoeff[13];
        };
        
        //-----------------------
        /*
         * Article: Terrissol M, Beaudre A (1990) Simulation of space and time
         *          evolution of radiolytic species induced by electrons in water.
         *          Radiat Prot Dosimetry 31:171–175
         */
        struct Terrisol1990A{
            static void GetPenetration(G4double energy,
                                       G4ThreeVector& displacement);
            static double GetRmean(double energy);
            static double Get3DStdDeviation(double energy);
            //-----
            // Terrisol, 1990
            static const double gEnergies_T1990[11];
            static const double gStdDev_T1990[11];
        };
        
        // --------
        /* Ballarini
         */
        struct Ballarini{
            static void GetPenetration(G4double energy,
                                       G4ThreeVector& displacement);
            static double GetRmean(double energy);
        };
    }
}

/**
 * When an electron reaches the highest energy domain of
 * G4DNAOneStepThermalizationModel,
 * it is then automatically converted into a solvated electron and displace
 * from its original position using a published thermalization statistic.
 */

template<typename MODEL=DNA::Penetration::Ballarini>
class TsTDNAOneStepThermalizationModelh : public G4VEmModel
{
public:
    typedef MODEL Model;
    TsTDNAOneStepThermalizationModelh(const G4ParticleDefinition* p = 0,
                                      const G4String& nam =
                                      "TsDNAOneStepThermalizationModel");
    virtual ~TsTDNAOneStepThermalizationModelh();
    
    virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
    
    virtual G4double CrossSectionPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* p,
                                           G4double ekin,
                                           G4double emin,
                                           G4double emax);
    
    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4DynamicParticle*,
                                   G4double tmin,
                                   G4double maxEnergy);
    
    inline void SetVerbose(int flag){
        fVerboseLevel = flag;
    }
    
    void GetPenetration(G4double energy,
                        G4ThreeVector& displacement);
    
    double GetRmean(double energy);
    
protected:
    const std::vector<G4double>* fpWaterDensity;
    
    G4ParticleChangeForGamma* fParticleChangeForGamma;
    G4bool fIsInitialised;
    G4int fVerboseLevel;
    G4Navigator* fNavigator;
    
private:
    TsTDNAOneStepThermalizationModelh&
    operator=(const TsTDNAOneStepThermalizationModelh &right);
    TsTDNAOneStepThermalizationModelh(const TsTDNAOneStepThermalizationModelh&);
};

#include "TsDNAOneStepThermalizationModelh.hh"

typedef TsTDNAOneStepThermalizationModelh<DNA::Penetration::Ballarini> TsDNAOneStepThermalizationModel;

//typedef TsTDNAOneStepThermalizationModelh<DNA::Penetration::Terrisol1990> TsDNAOneStepThermalizationModel;
// Note: if you use the above distribution, it would be
// better to follow the electrons down to 6 eV and only then apply
// the one step thermalization
#endif

