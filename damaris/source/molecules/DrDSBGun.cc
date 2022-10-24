// Extra Class for use by DaMaRiS
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
#include "DrDSBGun.hh"
#include "DrReaderSDD.hh"
#include "TsDamagePhaseSpaceStore.hh"
#include "DrDefinitions.hh"
#include "TsParameterManager.hh"
#include "DrDSBMoleculeManager.hh"
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

using namespace std;

DrDSBGun::DrDSBGun() {
    fMoleculeGun = new DrMoleculeGun();
    fBoundingRadius = DrDefinitions::Instance()->GetBoundingCellOrNucleusRadius();
    DrDefinitions::Instance()->fInitialBreakNumber = 0;
}

DrDSBGun::~DrDSBGun() { delete fMoleculeGun; }

//------------------------------------------------------------------------------
//----------------Externally/internally called functions------------------------
//------------------------------------------------------------------------------

void DrDSBGun::DefineAllTracks() { fMoleculeGun->DefineTracks(); }

void DrDSBGun::PlaceClockMolecule() {
    fMoleculeGun->AddMolecule("Clock", G4ThreeVector(0., 0., 0.), 0. * second);
}

//------------------------------------------------------------------------------
//----------------Placing custom distrobutions of breaks------------------------
//------------------------------------------------------------------------------

void DrDSBGun::AddMoleculeToGun(const G4ThreeVector& position, const G4String& name, const G4double& time){
    auto* breakMol = new DrBreakMolecule();
//    DrBreakTable::Instance()->NewBreakID(breakMol);

    auto DSBMoleculeManager = new DrDSBMoleculeManager();
    DSBMoleculeManager->NewBreakID(breakMol);

    breakMol->sBreakEndA->fDamageTypes[0] = 0;
    breakMol->sBreakEndA->fDamageTypes[1] = 0;
    breakMol->sBreakEndA->fDamageTypes[2] = 1;
    breakMol->sBreakEndA->fOriginalPosition = position;
    breakMol->sBreakEndA->fLesionTime = {time};
    fMoleculeGun->AddMolecule(breakMol, name);

    delete DSBMoleculeManager;
}

G4int DrDSBGun::BuildDSBEnd_Offset(G4String name, G4double offset){

    G4ThreeVector position = G4ThreeVector(offset*nm,0.0,0.0);
    G4double time = 0.0;
    AddMoleculeToGun(position, name, time);

    G4cout << "Placed a " << name << " molecule at " << position << G4endl;
    DrDefinitions::Instance()->fInitialBreakNumber = (0);

    if (!fMoleculeGun->GetMoleculeShoot().empty()) return 1;
    else return 0;
}

G4int DrDSBGun::BuildDSBEnd_Origin(G4int number, G4String name) {
    for (G4int i{0}; i < number; i++) {
        G4ThreeVector position = G4ThreeVector(0.0,0.0,0.0);
        G4double time = 0.0;
        AddMoleculeToGun(position, name, time);
    }

    DrDefinitions::Instance()->fInitialBreakNumber = G4int(floor(number / 2));
    G4cout << "Placed " << number << " " << name << "(s) at origin" << G4endl;

    if (!fMoleculeGun->GetMoleculeShoot().empty()) return 1;
    else return 0;
}

G4int DrDSBGun::BuildDSB_Column(G4int number, G4double radius_nm) {

    for (G4int i{0}; i < number; i++) {

        G4bool inNucleus{false};
        G4bool inCircle{false};
        G4double randx{DBL_MAX};
        G4double randy{DBL_MAX};
        G4double randz{DBL_MAX};

        while(!inNucleus){
            while(!inCircle){
                randx = 2 * (G4UniformRand() - 0.5) * radius_nm*nm;
                randy = 2 * (G4UniformRand() - 0.5) * radius_nm*nm;
                G4double distFromCentre = sqrt(randx*randx + randy*randy);
                if(distFromCentre <= radius_nm*nm) inCircle = true;
            }
            randz = 2 * (G4UniformRand() - 0.5) * fBoundingRadius;
            G4double distFromCentre = sqrt(randx*randx + randy*randy + randz*randz);
            if(distFromCentre <= fBoundingRadius) inNucleus = true;

            G4ThreeVector position = G4ThreeVector(randx,randy,randz);
            G4double time = 0.0;
            std::pair<G4int,G4int> dirt = {0,0};
            fMoleculeGun->AddMolecule(dirt, "DSBEnd", position, time);
            fMoleculeGun->AddMolecule(dirt, "DSBEnd", position, time);
        }
    }

    DrDefinitions::Instance()->fInitialBreakNumber = G4int(floor(number / 2));
    G4cout << "Placed " << number << " DSB ends in a column" << G4endl;

    if (!fMoleculeGun->GetMoleculeShoot().empty()) return 1;
    else return 0;
}

G4int DrDSBGun::BuildDSB_SepSpaceAndTime(G4double separation,
                                         G4int backbone,
                                         G4int base,
                                         G4int DSB_or_DSBEnd,
                                         G4double timeDelay) {

    G4double LHS = -(separation / 2.);
    G4double RHS = separation / 2.;

    if(DSB_or_DSBEnd == 0){
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", G4ThreeVector(LHS, 0., 0.), 0. * second);
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", G4ThreeVector(LHS, 0., 0.), 0. * second);

        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", G4ThreeVector(RHS, 0., 0.), timeDelay * second);
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", G4ThreeVector(RHS, 0., 0.), timeDelay * second);

        DrDefinitions::Instance()->fInitialBreakNumber = (2);
        G4cout << G4endl << "placed >2< double strand breaks "
               << G4BestUnit(separation, "Length") << " apart"<<G4endl <<"with a "
               << G4BestUnit(timeDelay, "Time") <<" delay between them."<< G4endl;
    }
    else if(DSB_or_DSBEnd == 1){
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", G4ThreeVector(LHS, 0., 0.), 0. * second);
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", G4ThreeVector(RHS, 0., 0.), timeDelay * second);
        DrDefinitions::Instance()->fInitialBreakNumber = (1);
        G4cout << G4endl << "placed >2< double strand break ends "
               << G4BestUnit(separation, "Length") << " apart." << G4endl;
    }
    if (!fMoleculeGun->GetMoleculeShoot().empty()) return 1;
    else return 0;
}

G4int DrDSBGun::BuildDSB_Pattern(G4int backbone, G4int base,
                                 vector<G4int> dsbPattern) {
    vector<G4ThreeVector> placeList;
    auto numberOfGroups = (G4int)dsbPattern.size();
    G4double xAxisDisplacement = ((G4double)(numberOfGroups - 1) / 2) * 100 * nm;
    for (G4int i = 0; i < numberOfGroups; i++) {
        for (G4int j = 0; j < dsbPattern[i]; j++)
            placeList.push_back(G4ThreeVector(xAxisDisplacement, 0, 0));
        xAxisDisplacement += 100 * nm;
    }

    for (const auto& placeThis : placeList) {
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", placeThis, 0. * second);
        fMoleculeGun->AddMolecule(make_pair(backbone, base), "DSBEnd", placeThis, 0. * second);
    }
    DrDefinitions::Instance()->fInitialBreakNumber = (G4int)placeList.size();
    G4cout << G4endl << "placed custom distribution of >" << placeList.size()
           << "< double strand break ends" << G4endl;
    return 1;
}

//------------------------------------------------------------------------------
//----------Placing from input files or other sources of damage-----------------
//------------------------------------------------------------------------------

G4int DrDSBGun::PlaceFromSTDFile(const G4String& fileName) {
    G4cout << "Placing breaks from STD format file " + fileName << G4endl;

    pair<G4int, vector<vector<vector<DrDamageEvent*> > > > R_STDInput;

    G4cout <<"Using SDDv1.0 damage format reader"<<G4endl;
    auto reader = DrReaderSDD();
    TsDamagePhaseSpaceStore* TsPSS = TsDamagePhaseSpaceStore::Instance();
    reader.ReadSTDInput(fileName,TsPSS->GetDamageHeaderPtr(),TsPSS->GetDamageEventStore());
    R_STDInput = make_pair(1,TsDamagePhaseSpaceStore::Instance()->GetDamageEventStore());

    if (!R_STDInput.second.empty()) {
        G4cout << "Processing " << R_STDInput.second.size() << " lines." << G4endl;
    } else {
        G4cerr << "ERROR: STDInput not successfully read,"
               << " event vector not populated" << G4endl;
        DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
    }

    vector< vector<DrDamageEvent*> > STDInput;
    if (R_STDInput.first == 1) {
        // exposure
        G4int selection = SelectRandomFromSTDInput(R_STDInput.second);
        STDInput = R_STDInput.second[selection];
    }
    else if (R_STDInput.first == 2) {
        // single event
        STDInput = R_STDInput.second[0];
    }
    return PlaceFromSTDInput(STDInput);
}

G4int DrDSBGun::PlaceFromSTDInput(vector<vector<DrDamageEvent*> > STDInput) {
    G4cout << "Placing breaks from STD format" << G4endl;
    return PlaceBreaks(STDInput);
}

DoubleStrandBreak DrDSBGun::ParseSTDInput(DrDamageEvent* damageEvent) {
    DoubleStrandBreak DSB;

        if (damageEvent->fDamageTypes[2] > 0) {
            if (!damageEvent->fFullBreakStructure.empty()){
                DSB = ConstructBreak(damageEvent->fFullBreakStructure);
            } else {
                G4int baseLesions = damageEvent->fDamageTypes[0];
                G4int singleStrandBreaks = damageEvent->fDamageTypes[1];
                G4int doubleStrandBreaks = damageEvent->fDamageTypes[2];

                if (doubleStrandBreaks == 1) {
                    auto RHS_BL = G4int(round(fabs(G4UniformRand()*(baseLesions + 1) - 0.5) - 0.01));
                    G4int LHS_BL = baseLesions - RHS_BL;

                    auto RHS_SSB = G4int(round(fabs(G4UniformRand()*(singleStrandBreaks + 1) - 0.5) - 0.01));
                    G4int LHS_SSB = singleStrandBreaks - RHS_SSB;

                    DSB.fBackboneBreak1 = RHS_SSB;
                    DSB.fBaseBreak1 = RHS_BL;
                    DSB.fBackboneBreak2 = LHS_SSB;
                    DSB.fBaseBreak2 = LHS_BL;
                } else {
                    G4cerr << "ERROR: Not capable of dealing with multiple DSBs"
                           << " on the same line" << G4endl;
                    DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
                }
            }
        }
    return DSB;
}

DoubleStrandBreak DrDSBGun::ConstructBreak(vector<vector<G4int>> breakVector) {
    //--------------------------------------------------------------------------
    // Finds the sites at which the cluster of damage will split into break ends
    vector<pair<G4int, G4int>> proposedBreakSites; // pair<strand1BP, strand2BP>
    G4int currentMinSeperation(INT_MAX);           // in terms of BP
    // scans the current cluster for strand 1 backbone breaks
    for (auto strand1 : breakVector) {
        if (strand1[0] == 1) {
            // scans the opposite strand and finds the closest strand 2 backbone break
            for (auto strand2 : breakVector) {
                if (strand2[0] == 4) {
                    // if smallest seperation clear old vector of break sites amd start new list with this
                    if (fabs(strand1[1] - strand2[1]) < currentMinSeperation) {
                        currentMinSeperation = G4int(fabs(strand1[1] - strand2[1]));
                        proposedBreakSites.clear();
                        proposedBreakSites.emplace_back(strand1[1], strand2[1]);
                        // if same seperation as smallest add to list of alternative break sites to randomly choose from later
                    } else if (fabs(strand1[1] - strand2[1]) == currentMinSeperation) {
                        proposedBreakSites.emplace_back(strand1[1], strand2[1]);
                    }
                }
            }
        }
    }
    // proposedBreakSites now contains the possible break sites, select one at random
    auto numberOfPossibleBreakSites = G4int(proposedBreakSites.size());
    pair<G4int, G4int> proposedBreakSite;
    G4int rand{0};
    if (numberOfPossibleBreakSites > 1) {
        rand = G4int(round(G4UniformRand() * (numberOfPossibleBreakSites)-0.5));
    }
    proposedBreakSite = proposedBreakSites[rand];
    //--------------------------------------------------------------------------
    // assign insults not involved in break site to either DSBend
    pair<G4int, G4int> breakRHS = {0, 0}; //[SSB,BL]
    pair<G4int, G4int> breakLHS = {0, 0}; //[SSB,BL]
    for (auto II : breakVector) {
        if ((II[0] == proposedBreakSite.first && II[1] == 1) ||
            (II[0] == proposedBreakSite.second && II[1] == 2)) {
        } // ie. contributing to break site
        else {
            if (II[1] == 1 && II[2] != 2) {
                if (II[0] >= proposedBreakSite.first)
                    breakRHS.first++; //+SSB to RHS
                else
                    breakLHS.first++; //+SSB to LHS
            } else if (II[1] == 2 && II[2] != 2) {
                if (II[0] >= proposedBreakSite.second)
                    breakRHS.first++; //+SSB to RHS
                else
                    breakLHS.first++; //+SSB to LHS
            } else if (II[1] == 3 && II[2] != 2) {
                if (II[0] >= proposedBreakSite.first)
                    breakRHS.second++; //+BL to RHS
                else
                    breakLHS.second++; //+BL to LHS
            } else if (II[1] == 4 && II[2] != 2) {
                if (II[0] >= proposedBreakSite.second)
                    breakRHS.second++; //+BL to RHS
                else
                    breakLHS.second++; //+BL to LHS
            }
        }
    }
    // create a DSB object to store RHS ad LHS break information
    DoubleStrandBreak DSB;
    DSB.fBackboneBreak1 = breakRHS.first;
    DSB.fBaseBreak1 = breakRHS.second;
    DSB.fBackboneBreak2 = breakLHS.first;
    DSB.fBaseBreak2 = breakLHS.second;

    return DSB;
}

G4int DrDSBGun::SelectRandomFromSTDInput(storeTypeDef input) {

    auto count = G4int(input.size());
    G4int selected{-1};

    G4cout << "Found " << count << " \"new damage events\", ";

    // in case random is not desired can
    if (DrDefinitions::Instance()->GetSelectFromExposure() != -1){
        selected = DrDefinitions::Instance()->GetSelectFromExposure();
        G4cout << " using user selection " << selected << G4endl;
    } else if (count == 0) {
        selected = 0;
        G4cout << " using " << selected << G4endl;
    } else if (count > 0) {
        selected = G4int(round(G4UniformRand() * (count)-0.5));
        G4cout << " randomly selected number " << selected << G4endl;
    } else {
        G4cerr << "ERROR: cannot select experiment from STD input! Input vector"
               << " size is: " << input.size() << G4endl;
        DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
    }
    return selected;
}

G4int DrDSBGun::PlaceBreaks(vector<vector<DrDamageEvent*> > list) {

    auto DSBMoleculeManager = new DrDSBMoleculeManager();

    G4int countDSB{0};

    for (const auto& level1 : list){
        for(auto DSB: level1){
            auto* breakMolLHS = new DrBreakMolecule();
            auto* breakMolRHS = new DrBreakMolecule();

            //@@@@ To sort out additional backbone and base lesions for each
            //@@@@ side of the double strand break.
            DoubleStrandBreak tempDSB = ParseSTDInput(DSB);
            breakMolLHS->sBreakEndA->fDamageTypes[0] = tempDSB.fBaseBreak1;
            breakMolLHS->sBreakEndA->fDamageTypes[1] = tempDSB.fBackboneBreak1;
            breakMolRHS->sBreakEndA->fDamageTypes[0] = tempDSB.fBaseBreak2;
            breakMolRHS->sBreakEndA->fDamageTypes[1] = tempDSB.fBackboneBreak2;

            SetBreakMolPropertiesFromEvent(DSB,breakMolLHS);
            SetBreakMolPropertiesFromEvent(DSB,breakMolRHS);

            //@@@@ Set correct partner IDs
            DSBMoleculeManager->NewBreakID(breakMolLHS);
            DSBMoleculeManager->NewBreakID(breakMolRHS);

            breakMolLHS->sBreakEndA->fCorrectPartnerBreakMoleculeID = breakMolRHS->sBreakEndA->fOriginalBreakMoleculeID;
            breakMolRHS->sBreakEndA->fCorrectPartnerBreakMoleculeID = breakMolLHS->sBreakEndA->fOriginalBreakMoleculeID;

            if(DrDefinitions::Instance()->GetTurnOffTime()){
                breakMolLHS->sBreakEndA->fLesionTime = {0.0};
                breakMolRHS->sBreakEndA->fLesionTime = {0.0};
            }

            if(breakMolLHS->sBreakEndA->fChromasomeID[3] == 0){
                //break is on short arm
                breakMolLHS->sBreakEndA->fTowardsCentromere = false;
                breakMolRHS->sBreakEndA->fTowardsCentromere = true;
            }
            else if(breakMolLHS->sBreakEndA->fChromasomeID[3] == 1){
                //break is in long arm
                breakMolLHS->sBreakEndA->fTowardsCentromere = true;
                breakMolRHS->sBreakEndA->fTowardsCentromere = false;
            }

            fMoleculeGun->AddMolecule(breakMolLHS,"DSBEnd");
            fMoleculeGun->AddMolecule(breakMolRHS,"DSBEnd");

            countDSB++;
        }
    }

    delete DSBMoleculeManager;

    G4cout<<"Placed "<<countDSB*2<<" double strand break ends."<<G4endl;

    DrDefinitions::Instance()->fInitialBreakNumber = (countDSB);

    if (!fMoleculeGun->GetMoleculeShoot().empty())
        return 1;
    else
        return 0;
}

void DrDSBGun::SetBreakMolPropertiesFromEvent(DrDamageEvent* damage, DrBreakMolecule* breakMol) {
    breakMol->sBreakEndA->fOriginalPosition = damage->fPosition[0]*um;
    breakMol->sBreakEndA->fChromasomeID = damage->fChromasomeID;
    breakMol->sBreakEndA->fChromasomePosition = damage->fChromasomePosition;
    breakMol->sBreakEndA->fCause = damage->fCause;
    breakMol->sBreakEndA->fDamageTypes = damage->fDamageTypes;
    breakMol->sBreakEndA->fFullBreakStructure = damage->fFullBreakStructure;
    breakMol->sBreakEndA->fDNASequence = damage->fDNASequence;
    if (!damage->fLesionTime.empty()){
        breakMol->sBreakEndA->fLesionTime = damage->fLesionTime;
    }
    else{
        breakMol->sBreakEndA->fLesionTime.emplace_back(0.0);
    }
}
