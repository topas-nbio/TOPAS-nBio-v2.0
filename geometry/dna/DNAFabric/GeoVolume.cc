// Extra Class for TsDNAFabric
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

#include "GeoVolume.hh"

#include "G4NistManager.hh"
#include "G4PVPlacement.hh"


struct DNAPlacementData
{
    G4ThreeVector posCenterDNA;
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;
};

GeoVolume::GeoVolume(G4int verbose, G4double factor) :
    fVerbose(verbose), fFactor(factor)
{
    // Water is defined from NIST material database
    fWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    // All the variables are set to -1
    // They must be defined throught the setter methods before running any Build method
    // A check based on the -1 value is runned in every Build method
    fSugarTHFRadiusWater = -1*fFactor*m;
    fSugarTMPRadiusWater = -1*fFactor*m;
    fSugarTHFRadius = -1*fFactor*m;
    fSugarTMPRadius = -1*fFactor*m;
    fBaseRadiusWater = -1*fFactor*m;
    fBaseRadius = -1*fFactor*m;
    fFiberPitch = -1*fFactor*m;
    fFiberNbNuclPerTurn = -1*fFactor*m;
    fNucleoNum = -1;
    fBpNum = -1;
    fHistoneHeight = -1*fFactor*m;
    fHistoneRadius = -1*fFactor*m;

    fpDnaMoleculePositions = new std::map<G4String, std::vector<std::vector<G4double> > >();

    //DefineMaterial();
}

GeoVolume::~GeoVolume()
{
    delete fpDnaMoleculePositions;
}

G4LogicalVolume* GeoVolume::BuildLogicFiber(std::vector<std::vector<DNAPlacementData> >* dnaVolPos,
                                            std::vector<G4ThreeVector>* posNucleo,
                                            std::map<G4ThreeVector, G4double>* posAndRadiusMap,
                                            G4bool isVisu)
{
    G4Tubs* solidFiber = new G4Tubs("Histone", 0., 17.*fFactor*nm, 80.*fFactor*nm, 0, 360);

    G4LogicalVolume* logicFiber = new G4LogicalVolume(solidFiber, fWater,"logic fiber");

    //***************************************
    // Check
    //***************************************

    if(fSugarTHFRadius==-1||fSugarTMPRadius==-1||fBaseRadius==-1||fFiberPitch==-1||fFiberNbNuclPerTurn==-1
            ||fNucleoNum==-1||fBpNum==-1||fHistoneHeight==-1||fHistoneRadius==-1)
    {
        G4cerr<<"FatalError: GeoVolume::BuildLogicFiber. A class parameter has not been initialized and its value is still negative"<<G4endl;
        std::exit(EXIT_FAILURE);
    }

    //***************************************
    // Create the histone volume
    //***************************************

    G4VisAttributes white(G4Colour(1.0, 1.0, 1.0) );

    G4Tubs* solidHistone = new G4Tubs("Histone", 0., fHistoneRadius, fHistoneHeight, 0, 360);
    G4LogicalVolume* logicHistone = new G4LogicalVolume(solidHistone,fWater,"logic histone");
    logicHistone->SetVisAttributes(white);

    //***************************************
    // Generate the cut solids
    //***************************************

    // For the positions, we only use the second (number 1) nucleosome.
    // Indeed, it is a "middle" nucleosome and, thus, the two extemities will be cutted.
    // If we took the first, we would have a first volume non cutted and, thus, some overlaps.
    std::vector<DNAPlacementData>* nuclVolPos = &dnaVolPos->at(1);

    // We create here all the DNA volumes around the histone based on the second nucleosome (number 1) positions.
    // Only the volumes of the second nucleosome are generated. By placing them several times we will build the fiber.
    // This is done to save memory and improve speed.
    // We saved the volumes in a map.
    std::map<G4String, std::vector<G4LogicalVolume*> >* volMap
            = CreateNucleosomeCuttedSolidsAndLogicals(nuclVolPos, posAndRadiusMap, isVisu);

    // Print the mean volume number of the DNA volumes
    if(fVerbose>2)
    {
        CalculateMeanVol(volMap);
    }

    //*********************************************************
    // Save the positions of the first nucleosome in vectors
    //*********************************************************

    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;
    std::vector<G4ThreeVector> posSugarTMP1Vect;
    std::vector<G4ThreeVector> posSugarTHF1Vect;
    std::vector<G4ThreeVector> posBase1Vect;
    std::vector<G4ThreeVector> posBase2Vect;
    std::vector<G4ThreeVector> posSugarTHF2Vect;
    std::vector<G4ThreeVector> posSugarTMP2Vect;

    // iterate on each base pair around one nucleosome
    for(int j=0;j<fBpNum;++j)
    {
        // Get the base volume positions of the first nucleosome
        posSugarTMP1 = dnaVolPos->at(0)[j].posSugarTMP1;
        posSugarTHF1 = dnaVolPos->at(0)[j].posSugarTHF1;
        posBase1 = dnaVolPos->at(0)[j].posBase1;
        posBase2 = dnaVolPos->at(0)[j].posBase2;
        posSugarTHF2 = dnaVolPos->at(0)[j].posSugarTHF2;
        posSugarTMP2 = dnaVolPos->at(0)[j].posSugarTMP2;

        // Save each of them in a vector
        posSugarTMP1Vect.push_back(posSugarTMP1);
        posSugarTHF1Vect.push_back(posSugarTHF1);
        posBase1Vect.push_back(posBase1);
        posBase2Vect.push_back(posBase2);
        posSugarTHF2Vect.push_back(posSugarTHF2);
        posSugarTMP2Vect.push_back(posSugarTMP2);
    }

    // Save the first histone position
    G4ThreeVector posHistone = posNucleo->at(0);

    //***************************************
    // Do the placements
    //***************************************

    // Distance value used to place the volumes at the beginning of the fiber
    G4ThreeVector minusForFiber = G4ThreeVector(0.,0.,-solidFiber->GetDz() + fHistoneHeight);

    // Build the nucleosome helix inside the fiber

    G4double zShift = fFiberPitch/fFiberNbNuclPerTurn;
    G4int count = 0;

    // iterate on each histone
    for(int i=0;i<fNucleoNum;++i)
    {
        G4RotationMatrix* rotObj = new G4RotationMatrix();
        rotObj->rotateZ(i*-fFiberDeltaAngle + fFiberDeltaAngle);

        // iterate on each bp
        for(int j=0;j<fBpNum;++j)
        {
            // Rotate position vectors after the the first nucleosome placement
            if(i!=0)
            {
                posSugarTMP1 = posSugarTMP1Vect[j].rotateZ(fFiberDeltaAngle);
                posSugarTHF1 = posSugarTHF1Vect[j].rotateZ(fFiberDeltaAngle);
                posBase1 = posBase1Vect[j].rotateZ(fFiberDeltaAngle);
                posBase2 = posBase2Vect[j].rotateZ(fFiberDeltaAngle);
                posSugarTHF2 = posSugarTHF2Vect[j].rotateZ(fFiberDeltaAngle);
                posSugarTMP2 = posSugarTMP2Vect[j].rotateZ(fFiberDeltaAngle);
            }
            // If it is the first nucleosome placement then just use the coordinates of the vectors (no rotation)
            else
            {
                posSugarTMP1 = posSugarTMP1Vect[j];
                posSugarTHF1 = posSugarTHF1Vect[j];
                posBase1 = posBase1Vect[j];
                posBase2 = posBase2Vect[j];
                posSugarTHF2 = posSugarTHF2Vect[j];
                posSugarTMP2 = posSugarTMP2Vect[j];
            }

            // Add the z shift to build the helix
            posSugarTMP1 += G4ThreeVector(0.,0.,i*zShift);
            posSugarTHF1 += G4ThreeVector(0.,0.,i*zShift);
            posBase1 += G4ThreeVector(0.,0.,i*zShift);
            posBase2 += G4ThreeVector(0.,0.,i*zShift);
            posSugarTHF2 += G4ThreeVector(0.,0.,i*zShift);
            posSugarTMP2 += G4ThreeVector(0.,0.,i*zShift);

            // Place into the fiber (minus to start at the beginning (negative coord))
            //
            posSugarTMP1 += minusForFiber;
            posSugarTHF1 += minusForFiber;
            posBase1 += minusForFiber;
            posBase2 += minusForFiber;
            posSugarTHF2 += minusForFiber;
            posSugarTMP2 += minusForFiber;

            // Placements
            new G4PVPlacement(0,G4ThreeVector(),volMap->at("sugarTMP1")[j],"backboneTMP1",volMap->at("sugarTMP1Water")[j],false,count);
            (*fpDnaMoleculePositions)["Phosphate"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getX());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getY());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP1.getZ());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(count);
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(1);
            new G4PVPlacement(0,G4ThreeVector(),volMap->at("sugarTHF1")[j],"backboneTHF1",volMap->at("sugarTHF1Water")[j],false,count);
            (*fpDnaMoleculePositions)["Desoxyribose"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getX());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getY());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF1.getZ());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(count);
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(1);

            if(j%2) // odd
            {
                new G4PVPlacement(0,G4ThreeVector(),volMap->at("base1")[j],"base_cytosine",volMap->at("base1Water")[j],false,count);
                (*fpDnaMoleculePositions)["Cytosine"].push_back(std::vector<double>());
                (*fpDnaMoleculePositions)["Cytosine"].back().push_back(posBase1.getX());
                (*fpDnaMoleculePositions)["Cytosine"].back().push_back(posBase1.getY());
                (*fpDnaMoleculePositions)["Cytosine"].back().push_back(posBase1.getZ());
                (*fpDnaMoleculePositions)["Cytosine"].back().push_back(count);
                (*fpDnaMoleculePositions)["Cytosine"].back().push_back(1);
                new G4PVPlacement(0,G4ThreeVector(),volMap->at("base2")[j],"base_guanine",volMap->at("base2Water")[j],false,count);
                (*fpDnaMoleculePositions)["Guanine"].push_back(std::vector<double>());
                (*fpDnaMoleculePositions)["Guanine"].back().push_back(posBase2.getX());
                (*fpDnaMoleculePositions)["Guanine"].back().push_back(posBase2.getY());
                (*fpDnaMoleculePositions)["Guanine"].back().push_back(posBase2.getZ());
                (*fpDnaMoleculePositions)["Guanine"].back().push_back(count);
                (*fpDnaMoleculePositions)["Guanine"].back().push_back(2);
            }
            else // even
            {
                new G4PVPlacement(0,G4ThreeVector(),volMap->at("base1")[j],"base_thymine",volMap->at("base1Water")[j],false,count);
                (*fpDnaMoleculePositions)["Thymine"].push_back(std::vector<double>());
                (*fpDnaMoleculePositions)["Thymine"].back().push_back(posBase1.getX());
                (*fpDnaMoleculePositions)["Thymine"].back().push_back(posBase1.getY());
                (*fpDnaMoleculePositions)["Thymine"].back().push_back(posBase1.getZ());
                (*fpDnaMoleculePositions)["Thymine"].back().push_back(count);
                (*fpDnaMoleculePositions)["Thymine"].back().push_back(1);
                new G4PVPlacement(0,G4ThreeVector(),volMap->at("base2")[j],"base_adenine",volMap->at("base2Water")[j],false,count);
                (*fpDnaMoleculePositions)["Adenine"].push_back(std::vector<double>());
                (*fpDnaMoleculePositions)["Adenine"].back().push_back(posBase2.getX());
                (*fpDnaMoleculePositions)["Adenine"].back().push_back(posBase2.getY());
                (*fpDnaMoleculePositions)["Adenine"].back().push_back(posBase2.getZ());
                (*fpDnaMoleculePositions)["Adenine"].back().push_back(count);
                (*fpDnaMoleculePositions)["Adenine"].back().push_back(2);
            }

            new G4PVPlacement(0,G4ThreeVector(),volMap->at("sugarTHF2")[j],"backboneTHF2",volMap->at("sugarTHF2Water")[j],false,count);
            (*fpDnaMoleculePositions)["Desoxyribose"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getX());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getY());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(posSugarTHF2.getZ());
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(count);
            (*fpDnaMoleculePositions)["Desoxyribose"].back().push_back(2);
            new G4PVPlacement(0,G4ThreeVector(),volMap->at("sugarTMP2")[j],"backboneTMP2",volMap->at("sugarTMP2Water")[j],false,count);
            (*fpDnaMoleculePositions)["Phosphate"].push_back(std::vector<double>());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getX());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getY());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(posSugarTMP2.getZ());
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(count);
            (*fpDnaMoleculePositions)["Phosphate"].back().push_back(2);

            new G4PVPlacement(0,posSugarTMP1,volMap->at("sugarTMP1Water")[j],"sugarTMP1Hydra",logicFiber,false,count);
            new G4PVPlacement(0,posSugarTHF1,volMap->at("sugarTHF1Water")[j],"sugarTHF1Hydra",logicFiber,false,count);
            new G4PVPlacement(0,posBase1,volMap->at("base1Water")[j],"base1Hydra",logicFiber,false,count);
            new G4PVPlacement(0,posBase2,volMap->at("base2Water")[j],"base2Hydra",logicFiber,false,count);
            new G4PVPlacement(0,posSugarTHF2,volMap->at("sugarTHF2Water")[j],"sugarTHF2Hydra",logicFiber,false,count);
            new G4PVPlacement(0,posSugarTMP2,volMap->at("sugarTMP2Water")[j],"sugarTMP2Hydra",logicFiber,false,count);

            ++count;
        }

        // Place the histone volume
        //
        // Rotate
        G4ThreeVector posHistoneForNucleo = posHistone;
        posHistoneForNucleo.rotateZ(i*fFiberDeltaAngle);
        // Apply the z shift
        posHistoneForNucleo += G4ThreeVector(0.,0.,i*zShift);
        // Place into the fiber (minus to start at the beginning (negative coord)
        posHistoneForNucleo += minusForFiber;
        // Place
        new G4PVPlacement(0,posHistoneForNucleo,logicHistone,"histone",logicFiber,true,i);
        (*fpDnaMoleculePositions)["Histone"].push_back(std::vector<double>());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getX());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getY());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(posHistoneForNucleo.getZ());
        (*fpDnaMoleculePositions)["Histone"].back().push_back(i);
        (*fpDnaMoleculePositions)["Histone"].back().push_back(0);
    }

    return logicFiber;
}


std::map<G4String, std::vector<G4LogicalVolume*> >* GeoVolume::CreateNucleosomeCuttedSolidsAndLogicals(std::vector<DNAPlacementData>* nucleosomeVolumePositions,
                                                                                                       std::map<G4ThreeVector, G4double>* posAndRadiusMap,
                                                                                                       G4bool isVisu)
{
    std::map<G4String, std::vector<G4LogicalVolume*> >* logicSolidsMap = new std::map<G4String, std::vector<G4LogicalVolume*> >;

    G4int basePairNum = nucleosomeVolumePositions->size();

    //***************************************
    // Create elementary solids
    //***************************************

    if(fSugarTHFRadius==-1 || fSugarTMPRadius==-1 || fBaseRadius==-1)
    {
        G4cerr<<"************************************************************"<<G4endl;
        G4cerr<<"GeoVolume::CreateNucleosomeCuttedSolidsAndLogicals: fSugarTHFRadius, fSugarTMPRadius or fBaseRadius were not set. Fatal error."<<G4endl;
        G4cerr<<"************************************************************"<<G4endl;
        std::exit(EXIT_FAILURE);
    }

    G4VisAttributes red(G4Colour(1.0, 0.0, 0.0) );
    G4VisAttributes blue(G4Colour(0.0, 0.0, 1.0) );
    G4VisAttributes green(G4Colour(0.0, 1.0, 0.0) );
    G4VisAttributes yellow(G4Colour(1.0, 1.0, 0.0) );

    G4Orb* solidSugarTHFWater = new G4Orb("solid_sugar_THF_Water", fSugarTHFRadiusWater);
    G4Orb* solidSugarTMPWater = new G4Orb("solid_sugar_TMP_Water", fSugarTMPRadiusWater);
    G4Orb* solidBaseWater = new G4Orb("solid_base_Water", fBaseRadiusWater);

    G4Orb* solidSugarTHF = new G4Orb("solid_sugar_THF", fSugarTHFRadius);
    G4Orb* solidSugarTMP = new G4Orb("solid_sugar_TMP", fSugarTMPRadius);
    G4Orb* solidBase = new G4Orb("solid_base", fBaseRadius);

    //***************************************
    // Generate the cutted solids
    //***************************************
    
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;
    
    // iterate on each base pair
    for(int j=0;j<basePairNum;++j)
    {
        // Get the position
        posSugarTMP1 = nucleosomeVolumePositions->at(j).posSugarTMP1;
        posSugarTHF1 = nucleosomeVolumePositions->at(j).posSugarTHF1;
        posBase1 = nucleosomeVolumePositions->at(j).posBase1;
        posBase2 = nucleosomeVolumePositions->at(j).posBase2;
        posSugarTHF2 = nucleosomeVolumePositions->at(j).posSugarTHF2;
        posSugarTMP2 = nucleosomeVolumePositions->at(j).posSugarTMP2;

        // Create the solid volumes

        // Water layer solids
        //
        G4VSolid* sugarTMP1Water;
        G4VSolid* sugarTHF1Water;
        G4VSolid* base1Water;
        G4VSolid* base2Water;
        G4VSolid* sugarTHF2Water;
        G4VSolid* sugarTMP2Water;

        // DNA solids
        //
        G4VSolid* sugarTMP1;
        G4VSolid* sugarTHF1;
        G4VSolid* base1;
        G4VSolid* base2;
        G4VSolid* sugarTHF2;
        G4VSolid* sugarTMP2;

        // if isVisu is false it means we will run calculations so we need the cutted volumes
        if(!isVisu)
        {
            // Generate the cutted solids

            // Water layer solids
            //
            sugarTMP1Water = CreateCutSolid(solidSugarTMPWater,posSugarTMP1,posAndRadiusMap);
            sugarTHF1Water = CreateCutSolid(solidSugarTHFWater,posSugarTHF1,posAndRadiusMap);
            base1Water = CreateCutSolid(solidBaseWater,posBase1,posAndRadiusMap);
            base2Water = CreateCutSolid(solidBaseWater,posBase2,posAndRadiusMap);
            sugarTHF2Water = CreateCutSolid(solidSugarTHFWater,posSugarTHF2,posAndRadiusMap);
            sugarTMP2Water = CreateCutSolid(solidSugarTMPWater,posSugarTMP2,posAndRadiusMap);

            // DNA solids
            //
            sugarTMP1 = CreateCutSolid(solidSugarTMP,posSugarTMP1,posAndRadiusMap, "sugarTMP", true);
            sugarTHF1 = CreateCutSolid(solidSugarTHF,posSugarTHF1,posAndRadiusMap, "sugarTHF", true);
            base1 = CreateCutSolid(solidBase,posBase1,posAndRadiusMap, "base", true);
            base2 = CreateCutSolid(solidBase,posBase2,posAndRadiusMap, "base", true);
            sugarTHF2 = CreateCutSolid(solidSugarTHF,posSugarTHF2,posAndRadiusMap, "sugarTHF", true);
            sugarTMP2 = CreateCutSolid(solidSugarTMP,posSugarTMP2,posAndRadiusMap, "sugarTMP", true);
        }
        // if isVisu is true it means we just want to visualize the geometry so we do not need the cutted volumes
        else
        {
            sugarTMP1Water = solidSugarTMPWater;
            sugarTHF1Water = solidSugarTHFWater;
            base1Water = solidBaseWater;
            base2Water = solidBaseWater;
            sugarTHF2Water = solidSugarTHFWater;
            sugarTMP2Water = solidSugarTMPWater;

            sugarTMP1 = solidSugarTMP;
            sugarTHF1 = solidSugarTHF;
            base1 = solidBase;
            base2 = solidBase;
            sugarTHF2 = solidSugarTHF;
            sugarTMP2 = solidSugarTMP;
        }

        // Create the logical volumes
        //
        G4LogicalVolume* logicSugarTMP1Water;
        G4LogicalVolume* logicSugarTHF1Water;
        G4LogicalVolume* logicBase1Water;
        G4LogicalVolume* logicBase2Water;
        G4LogicalVolume* logicSugarTHF2Water;
        G4LogicalVolume* logicSugarTMP2Water;
        G4LogicalVolume* logicSugarTHF1;
        G4LogicalVolume* logicSugarTMP1;
        G4LogicalVolume* logicBase1;
        G4LogicalVolume* logicBase2;
        G4LogicalVolume* logicSugarTHF2;
        G4LogicalVolume* logicSugarTMP2;

        // Create the logical volumes
        //
        logicSugarTMP1Water = new G4LogicalVolume(sugarTMP1Water,fWater,"logic_sugarTMP_1_hydra_Backbone1");
        logicSugarTHF1Water = new G4LogicalVolume(sugarTHF1Water,fWater,"logic_sugarTHF_1_hydra_Backbone1");
        logicBase1Water = new G4LogicalVolume(base1Water, fWater,"logic_base_1_hydra_Base1");
        logicBase2Water = new G4LogicalVolume(base2Water, fWater,"logic_base_2_hydra_Base2");
        logicSugarTHF2Water = new G4LogicalVolume(sugarTHF2Water,fWater,"logic_sugarTHF_2_hydra_Backbone2");
        logicSugarTMP2Water = new G4LogicalVolume(sugarTMP2Water,fWater,"logic_sugarTMP_2_hydra_Backbone2");

        // Only water here
        //
        logicSugarTMP1 = new G4LogicalVolume(sugarTMP1,fWater,"logic_sugar_TMP_1_Backbone1");
        logicSugarTHF1 = new G4LogicalVolume(sugarTHF1,fWater,"logic_sugar_THF_1_Backbone1");
        if(j%2) // odd
        {
            logicBase1 = new G4LogicalVolume(base1,fWater,"logic_base_cytosine_Base1"); // PY
            logicBase2 = new G4LogicalVolume(base2,fWater,"logic_base_guanine_Base2"); // PU
        }
        else // even
        {
            logicBase1 = new G4LogicalVolume(base1,fWater,"logic_base_thymine_Base1"); // PY
            logicBase2 = new G4LogicalVolume(base2,fWater,"logic_base_adenine_Base2"); // PU
        }
        logicSugarTHF2 = new G4LogicalVolume(sugarTHF2,fWater,"logic_sugar_THF_2_Backbone2");
        logicSugarTMP2 = new G4LogicalVolume(sugarTMP2,fWater,"logic_sugar_TMP_2_Backbone2");

        /*
        logicSugarTMP1 = new G4LogicalVolume(sugarTMP1,backbone_TMP,"logic_sugar_TMP_1");
        logicSugarTHF1 = new G4LogicalVolume(sugarTHF1,backbone_THF,"logic_sugar_THF_1");
        if(j%2) // odd
        {
            logicBase1 = new G4LogicalVolume(base1,cytosine_PY,"logic_base_cytosine"); // PY
            logicBase2 = new G4LogicalVolume(base2,guanine_PU,"logic_base_guanine"); // PU
        }
        else // even
        {
            logicBase1 = new G4LogicalVolume(base1,thymine_PY,"logic_base_thymine"); // PY
            logicBase2 = new G4LogicalVolume(base2,adenine_PU,"logic_base_adenine"); // PU
        }
        logicSugarTHF2 = new G4LogicalVolume(sugarTHF2,backbone_THF,"logic_sugar_THF_2");
        logicSugarTMP2 = new G4LogicalVolume(sugarTMP2,backbone_TMP,"logic_sugar_TMP_2");
        */

        // Set the colors
        //
        logicSugarTMP1Water->SetVisAttributes(blue);
        logicSugarTHF1Water->SetVisAttributes(blue);
        logicBase1Water->SetVisAttributes(blue);
        logicBase2Water->SetVisAttributes(blue);
        logicSugarTHF2Water->SetVisAttributes(blue);
        logicSugarTMP2Water->SetVisAttributes(blue);

        logicSugarTMP1->SetVisAttributes(green);
        logicSugarTHF1->SetVisAttributes(yellow);
        logicBase1->SetVisAttributes(red);
        logicBase2->SetVisAttributes(red);
        logicSugarTHF2->SetVisAttributes(yellow);
        logicSugarTMP2->SetVisAttributes(green);

        // Save the logical volumes in the output map
        //
        (*logicSolidsMap)["sugarTMP1WaterBackbone1"].push_back(logicSugarTMP1Water);
        (*logicSolidsMap)["sugarTHF1WaterBackbone1"].push_back(logicSugarTHF1Water);
        (*logicSolidsMap)["base1WaterBase1"].push_back(logicBase1Water);
        (*logicSolidsMap)["base2WaterBase2"].push_back(logicBase2Water);
        (*logicSolidsMap)["sugarTHF2WaterBackbone2"].push_back(logicSugarTHF2Water);
        (*logicSolidsMap)["sugarTMP2WaterBackbone2"].push_back(logicSugarTMP2Water);

        (*logicSolidsMap)["sugarTMP1Backbone1"].push_back(logicSugarTMP1);
        (*logicSolidsMap)["sugarTHF1Backbone1"].push_back(logicSugarTHF1);
        (*logicSolidsMap)["base1Base1"].push_back(logicBase1);
        (*logicSolidsMap)["base2Base2"].push_back(logicBase2);
        (*logicSolidsMap)["sugarTHF2Backbone2"].push_back(logicSugarTHF2);
        (*logicSolidsMap)["sugarTMP2Backbone2"].push_back(logicSugarTMP2);
    }

    return logicSolidsMap;
}

G4VSolid* GeoVolume::CreateCutSolid(G4Orb *solidOrbRef,
                                               G4ThreeVector& posRef,
                                               std::map<G4ThreeVector,G4double>* tarMap,
                                               G4String volName,
                                               G4bool in)
{
    G4SubtractionSolid* solidCut(NULL);

    bool isCutted = false;
    bool isOurVol = false;

    G4double radiusRef (0);

    if(volName=="base") radiusRef = fBaseRadiusWater;
    else if(volName=="sugarTHF") radiusRef = fSugarTHFRadiusWater;
    else if(volName=="sugarTMP") radiusRef = fSugarTMPRadiusWater;
    // Ref is the solid on which we focus
    else radiusRef = solidOrbRef->GetRadius();

    std::map<G4ThreeVector,G4double>::iterator it;
    std::map<G4ThreeVector,G4double>::iterator ite;
    // iterate on all the target volumes
    for(it=tarMap->begin(), ite=tarMap->end();it!=ite;++it)
    {
        G4double radiusTar = it->second;
        G4ThreeVector posTar = it->first;
        G4double distance = std::abs( (posRef-posTar).getR() );

        // Check if target volume = reference volume (can only happen once)
        if(distance == 0 && !isOurVol)
        {
            // target volume is the reference volume
            isOurVol = true;
        }
        else if(distance == 0 && isOurVol)
        {
            G4cerr<<"DetectorConstruction::CreateCutSolid, Fatal Error. Two volumes are placed at the same position."<<G4endl;
            exit(EXIT_FAILURE);
        }
        // Check if the volumes are close enough to allow the cut
        else if(distance <= radiusRef+radiusTar)
        {
            // we allow the cut

            G4Box* solidBox = new G4Box("solid box for cut", 2*radiusTar, 2*radiusTar, 2*radiusTar);

            // To calculate the position of the intersection center
            //
            // diff vector to from ref to tar
            G4ThreeVector diff = posTar - posRef;
            // Find the intersection point
            G4double d = (pow(radiusRef,2)-pow(radiusTar,2)+pow(distance,2) ) / (2*distance) + solidBox->GetZHalfLength() - 0.001*fFactor*nm;
            if(in) d -= 0.002*fFactor*nm;
            // "* ( diff/diff.getR() )" is necessary to get a vector in the right direction as output
            G4ThreeVector pos = d *( diff/diff.getR() );

            G4double phi = std::acos(pos.getZ()/pos.getR());
            G4double theta = std::acos( pos.getX() / ( pos.getR()*std::cos(M_PI/2.-phi) ) );

            if(pos.getY()<0) theta = -theta;

            G4ThreeVector rotAxisForPhi(1*fFactor*nm,0.,0.);
            rotAxisForPhi.rotateZ(theta+M_PI/2);
            G4RotationMatrix *rotMat = new G4RotationMatrix;
            rotMat->rotate(-phi, rotAxisForPhi);

            G4ThreeVector rotZAxis(0.,0.,1*fFactor*nm);
            rotMat->rotate(theta, rotZAxis);

            // first time
            if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, pos);
            // other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, pos);

            isCutted = true;
        }
    }

    if(isCutted) return solidCut;
    else return solidOrbRef;
}

void GeoVolume::CalculateMeanVol(std::map<G4String, std::vector<G4LogicalVolume*> >* logicSolidsMap)
{
    G4double sugarTMP1Vol = 0;
    G4double sugarTHF1Vol = 0;
    G4double base1Vol = 0;
    G4double base2Vol = 0;
    G4double sugarTHF2Vol = 0;
    G4double sugarTMP2Vol = 0;

    G4double sugarTMP1WaterVol = 0;
    G4double sugarTHF1WaterVol = 0;

    G4double base1WaterVol = 0;
    G4double base2WaterVol = 0;
    G4double sugarTHF2WaterVol = 0;
    G4double sugarTMP2WaterVol = 0;

    // iterate on each bp
    for(int j=0;j<fBpNum;++j)
    {
        if(j==0) G4cout<<"Volume calculations\nThe fFactor value is taken into account already. We assume fFactor=1e+9."<<G4endl;
        assert(fFactor==1e+9);
        sugarTMP1Vol += logicSolidsMap->at("sugarTMP1")[j]->GetSolid()->GetCubicVolume();
        sugarTHF1Vol += logicSolidsMap->at("sugarTHF1")[j]->GetSolid()->GetCubicVolume();
        base1Vol += logicSolidsMap->at("base1")[j]->GetSolid()->GetCubicVolume();
        base2Vol += logicSolidsMap->at("base2")[j]->GetSolid()->GetCubicVolume();
        sugarTHF2Vol += logicSolidsMap->at("sugarTHF2")[j]->GetSolid()->GetCubicVolume();
        sugarTMP2Vol += logicSolidsMap->at("sugarTMP2")[j]->GetSolid()->GetCubicVolume();

        sugarTMP1WaterVol += logicSolidsMap->at("sugarTMPWater1")[j]->GetSolid()->GetCubicVolume();
        sugarTHF1WaterVol += logicSolidsMap->at("sugarTHFWater1")[j]->GetSolid()->GetCubicVolume();
        base1WaterVol += logicSolidsMap->at("baseWater1")[j]->GetSolid()->GetCubicVolume();
        base2WaterVol += logicSolidsMap->at("baseWater2")[j]->GetSolid()->GetCubicVolume();
        sugarTHF2WaterVol += logicSolidsMap->at("sugarTHFWater2")[j]->GetSolid()->GetCubicVolume();
        sugarTMP2WaterVol += logicSolidsMap->at("sugarTMPWater2")[j]->GetSolid()->GetCubicVolume();
    }

    sugarTMP1Vol = sugarTMP1Vol/(fBpNum*pow(fFactor,3) );
    sugarTHF1Vol = sugarTHF1Vol/(fBpNum*pow(fFactor,3) );
    base1Vol = base1Vol/(fBpNum*pow(fFactor,3));
    base2Vol = base2Vol/(fBpNum*pow(fFactor,3));
    sugarTHF2Vol = sugarTHF2Vol/(fBpNum*pow(fFactor,3));
    sugarTMP2Vol = sugarTMP2Vol/(fBpNum*pow(fFactor,3));

    sugarTMP1WaterVol = sugarTMP1WaterVol/(fBpNum*pow(fFactor,3) ) - sugarTMP1Vol;
    sugarTHF1WaterVol = sugarTHF1WaterVol/(fBpNum*pow(fFactor,3) ) - sugarTHF1Vol;
    base1WaterVol = base1WaterVol/(fBpNum*pow(fFactor,3)) - base1Vol;
    base2WaterVol = base2WaterVol/(fBpNum*pow(fFactor,3)) - base2Vol;
    sugarTHF2WaterVol = sugarTHF2WaterVol/(fBpNum*pow(fFactor,3)) - sugarTHF2Vol;
    sugarTMP2WaterVol = sugarTMP2WaterVol/(fBpNum*pow(fFactor,3)) - sugarTMP2Vol;

    G4cout<<"\nsugarTMP1Vol="<<sugarTMP1Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"\nsugarTHF1Vol="<<sugarTHF1Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"base1Vol="<<base1Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"base2Vol="<<base2Vol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"sugarTHF2Vol="<<sugarTHF2Vol/m3*1e+27<<" nm3\n"<<G4endl;
    G4cout<<"sugarTMP2Vol="<<sugarTMP2Vol/m3*1e+27<<" nm3\n"<<G4endl;

    G4cout<<"\nsugarTMP1WaterVol="<<sugarTMP1WaterVol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"\nsugarTHF1WaterVol="<<sugarTHF1WaterVol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"base1WaterVol="<<base1WaterVol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"base2WaterVol="<<base2WaterVol/m3*1e+27<<" nm3"<<G4endl;
    G4cout<<"sugarTHF2WaterVol="<<sugarTHF2WaterVol/m3*1e+27<<" nm3\n"<<G4endl;
    G4cout<<"sugarTMP2WaterVol="<<sugarTMP2WaterVol/m3*1e+27<<" nm3\n"<<G4endl;
}


