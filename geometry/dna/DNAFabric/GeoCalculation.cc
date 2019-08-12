// Extra Class for TsFiber
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

#include "GeoCalculation.hh"

#include "G4ios.hh"

struct DNAPlacementData
{
    //G4RotationMatrix *rotPosMat;
    //G4RotationMatrix *rotSolidMat;

    G4ThreeVector posCenterDNA;
    G4ThreeVector posSugarTMP1;
    G4ThreeVector posSugarTHF1;
    G4ThreeVector posBase1;
    G4ThreeVector posBase2;
    G4ThreeVector posSugarTHF2;
    G4ThreeVector posSugarTMP2;
};

GeoCalculation::GeoCalculation(G4int verbose, G4double factor) :
    fVerbose(verbose), fFactor(factor)
{

}

GeoCalculation::~GeoCalculation()
{
    // Free the memory
    //
    delete fPosNucleo;
    delete fPosDNA;
    delete fPosAndRadiusMap;
}

void GeoCalculation::Initialize()
{
    //***************************************
    // Elementary informations
    //***************************************

    //fFactor 1 vol eq
    fSugarTMPRadius = (0.270e-9)*fFactor*m;
    fSugarTHFRadius = (0.290e-9)*fFactor*m;
    //G4double ratio = 1.7;
    fSugarTMPRadiusWater = fSugarTMPRadius*1.7;
    fSugarTHFRadiusWater = fSugarTHFRadius*1.7;

    fBaseRadius = (0.300e-9)*fFactor*m;
    fBaseRadiusWater = fBaseRadius*1.7;

    // Coming from PDB barycenter
    fPosSugarTMP1 = G4ThreeVector(0.863185*fFactor*nm,-0.209463*fFactor*nm,-0.211519*fFactor*nm);
    fPosSugarTHF1 = G4ThreeVector(0.687466*fFactor*nm,0.136434*fFactor*nm,-0.103885*fFactor*nm);
    fPosBase1 = G4ThreeVector(0.334746*fFactor*nm,-0.159769*fFactor*nm,-0.0371322*fFactor*nm);
    fPosBase2 = G4ThreeVector(-0.270308*fFactor*nm,-0.0308147*fFactor*nm,0.0272545*fFactor*nm);
    fPosSugarTHF2 = G4ThreeVector(-0.712485*fFactor*nm,0.235113*fFactor*nm,0.114808*fFactor*nm);
    fPosSugarTMP2 = G4ThreeVector(-0.944741*fFactor*nm,-0.0830581*fFactor*nm,0.218929*fFactor*nm);


    //***************************************
    //***************************************
    // Calculation parameters
    //***************************************
    //***************************************

    //***************************************
    // Fiber parameters
    //***************************************

    fHistoneNum = 3; // should be 3 because the others are generated through rotations (see Construct() )
    fHistoneRadius = 2.4*fFactor*nm; //2.860*fFactor*nm;
    fHistoneHeight = 2.860*fFactor*nm; //2.860*fFactor*nm;

    // Histone helix parameters
    //
    fFiberPitch = 10*fFactor*nm;//fFiberPitch = 7.740*fFactor*nm;//fFiberPitch = 4.740*fFactor*nm;
    fFiberCentralRadius = 10.460*fFactor*nm;
    fFiberNbNuclPerTurn = 6;
    fFiberDeltaAngle = 360./fFiberNbNuclPerTurn *deg;

    //***************************************
    // DNA around histone parameters
    //***************************************

    bpNumAroundHistone = 154; //const G4int bpNumAroundHistone = 154;

    angleBpAroundHistone = 36*deg;

    // Second helix parameters (big simple helix)
    secondHelixPitch = 2.370*fFactor*nm;
    centralRadius = 4.045*fFactor*nm;
    nbBasePairPerTurn = 77;
    deltaAngle = 360./nbBasePairPerTurn *deg;

    //***************************************
    // DNA linker parameters
    //***************************************
    // Schema:
    //
    //nucleosome ---(straight part)-----\(
    //                                   \curved part
    //                                    \)
    //                                     \----(straight part)---- nucleosome

    bpNumForLinker = 46; //const G4int bpNumForLinker = 46;

    linkerCentralRadius = 14.505*fFactor*nm;
    // 2*0.01*15 is: 2 straight parts * 0.01 heigth per base in these parts * 15 bp in one straight part.
    // 3.555 is the height per base in the curved part.
    // The sum is the height per bp in the curve part of the linker.
    linkerHeightPerBp = ( -(3.555)/15)*fFactor*nm; //( -(3.555+2*0.01*15)/15)*fFactor*nm;
    nbBasePairPerTurnForLinker = 46; // do not change
    deltaLinkerAngle = 60./nbBasePairPerTurnForLinker *deg;

    //***************************************
    // Generate the positions
    //***************************************

    fPosNucleo = CalculateNucleosomePosition(fHistoneNum);

    fPosDNA = CalculateDNAPosition(fHistoneNum,fPosSugarTMP1,fPosSugarTHF1,fPosBase1,fPosBase2,fPosSugarTHF2,fPosSugarTMP2);

    fPosAndRadiusMap = GenerateCoordAndRadiusMap(fPosDNA);

    //***************************************
    // Output informations
    //***************************************

    if(fVerbose>0)
    {
        G4cout<<"********************************"<<G4endl;
        G4cout<<"fFactor="<<fFactor<<G4endl;
        G4cout<<"If fFactor=1 then the unit are correct (nanometric order)."<<G4endl;
        G4cout<<"********************************"<<G4endl;

        G4cout<<"********************************"<<G4endl;
        G4cout<<"Elementary informations"<<G4endl;
        G4cout<<"********************************"<<G4endl;
        G4cout<<"fSugarTHFRadius="<<fSugarTHFRadius/nm<<" nm"<<G4endl;
        G4cout<<"fSugarTMPRadius="<<fSugarTMPRadius/nm<<" nm"<<G4endl;
        G4cout<<"fBaseRadius="<<fBaseRadius/nm<<" nm"<<G4endl;
        G4cout<<"fPosSugarTMP1="<<fPosSugarTMP1/nm<<" nm"<<G4endl;
        G4cout<<"fPosSugarTHF1="<<fPosSugarTHF1/nm<<" nm"<<G4endl;
        G4cout<<"fPosBase1="<<fPosBase1/nm<<" nm"<<G4endl;
        G4cout<<"fPosBase2="<<fPosBase2/nm<<" nm"<<G4endl;
        G4cout<<"fPosSugarTHF2="<<fPosSugarTHF2/nm<<" nm"<<G4endl;
        G4cout<<"fPosSugarTMP2="<<fPosSugarTMP2/nm<<" nm"<<G4endl;

        G4cout<<"********************************"<<G4endl;
        G4cout<<"Important numbers"<<G4endl;
        G4cout<<"********************************"<<G4endl;
        G4cout<<"fNucleoNum="<<fNucleoNum<<G4endl;
        G4cout<<"fBpNum="<<fBpNum<<G4endl;

        G4cout<<"********************************"<<G4endl;
        G4cout<<"Fiber parameters"<<G4endl;
        G4cout<<"********************************"<<G4endl;
        G4cout<<"fHistoneRadius="<<fHistoneRadius/nm<<" nm"<<G4endl;
        G4cout<<"fHistoneHeight="<<fHistoneHeight/nm<<" nm"<<G4endl;
        G4cout<<"fFiberPitch="<<fFiberPitch/nm<<" nm"<<G4endl;
        G4cout<<"fFiberCentralRadius="<<fFiberCentralRadius/nm<<" nm"<<G4endl;
        G4cout<<"fFiberNbNuclPerTurn="<<fFiberNbNuclPerTurn<<G4endl;
        G4cout<<"fFiberDeltaAngle="<<fFiberDeltaAngle/deg<<" deg"<<G4endl;

        G4cout<<"********************************"<<G4endl;
        G4cout<<"DNA around histone parameters"<<G4endl;
        G4cout<<"********************************"<<G4endl;
        G4cout<<"bpNumAroundHistone="<<bpNumAroundHistone<<G4endl;
        G4cout<<"angleBpAroundHistone="<<angleBpAroundHistone/deg<<" deg"<<G4endl;
        G4cout<<"secondHelixPitch="<<secondHelixPitch/nm<<" nm"<<G4endl;
        G4cout<<"centralRadius="<<centralRadius/nm<<" nm"<<G4endl;
        G4cout<<"nbBasePairPerTurn="<<nbBasePairPerTurn<<G4endl;
        G4cout<<"deltaAngle="<<deltaAngle/deg<<" deg"<<G4endl;

        G4cout<<"********************************"<<G4endl;
        G4cout<<"DNA linker parameters"<<G4endl;
        G4cout<<"********************************"<<G4endl;
        G4cout<<"bpNumForLinker="<<bpNumForLinker<<G4endl;
        G4cout<<"linkerCentralRadius="<<linkerCentralRadius/nm<<" nm"<<G4endl;
        G4cout<<"linkerHeightPerBp="<<linkerHeightPerBp/nm<<" nm"<<G4endl;
        G4cout<<"nbBasePairPerTurnForLinker="<<nbBasePairPerTurnForLinker/nm<<" nm"<<G4endl;
        G4cout<<"deltaLinkerAngle="<<deltaLinkerAngle/nm<<" nm"<<G4endl;

        G4cout<<"********************************"<<G4endl;
    }
}

DNAPosData* GeoCalculation::CalculateDNAPosition(G4int histoneNum,
                                                 G4ThreeVector& posSugarTMP1,
                                                 G4ThreeVector& posSugarTHF1,
                                                 G4ThreeVector& posBase1,
                                                 G4ThreeVector& posBase2,
                                                 G4ThreeVector& posSugarTHF2,
                                                 G4ThreeVector& posSugarTMP2)
{
    DNAPosData* data = new DNAPosData;

    G4ThreeVector posOfLastBp, posOfFirstBp, posOfSecondBp;

    // Get the initial radius in xy plan
    G4double sugarTMP1Rxy = sqrt(pow(posSugarTMP1.getX(),2)+pow(posSugarTMP1.getY(),2) );
    G4double sugarTHF1Rxy = sqrt(pow(posSugarTHF1.getX(),2)+pow(posSugarTHF1.getY(),2) );
    G4double base1Rxy = sqrt(pow(posBase1.getX(),2)+pow(posBase1.getY(),2) );
    G4double base2Rxy = sqrt(pow(posBase2.getX(),2)+pow(posBase2.getY(),2) );
    G4double sugarTHF2Rxy = sqrt(pow(posSugarTHF2.getX(),2)+pow(posSugarTHF2.getY(),2) );
    G4double sugarTMP2Rxy = sqrt(pow(posSugarTMP2.getX(),2)+pow(posSugarTMP2.getY(),2) );

    // Get the initial angle in xy plan
    G4double sugarTMP1IniAngle = GetAngleToXAxis(posSugarTMP1);
    G4double sugarTHF1IniAngle = GetAngleToXAxis(posSugarTHF1);
    G4double base1IniAngle = GetAngleToXAxis(posBase1);
    G4double base2IniAngle = GetAngleToXAxis(posBase2);
    G4double sugarTHF2IniAngle = GetAngleToXAxis(posSugarTHF2);
    G4double sugarTMP2IniAngle = GetAngleToXAxis(posSugarTMP2);

    // Get the initial z pos
    G4double sugarTMP1ZIni = posSugarTMP1.getZ();
    G4double sugarTHF1ZIni = posSugarTHF1.getZ();
    G4double base1ZIni = posBase1.getZ();
    G4double base2ZIni = posBase2.getZ();
    G4double sugarTHF2ZIni = posSugarTHF2.getZ();
    G4double sugarTMP2ZIni = posSugarTMP2.getZ();


    // iterate for each  nucleosome
    for(int n=0;n<histoneNum;++n)
    {
        data->push_back(std::vector<DNAPlacementData>());

        //****************************************
        // Base pairs around the histone
        //****************************************

        // iterate on each bp around the histone
        for(int i=0;i<bpNumAroundHistone;++i)
        {
            // First DNA helix (small simple helix)
            //
            //G4double dBetweenEl = 0.33*fFactor*nm;
            G4double firstHelixMat[7][3] = {{0.,0.,0}, // CenterDNA //i*dBetweenEl
                                            {sugarTMP1Rxy*cos(i*angleBpAroundHistone/rad+sugarTMP1IniAngle/rad),sugarTMP1ZIni,sugarTMP1Rxy*sin(i*angleBpAroundHistone/rad+sugarTMP1IniAngle/rad)}, // SugarTMP1
                                            {sugarTHF1Rxy*cos(i*angleBpAroundHistone/rad+sugarTHF1IniAngle/rad),sugarTHF1ZIni,sugarTHF1Rxy*sin(i*angleBpAroundHistone/rad+sugarTHF1IniAngle/rad)}, // SugarTHF1
                                            {base1Rxy*cos(i*angleBpAroundHistone/rad+base1IniAngle/rad),base1ZIni,base1Rxy*sin(i*angleBpAroundHistone/rad+base1IniAngle/rad)}, // Base1
                                            {base2Rxy*cos(i*angleBpAroundHistone/rad+base2IniAngle/rad),base2ZIni,base2Rxy*sin(i*angleBpAroundHistone/rad+base2IniAngle/rad)}, // Base2
                                            {sugarTHF2Rxy*cos(i*angleBpAroundHistone/rad+sugarTHF2IniAngle/rad),sugarTHF2ZIni,sugarTHF2Rxy*sin(i*angleBpAroundHistone/rad+sugarTHF2IniAngle/rad)}, // SugarTHF2
                                            {sugarTMP2Rxy*cos(i*angleBpAroundHistone/rad+sugarTMP2IniAngle/rad),sugarTMP2ZIni,sugarTMP2Rxy*sin(i*angleBpAroundHistone/rad+sugarTMP2IniAngle/rad)}}; // SugarTMP2

            // Second helix matrix
            G4double secondHelixMat[3] = {centralRadius*cos(i*deltaAngle/rad), centralRadius*sin(i*deltaAngle/rad), i*secondHelixPitch/nbBasePairPerTurn};

            // Rotate first helix to be always ortho to the path of the second one
            G4double rotMat[3][3]={{cos(i*-deltaAngle/rad),-sin(i*-deltaAngle/rad),0},
                                   {sin(i*-deltaAngle/rad),cos(i*-deltaAngle/rad),0},
                                   {0,0,1}};

            G4double firstHelixRotatedMat[7][3] = {};

            // iterate to perform matrix product
            //
            // iterate on each line of output matrix
            for(int line=0;line<7;line++)
            {
                // iterate on each column of output matrix
                for(int column=0;column<3;column++)
                {
                    // iterate to do the sum
                    for(int p=0;p<3;++p)
                    {
                        firstHelixRotatedMat[line][column] += firstHelixMat[line][p]*rotMat[p][column];
                    }
                }
            }

            // Create super helix
            //
            G4double superHelixMat[7][3];

            for(int line=0;line<7;++line)
            {
                for(int column=0;column<3;++column)
                {
                    superHelixMat[line][column] = firstHelixRotatedMat[line][column] + secondHelixMat[column];
                }
            }

            // save info for linker pos (only for the first histone+DNA group)
            //
            if(i==bpNumAroundHistone-1 && n==0)
            {
                posOfLastBp = G4ThreeVector(superHelixMat[0][0],superHelixMat[0][1],superHelixMat[0][2]);
            }
            if(i==0 && n==0)
            {
                posOfFirstBp = G4ThreeVector(superHelixMat[0][0],superHelixMat[0][1],superHelixMat[0][2]);
            }
            if(i==1 && n==0)
            {
                posOfSecondBp = G4ThreeVector(superHelixMat[0][0],superHelixMat[0][1],superHelixMat[0][2]);
            }

            // Fiber part

            G4double superHelixFMat[7][3] = {};

            // iterate to perform matrix product

            // iterate on each line of output matrix+
            for(int line=0;line<7;line++)
            {
                // iterate on each column of output matrix
                for(int column=0;column<3;column++)
                {
                    // iterate to do the sum
                    for(int p=0;p<3;++p)
                    {
                        superHelixFMat[line][column] += superHelixMat[line][p]*fRotFiberMatVect[n][p][column];
                    }
                }
            }

            for(int line=0;line<7;++line)
            {
                for(int column=0;column<3;++column)
                {
                    superHelixFMat[line][column] += fFiberHelixMatVect[n][column];
                }
            }

            DNAPlacementData p;
            data->at(n).push_back(p);

            // Save the positions
            //
            data->at(n)[i].posCenterDNA = G4ThreeVector(superHelixFMat[0][0],superHelixFMat[0][1],superHelixFMat[0][2]);
            data->at(n)[i].posSugarTMP1 = G4ThreeVector(superHelixFMat[1][0],superHelixFMat[1][1],superHelixFMat[1][2]);
            data->at(n)[i].posSugarTHF1 = G4ThreeVector(superHelixFMat[2][0],superHelixFMat[2][1],superHelixFMat[2][2]);
            data->at(n)[i].posBase1 = G4ThreeVector(superHelixFMat[3][0],superHelixFMat[3][1],superHelixFMat[3][2]);
            data->at(n)[i].posBase2 = G4ThreeVector(superHelixFMat[4][0],superHelixFMat[4][1],superHelixFMat[4][2]);
            data->at(n)[i].posSugarTHF2 = G4ThreeVector(superHelixFMat[5][0],superHelixFMat[5][1],superHelixFMat[5][2]);
            data->at(n)[i].posSugarTMP2 = G4ThreeVector(superHelixFMat[6][0],superHelixFMat[6][1],superHelixFMat[6][2]);
        }

        //*********************************
        // Linker bp
        //*********************************

        // iterate on each linker bp
        for(int i=0;i<bpNumForLinker;++i)
        {
            // First DNA helix (small simple helix)
            //
            //G4double dBetweenEl = 0.33*fFactor*nm;
            G4double firstHelixMat[7][3] = {{0.,0.,0}, // CenterDNA //i*dBetweenEl
                                            {sugarTMP1Rxy*cos(i*angleBpAroundHistone/rad+sugarTMP1IniAngle/rad),sugarTMP1ZIni,sugarTMP1Rxy*sin(i*angleBpAroundHistone/rad+sugarTMP1IniAngle/rad)}, // SugarTMP1
                                            {sugarTHF1Rxy*cos(i*angleBpAroundHistone/rad+sugarTHF1IniAngle/rad),sugarTHF1ZIni,sugarTHF1Rxy*sin(i*angleBpAroundHistone/rad+sugarTHF1IniAngle/rad)}, // SugarTHF1
                                            {base1Rxy*cos(i*angleBpAroundHistone/rad+base1IniAngle/rad),base1ZIni,base1Rxy*sin(i*angleBpAroundHistone/rad+base1IniAngle/rad)}, // Base1
                                            {base2Rxy*cos(i*angleBpAroundHistone/rad+base2IniAngle/rad),base2ZIni,base2Rxy*sin(i*angleBpAroundHistone/rad+base2IniAngle/rad)}, // Base2
                                            {sugarTHF2Rxy*cos(i*angleBpAroundHistone/rad+sugarTHF2IniAngle/rad),sugarTHF2ZIni,sugarTHF2Rxy*sin(i*angleBpAroundHistone/rad+sugarTHF2IniAngle/rad)}, // SugarTHF2
                                            {sugarTMP2Rxy*cos(i*angleBpAroundHistone/rad+sugarTMP2IniAngle/rad),sugarTMP2ZIni,sugarTMP2Rxy*sin(i*angleBpAroundHistone/rad+sugarTMP2IniAngle/rad)}}; // SugarTMP2

            // Rotate first helix to be always ortho to the path of the second one
            G4double rotMat[3][3]={{cos(i*-deltaLinkerAngle/rad),-sin(i*-deltaLinkerAngle/rad),0},
                                   {sin(i*-deltaLinkerAngle/rad),cos(i*-deltaLinkerAngle/rad),0},
                                   {0,0,1}};

            // iterate to perform matrix product

            G4double firstHelixRotatedMat[7][3] = {};

            // iterate on each line of output matrix
            for(int line=0;line<7;line++)
            {
                // iterate on each column of output matrix
                for(int column=0;column<3;column++)
                {
                    // iterate to do the sum
                    for(int p=0;p<3;++p)
                    {
                        firstHelixRotatedMat[line][column] += firstHelixMat[line][p]*rotMat[p][column];
                    }
                }
            }

            linkerArcCircleMat[0] = linkerCentralRadius*cos(i*deltaLinkerAngle/rad);
            linkerArcCircleMat[1] = linkerCentralRadius*sin(i*deltaLinkerAngle/rad);

            // Corrections to increase the space between the linker and the DNA around the histone
            //
            if(i>=0 && i<15)
            {
                // to quit the histone: increase the z coord
                linkerArcCircleMat[2] = i*0.01*fFactor*nm;
            }

            else if(i>=15 && i<30)
            {
                // to do the link
                linkerArcCircleMat[2] = (i-15)*linkerHeightPerBp;
            }
            else
            {
                // to join next histone: decrease the z coord
                linkerArcCircleMat[2] = (i-30)*-0.01*fFactor*nm+15*linkerHeightPerBp;
            }

            // Create final arc equation
            //
            G4double arcFinalMat[7][3];

            for(int line=0;line<7;++line)
            {
                for(int column=0;column<3;++column)
                {
                    arcFinalMat[line][column] = firstHelixRotatedMat[line][column] + linkerArcCircleMat[column];
                }
            }

            // Create positions matrix
            //
            G4double posLinkerMat[7][3];
            for(int line=0;line<7;++line)
            {
                for(int column=0;column<3;++column)
                {
                    posLinkerMat[line][column] = arcFinalMat[line][column];

                    // x axis
                    if(column==0)
                    {
                        // Remove r distance on the X axis to make the linker start at (0,0,0)
                        posLinkerMat[line][column] -= linkerCentralRadius;
                        // Add the x position of the last bp placed around the nucleosome
                        posLinkerMat[line][column] += posOfLastBp.getX();
                    }

                    // y axis
                    else if(column==1)
                    {
                        // Add a shift to avoid placing the first bp of the linker at the same position as the last of the nucleosome
                        posLinkerMat[line][column] += (posOfSecondBp.getY() - posOfFirstBp.getY() );
                    }

                    // z axis
                    else if(column==2)
                    {
                        // Add the z position of the last bp placed around the nucleosome
                        posLinkerMat[line][column] += posOfLastBp.getZ();
                    }
                }
            }

            // Fiber part

            G4double posLinkerFMat[7][3] = {};

            // iterate to perform matrix product

            // iterate on each line of output matrix
            for(int line=0;line<7;line++)
            {
                // iterate on each column of output matrix
                for(int column=0;column<3;column++)
                {
                    // iterate to do the sum
                    for(int p=0;p<3;++p)
                    {
                        posLinkerFMat[line][column] += posLinkerMat[line][p]*fRotFiberMatVect[n][p][column];
                    }
                }
            }

            for(int line=0;line<7;++line)
            {
                for(int column=0;column<3;++column)
                {
                    posLinkerFMat[line][column] += fFiberHelixMatVect[n][column];
                }
            }

            DNAPlacementData p;
            data->at(n).push_back(p);

            // Save the positions
            data->at(n)[i+bpNumAroundHistone].posCenterDNA = G4ThreeVector(posLinkerFMat[0][0],posLinkerFMat[0][1],posLinkerFMat[0][2]);
            data->at(n)[i+bpNumAroundHistone].posSugarTMP1 = G4ThreeVector(posLinkerFMat[1][0],posLinkerFMat[1][1],posLinkerFMat[1][2]);
            data->at(n)[i+bpNumAroundHistone].posSugarTHF1 = G4ThreeVector(posLinkerFMat[2][0],posLinkerFMat[2][1],posLinkerFMat[2][2]);
            data->at(n)[i+bpNumAroundHistone].posBase1 = G4ThreeVector(posLinkerFMat[3][0],posLinkerFMat[3][1],posLinkerFMat[3][2]);
            data->at(n)[i+bpNumAroundHistone].posBase2 = G4ThreeVector(posLinkerFMat[4][0],posLinkerFMat[4][1],posLinkerFMat[4][2]);
            data->at(n)[i+bpNumAroundHistone].posSugarTHF2 = G4ThreeVector(posLinkerFMat[5][0],posLinkerFMat[5][1],posLinkerFMat[5][2]);
            data->at(n)[i+bpNumAroundHistone].posSugarTMP2 = G4ThreeVector(posLinkerFMat[6][0],posLinkerFMat[6][1],posLinkerFMat[6][2]);
        }
    }

    return data;
}

std::vector<G4ThreeVector>* GeoCalculation::CalculateNucleosomePosition(G4int nucleoNum)
{
    // Vector to fill and to output
    std::vector<G4ThreeVector>* posNucleoData = new std::vector<G4ThreeVector>(nucleoNum);

    fFiberHelixMatVect.reserve(nucleoNum);
    fRotFiberMatVect.reserve(nucleoNum);

    // Histone start coord (add something in z because of the super helix not zero centered)
    G4double histoneMat[1][3] = {{0.,0.,2.370*fFactor*nm}};

    // iterate on each nucleosome
    for(int n=0;n<nucleoNum;++n)
    {
        G4double fiberHelixMat[3];
        fiberHelixMat[0] = fFiberCentralRadius*cos(n*fFiberDeltaAngle/rad); // x
        fiberHelixMat[1] = fFiberCentralRadius*sin(n*fFiberDeltaAngle/rad); // y
        fiberHelixMat[2] = n*fFiberPitch/fFiberNbNuclPerTurn; // z

        // Save the matrix
        fFiberHelixMatVect.push_back(std::vector<G4double>() );
        fFiberHelixMatVect[n].push_back(fiberHelixMat[0]);
        fFiberHelixMatVect[n].push_back(fiberHelixMat[1]);
        fFiberHelixMatVect[n].push_back(fiberHelixMat[2]);

        // Rotation matrix (Rz) for placement inside the fiber.
        // We must rotate the objects first and then place them.
        //
        G4double rotFiberMat[3][3];
        rotFiberMat[0][0] = cos(n*-fFiberDeltaAngle/rad);rotFiberMat[0][1]=-sin(n*-fFiberDeltaAngle/rad);rotFiberMat[0][2]=0; // x0, y0, z0
        rotFiberMat[1][0] = sin(n*-fFiberDeltaAngle/rad);rotFiberMat[1][1]=cos(n*-fFiberDeltaAngle/rad);rotFiberMat[1][2]=0; // x1, y1, z1
        rotFiberMat[2][0] = 0;rotFiberMat[2][1]=0;rotFiberMat[2][2]=1; // x2, y2, z2

        // Save the matrix
        fRotFiberMatVect.push_back(std::vector<std::vector<G4double> >() );
        for(int line=0;line<3;++line)
        {
            fRotFiberMatVect[n].push_back(std::vector<G4double>() );

            for(int column=0;column<3;++column)
            {
                fRotFiberMatVect[n][line].push_back(rotFiberMat[line][column]);
            }
        }


        // Fiber part

        G4double histoneHelixFMat[1][3] = {};

        // Matrix product to rotate
        //
        // iterate on each line of output matrix+
        for(int line=0;line<1;line++)
        {
            // iterate on each column of output matrix
            for(int column=0;column<3;column++)
            {
                // iterate to do the sum
                for(int p=0;p<3;++p)
                {
                    histoneHelixFMat[line][column] += histoneMat[line][p]*rotFiberMat[p][column];
                }
            }
        }

        // Add the coord of the helix
        //
        for(int line=0;line<1;++line)
        {
            for(int column=0;column<3;++column)
            {
                histoneHelixFMat[line][column] += fiberHelixMat[column];
            }
        }

        // fill the output vector
        posNucleoData->at(n) = G4ThreeVector(histoneHelixFMat[0][0],histoneHelixFMat[0][1],histoneHelixFMat[0][2]);
    }

    return posNucleoData;
}

std::map<G4ThreeVector,G4double>* GeoCalculation::GenerateCoordAndRadiusMap(DNAPosData* dnaPosData)
{
    // To build the coord map used by the cut algorithm

    std::map<G4ThreeVector,G4double>* outMap = new std::map<G4ThreeVector,G4double>;

    G4ThreeVector coord;
    G4double radius;

    // iterate on each histone
    for(int i=0, ei=dnaPosData->size();i<ei;++i)
    {
        // iterate on each bp
        for(int j=0, ej=dnaPosData->at(i).size();j<ej;++j)
        {
            // sugarTMP 1
            coord = dnaPosData->at(i)[j].posSugarTMP1;
            radius = fSugarTMPRadiusWater;
            (*outMap)[coord] = radius;

            // sugarTHF 1
            coord = dnaPosData->at(i)[j].posSugarTHF1;
            radius = fSugarTHFRadiusWater;
            (*outMap)[coord] = radius;

            // base 1
            coord = dnaPosData->at(i)[j].posBase1;
            radius = fBaseRadiusWater;
            (*outMap)[coord] = radius;

            // base 2
            coord = dnaPosData->at(i)[j].posBase2;
            radius = fBaseRadiusWater;
            (*outMap)[coord] = radius;

            // sugarTHF 2
            coord = dnaPosData->at(i)[j].posSugarTHF2;
            radius = fSugarTHFRadiusWater;
            (*outMap)[coord] = radius;

            // sugarTMP 2
            coord = dnaPosData->at(i)[j].posSugarTMP2;
            radius = fSugarTMPRadiusWater;
            (*outMap)[coord] = radius;
        }
    }

    return outMap;
}

G4double GeoCalculation::GetAngleToXAxis(G4ThreeVector t)
{
    // Angle to x axis of the projection of t on the XY plane

    G4double x = t.getX();
    G4double y = t.getY();

    G4double dxy = sqrt(pow(x,2) + pow(y,2));

    G4double res = std::acos(x/dxy)*rad;

    if(t.getY()<0) res = -res;

    return res;
}
