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

#include "PhysGeoImport.hh"

PhysGeoImport::PhysGeoImport() :
    fIsVisu(false),
    fFactor(1.),
    fGeoName("")
{
    DefineMaterial();
}

PhysGeoImport::PhysGeoImport(bool isVisu) :
    fIsVisu(isVisu),
    fFactor(1.),
    fGeoName("")
{
    DefineMaterial();
}

PhysGeoImport::~PhysGeoImport()
{

}

G4LogicalVolume* PhysGeoImport::CreateLogicVolumeDNA(const std::string& fileName)
{
    // The idea is to return a pointer to the mother volume of the input file as output.


    // Define the colors

    G4VisAttributes red(G4Colour(1.0, 0.0, 0.0) );
    G4VisAttributes blue(G4Colour(0.0, 0.0, 1.0) );
    G4VisAttributes green(G4Colour(0.0, 1.0, 0.0) );
    G4VisAttributes yellow(G4Colour(1.0, 1.0, 0.0) );

    // Parse the input file

    ParseFile(fileName);

    // Create the volumes

    std::string boxNameSolid = fGeoName+"_solid";
    G4Box* box_solid = new G4Box(boxNameSolid, fSize/2, fSize/2, fSize/2);

    std::string boxNameLogic = fGeoName+"_logic";
    G4LogicalVolume* box_logic = new G4LogicalVolume(box_solid, fpWater, boxNameLogic);

    // sort the molecules
    std::cout<<"Sort the molecule container"<<std::endl;
    std::sort(fMolecules.begin(), fMolecules.end() );
    std::cout<<"Molecule container sorted"<<std::endl;

    // Loop on all the parsed molecules
    for(int i=0, ie=fMolecules.size(); i<ie; ++i)
    {
        // Print the number of the molecule being processed once every 100 times
        if(i%100==0)
            std::cout<<"Start molecule "<<fMolecules[i].fName<<" "<<i<<"/"<<ie<<std::endl;

        // Retrieve general molecule informations
        //
        std::string name = fMolecules[i].fName;
        std::string materialName = "water";
        double radius = fMolecules[i].fRadius;
        double waterRadius = fMolecules[i].fRadiusWater;
        G4ThreeVector moleculePosition = fMolecules[i].fPosition;
        int copyNum = fMolecules[i].fCopyNumber;

        // Water hydration shell volume part

        G4Orb* moleculeWater_solid = 0;
        G4VSolid* moleculeWaterCut_solid = 0;
        G4LogicalVolume* moleculeWater_logic = 0;

        // If water radius != 0 then we have a water hydration shell
        double tol = 0.0001;
        if(waterRadius > (0 + tol)*nm)
        {
            G4Material* mat = fpWater;

            std::string nameWaterSolid = name+"_water_solid";
            moleculeWater_solid = new G4Orb(nameWaterSolid, waterRadius);
            if(!fIsVisu) moleculeWaterCut_solid = CreateCutSolid(moleculeWater_solid, fMolecules[i], fMolecules, false);

            std::string nameWaterLogic = name+"_water_logic";
            if(!fIsVisu) moleculeWater_logic = new G4LogicalVolume(moleculeWaterCut_solid, mat, nameWaterLogic);
            else moleculeWater_logic = new G4LogicalVolume(moleculeWater_solid, mat, nameWaterLogic);

            moleculeWater_logic->SetVisAttributes(blue);

            std::string nameWaterPhys = name+"_water_phys";
            new G4PVPlacement(0, moleculePosition, moleculeWater_logic, nameWaterPhys, box_logic, false, copyNum);
        }

        // Dna volume part

        G4Orb* molecule_solid = 0;
        G4VSolid* moleculeCut_solid = 0;
        G4LogicalVolume* molecule_logic = 0;

        std::string nameSolid = fMolecules[i].fName+"_solid";
        molecule_solid = new G4Orb(nameSolid, radius);
        if(!fIsVisu) moleculeCut_solid = CreateCutSolid(molecule_solid, fMolecules[i], fMolecules, true);

        std::string nameLogic = name+"_logic";
        if(!fIsVisu) molecule_logic = new G4LogicalVolume(moleculeCut_solid, fpWater, nameLogic);
        else molecule_logic = new G4LogicalVolume(molecule_solid,fpWater, nameLogic);

        // If there was a water hydration shell volume then the current dna volume is
        // placed within it and, thus, its relative coordinates are 0,0,0.
        if(waterRadius > (0 + tol)*nm)
        {
            G4ThreeVector position(0.,0.,0.);

            std::string namePhys = name+"_phys";
            new G4PVPlacement(0, position, molecule_logic, namePhys, moleculeWater_logic, false, copyNum);
        }
        // If not, coordinates are those of the molecule
        else
        {
            G4ThreeVector position = moleculePosition;

            std::string namePhys = name+"_phys";
            new G4PVPlacement(0, position, molecule_logic, namePhys, box_logic, false, copyNum);
        }
    }

    // Clear the containers
    fMolecules.clear();
    fRadiusMap.clear();
    fWaterRadiusMap.clear();

    return box_logic;
}

void PhysGeoImport::ParseFile(const std::string& fileName)
{
    std::cout<<"Start parsing of "<<fileName<<std::endl;

    // Clear the containers
    fMolecules.clear();
    fRadiusMap.clear();
    fWaterRadiusMap.clear();

    // Setup the input stream
    std::ifstream file(fileName.c_str());

    // Check if the file was correctly opened
    if(!file.is_open())
    {
        // Geant4 exception
        G4cerr << "TOPAS is exiting due to a serious error in Geometry setup." << G4endl;
        G4cerr << "Could not open file " << fileName << G4endl;
        exit(1);
    }

    // Define the line string variable
    std::string line;

    // Read the file line per line
    while(std::getline(file, line) )
    {
        // Check the line to determine if it is empty
        if(line.empty()) continue; // skip the line if it is empty

        // Data string stream
        std::istringstream issLine(line);

        // String to determine the first letter/word
        std::string firstItem;

        // Put the first letter/word within the string
        issLine >> firstItem;

        // Check first letter to determine if the line is data or comment
        if(firstItem=="#") continue; // skip the line if it is comment

        // Use the file
        else if(firstItem=="_Name")
        {
            std::string name;
            issLine >> name;

            fGeoName = name;
        }
        else if(firstItem=="_Size")
        {
            double size;
            issLine >> size;
            size *= fFactor*nm;

            fSize = size;
        }
        else if(firstItem == "_Version")
        {

        }
        else if(firstItem=="_Number")
        {

        }
        else if(firstItem=="_Radius")
        {
            std::string name;
            issLine >> name;

            double radius;
            issLine >> radius;
            radius *= fFactor*nm;

            double waterRadius;
            issLine >> waterRadius;
            waterRadius *= fFactor*nm;

            fRadiusMap[name] = radius;
            fWaterRadiusMap[name] = waterRadius;
        }
        else if(firstItem=="_pl")
        {
            std::string name;
            issLine >> name;

            std::string material;
            issLine >> material;

            int strand;
            issLine >> strand;

            int copyNumber;
            issLine >> copyNumber;

            double x;
            issLine >> x;
            x *= fFactor*nm;

            double y;
            issLine >> y;
            y *= fFactor*nm;

            double z;
            issLine >> z;
            z *= fFactor*nm;

            Molecule1 molecule(name, copyNumber, G4ThreeVector(x, y, z), fRadiusMap[name], fWaterRadiusMap[name], material, strand);

            fMolecules.push_back(molecule);
        }
        else
        {
            // Geant4 exception
            G4cerr << "TOPAS is exiting due to a serious error in Geometry setup." << G4endl;
            G4cerr << firstItem << " is not defined in the parser. Check the input file: " << fileName  << G4endl;
            exit(1);
        }
    }

    // Close the file once the reading is done
    file.close();

    std::cout<<"End parsing of "<<fileName<<std::endl;
}


G4VSolid* PhysGeoImport::CreateCutSolid(G4Orb *solidOrbRef,
                                               Molecule1 &molRef,
                                               std::vector<Molecule1> &molList,
                                               G4bool in)
{
    // The idea behing this method is to cut overlap volumes by selecting one of them (the reference) and checking all the other volumes (the targets).
    // If a reference and a target volumes are close enough to overlap they will be cut.
    // The reference is already selected when we enter this method.

    // Use the tiny space to differentiate the frontiers (may not be necessary)
    G4double tinySpace = 0.001*fFactor*nm;

    // Cutted solid to be returned
    G4SubtractionSolid* solidCut(NULL);

    // Some flags
    G4bool isCutted = false;
    G4bool isOurVol = false;

    // Radius of the molecule to cut
    G4double radiusRef;
    if(molRef.fRadiusWater==0) radiusRef = molRef.fRadius;
    else radiusRef = molRef.fRadiusWater;

    // Reference volume position
    G4ThreeVector posRef = molRef.fPosition;

    // Before looping on all the volumes we check if the current reference volume overlaps with its container voxel boundaries.

    if(std::abs(posRef.x() ) + radiusRef > fSize/2 // along x
            || std::abs(posRef.y() ) + radiusRef > fSize/2 // along y
            || std::abs(posRef.z() ) + radiusRef > fSize/2) // along z
    {
        // If we enter here, then the reference volume overlaps with the boundaries of the container voxel

        // Box used to cut
        G4Box* solidBox = new G4Box("solid_box_for_cut", fSize/2, fSize/2, fSize/2);
        G4ThreeVector posBox;

        // Create a dummy rotation matrix
        G4RotationMatrix *rotMat = new G4RotationMatrix;
        rotMat->rotate(0, G4ThreeVector(0,0,1) );

        // To choose the cut direction

        // Up
        if(std::abs( posRef.y() + radiusRef ) > fSize/2 )
        {
           posBox = -posRef +  G4ThreeVector(0,fSize,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Down
        if(std::abs( posRef.y() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,-fSize,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
            isCutted = true;
        }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Left
        if(std::abs( posRef.x() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(fSize,0,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Right
        if(std::abs( posRef.x() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(-fSize,0,0);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Forward
        if(std::abs( posRef.z() + radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,0,fSize);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }

        // Backward
        if(std::abs( posRef.z() - radiusRef ) > fSize/2 )
        {
            posBox = -posRef + G4ThreeVector(0,0,-fSize);

            // If the volume is cutted for the first time
            if(!isCutted)
            {
                solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, posBox);
                isCutted = true;
            }
            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, posBox);
        }
    }

    // Look the other volumes of the voxel

    // Loop on all the target volumes (other volumes with potential overlaps)
    for(int i=0, ie=molList.size(); i<ie; ++i)
    {
        G4ThreeVector posTar = molList[i].fPosition;

        G4double rTar = posRef.z();
        G4double zTar = posTar.z();

        if(zTar>rTar+20*fFactor*nm)
        {
            break;
        }
        else if(zTar<rTar-20*fFactor*nm)
        {
            continue;
        }


        // Retrieve current target sphere informations
        G4double radiusTar;
        if(molList[i].fRadiusWater==0) radiusTar = molList[i].fRadius;
        else radiusTar = molList[i].fRadiusWater;



        // Compute the distance reference-target
        G4double distance = std::abs( (posRef - posTar).getR() );

        // Use the distance to check if the current target is also the reference.
        // This can only happen once per loop.
        if(distance==0 && !isOurVol)
        {
            // Target volume is also reference volume.

            // Set the flag
            isOurVol = true;

            // Next iteration
            continue;
        }
        // If the condition is correct more than one time then there is a mistake somewhere.
        else if(distance == 0 && isOurVol)
        {
            G4cerr<<"********************* Fatal Error **************************"<<G4endl;
            G4cerr<<"DetectorConstruction::CreateCutSolid: Two volumes are placed at the same position."<<G4endl;
            exit(EXIT_FAILURE);
        }

        // If the volumes are differents then we want to know if they are
        // close enough to overlap and, thus, to intiate a cut.
        else if(distance <= radiusRef+radiusTar)
        {
            // Volumes are close enough, there will be a cut

            // Box used to cut
            G4Box* solidBox = new G4Box("solid_box_for_cut", 2*radiusTar, 2*radiusTar, 2*radiusTar);

            // This part is tricky.
            // The goal is to calculate the position of the intersection center

            // diff vector to from ref to tar
            G4ThreeVector diff = posTar - posRef;

            // Find the intersection point and add to it half the length of the box used to cut
            G4double d = (pow(radiusRef,2)-pow(radiusTar,2)+pow(distance,2) ) / (2*distance) + solidBox->GetZHalfLength() - tinySpace;

            // If we are in another volume we choose to double the tinySpace to differentiate without ambiguities the inner and outer volume frontiers.
            // (may not be necessary)
            if(in) d -= 2*tinySpace;

            // Position of the box in order to achieve the cut.
            // "* ( diff/diff.getR() )" is necessary to get a vector in the right direction as output
            G4ThreeVector pos = d *( diff/diff.getR() );

            // Get the rotation angles because the box used to cut needs to be rotated
            // to give the right "cut angle".
            G4double phi = std::acos(pos.getZ()/pos.getR());
            G4double theta = std::acos( pos.getX() / ( pos.getR()*std::cos(M_PI/2.-phi) ) );

            if(pos.getY()<0) theta = -theta;

            G4ThreeVector rotAxisForPhi(1*fFactor*nm,0.,0.);
            rotAxisForPhi.rotateZ(theta+M_PI/2);

            // Create the rotation matrix
            G4RotationMatrix *rotMat = new G4RotationMatrix;
            rotMat->rotate(-phi, rotAxisForPhi);

            // Rotate it again
            G4ThreeVector rotZAxis(0.,0.,1*fFactor*nm);
            rotMat->rotate(theta, rotZAxis);

            // If the volume is cutted for the first time
            if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, pos);

            // For the other times
            else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, pos);

            // Set the cut flag
            isCutted = true;
        }
    }

    // If there was at least one cut then we return the cutted volume
    if(isCutted) return solidCut;

    // Otherwise, we return the original volume
    else return solidOrbRef;
}

void PhysGeoImport::DefineMaterial()
{
    G4NistManager * man = G4NistManager::Instance();
    fpWater = man->FindOrBuildMaterial("G4_WATER");
    fMaterialVect.push_back(fpWater);
    fVacuum = man->FindOrBuildMaterial("G4_Galactic");
    fMaterialVect.push_back(fVacuum);


}
