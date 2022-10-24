// Component for TsSolenoidFiber
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
/*
 *
 *  Created on: 15 Feb 2016
 *      Author: Nick Henthorn
 */
// Geometry for a linear segment of the chromatin Fiber
// The Fiber is a solenoid conformation of cylinrical histones,
// wrapped by 1.65 left handed turns of the double helix
// The Histones are placed, recording position and rotation
// the histone pos and rotation is passed to the DNA path and a general path is made
// the general path is segmented in 0.34nm steps
// the DNA position and rotation is created and DNA placed

#include "TsSolenoidFiber.hh"

#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VisAttributes.hh"


using namespace CLHEP;
using namespace std;

TsSolenoidFiber::TsSolenoidFiber(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsSolenoidFiber::~TsSolenoidFiber()
{
}



G4VPhysicalVolume* TsSolenoidFiber::Construct()
{
    BeginConstruction();

    G4double ChromatinRadius = (37.0879/2.0) *nm;
    G4double ChromatinLength = (198.0/2.0) *nm;

    if (fPm->ParameterExists(GetFullParmName("FiberHalfLength"))){
      ChromatinLength = fPm->GetDoubleParameter(GetFullParmName("FiberHalfLength"),"Length");
    }
    if (fPm->ParameterExists(GetFullParmName("FiberRadius"))){
      ChromatinRadius = fPm->GetDoubleParameter(GetFullParmName("FiberRadius"),"Length");
    }


    //****************************************************************************
    //                   Cylinder containing Histone & DNA (envelope)           //
    //****************************************************************************

    G4Tubs* gCylinder = new G4Tubs("Chromatin",
                                   0,
                                   ChromatinRadius,
                                   ChromatinLength,
                                   180.0 * deg,
                                   360.0 * deg);

    fEnvelopeLog = CreateLogicalVolume(gCylinder);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    //vis
    G4VisAttributes* EnvelopeVis = new G4VisAttributes(G4Colour(0.,0.,1.));
    EnvelopeVis->SetVisibility(true);
    EnvelopeVis->SetForceWireframe(true);
    EnvelopeVis->SetForceAuxEdgeVisible(true);
    fEnvelopeLog->SetVisAttributes(EnvelopeVis);


    //****************************************************************************
    //                              Histones                                    //
    //****************************************************************************
    //Record the position and rotation of each Histone placed
    vector<pair<G4ThreeVector, G4RotationMatrix*>> HistoneDetails;
    BuildHistones(HistoneDetails, ChromatinRadius, ChromatinLength);


    //****************************************************************************
    //                                DNA                                       //
    //****************************************************************************
    //Builds DNA around and between the Histones
    BuildDNA(HistoneDetails);


    InstantiateChildren(fEnvelopePhys);


    return fEnvelopePhys;
}

//build the Histones
void TsSolenoidFiber::BuildHistones(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                            G4double ChromatinRadius,
                            G4double ChromatinLength)
{
    //****************************************************************************
    //                              Histones
    //****************************************************************************

    //place the histones in a left handed solenoid conformation

    G4int nb_HistPerTurn = 6;

    if (fPm->ParameterExists(GetFullParmName("nb_HistPerTurn"))){
        nb_HistPerTurn = fPm->GetIntegerParameter(GetFullParmName("nb_HistPerTurn"));
    }


    G4double HistoneRadius = 3.3*nm;
    G4double HistoneLength = (5.7 / 2.0)*nm;

    G4String SubComp = "Histone";

    G4Tubs* gHistone = new G4Tubs("Histone",
                                  0,
                                  HistoneRadius,
                                  HistoneLength,
                                  180*deg,
                                  360*deg);

    G4String Water="G4_WATER";

    // G4LogicalVolume* lHistone = CreateLogicalVolume(SubComp, Water, gHistone);
    G4LogicalVolume* lHistone = CreateLogicalVolume(SubComp, gHistone);

    HistoneRadius+=3.0*nm;  //extra spacing (3nm) for the double helix


    //Generate path for nucleosomes -- solenoid
    G4double x = ((ChromatinRadius-(HistoneRadius+(1.1*nm))) * sin(0.0));
    G4double y = ((ChromatinRadius-(HistoneRadius+(1.1*nm))) * cos(0.0));
    G4double z = (-ChromatinLength + HistoneRadius);

    G4int nNuc = (6.0*2.0*(-z)) / (2.0*HistoneRadius); //how many histones fit along chromatin z
    G4double zStep = (2.0*HistoneRadius)/nb_HistPerTurn;
    G4ThreeVector position (x,y,z);

    G4double theta = 0.0;
    G4double thetaStep = (2.0*pi)/(G4double)nb_HistPerTurn;

    //nNuc=1;

    G4int built=0;
    for (G4int i=0;i<nNuc;i++){
        theta+=thetaStep;
        position[0]=((ChromatinRadius-HistoneRadius) * sin(theta));
        position[1]=((ChromatinRadius-HistoneRadius) * cos(theta));
        position[2]=z;
        z+=zStep;


        if (position[2]+HistoneRadius>=ChromatinLength){
            nNuc=built;
            break;
        }

        G4RotationMatrix *HistoneRotation = new G4RotationMatrix();
        HistoneRotation->rotateZ((-120.0+((G4double)(i)*360.0/(G4double)nb_HistPerTurn))*deg);
        HistoneRotation->rotateY(90.0*deg);

        //record the position and rotation of the histone
        pair<G4ThreeVector,G4RotationMatrix*> Details;
        Details.first=position;
        Details.second=HistoneRotation;
        HistoneDetails.push_back(Details);

        //        G4Transform3D transformHistone = G4Transform3D(HistoneRotation,position);

        G4ThreeVector *pos = &position;
        CreatePhysicalVolume("Histone", i, true, lHistone, HistoneRotation, pos, fEnvelopePhys);

      //  new G4PVPlacement(HistoneRotation,
      //                    position,
      //                    lHistone,
      //                    "Histone",
      //                    fEnvelopeLog,
      //                    false,
      //                    0,
      //                    false);

        built++;
    }

    G4VisAttributes* HistoneVis = new G4VisAttributes(G4Colour(0.,0.,1.));
    HistoneVis->SetVisibility(true);
    HistoneVis->SetForceSolid(true);
    lHistone->SetVisAttributes(HistoneVis);
}

//build the DNA
void TsSolenoidFiber::BuildDNA(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails)
{
    //****************************************************************************
    //                                DNA
    //****************************************************************************

    G4bool BuildHalfCyl=false;
    G4bool BuildQuartCyl=false;
    G4bool BuildSphere=false;

    G4String DNAModel = "QuartCyl";
    if (fPm->ParameterExists(GetFullParmName("DNA_Model"))){
        DNAModel = fPm->GetStringParameter(GetFullParmName("DNA_Model"));
    }

    if (DNAModel=="HalfCyl"){BuildHalfCyl=true;}
    else if (DNAModel=="Sphere"){BuildSphere=true;}
    else if (DNAModel=="QuartCyl"){BuildQuartCyl=true;}

    SetDNAVolumes(BuildHalfCyl, BuildQuartCyl, BuildSphere);

    vector<G4ThreeVector> DNAPath;
    GenerateDNAPath(HistoneDetails, DNAPath);

    SegmentDNAPath(DNAPath);

    if (BuildSphere){
        PlaceDNASphere(DNAPath);

    } else {
        PlaceDNA(DNAPath);
    }
}

//set up DNA volumes
void TsSolenoidFiber::SetDNAVolumes(G4bool BuildHalfCyl,
                            G4bool ,
                            G4bool BuildSphere)
{
    // A choice to build the DNA volumes as Half Cylinders, Quarter Cylinders, or Spheres
    // Quarter cylinders by default

    //a denser water for DNA
    G4NistManager* man = G4NistManager::Instance();
    man -> FindOrBuildMaterial("G4_WATER");

    //G4String WatMod = "G4_WATER_MODIFIED";
    // if additional densities are wanted,use:
    //G4Material*H2O_mod=man -> BuildMaterialWithNewDensity(WatMod,
    //                                                      "G4_WATER",
    //                                                      1.407 *g/cm/cm/cm);

    //SubCompartment names
    G4String back1 = "Backbone1";
    G4String back2 = "Backbone2";
    G4String base1 = "Base1";
    G4String base2 = "Base2";

    //sphere DNA
    if (BuildSphere){
        G4Sphere* gDNA_base = new G4Sphere("DNA_base",
                                           0*nm,
                                           0.208*nm,
                                           0*deg,
                                           360*deg,
                                           0*deg,
                                           180*deg);

        lBase1 = CreateLogicalVolume(base1, gDNA_base);
        lBase2 = CreateLogicalVolume(base2, gDNA_base);


        G4Sphere* gDNA_backbone = new G4Sphere("DNA_backbone",
                                               0*nm,
                                               0.24*nm,
                                               0*deg,
                                               360*deg,
                                               0*deg,
                                               180*deg);

        lBack1 = CreateLogicalVolume(back1, gDNA_backbone);
        lBack2 = CreateLogicalVolume(back2, gDNA_backbone);
    }
    //half cylinder
    else if (BuildHalfCyl){
        G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        180.*deg,
                                        180.*deg);
        G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        0.*deg,
                                        180.*deg);

        lBase1 = CreateLogicalVolume(base1, gDNA_base1);
        lBase2 = CreateLogicalVolume(base2, gDNA_base2);


        G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_backbone1",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            180.*deg,
                                            180.*deg);
        G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_backbone2",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            0*deg,
                                            180.*deg);

        lBack1 = CreateLogicalVolume(back1, gDNA_backbone1);
        lBack2 = CreateLogicalVolume(back2, gDNA_backbone2);
    }
    //quarter cylinder
    else {
        G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        135.*deg,
                                        180.*deg);
        G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
                                        0*nm,
                                        0.5*nm,
                                        ((0.34-0.01)/2.0)*nm,
                                        -45.*deg,
                                        180.*deg);

        lBase1 = CreateLogicalVolume(base1, gDNA_base1);
        lBase2 = CreateLogicalVolume(base2, gDNA_base2);

        G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_backbone1",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            180.*deg,
                                            90.*deg);
        G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_backbone2",
                                            0.5*nm,
                                            1.15*nm,
                                            ((0.34-0.01)/2.0)*nm,
                                            0*deg,
                                            90*deg);

        lBack1 = CreateLogicalVolume(back1, gDNA_backbone1);
        lBack2 = CreateLogicalVolume(back2, gDNA_backbone2);
    }

    //        G4cout<<"MATERIAL NAME "<<lBack1->GetMaterial()->GetName()<<G4endl;
    //        G4cout<<"MATERIAL DENS "<<lBack1->GetMaterial()->GetDensity()/(g/cm/cm/cm)<<" g/cm3"<<G4endl;

    //vis
    G4VisAttributes* Back1Vis = new G4VisAttributes(G4Colour(1.0, 0.1, 0.02));
    Back1Vis->SetVisibility(true);
    Back1Vis->SetForceSolid(true);
    lBack1->SetVisAttributes(Back1Vis);

    G4VisAttributes * Back2Vis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.2));
    Back2Vis->SetVisibility(true);
    Back2Vis->SetForceSolid(true);
    lBack2 -> SetVisAttributes(Back2Vis);

    G4VisAttributes * Base1Vis = new G4VisAttributes(G4Colour(0.02, 0.92, 1.0));
    Base1Vis->SetVisibility(true);
    Base1Vis->SetForceSolid(true);
    lBase1 -> SetVisAttributes(Base1Vis);

    G4VisAttributes * Base2Vis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.8));
    Base2Vis->SetVisibility(true);
    Base2Vis->SetForceSolid(true);
    lBase2 -> SetVisAttributes(Base2Vis);

}

//generate a path for DNA around the histones
void TsSolenoidFiber::GenerateDNAPath(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                              vector<G4ThreeVector> &path)
{
    //generate a path of 3vector points that defines the centre of the dbl helix

    G4double TraceBP=15.0; //the number of BP to trace back for bezier curve (more = bendier or smoother)
    G4double nCoils = 1.65; //#coils around histone
    G4int nHistones=HistoneDetails.size();
    G4double offset = (pi/2.0) - (2.0*pi*nCoils);
    G4double zHist = (5.7 / 2.0)*nm;
    G4double Radius = 3.3*nm + 1.15*nm + 0.1*nm; //central radius of double helix around histone

    G4int nSteps=200; //go around the histone in 200 steps (segmented later)

    for (G4int i=0;i<nHistones;i++){
        G4RotationMatrix *rot=HistoneDetails[i].second;

        //histone DNA
        for (G4int j=0;j<nSteps;j++){
            G4double angle = ((G4double)j * 2.0 * pi * nCoils / (G4double)nSteps)+(pi/2.0)+offset;
            //take a pt around the histone
            G4ThreeVector pt( (Radius * sin(angle)),
                             (Radius * cos(angle)),
                             (-zHist + ((G4double)j/(G4double)nSteps)*2.0*(zHist)));

            //rotate the step to the histone
            G4ThreeVector RotPt;
            RotPt[0] = pt.x()*rot->xx() + pt.y()*rot->yx() + pt.z()*rot->zx();
            RotPt[1] = pt.x()*rot->xy() + pt.y()*rot->yy() + pt.z()*rot->zy();
            RotPt[2] = pt.x()*rot->xz() + pt.y()*rot->yz() + pt.z()*rot->zz();

            //translate to the histone position
            RotPt += HistoneDetails[i].first;

            path.push_back(RotPt);
        }


        //link to nextPt -- off histone (start from end of histone i -> start of histone i+1)
        if (i!=nHistones-1){
            G4ThreeVector start = path[path.size()-1]; //last pt on histone

            //calc 1st point on next histone
            G4double angle = (pi/2.0)+offset;
            G4ThreeVector pt( (Radius * sin(angle)),
                             (Radius * cos(angle)),
                             -zHist );
            G4RotationMatrix *nextRot = HistoneDetails[i+1].second;
            G4ThreeVector RotPt;
            RotPt[0] = pt.x()*nextRot->xx() + pt.y()*nextRot->yx() + pt.z()*nextRot->zx();
            RotPt[1] = pt.x()*nextRot->xy() + pt.y()*nextRot->yy() + pt.z()*nextRot->zy();
            RotPt[2] = pt.x()*nextRot->xz() + pt.y()*nextRot->yz() + pt.z()*nextRot->zz();
            G4ThreeVector end = RotPt + HistoneDetails[i+1].first; //1st pt on next histone

            //calc 2nd point on next histone
            angle= (1.0 * 2.0 * pi * nCoils / (G4double)nSteps)+(pi/2.0)+offset;
            G4ThreeVector secpt( (Radius * sin(angle)),
                                (Radius * cos(angle)),
                                (-zHist + (1.0/(G4double)nSteps)*2.0*(zHist)));
            RotPt[0] = secpt.x()*nextRot->xx() + secpt.y()*nextRot->yx() + secpt.z()*nextRot->zx();
            RotPt[1] = secpt.x()*nextRot->xy() + secpt.y()*nextRot->yy() + secpt.z()*nextRot->zy();
            RotPt[2] = secpt.x()*nextRot->xz() + secpt.y()*nextRot->yz() + secpt.z()*nextRot->zz();
            secpt = RotPt + HistoneDetails[i+1].first; //2nd pt on next histone

            G4ThreeVector StartDir = (start - path[path.size()-2]).unit(); //direction coming off histone
            G4ThreeVector MidPoint1 = ((TraceBP*0.3)*StartDir); //trace it forward x bp
            MidPoint1[0] *= nm;
            MidPoint1[1] *= nm;
            MidPoint1[2] *= nm;
            MidPoint1 += start;

            G4ThreeVector EndDir = (secpt-end).unit();
            G4ThreeVector MidPoint2 = ((TraceBP*0.3)*EndDir); //trace it forward x bp
            MidPoint2[0] *= nm;
            MidPoint2[1] *= nm;
            MidPoint2[2] *= nm;
            MidPoint2 = end - MidPoint2;

            //link the 4 points with Bezier in 100 steps
            Bezier(start, MidPoint1, MidPoint2, end, path, 100);
        }
    }

}

//Bezier curve to smoothly join 2 points (bending through 2 other points)
void TsSolenoidFiber::Bezier(G4ThreeVector &start,
                     G4ThreeVector &MidPoint1,
                     G4ThreeVector &MidPoint2,
                     G4ThreeVector &end,
                     vector<G4ThreeVector> &path,
                     G4int nSteps)
{
    G4double j=0.0;
    for (G4int k=0;k<nSteps-1;k++){
        j+=1.0/(G4double)nSteps;
        G4ThreeVector BezPt = pow(1.0-j,3.0)*start + 3.0*pow(1.0-j,2.0)*j*MidPoint1 +
        3.0*(1.0-j)*j*j*MidPoint2 + j*j*j*end;
        path.push_back(BezPt);
    }
}

//split the DNA path into 0.34nm steps
void TsSolenoidFiber::SegmentDNAPath(vector<G4ThreeVector> &path)
{
    vector<G4ThreeVector> newPath;
    double rise=0.34*nm;
    int nPts=path.size();
    int counter=0;
    newPath.push_back(path[0]);
    for (int i=0;i<nPts-1;i++){
        G4ThreeVector vector = (path[i+1]-newPath[counter]);
        double length = vector.mag();
        int nDiv = length/rise;
        G4ThreeVector unit = vector.unit();
        unit *= rise;
        for (int j=0;j<nDiv;j++){
            G4ThreeVector nuwe = newPath[counter] + unit;
            newPath.push_back(nuwe);
            counter++;
        }
    }

    newPath=path;
    //path.clear();

}

//use DNA path to place sphere DNA volumes
void TsSolenoidFiber::PlaceDNASphere(vector<G4ThreeVector> &newPath)
{
    G4double helixRadius = 1.0*nm;
    G4double rotPair = ((2.0*pi)/10.0);   //10bp per turn
    G4int nBP=newPath.size();
    G4double rBack=helixRadius - 0.24*nm;
    G4double rBase=rBack - 0.24*nm - 0.208*nm;

    for (int bp=0; bp<nBP-1; bp++){
        G4double angle1 = (G4double)bp * rotPair;
        G4double angle2 = angle1+pi;// + (120.0*pi/180.0); //offset for strand2 (major and minor groove)

        //temporary positions
        G4ThreeVector back1temp = G4ThreeVector((rBack*cos(angle1)), (rBack*sin(angle1)), 0.0);
        G4ThreeVector back2temp = G4ThreeVector((rBack*cos(angle2)), (rBack*sin(angle2)), 0.0);
        G4ThreeVector base1temp = G4ThreeVector((rBase*cos(angle1)), (rBase*sin(angle1)), 0.0);
        G4ThreeVector base2temp = G4ThreeVector((rBase*cos(angle2)), (rBase*sin(angle2)), 0.0);

        //Rotation to point to next plane
        G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
        G4ThreeVector norm (0.,0.,-1.); //the normal to the plane (G4 build planes facing -z)
        G4double DotProd = norm.dot(vecNext);
        G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
        G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm

        //set up new 3Vectors for rotated pos
        G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.);

        //Apply rotation
        G4RotationMatrix *rot = new G4RotationMatrix;
        rot->rotate(AngBetween, cross);

        ApplyRotation(back1, back1temp, rot);
        ApplyRotation(back2, back2temp, rot);
        ApplyRotation(base1, base1temp, rot);
        ApplyRotation(base2, base2temp, rot);
        delete rot;

        //Translate
        back1+=newPath[bp];
        back2+=newPath[bp];
        base1+=newPath[bp];
        base2+=newPath[bp];

        int bpID=bp+1;

        G4ThreeVector *posBack1 = &back1;
        CreatePhysicalVolume("Backbone1", bpID, true, lBack1, 0, posBack1, fEnvelopePhys);
        G4ThreeVector *posBase1 = &base1;
        CreatePhysicalVolume("Base1", bpID, true, lBase1, 0, posBase1, fEnvelopePhys);
        G4ThreeVector *posBack2 = &back2;
        CreatePhysicalVolume("Backbone2", bpID, true, lBack2, 0, posBack2, fEnvelopePhys);
        G4ThreeVector *posBase2 = &base2;
        CreatePhysicalVolume("Base2", bpID, true, lBase2, 0, posBase2, fEnvelopePhys);

      //  new G4PVPlacement(0,
      //                    back1,
      //                    lBack1,
      //                    "back1",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
      //  new G4PVPlacement(0,
      //                    base1,
      //                    lBase1,
      //                    "base1",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
       //
      //  new G4PVPlacement(0,
      //                    back2,
      //                    lBack2,
      //                    "back2",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
      //  new G4PVPlacement(0,
      //                    base2,
      //                    lBase2,
      //                    "base2",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);

    }

    G4cout<<G4endl<<"DNA BUILT: "<<(G4double)nBP/1000.0<<" kbp"<<G4endl<<G4endl;
}

//rotate a 3vector by a rot matrix and return new coordinates
void TsSolenoidFiber::ApplyRotation(G4ThreeVector &rotated, G4ThreeVector &vector, G4RotationMatrix *rot)
{
    rotated[0] = vector[0]*rot->xx() + vector[1]*rot->yx() + vector[2]*rot->zx();
    rotated[1] = vector[0]*rot->xy() + vector[1]*rot->yy() + vector[2]*rot->zy();
    rotated[2] = vector[0]*rot->xz() + vector[1]*rot->yz() + vector[2]*rot->zz();

}

//Use DNA path to place halfCyl or quartCyl DNA volumes
void TsSolenoidFiber::PlaceDNA(vector<G4ThreeVector> &newPath)
{
    //G4double helixRadius = 1.0*nm;
    G4double rotPair = ((2.0*pi)/10.0);   //10bp per turn
    G4int nBP=newPath.size();

    for (G4int bp=0; bp<nBP-1; bp++){

        //Position of base + back in xy
        //Definitely gives right handed coil (checked) -- left handed in -ve z?
        G4double angle1 = -(G4double)bp * rotPair;

        //Rotation to point to next plane
        G4ThreeVector vecNext = (newPath[bp]-newPath[bp+1]).unit(); //unit vec pointing to next
        G4ThreeVector norm (0.,0.,1.); //the normal to the plane (G4 build planes facing -z)
        G4double DotProd = norm.dot(vecNext);
        G4double AngBetween = acos(DotProd); //angle between this plane and next (rad)
        G4ThreeVector cross = (vecNext.cross(norm)).unit(); //vector perp to vecnext and norm

        //set up new 3Vectors for rotated pos
        G4ThreeVector back1(0.,0.,0.), back2(0.,0.,0.), base1(0.,0.,0.), base2(0.,0.,0.);

        //Apply rotation
        G4RotationMatrix *rot1 = new G4RotationMatrix();
        rot1->rotate(AngBetween, cross);
        rot1->rotateZ(angle1);

        //Translate
        back1+=newPath[bp];
        back2+=newPath[bp];
        base1+=newPath[bp];
        base2+=newPath[bp];

        int bpID=bp+1;

        //G4cout<<bp<<" "<<back1<<" "<<newPath[bp]<<G4endl;
        //        G4Transform3D transformStrand1 = G4Transform3D(rot1,newPath[bp]);
        //        G4Transform3D transformStrand2 = G4Transform3D(rot2,newPath[bp]);



        G4ThreeVector *posBack1 = &back1;
        CreatePhysicalVolume("Backbone1", bpID, true, lBack1, rot1, posBack1, fEnvelopePhys);
        G4ThreeVector *posBase1 = &base1;
        CreatePhysicalVolume("Base1", bpID, true, lBase1, rot1, posBase1, fEnvelopePhys);
        G4ThreeVector *posBack2 = &back2;
        CreatePhysicalVolume("Backbone2", bpID, true, lBack2, rot1, posBack2, fEnvelopePhys);
        G4ThreeVector *posBase2 = &base2;
        CreatePhysicalVolume("Base2", bpID, true, lBase2, rot1, posBase2, fEnvelopePhys);


      //  new G4PVPlacement(rot1,
      //                    newPath[bp],
      //                    lBack1,
      //                    "back1",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
      //  new G4PVPlacement(rot1,
      //                    newPath[bp],
      //                    lBase1,
      //                    "base1",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
       //
      //  new G4PVPlacement(rot1,
      //                    newPath[bp],
      //                    lBack2,
      //                    "back2",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
      //  new G4PVPlacement(rot1,
      //                    newPath[bp],
      //                    lBase2,
      //                    "base2",
      //                    fEnvelopeLog,
      //                    false,
      //                    bpID,
      //                    false);
    }

    G4cout<<G4endl<<"DNA BUILT: "<<(G4double)nBP/1000.0<<" kbp"<<G4endl<<G4endl;
}
