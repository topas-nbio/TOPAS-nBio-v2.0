// Component for TsPlasmidSupercoiled
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
 *  Developed by Nicholas Henthorn, PRECISE, University of Manchester
 *  nicholas.henthorn@manchester.ac.uk
 *  https://tinyurl.com/y7xuxw9d
 *  17/10/2018, LINK TO PUBLICATION
 */

#include "TsPlasmidSupercoiled.hh"

#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Orb.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VisAttributes.hh"

TsPlasmidSupercoiled::TsPlasmidSupercoiled(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsPlasmidSupercoiled::~TsPlasmidSupercoiled()
{
}



G4VPhysicalVolume* TsPlasmidSupercoiled::Construct()
{
    BeginConstruction();

	
    //****************************************************************************
    //                   Read Plasmid CoOrdinates                               //
    //****************************************************************************
	G4String fileName = fPm->GetStringParameter(GetFullParmName("FileName"));
	
    //build the dna
	std::vector<G4ThreeVector> path; //DNA path
	std::vector<DNA*> DNAPts; //ptrs to DNA
    
    // Read in plasmid vertices from file
	std::ifstream infile(fileName);
    if (!infile){
		G4cout<<"Can't open input file: " << fileName << G4endl;
		fPm->AbortSession(true);
    }
    
	G4double x, y, z;
	G4double xmin = 1.0*mm, xmax = 0.0*mm;
	G4double ymin = 1.0*mm, ymax = 0.0*mm;
	G4double zmin = 1.0*mm, zmax = 0.0*mm;

    while (1){
		infile >> x >> y >> z;
        if ( !infile.good() )
			break;
		x *= nm;
		y *= nm;
		z *= nm;
        G4ThreeVector PathPoint = G4ThreeVector(x,y,z);
        path.push_back(PathPoint);
		if ( xmin > x )
			xmin = x;
		if ( ymin > y )
			ymin = y;
		if ( zmin > z )
			zmin = z;
		if ( xmax < x )
			xmax = x;
		if ( ymax < y )
			ymax = y;
		if ( zmax < z )
			zmax = z;
    }
	//add the 1st point again to close the plasmid
	path.push_back(path[0]);
	

	G4double HLX = 0.5*(xmax-xmin);
	G4double HLY = 0.5*(ymax-ymin);
	G4double HLZ = 0.5*(zmax-zmin);
	G4ThreeVector offset = G4ThreeVector(0.5*(xmax+xmin),
										 0.5*(ymax+ymin),
										 0.5*(zmax+zmin));

	for ( size_t u = 0; u < path.size(); u++ )
		path[u] -= offset;
	
	G4Box* envelope = new G4Box(fName, HLX+0.5*3.4*nm, HLY+0.5*3.4*nm, HLZ+0.5*3.4*nm);
	fEnvelopeLog = CreateLogicalVolume(envelope);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	

	G4bool segment = false;
	if (fPm->ParameterExists(GetFullParmName("SegmentPlasmidPath")))
		segment = fPm->GetBooleanParameter(GetFullParmName("SegmentPlasmidPath"));
	
	vector<G4ThreeVector> SmoothedPath;

	if ( segment ) {
		//****************************************************************************
		//                   Smooth the path                                        //
		//****************************************************************************
		
		for (size_t i=1;i<path.size()-1;i++){
			G4ThreeVector OffDirect(0,0,0), OnDirect(0,0,0);
			OffDirect=(path[i]-path[i-1]).unit();
			OnDirect=(path[i]-path[i+1]).unit();
			
			//trace 2nm
			G4ThreeVector MidPoint1(0,0,0), MidPoint2(0,0,0);
			MidPoint1 = path[i-1]+(2.0*OffDirect*nm);
			MidPoint2 = path[i]+(2.0*OnDirect*nm);
			Bezier(path[i-1],
				   MidPoint1,
				   MidPoint2,
				   path[i],
				   SmoothedPath,
				   500);
		}
		
		//Join last pt to first pt
		G4ThreeVector OffDirect(0,0,0), OnDirect(0,0,0);
		OffDirect=(SmoothedPath[SmoothedPath.size()-1]-SmoothedPath[0]).unit();
		OnDirect=(SmoothedPath[0]-path[1]).unit();
		
		//trace 2nm
		G4ThreeVector MidPoint1(0,0,0), MidPoint2(0,0,0);
		MidPoint1 = SmoothedPath[SmoothedPath.size()-1]+(2.0*OffDirect*nm);
		MidPoint2 = SmoothedPath[0]+(2.0*OnDirect*nm);
		Bezier(SmoothedPath[SmoothedPath.size()-1],
			   MidPoint1,
			   MidPoint2,
			   SmoothedPath[0],
			   SmoothedPath,
			   500);
 
	} else {
		SmoothedPath = path;
	}
	
    G4bool BuildHalfCyl=false;
    G4bool BuildQuartCyl=false;
    G4bool BuildSphere=false;
    
    G4String DNAModel = "QuartCyl";
    if (fPm->ParameterExists(GetFullParmName("DNA_Model"))){
        DNAModel = fPm->GetStringParameter(GetFullParmName("DNA_Model"));
    }
    
	if (DNAModel=="HalfCyl") {
		BuildHalfCyl=true;
	} else if (DNAModel=="Sphere") {
		BuildSphere=true;
	} else if (DNAModel=="QuartCyl") {
		BuildQuartCyl=true;
	} else {
		G4cerr << "Error. DNA model " << DNAModel << " not found. " << G4endl;
		G4cerr << "Use: HalfCyl, QuartCyl or Sphere" << G4endl;
	}
	
    DNACoordinates * co = new DNACoordinates;
    co->Generate(SmoothedPath,
                 DNAPts,
                 BuildHalfCyl,
                 BuildQuartCyl,
                 BuildSphere,
				 segment);
    delete co;
    
    
    //****************************************************************************
    //                   Build DNA                                              //
    //****************************************************************************
    PlaceDNA(DNAPts,
             BuildHalfCyl,
             BuildQuartCyl,
			 BuildSphere);
    
    InstantiateChildren(fEnvelopePhys);
    return fEnvelopePhys;
}


//****************************************************************************
//                   Place DNA                                              //
//****************************************************************************
void TsPlasmidSupercoiled::PlaceDNA(vector<DNA *> &DNApt,
                         G4bool BuildHalfCyl,
                         G4bool BuildQuartCyl,
                         G4bool BuildSphere)
{
    G4LogicalVolume* lBase1 = 0;
    G4LogicalVolume* lBase2 = 0;
    G4LogicalVolume* lBack1 = 0;
    G4LogicalVolume* lBack2 = 0;
	
	//Spherical DNA
    if (BuildSphere){
		G4Orb* gDNA_base = new G4Orb("DNA_base",  0.208*nm);
        
        lBase1 = CreateLogicalVolume("Base1", gDNA_base);
        lBase2 = CreateLogicalVolume("Base2", gDNA_base);
        
		G4Orb* gDNA_backbone = new G4Orb("DNA_deoxyribose", 0.29*nm); 
        
        lBack1 = CreateLogicalVolume("Backbone1", gDNA_backbone);
        lBack2 = CreateLogicalVolume("Backbone2", gDNA_backbone);
    }
    //Cylinder DNA
    else {
        G4double SweptAngle=0.0;
        G4double StartAngle1=0.0;
        G4double StartAngle2=0.0;
        if (BuildQuartCyl){
            SweptAngle=90.0;
            StartAngle1=135.0;
            StartAngle2=-45.0;
        }
        else if (BuildHalfCyl){
            SweptAngle=180.0;
            StartAngle1=180.0;
            StartAngle2=0.0;
        }
        
        G4Tubs* gDNA_base1 = new G4Tubs("DNA_base1",
                                       0*nm, //pRmin
                                       0.5*nm,
                                       ((0.34-0.01)/2.0)*nm,
                                       StartAngle1*deg, //pSphi //180
                                       180*deg); //pDphi
        G4Tubs* gDNA_base2 = new G4Tubs("DNA_base2",
                                       0*nm, //pRmin
                                       0.5*nm,
                                       ((0.34-0.01)/2.0)*nm,
                                       StartAngle2*deg, //pSphi //0
                                       180*deg); //pDphi
        
        
        lBase1 = CreateLogicalVolume("Base1", gDNA_base1);
        lBase2 = CreateLogicalVolume("Base2", gDNA_base2);
        
        
        G4Tubs* gDNA_backbone1 = new G4Tubs("DNA_deoxyribose1",
                                           0.5*nm,
                                           1.15*nm,
                                           ((0.34-0.01)/2.0)*nm,
                                           180*deg, //pSphi
                                           SweptAngle*deg); //pDphi
        G4Tubs* gDNA_backbone2 = new G4Tubs("DNA_deoxyribose2",
                                           0.5*nm,
                                           1.15*nm,
                                           ((0.34-0.01)/2.0)*nm,
                                           0*deg, //pSphi
                                           SweptAngle*deg); //pDphi
        
        lBack1 = CreateLogicalVolume("Backbone1", gDNA_backbone1);
        lBack2 = CreateLogicalVolume("Backbone2", gDNA_backbone2);
    }
    
    //place the DNA
    for (size_t bp=0; bp<DNApt.size(); bp++){
        G4int bpID = DNApt[bp]->GetBP();
        G4int Stra = DNApt[bp]->GetStrand();
        G4bool IsBac = DNApt[bp]->GetIsBack();
        G4bool IsBas = DNApt[bp]->GetIsBase();
        
		G4RotationMatrix* rotation = DNApt[bp]->GetRot();
		G4ThreeVector* position = new G4ThreeVector(DNApt[bp]->GetPos());
		
        if (IsBac){
            if (Stra==1){
                CreatePhysicalVolume("Backbone1", bpID, true, lBack1, rotation, position, fEnvelopePhys);
            }
            else if (Stra==2){
                CreatePhysicalVolume("Backbone2", bpID, true, lBack2, rotation, position, fEnvelopePhys);
            }
        } else if (IsBas) {
            if (Stra==1){
                CreatePhysicalVolume("Base1", bpID, true, lBase1, rotation, position, fEnvelopePhys);
            }
            else if (Stra==2){
                CreatePhysicalVolume("Base2", bpID, true, lBase2, rotation, position, fEnvelopePhys);
            }
        }
        
    }
    G4cout<<G4endl<<"base pairs: "<<(0.25*DNApt.size()/1000.0)<<" kbp"<<G4endl;
}


// -------------------------------------------------------------------------
//          Bezier Curve
// -------------------------------------------------------------------------
void TsPlasmidSupercoiled::Bezier(G4ThreeVector &start,
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
