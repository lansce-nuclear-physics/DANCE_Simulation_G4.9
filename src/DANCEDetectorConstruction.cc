
////////////////////////////////////////////////////////////////////////
//                                                                    //
//   Software Name: DANCE Data Acquisition and Analysis Package       //
//     Subpackage: DANCE_Simulation_G4.9                              //
//   Identifying Number: C18105                                       // 
//                                                                    //
////////////////////////////////////////////////////////////////////////
//                                                                    //
//                                                                    //
// Copyright 2019.                                                    //
// Triad National Security, LLC. All rights reserved.                 //
//                                                                    //
//                                                                    //
//                                                                    //
// This program was produced under U.S. Government contract           //
// 89233218CNA000001 for Los Alamos National Laboratory               //
// (LANL), which is operated by Triad National Security, LLC          //
// for the U.S. Department of Energy/National Nuclear Security        //
// Administration. All rights in the program are reserved by          //
// Triad National Security, LLC, and the U.S. Department of           //
// Energy/National Nuclear Security Administration. The Government    //
// is granted for itself and others acting on its behalf a            //
// nonexclusive, paid-up, irrevocable worldwide license in this       //
// material to reproduce, prepare derivative works, distribute        //
// copies to the public, perform publicly and display publicly,       //
// and to permit others to do so.                                     //
//                                                                    //
// This is open source software; you can redistribute it and/or       //
// modify it under the terms of the GPLv2 License. If software        //
// is modified to produce derivative works, such modified             //
// software should be clearly marked, so as not to confuse it         //
// with the version available from LANL. Full text of the GPLv2       //
// License can be found in the License file of the repository         //
// (GPLv2.0_License.txt).                                             //
//                                                                    //
////////////////////////////////////////////////////////////////////////



// DANCEDetector Construction   -  M. Jandel - 2006 - May 9
//  
// Notes: 
//	The DANCE full array is constructed - 160 crystals of A,B,C,D type,
//	wrapping material, Aluminium holders of all crystals, the supporting Al ball,
//	photomultiplier lense (Quarz), beampipe, LiH shell	
//
//	to visualize the DANCE one should use OPENGL or DAWN and torn on the
//	corresponding GEANT4 environment variables
//
//	The mapping of the DANCE crystal might seem to be puzzling but actually it's not
//	How it is done:
//		If GEANT4 registers a hit it is the Logical Volume or SensitiveDetector that its being referenced to
//		That is why one has to map the logical volumes or sensitive detectors 
//		to the real DANCE crystal_ID labels.
//		In this version half of the DANCE is build using an input file angles.txt
//		and half (80) detectors are constructed using the angles given in the file.
//		For the other half the 180 reflection symmetry is used. The crystals are being added in pairs
//		first one positioned according to the angles.txt file and the second one is reflected around 
//		the symmetry axis.
//		Labeling goes as followin: Arrays of ID2D[CrysType][CrysNum]=RealCrystalID
//			(This is done while building the array crystal by crystal,
//				DetType: A=0,B=1,C=2,D=3)
//		When the Logical Volumes are registered it is important to keep the track of the naming conventions
//		and of course mapping the crystals.
//		This is achieved by creating the second map GeantID2D which maps logical volumes and CrysType,CrysNum.
//		Any hit in the DANCEEventAction coming from Logical volume i is then easily converted
//		using the maps ID2D and GeantID2D matrices. More is obvious from implementation of these maps
//		in this class and DANCEEventAction class. Just to mention,  the maps ID2D and GeantID2D,
//		are defined in the global class DANCEGlobal.hh
//
//		Naming convention for Sensitive detectors is : BaF_A_0  (BaF_"DetType"_DetNo)

#include "DANCEDetectorConstruction.hh"

#include "G4RunManager.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Vector3D.hh"

#include "G4ReflectionFactory.hh"
#include "G4NistManager.hh"

#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4VSDFilter.hh"
#include "G4SDParticleFilter.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSMinKinEAtGeneration.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4SDManager.hh"

#include "G4Point3D.hh"

#include "DANCEPrimaryGeneratorAction.hh"

#include "C_crystal.hh"
#include "B_crystal.hh"

#include "DANCEDetector.hh"
#include "DANCEGlobal.hh"

using namespace std;

extern int A_counter,B_counter,C_counter,D_counter;


DANCEDetectorConstruction::DANCEDetectorConstruction()
 : experimentalHall_log(0), experimentalHall_phys(0)
 {
 
 // define names for the detectors
 G4String Num;
    
 for(size_t j=0;j<15;j++){
 BaF_A_crystal_phys[j]=0;
 BaF_A_crystal_log[j]=0;

 
   ostringstream oss;
   oss<<j;
   Num=oss.str();
  
 detNameA[j]="BaF_A_"+Num;
 
// G4cout << detNameA[j] << G4endl;
 }

for(size_t j=0;j<40;j++){

 BaF_D_crystal_phys[j]=0;
 BaF_D_crystal_log[j]=0;
 
   ostringstream oss;
   oss<<j;
   Num=oss.str();
  
 detNameD[j]="BaF_D_"+Num;
 
// G4cout << detNameD[j] << G4endl;
}
 
  for(size_t j=0;j<60;j++){
	 BaF_B_crystal_phys[j]=0;
	 BaF_B_crystal_log[j]=0;
 
	   ostringstream oss1;
	   oss1<<j;
	   Num=oss1.str();
  
	 detNameB[j]="BaF_B_"+Num;
 
//	 G4cout << detNameB[j] << G4endl; 
 }
 
   for(size_t j=0;j<60;j++){
 BaF_C_crystal_phys[j]=0;
 BaF_C_crystal_log[j]=0;

   ostringstream oss1;
   oss1<<j;
   Num=oss1.str();
  
 detNameC[j]="BaF_C_"+Num;
 
 // G4cout << detNameC[j] << G4endl;
 
 }
// detName[i]="BaF_A_"+i;
 
}


DANCEDetectorConstruction::~DANCEDetectorConstruction()
{
}


G4VPhysicalVolume* DANCEDetectorConstruction::Construct()
{


  DANCEGlobal* GPP=DANCEGlobal::Instance();

  G4NistManager *man=G4NistManager::Instance();
  man->SetVerbose(0);

//______________________________________________________________________//
//									//
// 			Materials, Elements, Isotopes definitions	//
//									//
//______________________________________________________________________//

  G4Material* BaF = man->FindOrBuildMaterial("G4_BARIUM_FLUORIDE");
  G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");
  G4Material* Teflon = man->FindOrBuildMaterial("G4_TEFLON");
  G4Material* Al = man->FindOrBuildMaterial("G4_Al");
  G4Material* Ti = man->FindOrBuildMaterial("G4_Ti");
  G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");
  G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
  G4Material* Au = man->FindOrBuildMaterial("G4_Au");
  G4Material* Gd = man->FindOrBuildMaterial("G4_Gd");

  // stainless steel
	G4Element* Fe = new G4Element("Iron"  ,"Fe" , 26., 55.85*g/mole);  
	G4Element* C  = new G4Element( "Carbon", "C",   6. , 12.011*g/mole);
	G4Element* Co = new G4Element( "Cobalt", "Co", 27. , 58.9332*g/mole);
  G4Material* ssteel = new G4Material ("Steel", 7.7*g/cm3, 3);
	ssteel->AddElement(C, 0.04);
	ssteel->AddElement(Fe, 0.88);
	ssteel->AddElement(Co, 0.08);

  G4Material* TargetMat=new G4Material("Ysource", 1.*g/cm3, 1);
  	G4Isotope* aIsotope=new G4Isotope("Sr",38,88,88*g/mole); 


// Lanthanum Bromide

  G4Element* La=new G4Element("Lanthanum","La",57.,139.9055*g/mole);
  G4Element* Br=new G4Element("Bromine","Br",35.,79.904*g/mole);
  
  G4Material* LaBr3 = new G4Material("LaBr3",5.08 *g/cm3,2);
  LaBr3->AddElement(La, 0.36688);
  LaBr3->AddElement(Br, 0.63312);

//aIsotope->GetIsotopeTable();	


//aIsotope->GetIsotopeTable();	

	G4Element* aElement=new G4Element("Sr",38,1);
	aElement->AddIsotope(aIsotope,1);

	aIsotope=new G4Isotope("Y",39,88,88*g/mole);
	aElement=new G4Element("Y",39,1);
	aElement->AddIsotope(aIsotope,1);

	aIsotope=new G4Isotope("Cs",55,137,137*g/mole);
	aElement=new G4Element("Cs",55,1);
	aElement->AddIsotope(aIsotope,1);
 
	  G4Element* elO  = new G4Element("Oxygen",   "O", 8., 16.00*g/mole);
	  G4Element* elSi = new G4Element("Silicon", "Si", 14., 28.09*g/mole);

  
  	G4Element* elLi = new G4Element("Lithium","Li" , 3., 6.015*g/mole);
	G4Element* elH = new G4Element("Hydrogen" ,"H" , 1., 1.01*g/mole);

//  G4Material* LiH = new G4Material("LiH",0.85*g/cm3,2);
//	LiH->AddElement(elLi, 0.873);
//	LiH->AddElement(elH, 0.127);

	G4Isotope* Li6 = new G4Isotope("Li6", 3, 6, 6.015*g/mole);
	G4Isotope* Li7 = new G4Isotope("Li7", 3, 7, 7.016*g/mole);

  	G4Element* elLi_compound = new G4Element("Lithium_compound","Li" , 2);
	elLi_compound->AddIsotope(Li6,96.*perCent);
	elLi_compound->AddIsotope(Li7,4.*perCent);
	
//     	G4Material* LiH = new G4Material("LiH",0.7288*g/cm3,2);
     	G4Material* LiH = new G4Material("LiH",GPP->LiHDensity*g/cm3,2);
	LiH->AddElement(elLi_compound, 0.8572);
	LiH->AddElement(elH, 0.1427);

// Photomultiplier window is a spherical section made of quartz

// Quartz
// -------
//  density = 2.200*g/cm3; // fused quartz 
  G4double   density = 2.64*g/cm3;  // crystalline quartz (c.f. PDG) 
  G4Material *Quartz = new G4Material("Quartz",density, 2);
  Quartz->AddElement(elSi, 1) ;
  Quartz->AddElement(elO , 2) ;



// PVC
//
//      DATA APVC/12.011,1.008,35.453/
//      DATA ZPVC/6.,1.,17./
//      DATA WPVC/2.,3.,1./
//      DATA RHOPVC/1.406/

  G4Element* elCl = new G4Element("Chlorine","Cl" , 17., 35.453*g/mole);
      
  G4Material* PVC = new G4Material ("PVC", 1.406*g/cm3, 3);
	PVC->AddElement(C, 2);
	PVC->AddElement(elH, 3);
	PVC->AddElement(elCl, 1);


// print out tables

//	G4cout << *(G4Isotope::GetIsotopeTable());
//	G4cout << *(G4Element::GetElementTable());
//	G4cout << *(G4Material::GetMaterialTable());

 /*
  G4Element* elAm241 = new G4Element("elAm241", "Am", 95.,241.01*g/mole);
//  elAm241->AddIsotope(isoAm241,1.);
  
  G4double density=241.01*g/mole;
  
  G4Material *Am241=new G4Material("Am241",density,1);
  Am241->AddElement(elAm241,1.);
   */



//______________________________________________________________________//
//									//
// Define Visibility attributes here - good for later use
//									//
//______________________________________________________________________//

// beampipe
  G4VisAttributes* BeamPipeVA= new G4VisAttributes(G4Colour(0.1,0.9,0.9));
  BeamPipeVA->SetVisibility(true);

// LiH shell
  G4VisAttributes* LiHVA= new G4VisAttributes(G4Colour(0.2,0.2,0.3));
  LiHVA->SetVisibility(true);

// target
  G4VisAttributes* TargetVA= new G4VisAttributes(G4Colour(0.1,0.9,0.9));
  TargetVA->SetVisibility(true);

// crystals
  G4VisAttributes* Crystal_A_VA= new G4VisAttributes(G4Colour(0.5,1.0,1.0));
  Crystal_A_VA->SetVisibility(true);

  G4VisAttributes* Crystal_B_VA= new G4VisAttributes(G4Colour(1.0,0.1,0.1));
  Crystal_B_VA->SetVisibility(true);

  G4VisAttributes* Crystal_C_VA= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  Crystal_C_VA->SetVisibility(true);

  G4VisAttributes* Crystal_D_VA= new G4VisAttributes(G4Colour(0.1,1.0,0.1));
  Crystal_D_VA->SetVisibility(true);

// Aluminium
  G4VisAttributes* sAl_VA= new G4VisAttributes(G4Colour(0.9,1.0,1.0));
  sAl_VA->SetVisibility(true);

// you can hide things, but they will be there
  G4VisAttributes* hideVA= new G4VisAttributes(G4Colour(1.0,0.1,0.1));
  hideVA->SetVisibility(false);

  G4VisAttributes* sVA_Wrap= new G4VisAttributes(G4Colour(0.0,0.0,0.1));
  sVA_Wrap->SetVisibility(true);
//  sVA_Wrap->SetVisibility(false);

//______________________________________________________________________//
//									//
// 			Volumes definition				//
//									//
//______________________________________________________________________//

//------------------------------ experimental hall (world volume)- Air
//------------------------------ beam line along x axis

   
  G4double expHall_x = 350*cm;
  G4double expHall_y = 350*cm;
  G4double expHall_z = 350*cm;

  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Air,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log,"expHall",0,false,0);
  experimentalHall_log->SetVisAttributes (G4VisAttributes::Invisible);



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DANCE CRYSTALS
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "DANCEcrystals.inc"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
///////////////////
// CRYSTAL B

// constructor of class B_crystal
// B_crystal(G4double a1,G4double a2,G4double a3,G4double a4,G4double a5,G4double a6,G4double a7,G4double a8){
//
// a1 - length
// a2-a6 
// a7 - front distance

//	pDz_B1=a1/2.;
//	pDy1_B1=a2;
//	pDx1_B1=a3;
//	pDx2_B1=a4; 
//	BaF_B_ratio=(1.+a1/a7);
//	pDy1_B2=a5; //0.1140406/2.*RI;
//	pDx2_B2=a6; 

RI=17.*cm;

B_crystal *BaF_B=new B_crystal(
	15.*cm,
	0.144641/2.*RI,
	0.21428*RI/2.,
	0.20783*RI,
	0.5*RI*0.1663*sin(43.295*deg),
	0.17358*RI/2.,
	17.*cm);

// let there be crystal B
G4UnionSolid *BaF_B_crystal=BaF_B->Build_B_crystal();


// now wrap material

B_crystal *BaF_B_thicker=new B_crystal(
	15.1*cm,
	0.144641/2.*RI+0.1*cm,
	0.21428*RI/2.+0.1*cm,
	0.20783*RI+0.1*cm,
	0.5*RI*0.1663*sin(43.295*deg)+0.1*cm,
	0.17358*RI/2.+0.1*cm,
	16.9*cm);

G4UnionSolid *BaF_B_crystal_thicker=BaF_B_thicker->Build_B_crystal();
//G4SubtractionSolid *BaF_B_wrap=new G4SubtractionSolid("BaF_B_wrap",BaF_B_crystal_thicker,BaF_B_crystal,0,G4Vector3D(0,0,0.1*cm));
G4SubtractionSolid *BaF_B_wrap=new G4SubtractionSolid("BaF_B_wrap",BaF_B_crystal_thicker,BaF_B_crystal,0,G4Vector3D(0,0,0));


// here holder is being born

B_crystal *BaF_B_holder=new B_crystal(
	6.0*inch,
	0.144641/2.*(RI+BaF_B->BaF_B_length),
	0.21428*(RI+BaF_B->BaF_B_length)/2.,
	0.20783*(RI+BaF_B->BaF_B_length),
	0.5*(RI+BaF_B->BaF_B_length)*0.1663*sin(43.295*deg),
	0.17358*(RI+BaF_B->BaF_B_length)/2.,
	32.*cm);


//  G4UnionSolid *BaF_B_crystal_holder_1=BaF_B_holder->Build_B_crystal();

//	G4cout << "Constructing B_holder" << G4endl;
	
  G4UnionSolid *BaF_B_crystal_holder_1=BaF_B_holder->Build_B_crystal();

//	G4cout << "Finished Constructing B_holder" << G4endl;
	

  G4Cons *BaF_B_crystal_holder_part1=new G4Cons("BaF_B_crystal_holder_part1",0*inch,1.66*inch,0*inch,1.8*inch,3.1*inch,0.*deg,360.*deg);
  G4Cons *BaF_B_crystal_holder_part2=new G4Cons("BaF_B_crystal_holder_part2",2.375*inch,4*inch,2.375*inch,4*inch,3.1*inch,0.*deg,360.*deg);


  G4SubtractionSolid *BaF_B_crystal_holder1=new G4SubtractionSolid("BaF_B_crystal_holder1",BaF_B_crystal_holder_1,BaF_B_crystal_holder_part1,0,G4Vector3D(0,1.9*cm,0));
  G4SubtractionSolid *BaF_B_crystal_holder=new G4SubtractionSolid("BaF_B_crystal_holder",BaF_B_crystal_holder1,BaF_B_crystal_holder_part2,0,G4Vector3D(0,1.9*cm,0));

  
//  BaF_B_crystal_holder_log=new G4LogicalVolume(BaF_B_crystal_holder,Al,"BaF_B_crystal_holder_log",0,0,0);


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

///////////////////
// CRYSTAL C
// creating the crystals new method

// constructor of class C_
//C_crystal(G4double a1,G4double a2,G4double a3,G4double a4,G4double a5,G4double a6,G4double a7,G4double a8){
//	BaF_C_length=a1;
//	BaF_C_front_th1=a2;
//	BaF_C_front_th2=a3;
//	BaF_C_front_th3=a4;
//	BaF_C_front_height1=a5;
//	BaF_C_front_height2=a6;
//	BaF_C_frontdist=a7;
//	shift_C=a8;

// make sure RI is defined
// perhaps for crystals one should make it constant and fixed
// but for flexibility I chose have it as a variable


  C_crystal *BaF_C=new C_crystal(
  	7.5*cm,
	  0.087482477*RI,
	  0.5*0.19423*cos(59.7*deg)*RI,
	  0.5*0.16622*sin(9.7*deg)*RI,
	  0.19423*sin(59.7*deg)*RI,
	  0.16622*cos(9.7*deg)*RI,
	  17.*cm,
	  1.0*cm);

// let there be crystal C
  G4UnionSolid *BaF_C_crystal=BaF_C->Build_C_crystal();
  
// preparing for wrapping
  C_crystal *BaF_C_thicker=new C_crystal(
	  7.5*cm,
	  0.087482477*RI+1*mm,
	  0.5*0.19423*cos(59.7*deg)*RI+1*mm,
	  0.5*0.16622*sin(9.7*deg)*RI+1*mm,
	  0.19423*sin(59.7*deg)*RI+1*mm,
	  0.16622*cos(9.7*deg)*RI+1*mm,
	  18.5*cm,
	  0.9*cm);

  G4UnionSolid *BaF_C_crystal_thicker=BaF_C_thicker->Build_C_crystal();

  G4SubtractionSolid *BaF_C_wrap=new G4SubtractionSolid("C_wrap",BaF_C_crystal_thicker,BaF_C_crystal,0,G4Vector3D(0,0,0));

// here holder is being born
  C_crystal *BaF_C_holder=new C_crystal(
	  3.0*inch,
	  0.087482477*(RI+2*BaF_C->BaF_C_length),
	  0.5*0.19423*cos(59.7*deg)*(RI+2*BaF_C->BaF_C_length),
	  0.5*0.16622*sin(9.7*deg)*(RI+2*BaF_C->BaF_C_length),
	  0.19423*sin(59.7*deg)*(RI+2*BaF_C->BaF_C_length),
	  0.16622*cos(9.7*deg)*(RI+2*BaF_C->BaF_C_length),
	  32.*cm,
	  1.0*cm);

  G4UnionSolid *BaF_C_crystal_holder_1=BaF_C_holder->Build_C_crystal();

  G4Cons *BaF_C_crystal_holder_part1=new G4Cons("BaF_C_crystal_holder_part1",0*inch,1.70*inch,0*inch,1.9*inch,3.1*inch,0.*deg,360.*deg);
  G4Cons *BaF_C_crystal_holder_part2=new G4Cons("BaF_C_crystal_holder_part2",2.375*inch,4*inch,2.375*inch,4*inch,3.1*inch,0.*deg,360.*deg);


  G4SubtractionSolid *BaF_C_crystal_holder1=new G4SubtractionSolid("BaF_C_crystal_holder1",BaF_C_crystal_holder_1,BaF_C_crystal_holder_part1,0,G4Vector3D(0,-1.0*cm,0));
  G4SubtractionSolid *BaF_C_crystal_holder=new G4SubtractionSolid("BaF_C_crystal_holder",BaF_C_crystal_holder1,BaF_C_crystal_holder_part2,0,G4Vector3D(0,-1.0*cm,0));

  
  BaF_C_crystal_holder_log=new G4LogicalVolume(BaF_C_crystal_holder,Al,"BaF_C_crystal_holder_log",0,0,0);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  DANCE Detector - crystal holders

BaF_A_crystal_holder_log=new G4LogicalVolume(BaF_A_crystal_holder,Al,"BaF_A_crystal_holder_log",0,0,0);
BaF_D_crystal_holder_log=new G4LogicalVolume(BaF_D_crystal_holder,Al,"BaF_D_crystal_holder_log",0,0,0);
BaF_C_crystal_holder_log=new G4LogicalVolume(BaF_C_crystal_holder,Al,"BaF_C_crystal_holder_log",0,0,0);
BaF_B_crystal_holder_log=new G4LogicalVolume(BaF_B_crystal_holder,Al,"BaF_B_crystal_holder_log",0,0,0);

BaF_C_crystal_holder_log->SetVisAttributes(sAl_VA);
BaF_A_crystal_holder_log->SetVisAttributes(sAl_VA);
BaF_D_crystal_holder_log->SetVisAttributes(sAl_VA);
BaF_B_crystal_holder_log->SetVisAttributes(sAl_VA);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  DANCE Detector - wraping material

BaF_C_wrap_log=new G4LogicalVolume(BaF_C_wrap,PVC,"BaF_C_wrap_log",0,0,0);
BaF_C_wrap_log->SetVisAttributes(sVA_Wrap);

BaF_B_wrap_log=new G4LogicalVolume(BaF_B_wrap,PVC,"BaF_B_wrap_log",0,0,0);
BaF_B_wrap_log->SetVisAttributes(sVA_Wrap);

BaF_A_wrap_log=new G4LogicalVolume(BaF_A_wrap,PVC,"BaF_A_wrap_log",0,0,0);
BaF_A_wrap_log->SetVisAttributes(sVA_Wrap);

BaF_D_wrap_log=new G4LogicalVolume(BaF_D_wrap,PVC,"BaF_D_wrap_log",0,0,0);
BaF_D_wrap_log->SetVisAttributes(sVA_Wrap);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// beam pipe
  	G4RotationMatrix *BeamPipeRM=new G4RotationMatrix();
 	BeamPipeRM->rotateY(90.*deg);

  BeamPipe_log= new G4LogicalVolume(BeamPipe,Al,"BeamPipe_log",0,0,0);
  BeamPipe_phys = new G4PVPlacement(BeamPipeRM,G4ThreeVector(0,0,0),BeamPipe_log,"Al",experimentalHall_log,false,0);

  BeamPipe_log->SetVisAttributes(BeamPipeVA);
//  BeamPipe_log->SetVisAttributes(hideVA);

// vacuum in pipe

// HERE BEAMPIPE - Aluminium 2.1825 - 2.5 cm

//G4double pRMin=0.5*1.75*inch; //2.1825*cm;
//G4double pRMax=0.5*1.83*inch; //2.5*cm;
pDz=175*cm;
//G4double pSPhi=0.*deg;
//G4double pDPhi=360.*deg;

  G4Tubs *BeamPipeVacuum=new G4Tubs("BeamPipe",0.,pRMin,0.5*350.*cm,pSPhi,pDPhi);

if(GPP->BeamPipeVacuum){
  BeamPipeVacuum_log= new G4LogicalVolume(BeamPipeVacuum,Vacuum,"BeamPipe_log",0,0,0);
  BeamPipeVacuum_phys = new G4PVPlacement(BeamPipeRM,G4ThreeVector(0,0,0),BeamPipeVacuum_log,"BeamPipeVacuum",experimentalHall_log,false,0);
}
//  BeamPipeVacuum_log->SetVisAttributes(BeamPipeVA);

//RTH

  RTH_log= new G4LogicalVolume(RTH,Al,"RTH_log",0,0,0);
  RTH_phys = new G4PVPlacement(BeamPipeRM,G4ThreeVector(0,0,0),RTH_log,"RTH",experimentalHall_log,false,0);


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lithium hydrate sphere

if(GPP->LiHShell[0]){
  LiHSphere_log 	= 	new G4LogicalVolume(LiHSphere2,LiH,"LiHSphere_log",0,0,0);
  LiHSphere_phys 	= 	new G4PVPlacement(BeamPipeRM,G4ThreeVector(),LiHSphere_log,"LiH",experimentalHall_log,false,0);


  if(GPP->LiHShell[1]) LiHSphere_log->SetVisAttributes(LiHVA);
  else LiHSphere_log->SetVisAttributes(hideVA);
}

if(GPP->PVCShell){
	G4Sphere* PVCShell	=	new G4Sphere("PVCSphere",16.6*cm,16.7*cm,0.*deg,360.*deg,0.*deg,180.*deg);
	PVCShell_log 		= 	new G4LogicalVolume(PVCShell,PVC,"PVCShell_log",0,0,0);
  	PVCShell_phys 		= 	new G4PVPlacement(0,G4ThreeVector(),PVCShell_log,"PVCShell",experimentalHall_log,false,0);
	PVCShell_log->SetVisAttributes(hideVA);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//SuppSphere_log 	= 	new G4LogicalVolume(SuppSphere,Al,"SuppSphere_log",0,0,0);
//SuppSphere_phys = 	new G4PVPlacement(0,G4ThreeVector(),SuppSphere_log,"Al",experimentalHall_log,false,0);


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Target in this case is just a steel disk 0.2 mm

	G4RotationMatrix *TargetRM=new G4RotationMatrix();
	TargetRM->rotateY(90.*deg);

	pRMin=0*cm;
	pRMax=0.7*cm;
	pDz=10*cm;
	pSPhi=0.*deg;
	pDPhi=360.*deg;

	G4Tubs *Target=new G4Tubs("Target",pRMin,pRMax,pDz,pSPhi,pDPhi);

//  Target_log= new G4LogicalVolume(Target,ssteel,"Target_log",0,0,0);


  Target_log= new G4LogicalVolume(Target,Au,"Target_log",0,0,0);
if(GPP->BeamPipeVacuum){
//  Target_phys = new G4PVPlacement(TargetRM,G4ThreeVector(0,0,0),Target_log,"TargetMat",BeamPipeVacuum_log,false,0);
}
else{
//  Target_phys = new G4PVPlacement(TargetRM,G4ThreeVector(0,0,0),Target_log,"TargetMat",experimentalHall_log,false,0);
} 
 
//  Target_log->SetVisAttributes(TargetVA);


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Photomultiplier window

  	G4double PMT_thick   =   1.0*mm ; // Thickness of PMT window
  	G4double PMT_curv    =  65.5*mm ; // Radius of curvature of PMT window
//	G4double StartTheta  = (180.0-31.2)*pi/180. ;
//  	G4double EndTheta    = 31.2*pi/180. ;

//G4Sphere *PMTwindow ;
//G4Sphere *PMTwindow = new G4Sphere("PMT_solid",PMT_curv-PMT_thick,PMT_curv,0.0,twopi,StartTheta,EndTheta);

  	G4Tubs *PMTwindow=new G4Tubs("PMTwindow",0.,1.65*2.54*cm,2*mm,0.*deg,360*deg);

if(GPP->PMLense[0]){
	PMT_log = new G4LogicalVolume(PMTwindow,Quartz,"PMT_log",0,0,0);
	if(GPP->PMLense[1]) PMT_log->SetVisAttributes(Crystal_A_VA);
	else PMT_log->SetVisAttributes(hideVA);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Read in angles from input file "angles.txt"


	G4double BaF_theta[200],BaF_phi[200],BaF_psi[200],BaF_R[200];
	G4int Cr_Type[200],Cr_ID[200];

//	Cr_Type[]=1...A 2...B 3...C 4...D

	ifstream in("config/angles.txt");

	for(int i=0;i<162;i++){

		in >> Cr_Type[i] >> Cr_ID[i] >> BaF_theta[i] >> BaF_phi[i] >> BaF_psi[i] >> BaF_R[i];
		
		GPP->SetTheta(Cr_ID[i],BaF_theta[i]);
		GPP->SetPhi(Cr_ID[i],BaF_phi[i]);
		
		BaF_theta[i]=(BaF_theta[i])*deg;
		BaF_phi[i]=(BaF_phi[i]-180.)*deg;
		BaF_psi[i]=BaF_psi[i]*1.*deg;
		//else BaF_psi[i]=(BaF_psi[159-i]*1.)*deg;

		BaF_R[i]=GPP->CrystalDistance*cm; //17.5*cm; //BaF_R[i]*cm;

//		GPP->SetGeantID2D(Cr_Type[i],Cr_ID[i],i);
		
		//G4cout << Cr_Type[i] << "	" << Cr_ID[i] << "	" <<BaF_theta[i] << "	" << BaF_phi[i] << G4endl;
	}
	in.close();

//



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// defining variables for geometry calculations

	int i=0;

	G4Vector3D tC,tD,tA,tB;
	G4Point3D a3D,b3D,c3D,d3D;
	G4double angle;
	G4ThreeVector *base_trans[200];
	G4RotationMatrix* BaF_C_rot[200];
	G4Transform3D BaF_C_rotation[200];
	G4Vector3D yz(1,0.,0.); // beacause of o turn 90 degrees
	G4Vector3D rot_yz;

	double BaFPos_x ;
	double BaFPos_y ;
	double BaFPos_z ;
	int A_c,B_c,C_c,D_c;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// LOOP THROUGH ANGLES AND BUILD DANCE !!!

int holes_counter=0;
int pmt_counter=0;

A_c=B_c=C_c=D_c=0;  // counters for each tyupe of the crystal

//ofstream internalmap("InternalMap.dat");

for(i=0;i<80;i++){

//BaF_phi[i]=0*deg+i*0.*deg;;
//BaF_theta[i]=0.*deg+i*15*deg;
//BaF_psi[i]=0.*1.*deg;
//BaF_frontdist=BaF_R[i];

// case A-crystal
if(Cr_Type[i]==1 || Cr_Type[i]==4){
BaF_frontdist=BaF_R[i];
}
if(Cr_Type[i]==2 || Cr_Type[i]==3){
BaF_frontdist=(BaF_R[i]+BaF_length/2.);
}

BaFPos_x = BaF_frontdist*cos(BaF_theta[i]);
BaFPos_y = BaF_frontdist*sin(BaF_theta[i])*sin(BaF_phi[i]);
BaFPos_z = BaF_frontdist*sin(BaF_theta[i])*cos(BaF_phi[i]);

// here we will sit at the detector frame in its center
// and calculate the cross product of vector parralel 
// to beam(internal vector of detector)xvector towards the focus(0,0,0)
// this way we obtained the axis of final rotation 
// angle of rotation is given by the angle between the parallel to beam
// and vector toward the focus

// vectors
// initial

tA.set(1,0,0);
// desired
tB.set(-BaFPos_x,-BaFPos_y,-BaFPos_z);

tC=tA.cross(tB);
angle=tA.angle(tB);

//G4double twist=tC.angle(G4Vector3D(0.,0.,BaFPos_z/fabs(BaFPos_z)))*deg;


a3D.set(BaFPos_x,BaFPos_y,BaFPos_z);
b3D.set(BaFPos_x+tC.x(),BaFPos_y+tC.y(),BaFPos_z+tC.z());


//c3D.set(0,0,0);
//d3D.set(BaFPos_x,BaFPos_y,BaFPos_z);

G4Point3D origin(0,0,0);
G4Point3D phiAxis(0,0,1);
G4Point3D thetaAxis(0.,1.,0.);
G4Point3D psiAxis(1,0,0);

base_trans[i]=new G4ThreeVector(BaFPos_x, BaFPos_y, BaFPos_z);
BaF_C_rot[i]=new G4RotationMatrix();


BaF_C_rot[i]->rotateY(-90.*deg);
BaF_C_rot[i]->rotateX(-BaF_phi[i]+BaF_psi[i]);

BaF_C_rot[159-i]=new G4RotationMatrix();
BaF_C_rot[159-i]->rotateY(90.*deg);
BaF_C_rot[159-i]->rotateX(-BaF_phi[i]+BaF_psi[i]);

G4Transform3D posB1_1(*(BaF_C_rot[i]),G4ThreeVector(BaFPos_x,BaFPos_y,BaFPos_z));
G4Transform3D final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,0));

G4Transform3D final_1_rot;
G4Transform3D posB1_1_rot;
// here important correction

G4double perBaFPos_x = 0.;
G4double perBaFPos_y = 1;
G4double perBaFPos_z = -BaFPos_z/BaFPos_y;

tD=tB.cross(G4Vector3D(perBaFPos_x,perBaFPos_y,perBaFPos_z));  // vector perpendicular to position R or tB(here)
G4Vector3D tE=G4Vector3D(0.,0.,1.);


rot_yz=final_1*tE;
G4double twist_test=rot_yz.angle(tD);  // cross product of the tB and vector tanget to the circle at (x,0,0) with R=(0,perBaFPos_y,perBaFPos_z);

// apply correction and recalculate transforms

BaF_C_rot[i]->rotateX(-twist_test);  // rotation matrix with the proper twist
BaF_C_rot[159-i]->rotateX(-twist_test+180.*deg);
if(Cr_Type[i]==1) BaF_C_rot[159-i]->rotateX(36*deg);

posB1_1=G4Transform3D(*(BaF_C_rot[i]),G4ThreeVector(BaFPos_x,BaFPos_y,BaFPos_z));
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,0));

posB1_1_rot=G4Transform3D(*(BaF_C_rot[159-i]),G4ThreeVector(-BaFPos_x,-BaFPos_y,-BaFPos_z));
final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0,0));


// here the support structuire is created 

/*
if(Cr_Type[i]<5){
tube_subtract=new G4Tubs("Hole",0.,0.5*4.813*inch,50.*cm,pSPhi,pDPhi);
SuppSphereHoles[holes_counter+1]= new G4SubtractionSolid("SuppHoles", SuppSphereHoles[holes_counter], tube_subtract,final_1);
holes_counter++;

G4Tubs *tube_subtract2=new G4Tubs("Hole",0.,0.5*4.813*inch,50.*cm,pSPhi,pDPhi);
SuppSphereHoles[holes_counter+1]= new G4SubtractionSolid("SuppHoles", SuppSphereHoles[holes_counter], tube_subtract2,final_1_rot);
holes_counter++;
}
*/

/*
else{
hexa_subtract=new G4Polyhedra("Hole2b",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlane_Dhole,rInner_Dhole,rOuter_Dhole);
SuppSphereHoles[holes_counter+1]= new G4SubtractionSolid("SuppHoles", SuppSphereHoles[holes_counter], hexa_subtract,final_1);

holes_counter++;
G4Polyhedra *hexa_subtract2=new G4Polyhedra("Hole2b",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlane_Dhole,rInner_Dhole,rOuter_Dhole);
SuppSphereHoles[holes_counter+1]= new G4SubtractionSolid("SuppHoles", SuppSphereHoles[holes_counter], hexa_subtract2,final_1_rot);
holes_counter++;
}

*/


if(Cr_Type[i]==1){
//posB1_1=G4Transform3D(*(BaF_C_rot[i]),G4ThreeVector(BaFPos_x,BaFPos_y,BaFPos_z));
//final_1=G4Rotate3D(angle,sa3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,0));

//crystal
BaF_A_crystal_log[A_c]= new G4LogicalVolume(BaF_A_crystal,BaF,"BaF_A_crystal_log",0,0,0);
//BaF_A_crystal_log[A_c]= new G4LogicalVolume(BaF_A_crystal,LaBr3,"BaF_A_crystal_log",0,0,0);
BaF_A_crystal_phys[A_c] = new G4PVPlacement(final_1,BaF_A_crystal_log[A_c],"BaF",experimentalHall_log,false,A_c);

// wrapping
BaF_A_wrap_phys=new G4PVPlacement(final_1,BaF_A_wrap_log,"PVC_A",experimentalHall_log,false,A_c);

//holder
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,BaF_length+3*inch+1*mm));
BaF_A_crystal_holder_phys = new G4PVPlacement(final_1,BaF_A_crystal_holder_log,"Al_holder_A",experimentalHall_log,false,A_c);

if(Cr_ID[i]>=GPP->ShowCrystals[0] && Cr_ID[i]<=GPP->ShowCrystals[1]){ 
	BaF_A_crystal_log[A_c]->SetVisAttributes(Crystal_A_VA);
	BaF_A_crystal_holder_log->SetVisAttributes(sAl_VA);
}
else {
	BaF_A_crystal_log[A_c]->SetVisAttributes(hideVA);
	BaF_A_crystal_holder_log->SetVisAttributes(hideVA);
}

//PMT window
if(GPP->PMLense[0]){
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,BaF_length+2*mm));
PMT_phys=new G4PVPlacement(final_1,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
pmt_counter++;
}

GPP->SetID2D(Cr_Type[i]-1,A_c,Cr_ID[i]);

A_c++;

// tranform 180 degrees

// crystal
BaF_A_crystal_log[A_c]= new G4LogicalVolume(BaF_A_crystal,BaF,"BaF_A_crystal_log",0,0,0);
//BaF_A_crystal_log[A_c]= new G4LogicalVolume(BaF_A_crystal,LaBr3,"BaF_A_crystal_log",0,0,0);
BaF_A_crystal_phys[A_c] = new G4PVPlacement(final_1_rot,BaF_A_crystal_log[A_c],"BaF",experimentalHall_log,false,A_c);

// wrapping
BaF_A_wrap_phys=new G4PVPlacement(final_1_rot,BaF_A_wrap_log,"PVC_A",experimentalHall_log,false,A_c);

// holder
final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0,BaF_length+3.*inch+1*mm));
BaF_A_crystal_holder_phys = new G4PVPlacement(final_1_rot,BaF_A_crystal_holder_log,"Al_holder_A",experimentalHall_log,false,A_c);

if(GPP->GetSymID(Cr_ID[i])>=GPP->ShowCrystals[0] && GPP->GetSymID(Cr_ID[i])<=GPP->ShowCrystals[1]){ 
	BaF_A_crystal_log[A_c]->SetVisAttributes(Crystal_A_VA);
	BaF_A_crystal_holder_log->SetVisAttributes(sAl_VA);
}
else {
	BaF_A_crystal_log[A_c]->SetVisAttributes(hideVA);
	BaF_A_crystal_holder_log->SetVisAttributes(hideVA);
}

//PMT window
if(GPP->PMLense[0]){
	final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0,BaF_length+2*mm));
	PMT_phys=new G4PVPlacement(final_1_rot,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

GPP->SetID2D(Cr_Type[i]-1,A_c,GPP->GetSymID(Cr_ID[i]));
A_c++;

}

if(Cr_Type[i]==2){

final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,BaF_B->trans,0));

// crystal
BaF_B_crystal_log[B_c]= new G4LogicalVolume(BaF_B_crystal,BaF,"BaF_B_crystal_log",0,0,0);
//BaF_B_crystal_log[B_c]= new G4LogicalVolume(BaF_B_crystal,LaBr3,"BaF_B_crystal_log",0,0,0);
BaF_B_crystal_phys[B_c] = new G4PVPlacement(final_1,BaF_B_crystal_log[B_c],"BaF",experimentalHall_log,false,B_c);


// wrapping
BaF_B_wrap_phys=new G4PVPlacement(final_1,BaF_B_wrap_log,"PVC_B",experimentalHall_log,false,C_c);


//PMT window
if(GPP->PMLense[0]){
	final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,-0.5*cm,0.5*BaF_length+2.5*mm));
	PMT_phys=new G4PVPlacement(final_1,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

//holder
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,BaF_B_holder->trans,0.5*BaF_length+0.5*BaF_B_holder->BaF_B_length+0.7*mm));
BaF_B_crystal_holder_phys = new G4PVPlacement(final_1,BaF_B_crystal_holder_log,"Al_holder_B1",experimentalHall_log,false,A_c);

if(Cr_ID[i]>=GPP->ShowCrystals[0] && Cr_ID[i]<=GPP->ShowCrystals[1]){ 
	BaF_B_crystal_holder_log->SetVisAttributes(sAl_VA);
	BaF_B_crystal_log[B_c]->SetVisAttributes(Crystal_B_VA);
}
else{
	BaF_B_crystal_holder_log->SetVisAttributes(hideVA);
	BaF_B_crystal_log[B_c]->SetVisAttributes(hideVA);
}

GPP->SetID2D(Cr_Type[i]-1,B_c,Cr_ID[i]);

B_c++;

final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,BaF_B->trans,0));

// crystal
BaF_B_crystal_log[B_c]= new G4LogicalVolume(BaF_B_crystal,BaF,"BaF_B_crystal_log",0,0,0);
//BaF_B_crystal_log[B_c]= new G4LogicalVolume(BaF_B_crystal,LaBr3,"BaF_B_crystal_log",0,0,0);
BaF_B_crystal_phys[B_c] = new G4PVPlacement(final_1_rot,BaF_B_crystal_log[B_c],"BaF",experimentalHall_log,false,B_c);

// wrapping
BaF_B_wrap_phys=new G4PVPlacement(final_1_rot,BaF_B_wrap_log,"PVC_B",experimentalHall_log,false,C_c);

//PMT window
if(GPP->PMLense[0]){
	final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,-0.5*cm,0.5*BaF_length+2.5*mm));
	PMT_phys=new G4PVPlacement(final_1_rot,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

//holder
final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,BaF_B_holder->trans,0.5*BaF_length+0.5*BaF_B_holder->BaF_B_length+0.7*mm));
BaF_B_crystal_holder_phys = new G4PVPlacement(final_1_rot,BaF_B_crystal_holder_log,"Al_holder_B2",experimentalHall_log,false,A_c);

if(GPP->GetSymID(Cr_ID[i])>=GPP->ShowCrystals[0] && GPP->GetSymID(Cr_ID[i])<=GPP->ShowCrystals[1]){ 
	BaF_B_crystal_holder_log->SetVisAttributes(sAl_VA);
	BaF_B_crystal_log[B_c]->SetVisAttributes(Crystal_B_VA);
}
else{
	BaF_B_crystal_holder_log->SetVisAttributes(hideVA);
	BaF_B_crystal_log[B_c]->SetVisAttributes(hideVA);
}

GPP->SetID2D(Cr_Type[i]-1,B_c,GPP->GetSymID(Cr_ID[i]));

B_c++;


}


if(Cr_Type[i]==3){
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0-BaF_C->mid_trans,0));
//G4Transform3D final_2=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,left_trans,0));
//G4Transform3D final_3=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,right_trans,0));

BaF_C_crystal_log[C_c]= new G4LogicalVolume(BaF_C_crystal,BaF,"BaF_C_crystal_log",0,0,0);
//BaF_C_crystal_log[C_c]= new G4LogicalVolume(BaF_C_crystal,LaBr3,"BaF_C_crystal_log",0,0,0);
BaF_C_crystal_phys[C_c] = new G4PVPlacement(final_1,BaF_C_crystal_log[C_c],"BaF",experimentalHall_log,false,C_c);
//BaF_C_crystal_log[C_c]->SetVisAttributes(Crystal_C_VA);

// wrapping
BaF_C_wrap_phys=new G4PVPlacement(final_1,BaF_C_wrap_log,"PVC",experimentalHall_log,false,C_c);

// holder C
//posB1_1=G4Transform3D(*(BaF_C_rot[i]),G4ThreeVector(32./17.*BaFPos_x,32./17.*BaFPos_y,32./17.*BaFPos_z));
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0-BaF_C_holder->mid_trans+1.1*cm,0.5*(2*BaF_C->BaF_C_length+2*BaF_C_holder->BaF_C_length)));
BaF_C_crystal_holder_phys = new G4PVPlacement(final_1,BaF_C_crystal_holder_log,"Al_holder_C",experimentalHall_log,false,C_c);

if(Cr_ID[i]>=GPP->ShowCrystals[0] && Cr_ID[i]<=GPP->ShowCrystals[1]){ 
	BaF_C_crystal_holder_log->SetVisAttributes(sAl_VA);
	BaF_C_crystal_log[C_c]->SetVisAttributes(Crystal_C_VA);
}
else{
	BaF_C_crystal_holder_log->SetVisAttributes(hideVA);
	BaF_C_crystal_log[C_c]->SetVisAttributes(hideVA);
}

// PMT window
if(GPP->PMLense[0]){
	final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0.5*cm,0.5*(2*BaF_C->BaF_C_length+2.5*mm)));
	PMT_phys=new G4PVPlacement(final_1,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

GPP->SetID2D(Cr_Type[i]-1,C_c,Cr_ID[i]);

C_c++;

final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,-BaF_C->mid_trans,0));
//G4Transform3D final_2_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,left_trans,0));
//G4Transform3D final_3_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,right_trans,0));


BaF_C_crystal_log[C_c]= new G4LogicalVolume(BaF_C_crystal,BaF,"BaF_C_crystal_log",0,0,0);
//BaF_C_crystal_log[C_c]= new G4LogicalVolume(BaF_C_crystal,LaBr3,"BaF_C_crystal_log",0,0,0);
BaF_C_crystal_phys[C_c] = new G4PVPlacement(final_1_rot,BaF_C_crystal_log[C_c],"BaF",experimentalHall_log,false,C_c);
//BaF_C_crystal_log[C_c]->SetVisAttributes(Crystal_C_VA);

// wrapping
BaF_C_wrap_phys=new G4PVPlacement(final_1_rot,BaF_C_wrap_log,"PVC",experimentalHall_log,false,C_c);

// holder
final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0-BaF_C_holder->mid_trans+1.1*cm,0.5*(2*BaF_C->BaF_C_length+2*BaF_C_holder->BaF_C_length)));
BaF_C_crystal_holder_phys = new G4PVPlacement(final_1_rot,BaF_C_crystal_holder_log,"Al_holder_C",experimentalHall_log,false,C_c);

if(GPP->GetSymID(Cr_ID[i])>=GPP->ShowCrystals[0] && GPP->GetSymID(Cr_ID[i])<=GPP->ShowCrystals[1]){ 
	BaF_C_crystal_holder_log->SetVisAttributes(sAl_VA);
	BaF_C_crystal_log[C_c]->SetVisAttributes(Crystal_C_VA);
}
else{
	BaF_C_crystal_holder_log->SetVisAttributes(hideVA);
	BaF_C_crystal_log[C_c]->SetVisAttributes(hideVA);
}

// PMT window
if(GPP->PMLense[0]){
	final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0.5*cm,0.5*(2*BaF_C->BaF_C_length+2.5*mm)));
	PMT_phys=new G4PVPlacement(final_1_rot,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

GPP->SetID2D(Cr_Type[i]-1,C_c,GPP->GetSymID(Cr_ID[i]));

C_c++;


}

if(Cr_Type[i]==4){

BaF_D_crystal_log[D_c]= new G4LogicalVolume(BaF_D_crystal,BaF,"BaF_D_crystal_log",0,0,0);
//BaF_D_crystal_log[D_c]= new G4LogicalVolume(BaF_D_crystal,LaBr3,"BaF_D_crystal_log",0,0,0);
BaF_D_crystal_phys[D_c] = new G4PVPlacement(final_1,BaF_D_crystal_log[D_c],"BaF",experimentalHall_log,false,D_c);

//BaF_D_crystal_log[D_c]->SetVisAttributes(Crystal_D_VA);

// wrapping
BaF_D_wrap_phys=new G4PVPlacement(final_1,BaF_D_wrap_log,"PVC_A",experimentalHall_log,false,D_c);

// here holder
final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,BaF_length+3.*inch));
BaF_D_crystal_holder_phys = new G4PVPlacement(final_1,BaF_D_crystal_holder_log,"Al",experimentalHall_log,false,D_c);

if(Cr_ID[i]>=GPP->ShowCrystals[0] && Cr_ID[i]<=GPP->ShowCrystals[1]){ 
	BaF_D_crystal_holder_log->SetVisAttributes(sAl_VA);
	BaF_D_crystal_log[D_c]->SetVisAttributes(Crystal_D_VA);
}
else{
	BaF_D_crystal_holder_log->SetVisAttributes(hideVA);
	BaF_D_crystal_log[D_c]->SetVisAttributes(hideVA);
}

// PMT window
if(GPP->PMLense[0]){
	final_1=G4Rotate3D(angle,a3D,b3D)*posB1_1*G4Translate3D(G4Vector3D(0,0,BaF_length+2.5*mm));
	PMT_phys=new G4PVPlacement(final_1,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

GPP->SetID2D(Cr_Type[i]-1,D_c,Cr_ID[i]);

D_c++;

// G4cout << "D:	" << D_c << G4endl;

//final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0,0));

BaF_D_crystal_log[D_c]= new G4LogicalVolume(BaF_D_crystal,BaF,"BaF_D_crystal_log",0,0,0);
//BaF_D_crystal_log[D_c]= new G4LogicalVolume(BaF_D_crystal,LaBr3,"BaF_D_crystal_log",0,0,0);
BaF_D_crystal_phys[D_c] = new G4PVPlacement(final_1_rot,BaF_D_crystal_log[D_c],"BaF",experimentalHall_log,false,D_c);
//BaF_D_crystal_log[D_c]->SetVisAttributes(Crystal_D_VA);

// wrapping
BaF_D_wrap_phys=new G4PVPlacement(final_1_rot,BaF_D_wrap_log,"PVC_A",experimentalHall_log,false,D_c);

final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0,BaF_length+3*inch));
BaF_D_crystal_holder_phys = new G4PVPlacement(final_1_rot,BaF_D_crystal_holder_log,"Al",experimentalHall_log,false,D_c);

if(GPP->GetSymID(Cr_ID[i])>=GPP->ShowCrystals[0] && GPP->GetSymID(Cr_ID[i])<=GPP->ShowCrystals[1]){ 
	BaF_D_crystal_holder_log->SetVisAttributes(sAl_VA);
	BaF_D_crystal_log[D_c]->SetVisAttributes(Crystal_D_VA);
}
else{
	BaF_D_crystal_holder_log->SetVisAttributes(hideVA);
	BaF_D_crystal_log[D_c]->SetVisAttributes(hideVA);
}
// PMT window
if(GPP->PMLense[0]){
	final_1_rot=G4Rotate3D(-angle,-a3D,-b3D)*posB1_1_rot*G4Translate3D(G4Vector3D(0,0,BaF_length+2.5*mm));
	PMT_phys=new G4PVPlacement(final_1_rot,PMT_log,"Quartz",experimentalHall_log,false,pmt_counter);
	pmt_counter++;
}

GPP->SetID2D(Cr_Type[i]-1,D_c,GPP->GetSymID(Cr_ID[i]));

D_c++;
// G4cout << "D:	" << D_c << G4endl;

}

}

if(GPP->SupportingSphere[0]){
  SuppSphere_log = 	new G4LogicalVolume(SuppSphereHoles[holes_counter],Al,"SuppSphere_log",0,0,0);
  SuppSphere_phys = 	new G4PVPlacement(0,G4ThreeVector(),SuppSphere_log,"Al",experimentalHall_log,false,0);

//  if(GPP->SupportingSphere) SuppSphere_log->SetVisAttributes(sAl_VA);
  SuppSphere_log->SetVisAttributes(hideVA);
}

A_counter=A_c;
B_counter=B_c;
C_counter=C_c;
D_counter=D_c;

G4SDManager* SDman = G4SDManager::GetSDMpointer();

int counter=0;

for(i=0;i<A_c;i++){

//    	BaF_A_crystal_log[i]->SetRegion(aRegion);
    	G4String tName = "/"+detNameA[i];

    	G4VSensitiveDetector* det = new DANCEDetector(tName);
        SDman->AddNewDetector(det);

    	BaF_A_crystal_log[i]->SetSensitiveDetector(det);

	GPP->SetGeantID2D(0,i,counter);
	counter++;
	
//    	G4cout << "Setting the multidetector name to:" << tName << G4endl;
}  


for(i=0;i<B_c;i++){

//    	BaF_B_crystal_log[i]->SetRegion(bRegion);
    	G4String tName = "/"+detNameB[i];

    	G4VSensitiveDetector* det = new DANCEDetector(tName);
        SDman->AddNewDetector(det);

    	BaF_B_crystal_log[i]->SetSensitiveDetector(det);

	GPP->SetGeantID2D(1,i,counter);
	counter++;
	
//    	G4cout << "Setting the multidetector name to:" << tName << G4endl;
}  


for(i=0;i<C_c;i++){
//    	BaF_C_crystal_log[i]->SetRegion(cRegion);
    	G4String tName = "/"+detNameC[i];

    	G4VSensitiveDetector* det = new DANCEDetector(tName);
        SDman->AddNewDetector(det);

    	BaF_C_crystal_log[i]->SetSensitiveDetector(det);

	GPP->SetGeantID2D(2,i,counter);
	counter++;
	
//    	G4cout << "Setting the multidetector name to:" << tName << G4endl;
}  


for(i=0;i<D_c;i++){
//    	BaF_D_crystal_log[i]->SetRegion(dRegion);
    	G4String tName = "/"+detNameD[i];

    	G4VSensitiveDetector* det = new DANCEDetector(tName);
        SDman->AddNewDetector(det);

    	BaF_D_crystal_log[i]->SetSensitiveDetector(det);

	GPP->SetGeantID2D(3,i,counter);
	counter++;
	
//    	G4cout << "Setting the multidetector name to:" << tName << G4endl;


}  

// G4cout << "The final Map for Crystal IDing: " << G4endl;

/*
for(i=0;i<counter;i++){

G4cout << "Logical volume: " << i << "	" << "DetType: " << GPP->GetGeantID2D_DetType(i) << " DetNo: " <<
GPP->GetGeantID2D_DetNo(i) << "  CrystalID: " << GPP->GetID2D(GPP->GetGeantID2D_DetType(i),GPP->GetGeantID2D_DetNo(i)) <<
G4endl;

}
*/
   return experimentalHall_phys;
}

