
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



#ifndef DANCEDetectorConstruction_H
#define DANCEDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DANCEDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DANCEDetectorConstruction();
    ~DANCEDetectorConstruction();

    G4VPhysicalVolume* Construct();

//	G4SensitiveDetector *BaFSD;
	
  private:
    
    G4bool constructed;
    // Logical volumes
    //
    G4LogicalVolume* experimentalHall_log;

    G4LogicalVolume* BaF_A_crystal_log[20];
    G4LogicalVolume* BaF_B_crystal_log[70];
    G4LogicalVolume* BaF_C_crystal_log[70];
    G4LogicalVolume* BaF_D_crystal_log[40];
	G4LogicalVolume* BeamPipe_log;
	G4LogicalVolume* BeamPipeVacuum_log;
	G4LogicalVolume* Target_log;
	G4LogicalVolume* BaF_A_crystal_holder_log;
	G4LogicalVolume* BaF_B_crystal_holder_log;
	G4LogicalVolume* BaF_C_crystal_holder_log;
	G4LogicalVolume* BaF_D_crystal_holder_log;
	G4LogicalVolume* LiHSphere_log;
	G4LogicalVolume* SuppSphere_log;
	G4LogicalVolume* PMT_log;
	G4LogicalVolume* RTH_log;
	G4LogicalVolume* PVCShell_log;

	G4LogicalVolume* BaF_C_wrap_log;
	G4LogicalVolume* BaF_A_wrap_log;
	G4LogicalVolume* BaF_B_wrap_log;
	G4LogicalVolume* BaF_D_wrap_log;


    // Physical volumes
    //
    G4VPhysicalVolume* experimentalHall_phys;   

    G4VPhysicalVolume* BaF_A_crystal_phys[20];
    G4VPhysicalVolume* BaF_B_crystal_phys[70];
    G4VPhysicalVolume* BaF_C_crystal_phys[70];
    G4VPhysicalVolume* BaF_D_crystal_phys[40];
    G4VPhysicalVolume* BeamPipe_phys;
    G4VPhysicalVolume* BeamPipeVacuum_phys;
    G4VPhysicalVolume* Target_phys;
    G4VPhysicalVolume* BaF_A_crystal_holder_phys;
    G4VPhysicalVolume* BaF_B_crystal_holder_phys;
    G4VPhysicalVolume* BaF_C_crystal_holder_phys;
    G4VPhysicalVolume* BaF_D_crystal_holder_phys;
    G4VPhysicalVolume* PMT_phys;
    G4VPhysicalVolume* RTH_phys;

    G4VPhysicalVolume* BaF_C_wrap_phys;
    G4VPhysicalVolume* BaF_A_wrap_phys;
    G4VPhysicalVolume* BaF_B_wrap_phys;
    G4VPhysicalVolume* BaF_D_wrap_phys;

    
    G4VPhysicalVolume*  LiHSphere_phys; 
    G4VPhysicalVolume*  SuppSphere_phys; 
    G4VPhysicalVolume*  PVCShell_phys; 

    G4String detNameA[20];
    G4String detNameB[60];
    G4String detNameC[60];
    G4String detNameD[40];
 

};

#endif

