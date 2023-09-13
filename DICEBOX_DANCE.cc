
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


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "Randomize.hh"


#include "DANCEDetectorConstruction.hh"
#include "DANCEPhysicsList.hh"
#include "DANCEPrimaryGeneratorAction.hh"
#include "DANCEEventAction.hh"

#include "DANCENNNearestNeighbors.hh"

#include "DANCEGlobal.hh"


int A_counter,B_counter,C_counter,D_counter;
int AutoSave;

int main(int argc,char** argv) 
{

	G4cout << "-----------------------------------------------------------------------------" <<G4endl;
	G4cout << "|		  Start of the DANCE - GEANT4 simulation			|" <<G4endl;
	G4cout << "|										|" << G4endl;
	G4cout << "|		  Author:	Marian Jandel					|" <<G4endl;
	G4cout << "|		  		C-INC, LANL					|" <<G4endl;
	G4cout << "|		  comments & questions to: mjandel@lanl.gov			|" <<G4endl;
	G4cout << "-----------------------------------------------------------------------------" <<G4endl;
	G4cout << G4endl;
	G4cout << G4endl;

	//G4cout << "Segmented Mode is set to discriminate Crystal Mult=1 !!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;

	
// define the global instance of DANCEGlobal

DANCEGlobal* GlobalParams=DANCEGlobal::Instance();


if (argc==1)
{

	GlobalParams->MacroInputFile="vis.mac";
	GlobalParams->MasterInputFile="MasterInput.txt";

	G4cout << " *** An interactive session selected *** " << G4endl;
	G4cout << G4endl;

}  
else   // Define UI session for interactive mode.
{
	G4cout << "*** Batch mode selected ***" << G4endl;
	G4cout << G4endl;

	GlobalParams->MasterInputFile=argv[1];
	GlobalParams->MacroInputFile="To be found in MasterInputFile";
}  

G4cout << "Default macro file expected: " << GlobalParams->MacroInputFile << G4endl;
G4cout << "Default MasterInput file expected: " << GlobalParams->MasterInputFile << G4endl;
  
GlobalParams->Init();


NNNearestNeighbors *DANCENeiObj=NNNearestNeighbors::Instance();
DANCENeiObj->NNInitialize("config/DetectorMatrix.txt");

  // Construct the default run manager

 G4RunManager* runManager = new G4RunManager;


 // set mandatory initialization classes
 runManager->SetUserInitialization(new DANCEDetectorConstruction);
 runManager->SetUserInitialization(new DANCEPhysicsList);

 // Initialize G4 kernel
 runManager->Initialize();


 // set user action classes

 runManager->SetUserAction(new DANCEPrimaryGeneratorAction);
 runManager->SetUserAction(new DANCEEventAction);



#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

 // get the pointer to the User Interface manager 
 G4UImanager* UI = G4UImanager::GetUIpointer();  

// here we set if we want to read from external input or not

 if(GlobalParams->ExternalInput){
	while(GlobalParams->EndOfExternalInput==false){
 		if(GlobalParams->EndOfExternalInput==false) runManager->BeamOn(1);
	}
//	runManager->BeamOn(1); 
}

 else{
 if (argc==1)   // Define UI session for interactive mode.
 {
   G4UIsession* session=0;
  
   	// G4UIterminal is a (dumb) terminal.
	#ifdef G4UI_USE_TCSH
  	 session = new G4UIterminal(new G4UItcsh);      
	#else
   	session = new G4UIterminal();
	#endif    
      
   	UI->ApplyCommand("/control/execute vis.mac");    
   	session->SessionStart();

   delete session;
 }
 else           // Batch mode
 { 
   G4String command = "/control/execute ";
//   G4String fileName = argv[1];
   G4String fileName = GlobalParams->MacroInputFile;
   UI->ApplyCommand(command+fileName);
 }
 }

// GlobalParams->EndOfRun();
//  runManager->EndOfRun();
 
 // job termination
#ifdef G4VIS_USE
 delete visManager;
#endif
 delete runManager;

 return 0;
}


