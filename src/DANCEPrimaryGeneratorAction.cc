
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



#include "DANCEPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "CLHEP/Random/RandFlat.h"

#include "DANCEGlobal.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"

#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"

using namespace std;

DANCEPrimaryGeneratorAction::DANCEPrimaryGeneratorAction()
{

//  particleGun = new G4ParticleGun ();

  DANCEGlobal* GPP=DANCEGlobal::Instance();
   
  if(GPP->ExternalInput){
 	
	G4String file_name=G4String(GPP->ExternalInputFile);
	
//	G4cout << "Opening File: " << GPP->ExternalInputFile.c_str() << G4endl;
	G4cout << "Opening File: " << file_name.c_str() << G4endl;
  	
//	input_stream= new ifstream(GPP->ExternalInputFile.c_str(),ios::in);
	input_stream= new ifstream(file_name.c_str(),ios::in);
	 
	if(!input_stream->is_open())   {
	
//		G4cout << "File: " << GPP->ExternalInputFile.c_str() << " not found ! Exiting " << G4endl;
		G4cout << "File: " << file_name.c_str() << " not found ! Exiting " << G4endl;
		exit(0);
	}
  } 
 
  else{
  	sourceGun = new G4GeneralParticleSource ();
	//particleGun = new G4ParticleGun ();
  }



}

DANCEPrimaryGeneratorAction::~DANCEPrimaryGeneratorAction()
{
//  delete particleGun;
    DANCEGlobal* GPP=DANCEGlobal::Instance();
    if(GPP->ExternalInput){
	delete input_stream;
    }
	
    else {
	delete sourceGun;
    }
}

void DANCEPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  DANCEGlobal* GPP=DANCEGlobal::Instance();

//  GPP->MultiplicityOfPrimaryEvent=10;

  GPP->Prim_Mult=0;
  GPP->PrimaryTrackID[0]=-1;
  
  if(GPP->ExternalInput){
		
		// read the file
		
		(*input_stream) >> GPP->MultiplicityOfPrimaryEvent;
		//G4cout << "Multiplicity: " << GPP->MultiplicityOfPrimaryEvent << G4endl;
		if(input_stream->eof()) {
			GPP->EndOfExternalInput=true;
			//G4cout << " End of file reached in Generate Primaries" << G4endl;
			
			return;
		}
		
		else{

//		G4PrimaryVertex* g4vtx=new G4PrimaryVertex(4.*cm,0.,0.,0.);   // Displacement along the beam line
		G4PrimaryVertex* g4vtx=new G4PrimaryVertex(0.,0.,0.,0.);
		for(int i=0;i<GPP->MultiplicityOfPrimaryEvent;i++){
	

			// isotropic distribution
			
    			G4double phi = CLHEP::RandFlat::shoot(0.*deg,360.*deg);
    			G4double costheta = CLHEP::RandFlat::shoot(-1.,1.);
    			G4double sintheta = std::sqrt(1.-costheta*costheta);
    						
			
			G4double randomX = sintheta*std::cos(phi);
    			G4double randomY = sintheta*std::sin(phi);
    			G4double randomZ = costheta;
		
			G4ThreeVector v(randomX,randomY,randomZ);
    		
			G4double Particle_Energy;
		
			(*input_stream) >> Particle_Energy;
		
			//G4cout << Particle_Energy << "	";
			
			G4PrimaryParticle* g4prim= new G4PrimaryParticle(i,Particle_Energy*randomX*MeV,Particle_Energy*randomY*MeV,Particle_Energy*randomZ*MeV);

			G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  		G4String particleName;
	  		g4prim->SetG4code(particleTable->FindParticle(particleName="gamma"));

			g4vtx->SetPrimary(g4prim);		
		
		}
		anEvent->AddPrimaryVertex(g4vtx);
		//G4cout << G4endl;
		}
  }
  	
  else{
  	sourceGun->GeneratePrimaryVertex(anEvent);
  	//particleGun->GeneratePrimaryVertex(anEvent);

  }
}


