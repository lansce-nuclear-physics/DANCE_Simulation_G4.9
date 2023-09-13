
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


//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: DANCEDetector.cc,v 1.6 2003/08/20 16:32:50 duns Exp $
// --------------------------------------------------------------
//
#include "DANCEDetector.hh"
#include "DANCEDetectorHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"


#include "DANCEGlobal.hh"

DANCEDetector::DANCEDetector(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="hits");
  HCID = -1;
}

DANCEDetector::~DANCEDetector(){;}

void DANCEDetector::Initialize(G4HCofThisEvent*HCE)
{
  hitsCollection = new DANCEDetectorHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);
 
  
}


G4bool DANCEDetector::ProcessHits(G4Step*aStep,G4TouchableHistory* /*ROhist*/)
{

/*
  G4Track *aTrack=aStep->GetTrack();
  int track_ID=aTrack->GetTrackID();
  G4ParticleDefinition *particle_def_ID=aTrack->GetDefinition();
  
  bool Is_Stable=particle_def_ID->GetPDGStable();
  bool Is_ShortLived=particle_def_ID->IsShortLived();
  
  G4String particle_Name=particle_def_ID->GetParticleName();
  G4String cre_process="Target";
  
  if(track_ID>1){  cre_process=aTrack->GetCreatorProcess()->GetProcessName();}
  
  G4cout << track_ID << "	" << particle_Name << G4endl << "Is_Stable=" << Is_Stable << 
  "	Is_ShortLived=" << Is_ShortLived<< "	LifeTime=" << particle_def_ID->GetPDGLifeTime() << G4endl ;
  
  G4cout << track_ID << "	process: " << cre_process << G4endl;


  // here we fill the gpp primary event
  
  DANCEGlobal* GPP=DANCEGlobal::Instance();
  G4int new_prim_particle=1;

  G4double edep = aStep->GetTotalEnergyDeposit();
  
  if(cre_process=="RadioactiveDecay" && particle_Name=="gamma"){
  	
	for(int i=0;i<GPP->Prim_Mult;i++){
		if(track_ID==GPP->PrimaryTrackID[i]) new_prim_particle=0;
	}
	
	if(new_prim_particle==1) {
	
		G4cout << "-----------------------------------------" << G4endl;
	
		GPP->PrimaryTrackID[GPP->Prim_Mult]=track_ID;
		GPP->PrimaryEnergy[GPP->Prim_Mult]=aTrack->GetDynamicParticle()->GetTotalEnergy();

		G4cout << "Adding particle No: " << GPP->Prim_Mult << " to primaries" << G4endl;
		G4cout << "With kinetic energy: " << GPP->PrimaryEnergy[GPP->Prim_Mult] <<
		G4cout << "With deposited energy: " << edep <<	G4endl;
		G4cout << "making the total energy: " << edep + GPP->PrimaryEnergy[GPP->Prim_Mult] <<	G4endl;

		GPP->Prim_Mult++;

		G4cout << "-----------------------------------------" << G4endl;
	
	}
  	
  }
 */

  G4double edep = aStep->GetTotalEnergyDeposit();
  
  if(edep<0 || edep>1e10) G4cout << "This edep is out of bounds :  " << edep << "	it sucks !!!" <<
  G4endl;  
  if(edep==0.) return true;


  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  G4int copyNo = theTouchable->GetVolume()->GetCopyNo();
  G4double hitTime = preStepPoint->GetGlobalTime();
 
  // check if this finger already has a hit
  G4int ix = -1;
  for(G4int i=0;i<hitsCollection->entries();i++)
  {
    if((*hitsCollection)[i]->GetID()==copyNo)
    {
      ix = i;
//      G4cout << "Hit in Det No: " << copyNo << G4endl;
      break;
    }
  }
  // if it has, then take the earlier time
  if(ix>=0)
  {
  
 /* 
    if((*hitsCollection)[ix]->GetTime()>hitTime)
    { 
//    	G4cout << "Looks Like Multiple Hits at Time: " << (*hitsCollection)[ix]->GetTime() << G4endl;
    	(*hitsCollection)[ix]->SetTime(hitTime); 
    }
    */

//	G4cout << "Adding edep: " << edep << G4endl;
	
    	(*hitsCollection)[ix]->AddEdep(edep);
//    	G4cout << "Total Edep: " << (*hitsCollection)[ix]->GetEdep() << G4endl;
  }
  else
  // if not, create a new hit and set it to the collection
  {
//    G4cout << "Creating Hit in: " << copyNo << "	 with Edep:" << edep << G4endl;
    
    DANCEDetectorHit* aHit = new DANCEDetectorHit(copyNo,hitTime);
    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
    aHit->SetLogV(thePhysical->GetLogicalVolume());
    G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
    aTrans.Invert();
    aHit->SetRot(aTrans.NetRotation());
    aHit->SetPos(aTrans.NetTranslation());
    hitsCollection->insert( aHit );
    aHit->SetEdep(edep);
      }

  return true;
}

void DANCEDetector::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{;}

