
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
// $Id: DANCEPhysicsList.cc,v 1.7 2004/11/23 04:02:03 tkoi Exp $
// --------------------------------------------------------------
//
// 28-Jan-04 Add QGSP_BERT and QGSP_BIC for hadronic lists. T. Koi
// 22-Nov-04 Comment out QGSP_BERT and QGSP_BIC
//           Output Notificaiton message             
//           All Particles are created in GeneralPhysics 

#include "DANCEPhysicsList.hh"
#include "DANCEPhysListParticles.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>

#include "DANCEGeneralPhysics.hh"
#include "DANCEEMPhysics.hh"
//#include "DANCEMuonPhysics.hh"
#include "DANCEHadronPhysics.hh"
//#include "DANCEIonPhysics.hh"

//#include "HadronPhysicsQGSP_BERT.hh"
//#include "HadronPhysicsQGSP_BIC.hh"


DANCEPhysicsList::DANCEPhysicsList():  G4VModularPhysicsList()
{


  // default cut value  (1.0mm)
//  defaultCutValue = 0.005*mm;
  defaultCutValue = 0.005*mm;
  SetVerboseLevel(0);
 
 
  RegisterPhysics(  new DANCEPhysListParticles("particles"));

  // General Physics ( Create ALL Particle and apply Decay )
  RegisterPhysics( new DANCEGeneralPhysics("general") );

  // EM Physics ( Apply related Processes to gamma and e-/+)
  RegisterPhysics( new DANCEEMPhysics("standard EM"));

  // Muon Physics ( Apply related processes to mu and tau
//  RegisterPhysics(  new DANCEMuonPhysics("muon"));

   // Hadron Physics ( Apply related processes to hadrons )
  RegisterPhysics(  new DANCEHadronPhysics("hadron"));
// We do not use hadronic lists since v7.
  //RegisterPhysics(  new HadronPhysicsQGSP_BERT("hadron"));
  //RegisterPhysics(  new HadronPhysicsQGSP_BIC("hadron"));

  // Ion Physics ( Apply related processes to ions )
//  RegisterPhysics( new DANCEIonPhysics("ion"));

}

DANCEPhysicsList::~DANCEPhysicsList()
{
}

void DANCEPhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
//  SetCutsWithDefault();

  if (verboseLevel >1)
    G4cout << "DMXPhysicsList::SetCuts:";
  
  if (verboseLevel>0){
    G4cout << "DMXPhysicsList::SetCuts:";
    G4cout << "CutLength : " 
	   << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  //special for low energy physics
  G4double lowlimit=100*eV;  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,40.*MeV);

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "neutron");
  
  if (verboseLevel>0) DumpCutValuesTable();

}



