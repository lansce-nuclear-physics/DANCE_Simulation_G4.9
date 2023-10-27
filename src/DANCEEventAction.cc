
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



// DANCEEventAction.cc   -  M. Jandel - 2006 - May 9
//  
// Everything important takes place here
// Events are being processed and corresponding output is dumped
// available outputs are being set in Master_Input.txt file
// Options include: root histograms and binary files
// Guidelines for adding new root histogram:
// 	1. Define the pointer in DANCEEventAction.hh
// 	2. Declare the histogram at the end of DANCEEventAction::DANCEEventAction()
// 	3. Fill it appropriately in DANCEEventAction::EndOfEventAction()
// 	4. Dont forget to add it to a list of histograms for writing in DANCEEventAction::EndOfEventAction()
//		


#include "DANCEEventAction.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4THitsMap.hh"

#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VHitsCollection.hh"

#include "G4RunManager.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TMath.h"
#include "TRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "DANCEDetectorHit.hh"

#include "DANCEGlobal.hh"
#include "DANCENNNearestNeighbors.hh"

#ifdef G4ANALYSIS_USE
#include "DANCEAnalysisManager.hh"
#endif // G4ANALYSIS_USE


extern int A_counter,B_counter,C_counter,D_counter;
extern int AutoSave;

using namespace std;

DANCEEventAction::DANCEEventAction()
{
 
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  DANCEGlobal* GPP=DANCEGlobal::Instance();
  
	G4String Num;
   	numberOfEvent=0;
	gatedEvents=0;
		
	G4int leadcounter=0;
	
  for(size_t i=0;i<A_counter;i++)
  {
	   ostringstream oss;
	   oss<<i;
	   Num=oss.str();

	   fullName[leadcounter] = "BaF_A_"+Num+"/hits";
	   DHID[leadcounter] = SDman->GetCollectionID(fullName[leadcounter]);
//	   cout << "Pointer:   " << DHID[leadcounter] << endl;
	   leadcounter++;
      
  };	

  for(size_t i=0;i<B_counter;i++)
  {
	   ostringstream oss;
	   oss<<i;
	   Num=oss.str();
      	   fullName[leadcounter] = "BaF_B_"+Num+"/hits";
	   DHID[leadcounter] = SDman->GetCollectionID(fullName[leadcounter]);
//	   cout << "Pointer:   " << DHID[leadcounter] << endl;
      	   leadcounter++;

  };	

  for(size_t i=0;i<C_counter;i++)
  {
	   ostringstream oss;
	   oss<<i;
	   Num=oss.str();

           fullName[leadcounter] = "BaF_C_"+Num+"/hits";
	   DHID[leadcounter] = SDman->GetCollectionID(fullName[leadcounter]);
//	   cout << "Pointer:   " << DHID[leadcounter] << endl;
          leadcounter++;
  };
  
    for(size_t i=0;i<D_counter;i++)
  {
	   ostringstream oss;
	   oss<<i;
	   Num=oss.str();

     	   fullName[leadcounter] = "BaF_D_"+Num+"/hits";
           DHID[leadcounter] = SDman->GetCollectionID(fullName[leadcounter]);
//	   cout << "Pointer:   " << DHID[leadcounter] << endl;
           leadcounter++;
  };	
  
//  for(size_t i=0;i<leadcounter;i++){  
//  	G4cout <<  i << "	" << fullName[i] << G4endl; 
//  }
  	  

//G4cout << "A_counter = " << A_counter << G4endl;
//G4cout << "B_counter = " << B_counter << G4endl;
//G4cout << "C_counter = " << C_counter << G4endl;
//G4cout << "D_counter = " << D_counter << G4endl;


#ifdef G4ANALYSIS_USE
  plotter = 0;
  dc1Hits = dc2Hits = multHits = IDdet = 0;
 // dc1XY = dc2XY = evstof = 0;

  // Do some analysis

  DANCEAnalysisManager* analysisManager = DANCEAnalysisManager::getInstance();
  IHistogramFactory* hFactory = analysisManager->getHistogramFactory();

  if (hFactory)
  {
    // Create some AIDA histograms
    
    	dc1Hits = hFactory->createHistogram1D("Gamma Ray Energy [MeV]",1000,0,15);
	dc2Hits = hFactory->createHistogram1D("Total Gamma Ray Energy [MeV]",1000,0,15);

	multHits = hFactory->createHistogram1D("Crystal Multiplicity",10,0,10);
	IDdet = hFactory->createHistogram1D("Time spectrum",1000,0,1000);


    plotter = analysisManager->createPlotter();
    if (plotter)
    {
       plotter->createRegions(2,2);
       plotter->region(0)->plot(*dc1Hits);
       plotter->region(1)->plot(*dc2Hits);
       plotter->region(2)->plot(*multHits);
       plotter->region(3)->plot(*IDdet);

       plotter->show();
     }
  }
#endif // G4ANALYSIS_USE

//---------------------------------------------------------
// define ROOT HISTOGRAMS

// basic histograms
f1=new TH1F("Etot","Etot",1000,0,16);
f2=new TH1F("EgammaCrystal","EgammaCrystal",1000,0,16);

f3=new TH1F("MultCrystal","MultCrystal",170,-0.5,169.5);

f4=new TH1F("EgammaCluster","EgammaCluster",1000,0,16);

f5=new TH1F("TimeSpectrumCrystal","TimeSpectrumCrystal",1000,0,200);
f6=new TH1F("MultCluster","MultCluster",170,-0.5,169.5);

f7=new TH1F("ThetaCluster","ThetaCluster",360,0.,180.);
f8=new TH1F("PhiCluster","PhiCluster",720,0.,360.);


f9=new TH1F("ThetaCrystal","ThetaCrystal",360,0.,180.);
f10=new TH1F("PhiCrystal","PhiCrystal",720,0.,360.);

f11=new TH1F("TimeSpectrumCluster","TimeSpectrumCluster",1000,0,200);
f12=new TH1F("Crystal_ID","Crystal_ID",170,-0.5,169.5);

f13=new TH1F("PrimaryEgamma","PrimaryEgamma",1000,0,16);

f14=new TH2F("Mult_PrimVsCrystal","Mult_PrimVsCrystal",50,-0.5,49.5,50,-0.5,49.5);
f15=new TH2F("Mult_PrimVsCluster","Mult_PrimVsCluster",50,-0.5,49.5,50,-0.5,49.5);

f16=new TH1F("PrimaryMult","PrimaryMult",170,-0.5,169.5);

f17=new TH1F("Scaler","Scaler",10,0,10);

ECr_ID = new TH2F( "ECr_ID", "ECr_ID", 2000, 0, 20, 170, -0.5, 169.5);
/*
f20=new TH1F("EgammaCluster_M2","EgammaCluster_M2",256,0,16);
f21=new TH1F("EgammaCluster_M3","EgammaCluster_M3",256,0,16);
f22=new TH1F("EgammaCluster_M4","EgammaCluster_M4",256,0,16);
f23=new TH1F("EgammaCluster_M5","EgammaCluster_M5",256,0,16);
f24=new TH1F("EgammaCluster_M6","EgammaCluster_M6",256,0,16);

f29=new TH1F("EgammaCluster_M1_Qgated","EgammaCluster_M1_Qgated",256,0,16);
f30=new TH1F("EgammaCluster_M2_Qgated","EgammaCluster_M2_Qgated",256,0,16);
f31=new TH1F("EgammaCluster_M3_Qgated","EgammaCluster_M3_Qgated",256,0,16);
f32=new TH1F("EgammaCluster_M4_Qgated","EgammaCluster_M4_Qgated",256,0,16);
f33=new TH1F("EgammaCluster_M5_Qgated","EgammaCluster_M5_Qgated",256,0,16);
f34=new TH1F("EgammaCluster_M6_Qgated","EgammaCluster_M6_Qgated",256,0,16);

//Etot gated on cluster multiplicity
if(GPP->ClusterMGatedHists){
	f1_cl1=new TH1F("Etot_clM1","Etot_clM1",512,0,16);
	f1_cl2=new TH1F("Etot_clM2","Etot_clM2",512,0,16);
	f1_cl3=new TH1F("Etot_clM3","Etot_clM3",512,0,16);
	f1_cl4=new TH1F("Etot_clM4","Etot_clM4",512,0,16);
	f1_cl5=new TH1F("Etot_clM5","Etot_clM5",512,0,16);
	f1_cl6=new TH1F("Etot_clM6","Etot_clM6",512,0,16);
	f1_cl7=new TH1F("Etot_clM7","Etot_clM7",512,0,16);
	f1_cl8=new TH1F("Etot_clM8","Etot_clM8",512,0,16);
}

//Etot gated on cluster multiplicity
if(GPP->CrystalMGatedHists){
	f1_cr1=new TH1F("Etot_crM1","Etot_crM1",512,0,16);
	f1_cr2=new TH1F("Etot_crM2","Etot_crM2",512,0,16);
	f1_cr3=new TH1F("Etot_crM3","Etot_crM3",512,0,16);
	f1_cr4=new TH1F("Etot_crM4","Etot_crM4",512,0,16);
	f1_cr5=new TH1F("Etot_crM5","Etot_crM5",512,0,16);
	f1_cr6=new TH1F("Etot_crM6","Etot_crM6",512,0,16);
	f1_cr7=new TH1F("Etot_crM7","Etot_crM7",512,0,16);
	f1_cr8=new TH1F("Etot_crM8","Etot_crM8",512,0,16);
}
*/

///////////////////////////////////////////////////////////////
// these histograms match geant3 output

f20=new TH1F("EgammaCluster_M2","EgammaCluster_M2",1500,0,15);
f21=new TH1F("EgammaCluster_M3","EgammaCluster_M3",1500,0,15);
f22=new TH1F("EgammaCluster_M4","EgammaCluster_M4",1500,0,15);
f23=new TH1F("EgammaCluster_M5","EgammaCluster_M5",1500,0,15);
f24=new TH1F("EgammaCluster_M6","EgammaCluster_M6",1500,0,15);

f29=new TH1F("EgammaCluster_M1_Qgated","EgammaCluster_M1_Qgated",1500,0,15);
f30=new TH1F("EgammaCluster_M2_Qgated","EgammaCluster_M2_Qgated",1500,0,15);
f31=new TH1F("EgammaCluster_M3_Qgated","EgammaCluster_M3_Qgated",1500,0,15);
f32=new TH1F("EgammaCluster_M4_Qgated","EgammaCluster_M4_Qgated",1500,0,15);
f33=new TH1F("EgammaCluster_M5_Qgated","EgammaCluster_M5_Qgated",1500,0,15);
f34=new TH1F("EgammaCluster_M6_Qgated","EgammaCluster_M6_Qgated",1500,0,15);

//Etot gated on cluster multiplicity
if(GPP->ClusterMGatedHists){
	f1_cl1=new TH1F("Etot_clM1","Etot_clM1",1500,0,15);
	f1_cl2=new TH1F("Etot_clM2","Etot_clM2",1500,0,15);
	f1_cl3=new TH1F("Etot_clM3","Etot_clM3",1500,0,15);
	f1_cl4=new TH1F("Etot_clM4","Etot_clM4",1500,0,15);
	f1_cl5=new TH1F("Etot_clM5","Etot_clM5",1500,0,15);
	f1_cl6=new TH1F("Etot_clM6","Etot_clM6",1500,0,15);
	f1_cl7=new TH1F("Etot_clM7","Etot_clM7",1500,0,15);
	f1_cl8=new TH1F("Etot_clM8","Etot_clM8",1500,0,15);
}

//Etot gated on cluster multiplicity
if(GPP->CrystalMGatedHists){
	f1_cr1=new TH1F("Etot_crM1","Etot_crM1",1500,0,15);
	f1_cr2=new TH1F("Etot_crM2","Etot_crM2",1500,0,15);
	f1_cr3=new TH1F("Etot_crM3","Etot_crM3",1500,0,15);
	f1_cr4=new TH1F("Etot_crM4","Etot_crM4",1500,0,15);
	f1_cr5=new TH1F("Etot_crM5","Etot_crM5",1500,0,15);
	f1_cr6=new TH1F("Etot_crM6","Etot_crM6",1500,0,15);
	f1_cr7=new TH1F("Etot_crM7","Etot_crM7",1500,0,15);
	f1_cr8=new TH1F("Etot_crM8","Etot_crM8",1500,0,15);
}

////////////////////////////////////////////////////////////

// create spectra gated on crystal energy
if(GPP->SingleGateECrystal_E[0]>0){
	ostringstream oss;
	oss<<"GateOnE=" << GPP->SingleGateECrystal_E[0] << "MeV_dE=" << GPP->SingleGateECrystal_dE[0] <<"MeV";
	Num=oss.str();
 	f2_crgate1=new TH1F("EgammaCr_gate1",Num,1000,0,16);
}
if(GPP->SingleGateECrystal_E[1]>0){
	ostringstream oss;
	oss<<"GateOnE=" << GPP->SingleGateECrystal_E[1] << "MeV_dE=" << GPP->SingleGateECrystal_dE[1] <<"MeV";
	Num=oss.str();
	f2_crgate2=new TH1F("EgammaCr_gate2",Num,1000,0,16);

}
if(GPP->SingleGateECrystal_E[2]>0){
	ostringstream oss;
	oss<<"GateOnE=" << GPP->SingleGateECrystal_E[2] << "MeV_dE=" << GPP->SingleGateECrystal_dE[2] <<"MeV";
	Num=oss.str();
	f2_crgate3=new TH1F("EgammaCr_gate3",Num,1000,0,16);
}

// create spectra gated on cluster energy
if(GPP->SingleGateECluster_E[0]>0){
	ostringstream oss;
	oss<<"GateOnE=" << GPP->SingleGateECluster_E[0] << "MeV_dE=" << GPP->SingleGateECluster_dE[0] <<"MeV";
	Num=oss.str();
 	f2_clgate1=new TH1F("EgammaCl_gate1",Num,1000,0,16);
}
if(GPP->SingleGateECluster_E[1]>0){
	ostringstream oss;
	oss<<"GateOnE=" << GPP->SingleGateECluster_E[1] << "MeV_dE=" << GPP->SingleGateECluster_dE[1] <<"MeV";
	Num=oss.str();
	f2_clgate2=new TH1F("EgammaCl_gate2",Num,1000,0,16);

}
if(GPP->SingleGateECluster_E[2]>0){
	ostringstream oss;
	oss<<"GateOnE=" << GPP->SingleGateECluster_E[2] << "MeV_dE=" << GPP->SingleGateECluster_dE[2] <<"MeV";
	Num=oss.str();
	f2_clgate3=new TH1F("EgammaCl_gate3",Num,1000,0,16);
}

// --------- Begin Components for simulation binaries -------------------

  //Make binary output structure
  std::stringstream outfilename;
  outfilename.str();
  outfilename << "Test2.bin";
 
  outputbinfile.open(outfilename.str().c_str(), std::ios::out | std::ios::binary);

  global_timestamp=6000;
  dance_T0_counter=0;

// --------- End Components for simulation binaries ---------------------

}

DANCEEventAction::~DANCEEventAction()
{

// Save ROOT Histograms
  DANCEGlobal* GPP=DANCEGlobal::Instance();
  if(GPP->RootFile) RootAutoSave();

#ifdef G4ANALYSIS_USE
  DANCEAnalysisManager::dispose();
#endif // G4ANALYSIS_USE


// --------- Begin Components for simulation binaries -------------------
  outputbinfile.close();
// --------- End Components for simulation binaries ---------------------

}

void DANCEEventAction::BeginOfEventAction(const G4Event* evt)
{
/*
  DANCEGlobal* GPP=DANCEGlobal::Instance();

  if(GPP->EndOfExternalInput){

 	G4RunManager* runManager = G4RunManager::GetRunManager();
 
  	runManager->AbortRun(true);
//  	SetEventAborted();
	
	return;
  }
*/
}
void DANCEEventAction::EndOfEventAction(const G4Event* evt)
{
  
  DANCEGlobal* GPP=DANCEGlobal::Instance();
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();

  NNNearestNeighbors *NNP=NNNearestNeighbors::Instance();
  
  if(GPP->ExternalInput){
  	if(GPP->EndOfExternalInput) {
		
		G4cout << " End of file reached " << G4endl;
  		if(GPP->RootFile) RootAutoSave();
		
		return;
		
		}
  } 

  if(!HCE) return;

//  G4cout << "Primary particles --------------------------------------------------------------" << G4endl;

  G4int n_vertex = evt->GetNumberOfPrimaryVertex();

//  G4cout << "No of vertices: " << n_vertex << G4endl;

  G4int Prim_Mult=0;
  G4double EkinTot;
 
  if(GPP->ExternalInput){
  
  for(G4int iv=0;iv<n_vertex;iv++)
  {
    	G4PrimaryVertex* pv = evt->GetPrimaryVertex(iv);
    	G4PrimaryParticle* pp = pv->GetPrimary();
 	G4ThreeVector M=pp->GetMomentum();
	G4double Mass=pp->GetMass();
	
//	EkinTot=(sqrt(M.x()*M.x()+M.y()*M.y()+M.z()*M.z()+Mass*Mass)-Mass)/MeV;

// for gammas
	EkinTot=sqrt(M.x()*M.x()+M.y()*M.y()+M.z()*M.z())/MeV;

//	G4cout << "Energy: " << EkinTot << G4endl;
//	G4cout << "Energy: " << EkinTot << "	Mass: " << Mass << G4endl;


/*
	if(pp->GetPDGcode()==0){
		G4PrimaryParticle* daughter=pp->GetDaughter();	
		while(daughter){
			PrintPrimary(daughter,0);
			daughter->GetNext();
		}
	}
*/
	
	f13->Fill(EkinTot);
	Prim_Mult++;

   	while(pp)
    	{
      		//PrintPrimary(pp,0);
        	pp = pp->GetNext();
		if(pp){
 		M=pp->GetMomentum();
		Mass=pp->GetMass();
		EkinTot=sqrt(M.x()*M.x()+M.y()*M.y()+M.z()*M.z())/MeV;
//		G4cout << "Energy: " << EkinTot << "	Mass: " << Mass << G4endl;
		f13->Fill(EkinTot);
		Prim_Mult++;
	}
    }
  }  
  }
  
  else{
  
    for(G4int iv=0;iv<n_vertex;iv++)
  	{
    	G4PrimaryVertex* pv = evt->GetPrimaryVertex(iv);
    	G4PrimaryParticle* pp = pv->GetPrimary();
 	G4ThreeVector M=pp->GetMomentum();
	G4double Mass=pp->GetMass();
	
//	EkinTot=(sqrt(M.x()*M.x()+M.y()*M.y()+M.z()*M.z()+Mass*Mass)-Mass)/MeV;

// for gammas
	EkinTot=sqrt(M.x()*M.x()+M.y()*M.y()+M.z()*M.z())/MeV;

//	G4cout << "Energy: " << EkinTot << G4endl;
//	G4cout << "Energy: " << EkinTot << "	Mass: " << Mass << G4endl;


/*
	if(pp->GetPDGcode()==0){
		G4PrimaryParticle* daughter=pp->GetDaughter();	
		while(daughter){
			PrintPrimary(daughter,0);
			daughter->GetNext();
		}
	}
*/
	
	f13->Fill(EkinTot);
	Prim_Mult++;

   		while(pp)
    		{
      			//PrintPrimary(pp,0);
        		pp = pp->GetNext();
			if(pp){
 				M=pp->GetMomentum();
				Mass=pp->GetMass();
				EkinTot=sqrt(M.x()*M.x()+M.y()*M.y()+M.z()*M.z())/MeV;
//				G4cout << "Energy: " << EkinTot << "	Mass: " << Mass << G4endl;
				f13->Fill(EkinTot);
				Prim_Mult++;
			}
    		}
  	}  
  }
  
  numberOfEvent++;
  f17->Fill(0);
  
  double E_ball=0.;
  double min_time=1e40;
  
  double E_crystal_energy[160];
  double E_crystal_time[160];
  double E_cluster_energy[160];
  
  int E_crystal_ID[160];
  int Crystal_Mult=0;
  int Cluster_Mult=0;
  int counter=0;
 
  int CrystalGate_Flag[30];
  int CrystalGate_ID[30];
  
  int ClusterGate_Flag[30];
  int ClusterGate_ID[30];
  
  //cout << "min_time: " << min_time << endl;
  
//  G4cout <<"EventAction ----------------------------------------------------------------------------"<< G4endl;


// if we have a correct event retrieve the data
  double totE_noSmearing_thr=0.;
  double noSmearing=0;
  double totE_noSmearing=0.;

  if(HCE){
  
  for(int i=0;i<3;i++){
 	CrystalGate_Flag[i]=0;
 	ClusterGate_Flag[i]=0;
   }

  for(int i=0;i<160;i++){
  	DANCEDetectorHitsCollection* DHC=0;  
  	
	DHC=(DANCEDetectorHitsCollection*) (HCE->GetHC(DHID[i]));

	double time;
  	double E_crystal_tot=0.;
	
  	// if there is a Collection
	if(DHC!=0){
		int n_hit=DHC->entries();
		//double totE_noSmearing=0.;
		if(n_hit>0) 
		{
//			cout << "Det No:" << i << "	Hits:" << n_hit << endl;
			
			time=0;;
			
			for(int j=0;j<n_hit;j++){
			
				DANCEDetectorHit* Hit=(*DHC)[j];
				double eDep=Hit->GetEdep();
				time=time + Hit->GetTime();
				
				E_crystal_tot+=eDep;
				
//				#ifdef G4ANALYSIS_USE
//					if(time>=0 && time<1000. ) IDdet->fill(time);
//				#endif
				
//				G4cout << "Energy Deposited: " << eDep << "	at time: " << time <<  G4endl;
	
			}
			
		if(time<min_time) min_time=time;
			
		//    G4double E_smeared = CLHEP::RandFlat::shoot(Energy,0.1*Energy/2.35);

		// Resolution function
		// old approximate linear one
		//G4double sig=0.030915+(0.031302*E_crystal_tot);

		// new de/E=a1+a2/sqrt(E) ---> sigma=1/2.35 * E*a1+a2*sqrt(E)
		
		//G4double sig=1./2.354820045*(1.089e-2*E_crystal_tot+(0.136419*sqrt(E_crystal_tot)));
		G4double sig=1./2.354820045*(1.089e-2*E_crystal_tot+(0.146*sqrt(E_crystal_tot)));   //  The resolution for BaF2
//		G4double sig=1./2.354820045*(0.0*E_crystal_tot+(0.024128179*sqrt(E_crystal_tot)));  // The resolution for BrilLanCe
//		FWHM = 2.35*sig
//		FWHM/E = 76.3/sqrt(1000*E)

		
		// smearing function
	        G4double E_smeared = CLHEP::RandGauss::shoot(E_crystal_tot,sig);
//	        G4double E_smeared = E_crystal_tot;
//		double threshold= CLHEP::RandFlat::shoot(GPP->E_threshold,0.25);

		double threshold= GPP->E_threshold; //CLHEP::RandFlat::shoot(GPP->E_threshold,0.25);
		  
		 if(E_crystal_tot>=threshold) f17->Fill(1);
		 if(E_crystal_tot<threshold) f17->Fill(2);
		 if(fabs(E_crystal_tot-EkinTot)<0.05) f17->Fill(3);

		noSmearing=E_crystal_tot;
		totE_noSmearing=totE_noSmearing+noSmearing;
		
		
		E_crystal_tot=E_smeared;

		// here we have crystal energy and we are smearing it using above formula
		
		if(E_crystal_tot>=threshold){
			
			totE_noSmearing_thr=totE_noSmearing_thr+noSmearing;
			
			E_crystal_energy[counter]=E_crystal_tot;
			E_crystal_time[counter]=time;
		
			// first get DetType and DetNo
		
			int DetType=GPP->GetGeantID2D_DetType(i);
			int DetNo=GPP->GetGeantID2D_DetNo(i);

			E_crystal_ID[counter]=GPP->GetID2D(DetType,DetNo);

			// now crystal gate flags
			for(int gatei=0;gatei<3;gatei++){
						
				if(GPP->SingleGateECrystal_E[gatei]>0){
				
					if(E_crystal_tot>=(GPP->SingleGateECrystal_E[gatei]-GPP->SingleGateECrystal_dE[gatei]) &&
					E_crystal_tot<=(GPP->SingleGateECrystal_E[gatei]+GPP->SingleGateECrystal_dE[gatei])){
					
						CrystalGate_Flag[gatei]=100;
						CrystalGate_ID[gatei]=counter;
					
					}
				
				}
			
			}
		
		
			counter++;
			
		}
		
		
		}
	}
  }

// check if multiplicity is 0 if yes return


  if(totE_noSmearing_thr>0) f17->Fill(4);
  if(fabs(totE_noSmearing_thr-EkinTot)<0.05) f17->Fill(5);		

  if(totE_noSmearing>0) f17->Fill(6);
  if(fabs(totE_noSmearing-EkinTot)<0.05) f17->Fill(7);		

  if(counter==0) return;
  
  
// --------- Begin Components for simulation binaries -------------------
// We have an event with M>0, so we check it we need to increment T0s for SimBin
    //Make the timestamp for this DANCE event
    if(fSetTau) 
    {
      devent_tau = fSetTau;
    }
    else 
    {
    	// std::cout<<"Using the default tau value from 60Co source run"<<std::endl;
    	//devent_tau = 5818.4;
    	std::cout<<"Using a nominal tau value of 500 us, which should mean no pileup"<<std::endl;
    	devent_tau = 5.e5;
    	fSetTau = devent_tau;
    }
    
    global_timestamp += gRandom->Exp(devent_tau);

    //Put a T0 in the data stream every 50 ms
    if(global_timestamp > dance_T0_counter*50000000.0 + 6000.0) 
    {
      //Make a T0 event
      devt_out.timestamp = dance_T0_counter*50000000.0 + 6000.0; //6000 is an arbitrary offset
      devt_out.wfintegral = 1.;  //  This is currently arbitrary.  Need to check
      devt_out.Ifast = 1000;
      devt_out.Islow = 1000;
      devt_out.ID = 200;
      outputbinfile.write(reinterpret_cast<char*>(&devt_out),sizeof(DEVT_STAGE1_WF));
      
      //T0 is 50 ms apart 
      dance_T0_counter++;
      if ( ( dance_T0_counter % 100 ) == 0) 
         std::cout<<  " ++++ Writing T0 ++++ for T0 " << dance_T0_counter << std::endl;
    }

// --------- End Components for simulation binaries ---------------------

  Crystal_Mult=counter;


  // for segmented mode we require two crystals to fire up
  //if(Crystal_Mult<2) return;

  f3->Fill(Crystal_Mult);


// primary multiplicieties filled
// if there is no hit in DANCE primaries will not fill in as this method is exited in that case in the beginning

  f14->Fill(Prim_Mult,Crystal_Mult);
  f16->Fill(Prim_Mult);
  
// loop through crystals and fill crystal related histograms
// G4cout << "loop through crystals and fill crystal related histograms" << G4endl;

  int Cluster_ID[100];
  
  for(int i=0;i<counter;i++)
  {
    Cluster_ID[i]=i+1; // this is for preparation for clusterization routine
		
    E_crystal_time[i]=E_crystal_time[i]-min_time;
		
    #ifdef G4ANALYSIS_USE
      if(E_crystal_time[i]>=0 && E_crystal_time[i]<1000.) IDdet->fill(E_crystal_time[i]);
      if(E_crystal_energy[i]<15. && E_crystal_time[i]<200) dc1Hits->fill(E_crystal_energy[i]);
    #endif
    f5->Fill(E_crystal_time[i]);
    f2->Fill(E_crystal_energy[i]);
    //f3->Fill(E_crystal_ID[i]);

    ECr_ID->Fill(E_crystal_energy[i], E_crystal_ID[i]);
			
    f9->Fill(GPP->GetTheta(E_crystal_ID[i]));
    f10->Fill(GPP->GetPhi(E_crystal_ID[i]));
    f12->Fill(E_crystal_ID[i]);

    E_ball+=E_crystal_energy[i];


// --------- Begin Components for simulation binaries -------------------

    devt_out.timestamp = global_timestamp+gRandom->Gaus(0,1.5);   //account for the width of a dance devent 
    devt_out.wfintegral = 0.12;
    devt_out.Ifast = E_crystal_energy[i]*1000.0;
    devt_out.Islow = E_crystal_energy[i]*1000.0;
    devt_out.ID = E_crystal_ID[i];
    outputbinfile.write(reinterpret_cast<char*>(&devt_out),sizeof(DEVT_STAGE1));
 

// --------- End Components for simulation binaries ---------------------

  }
  	if(CrystalGate_Flag[0]==100){
		f2_crgate1->Fill(E_ball-E_crystal_energy[CrystalGate_ID[0]]);
  	}
	if(CrystalGate_Flag[1]==100){
		f2_crgate2->Fill(E_ball-E_crystal_energy[CrystalGate_ID[1]]);
	}
	if(CrystalGate_Flag[2]==100){
		f2_crgate3->Fill(E_ball-E_crystal_energy[CrystalGate_ID[2]]);
	}

	if(GPP->CrystalMGatedHists){
						
		if(Crystal_Mult==1)	f1_cr1->Fill(E_ball);
		if(Crystal_Mult==2)	f1_cr2->Fill(E_ball);
		if(Crystal_Mult==3)	f1_cr3->Fill(E_ball);
		if(Crystal_Mult==4)	f1_cr4->Fill(E_ball);
		if(Crystal_Mult==5)	f1_cr5->Fill(E_ball);
		if(Crystal_Mult==6)	f1_cr6->Fill(E_ball);
		if(Crystal_Mult==7)	f1_cr7->Fill(E_ball);
		if(Crystal_Mult>=8) 	f1_cr8->Fill(E_ball);	
	}
  
  #ifdef G4ANALYSIS_USE
	if(E_ball>0 && E_ball<15.) dc2Hits->fill(E_ball);
  #endif
	f1->Fill(E_ball);

// --------------------------------------------------------------------------------------------------------
// here we will clusterize
// result will be in Cluster_ID[i] where i runs through all the crystals
// minimum of Cluster_ID is 1 so don't forget to lower it by 1
// it is somewhat different then Jan Wouters routine but should lead to same results
// if used stand alone, set Cluster_ID[i]=i+1 for i=0,Crystal_Mult
// G4cout << "here we will clusterize" << G4endl;

	int i=0;
	int j=0;
	int Label=0;
	int intcounter;
	int Hits;
	int Internal_ID[100];
	// Cluster_ID[i]=i+1 was set before

	for(int loop=0;loop<Crystal_Mult;loop++){

	if(loop==0 || Cluster_ID[loop]>Label){

		Label++;
		intcounter=0;
		Hits=1;
		Internal_ID[0]=E_crystal_ID[loop];
	
		for(int l=0;l<Hits;l++){
			for(j=0;j<Crystal_Mult;j++){
				if(Internal_ID[l]!=E_crystal_ID[j]){ // no need to check the same crystals
			
//				G4cout << Internal_ID[l] << "	" << E_crystal_ID[j] << "	-----	" <<
//				NNP->NNIsNeighbor(Internal_ID[l],E_crystal_ID[j]) << G4endl;
		
					if(NNP->NNIsNeighbor(Internal_ID[l],E_crystal_ID[j])==1){	
					// if neighbor is found label it and add it to internal list (Internal_ID[]) for checking
					// further in l loop (Hits is incremented by one)
					// this way a new participant of the cluster is being chacked agains the other crystals
				
						if(Cluster_ID[j]>=0 && Cluster_ID[j]>Label) {
						
						// here if the crystal is already labeled we skip
						// so that we do not have to repeat already labeled crystals
						// more speed to it
											
							Cluster_ID[j]=Label; // this one is important- Adding the
									     //	crystals to the cluster
					
							Hits++;
							intcounter++;
							Internal_ID[intcounter]=E_crystal_ID[j];	
						}
					}
			
				}
			}
		}
	if(Hits==1) Cluster_ID[loop]=Label;
	}
	}
// end of cluster analysis
	
// --------------------------------------------------------------------------------------------------------
	// create clusters energy
// G4cout << "create clusters energy" << G4endl;

	Cluster_Mult=Label;

//	G4cout << Crystal_Mult << "	" << Cluster_Mult << G4endl;

	f6->Fill(Cluster_Mult);
  	f15->Fill(Prim_Mult,Cluster_Mult);

	double Cluster_time[160];
	double Cluster_theta[160];
	double Cluster_phi[160];
	int Crystals_in_Cluster[160];
	
	// zero the array
	for(i=0;i<Cluster_Mult;i++){
		Crystals_in_Cluster[i]=0;
		Cluster_theta[i]=0.;
		Cluster_phi[i]=0.;
		Cluster_time[i]=0.;		
		E_cluster_energy[i]=0;

	}
	// sum up the energy in the crystals of the same cluster
// G4cout << "sum up the energy in the crystals of the same cluster" << G4endl;

		
// write out crystal binary file the same way as in DANCE analysis
	unsigned char buffer[3];
	if(GPP->BinaryCrystalFile){
		
		buffer[0]=(unsigned char)(int (0));
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//1
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//2
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//3
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//4
	
		buffer[0]=(unsigned char)(int (2.0)); // PPAC_mult
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//5
		
		buffer[0]=(unsigned char)(int (Crystal_Mult)); // Crystal_mult
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//6
		
		buffer[0]=(unsigned char)(int (Cluster_Mult)); // Cluster_mult
		GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//7
		
	}
	
	for(i=0;i<Crystal_Mult;i++){
	
	
		if((Cluster_ID[i]-1)<0) G4cout << "(Cluster_ID[i]-1)<0 !!!! " << G4endl;
		
		Crystals_in_Cluster[Cluster_ID[i]-1]++;
		Cluster_theta[Cluster_ID[i]-1]+=GPP->GetTheta(E_crystal_ID[i]);
		Cluster_phi[Cluster_ID[i]-1]+=GPP->GetPhi(E_crystal_ID[i]);
		Cluster_time[Cluster_ID[i]-1]+=E_crystal_time[i];		
		E_cluster_energy[Cluster_ID[i]-1]+=E_crystal_energy[i];
		
		if(GPP->BinaryCrystalFile){
			buffer[0]=(unsigned char)(E_crystal_ID[i]); // crystal_ID
			GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//7

			buffer[0]=(unsigned char)(Cluster_ID[i]-1); // cluster ID			
			GPP->binary_crystal_file->write((char *) &buffer[0],sizeof(unsigned char));		//7
			
			int intEcrystal=int(E_crystal_energy[i]*65535./16.);
			buffer[1]=(unsigned char)(int (1.*intEcrystal/256.));
			intEcrystal=intEcrystal-(int)(buffer[1]*256.);		
			buffer[2]=(unsigned char)(intEcrystal);
			
			
			GPP->binary_crystal_file->write((char *) &buffer[1],sizeof(unsigned char));
			GPP->binary_crystal_file->write((char *) &buffer[2],sizeof(unsigned char));		
		}
	}
	
// Fill the cluster related histograms
// G4cout << "Fill the cluster related histograms" << G4endl;
	
	for(i=0;i<Cluster_Mult;i++){
		if(Crystals_in_Cluster[i]<=0) G4cout << "Crystals_in_Cluster[i]=" << Crystals_in_Cluster[i] << G4endl;
		Cluster_theta[i]=Cluster_theta[i]/(1.*Crystals_in_Cluster[i]);
		Cluster_phi[i]=Cluster_phi[i]/(1.*Crystals_in_Cluster[i]);
		Cluster_time[i]=Cluster_time[i]/(1.*Crystals_in_Cluster[i]);

		f7->Fill(Cluster_theta[i]);
		f8->Fill(Cluster_phi[i]);
		f11->Fill(Cluster_time[i]);	
		f4->Fill(E_cluster_energy[i]);
		
		// now fill in cluster energy
		
		if(Cluster_Mult==2) f20->Fill(E_cluster_energy[i]);	
		if(Cluster_Mult==3) f21->Fill(E_cluster_energy[i]);	
		if(Cluster_Mult==4) f22->Fill(E_cluster_energy[i]);	
		if(Cluster_Mult==5) f23->Fill(E_cluster_energy[i]);	
		if(Cluster_Mult==6) f24->Fill(E_cluster_energy[i]);	
		
		
		if(E_ball>=(GPP->Q_energy-GPP->dQ_energy) && E_ball<(GPP->Q_energy+GPP->dQ_energy)){
			if(Cluster_Mult==1) f29->Fill(E_cluster_energy[i]);
			if(Cluster_Mult==2) f30->Fill(E_cluster_energy[i]);	
			if(Cluster_Mult==3) f31->Fill(E_cluster_energy[i]);	
			if(Cluster_Mult==4) f32->Fill(E_cluster_energy[i]);	
			if(Cluster_Mult==5) f33->Fill(E_cluster_energy[i]);	
			if(Cluster_Mult==6) f34->Fill(E_cluster_energy[i]);	
		}
			// now cluster gate flags
			for(int gatei=0;gatei<3;gatei++){
						
				if(GPP->SingleGateECluster_E[gatei]>0){
				
					if(E_cluster_energy[i]>=(GPP->SingleGateECluster_E[gatei]-GPP->SingleGateECluster_dE[gatei]) &&
					E_cluster_energy[i]<=(GPP->SingleGateECluster_E[gatei]+GPP->SingleGateECluster_dE[gatei])){
					
						ClusterGate_Flag[gatei]=100;
						ClusterGate_ID[gatei]=i;
					
					}
				
				}
			
			}
	}
	
//----------------------------------------
// fill cluster energy gated histograms	
//------------------------------------------

  	if(ClusterGate_Flag[0]==100){
		f2_clgate1->Fill(E_ball-E_cluster_energy[ClusterGate_ID[0]]);
  	}
	if(ClusterGate_Flag[1]==100){
		f2_clgate2->Fill(E_ball-E_cluster_energy[ClusterGate_ID[1]]);
	}
	if(ClusterGate_Flag[2]==100){
		f2_clgate3->Fill(E_ball-E_cluster_energy[ClusterGate_ID[2]]);
	}


// cluster multiplicity gated

	if(GPP->ClusterMGatedHists){
						
		if(Cluster_Mult==1)	f1_cl1->Fill(E_ball);
		if(Cluster_Mult==2)	f1_cl2->Fill(E_ball);
		if(Cluster_Mult==3)	f1_cl3->Fill(E_ball);
		if(Cluster_Mult==4)	f1_cl4->Fill(E_ball);
		if(Cluster_Mult==5)	f1_cl5->Fill(E_ball);
		if(Cluster_Mult==6)	f1_cl6->Fill(E_ball);
		if(Cluster_Mult==7)	f1_cl7->Fill(E_ball);
		if(Cluster_Mult>=8) 	f1_cl8->Fill(E_ball);	
	}
	
/*
	G4cout << "Clusters Multiplicity: " << Cluster_Mult << G4endl;
	for(int loop=0;loop<Crystal_Mult;loop++){
		G4cout << "Crys No: " << E_crystal_ID[loop] << "	Cluster:	" << Cluster_ID[loop] << G4endl;
	}
*/

} // end of event analysis

// --------------------------------------------------------------------------------------------------------
// Autosave of ROOT histograms
	if(GPP->RootFile && (1.*int(numberOfEvent/(1.*GPP->AutoSave)))==(1.*numberOfEvent/(1.*GPP->AutoSave))) RootAutoSave();
// --------------------------------------------------------------------------------------------------------

// G4cout << "End of Event " << G4endl;
}

G4double DANCEEventAction::GetTotal(const G4THitsMap<G4double> &map) const
{

  G4double tot = 0.;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) 
//  { tot += *(itr->second); }
  { tot += *(itr->second); }
  
  return tot;
}

G4double DANCEEventAction::FindMinimum(const G4THitsMap<G4double> &map) const
{
/*
  G4double val = DBL_MAX;
  std::map<G4int,G4double*>::iterator itr = map.GetMap()->begin();
  for(; itr != map.GetMap()->end(); itr++) 
  { if(val>*(itr->second)) val = *(itr->second); }
  return val;
*/
return 0;
}

void DANCEEventAction::PrintPrimary(G4PrimaryParticle* pp,G4int ind)
{
  for(G4int ii=0;ii<=ind;ii++)
  { G4cout << "  "; }
  G4cout << "==PDGcode " << pp->GetPDGcode() << " ";
  if(pp->GetG4code()!=0)
  { G4cout << "(" << pp->GetG4code()->GetParticleName() << ")"; }
  else
  { G4cout << "is not defined in G4"; }
  G4cout << " " << pp->GetMomentum()/GeV << " [GeV] ";
  if(pp->GetTrackID()<0)
  { G4cout << G4endl; }
  else
  { G4cout << ">>> G4Track ID " << pp->GetTrackID() << G4endl; }

 // G4PrimaryParticle* daughter = pp->GetNext();
  G4PrimaryParticle* daughter = pp->GetDaughter();
  while(daughter)
  {
    G4cout << "daughter";
    PrintPrimary(daughter,ind+1);
    daughter = daughter->GetNext();
  }
}


void DANCEEventAction::RootAutoSave(){

	DANCEGlobal* GPP=DANCEGlobal::Instance();
	
	G4cout << "AutoSaving ROOT histograms in file: " << GPP->RootFileName << G4endl;
	printf("Events processed: %d \n",numberOfEvent);

	rootfileout=new TFile(GPP->RootFileName,"RECREATE");
	rootfileout->cd();

	f1->Write();
	f2->Write();
	ECr_ID->Write();
	f3->Write();
	f4->Write();
	f5->Write();
	f6->Write();
	f7->Write();
	f8->Write();
	f9->Write();
	f10->Write();
	f11->Write();
	f12->Write();
	f13->Write();
	f14->Write();
	f15->Write();
	f16->Write();
	f17->Write();

	f20->Write();
	f21->Write();
	f22->Write();
	f23->Write();
	f24->Write();
	
	f29->Write();
	f30->Write();
	f31->Write();
	f32->Write();
	f33->Write();
	f34->Write();
	
	if(GPP->ClusterMGatedHists){
		f1_cl1->Write();
		f1_cl2->Write();
		f1_cl3->Write();
		f1_cl4->Write();
		f1_cl5->Write();
		f1_cl6->Write();
		f1_cl7->Write();
		f1_cl8->Write();
	}

	if(GPP->CrystalMGatedHists){
		f1_cr1->Write();
		f1_cr2->Write();
		f1_cr3->Write();
		f1_cr4->Write();
		f1_cr5->Write();
		f1_cr6->Write();
		f1_cr7->Write();
		f1_cr8->Write();
	}

	if(GPP->SingleGateECrystal_E[0]>0){	
		f2_crgate1->Write();
	}
	if(GPP->SingleGateECrystal_E[1]>0){	
		f2_crgate2->Write();
	}
	if(GPP->SingleGateECrystal_E[2]>0){	
		f2_crgate3->Write();
	}

	if(GPP->SingleGateECluster_E[0]>0){	
		f2_clgate1->Write();
	}
	if(GPP->SingleGateECluster_E[1]>0){	
		f2_clgate2->Write();
	}
	if(GPP->SingleGateECluster_E[2]>0){	
		f2_clgate3->Write();
	}


	rootfileout->ls();
	rootfileout->Close();


	G4cout << "RootAutoSave completed" << G4endl;

}
  
