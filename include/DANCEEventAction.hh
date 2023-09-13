
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


#ifndef DANCEEventAction_h
#define DANCEEventAction_h 1

#include "globals.hh"
#include "G4Event.hh"

#include "G4THitsMap.hh"
#include "G4UserEventAction.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
using namespace AIDA;
#endif // G4ANALYSIS_USE

#include "TROOT.h"
#include "TObject.h"

#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttPad.h"

#include "TVirtualPad.h"

#include "TNamed.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TH1D.h"
#include "TCanvas.h"

class G4Event;

class DANCEEventAction : public G4UserEventAction
{
  public:
    DANCEEventAction();
    virtual ~DANCEEventAction();

  public:
    	virtual void BeginOfEventAction(const G4Event*);
    	virtual void EndOfEventAction(const G4Event*);
	void RootAutoSave();
	
  private:
    	G4double GetTotal(const G4THitsMap<G4double> &map) const;
    	G4double FindMinimum(const G4THitsMap<G4double> &map) const;
	void PrintPrimary(G4PrimaryParticle* pp,G4int ind);

    	G4int DHID[170];
    
  private:

//----------------------------------------------------
// Reserve the pointers for ROOT Histograms
	TFile *rootfileout;    
  
	TH1F *f1;
	TH1F *f2;
	TH2F *ECr_ID;
	TH1F *f4;
	TH1F *f5;
	
	TH1F *f6;
	TH1F *f3;

	TH1F *f7;
	TH1F *f8;
	TH1F *f9;
	TH1F *f10;
	TH1F *f11;
	TH1F *f13;

	TH1F *f12;
	TH2F *f14;
	TH2F *f15;
	TH1F *f16;
	TH1F *f17;

	TH1F *f20;
	TH1F *f21;
	TH1F *f22;
	TH1F *f23;
	TH1F *f24;

	TH1F *f29;
	TH1F *f30;
	TH1F *f31;
	TH1F *f32;
	TH1F *f33;
	TH1F *f34;

	TH1F *f1_cl1;
	TH1F *f1_cl2;
	TH1F *f1_cl3;
	TH1F *f1_cl4;
	TH1F *f1_cl5;
	TH1F *f1_cl6;
	TH1F *f1_cl7;
	TH1F *f1_cl8;

	TH1F *f1_cr1;
	TH1F *f1_cr2;
	TH1F *f1_cr3;
	TH1F *f1_cr4;
	TH1F *f1_cr5;
	TH1F *f1_cr6;
	TH1F *f1_cr7;
	TH1F *f1_cr8;

	TH1F *f2_clgate1;
	TH1F *f2_clgate2;
	TH1F *f2_clgate3;
	
	TH1F *f2_crgate1;
	TH1F *f2_crgate2;
	TH1F *f2_crgate3;
	
	G4double Total,Peak,Total_withthreshold;


//----------------------------------------------------


    	G4THitsMap<G4double> mapSum[1500][4];
    	G4int colIDSum[1500][4];

    	G4THitsMap<G4double> mapMin[20][4];
    	G4int colIDMin[20][4];
    	G4int numberOfEvent,gatedEvents;

	G4String detName[80];
	G4String primNameSum[4];
	G4String fullName[1500];
	
	
#ifdef G4ANALYSIS_USE
        IHistogram1D* dc1Hits;
	IHistogram1D* dc2Hits;

	IHistogram1D* multHits;
	IHistogram1D* IDdet;


//	ICloud2D* dc1XY;
//	ICloud2D* dc2XY;
//	ICloud2D* evstof;

	ITuple* tuple;
	IPlotter* plotter;
#endif // G4ANALYSIS_USE
	
  public:
    inline G4double GetTotalE(G4int i) const
    { return GetTotal(mapSum[i][0]); }
    inline G4double GetNGamma(G4int i) const
    { return GetTotal(mapSum[i][1]); }
    inline G4double GetNElectron(G4int i) const
    { return GetTotal(mapSum[i][2]); }
    inline G4double GetNPositron(G4int i) const
    { return GetTotal(mapSum[i][3]); }
    inline G4double GetTotalL(G4int i) const
    { return GetTotal(mapSum[i][4]); }
    inline G4double GetNStep(G4int i) const
    { return GetTotal(mapSum[i][5]); }
    inline G4double GetEMinGamma(G4int i) const
    { return FindMinimum(mapMin[i][0]); }
    inline G4double GetEMinElectron(G4int i) const
    { return FindMinimum(mapMin[i][1]); }
    inline G4double GetEMinPositron(G4int i) const
    { return FindMinimum(mapMin[i][2]); }
};

#endif

