
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


//#include "DANCENNNearestNeighbors.hh"
#include <fstream>
#include <TNamed.h>
#include <iostream>
#include "G4String.hh"
//#include "TString.h"

class DANCEGlobal : public TNamed
{
public:
	static DANCEGlobal* Instance();

	int GetID(){ return ID;};
	
	void Init(void );
	void Initialize(void );
	
	void SetID(int new_ID){ID=new_ID;};
	
	void SetID1(int No,int new_ID){mapID1[No]=new_ID;};
	void SetID2(int No,int new_ID){mapID2[No]=new_ID;};
	void SetID3(int No,int new_ID){mapID3[No]=new_ID;};
	void SetID4(int No,int new_ID){mapID4[No]=new_ID;};

// here we fill real crystal ID for all dettypes and detNo	

	void SetID2D(int DetType,int DetNo,int new_ID){    
		mapID2D[DetType][DetNo]=new_ID;
		mapID2D_DetType[new_ID]=DetType;	
		mapID2D_DetNo[new_ID]=DetNo;	
	};

	void SetGeantID2D(int DetType,int DetNo,int new_ID){
		GeantID2D[DetType][DetNo]=new_ID;
		GeantID2D_DetType[new_ID]=DetType;	
		GeantID2D_DetNo[new_ID]=DetNo;	
	};
	
	int GetSymID(int inID){
		return SymID[inID];
	};

	int GetID1(int No){return mapID1[No];};
	int GetID2(int No){return mapID2[No];};
	int GetID3(int No){return mapID3[No];};
	int GetID4(int No){return mapID4[No];};

	int GetID2D(int DetType,int DetNo){return mapID2D[DetType][DetNo];};
	int GetID2D_DetType(int old_ID){return mapID2D_DetType[old_ID];};
	int GetID2D_DetNo(int old_ID){return mapID2D_DetNo[old_ID];};

	int GetGeantID2D(int DetType,int DetNo){return GeantID2D[DetType][DetNo];};
	int GetGeantID2D_DetType(int old_ID){return GeantID2D_DetType[old_ID];};
	int GetGeantID2D_DetNo(int old_ID){return GeantID2D_DetNo[old_ID];};

	void SetTheta(int temp_ID,double temp_theta){
		Theta[temp_ID]=temp_theta;
	}
	void SetPhi(int temp_ID,double temp_phi){
		Phi[temp_ID]=temp_phi;
	}

	double GetTheta(int temp_ID) {return Theta[temp_ID];};	
	double GetPhi(int temp_ID) {return Phi[temp_ID];};
	
// reading the Master input helper method

	bool CheckForFlags(ifstream *in);
	
	void EndOfRun(){
	
		binary_crystal_file->close();
	
	};	
		
private:
	DANCEGlobal();
	~DANCEGlobal();
	
	int ID;
	int mapID1[162];
	int mapID2[162];	
	int mapID3[162];	
	int mapID4[162];	

	int SymID[80];	

	int GeantID2D[4][162];
	int GeantID2D_DetType[162];
	int GeantID2D_DetNo[162];

	int RealID2D[4][61];

	int mapID2D[4][61];
	int mapID2D_DetType[162];
	int mapID2D_DetNo[162];
	
	double Theta[162],Phi[162];

public:	
	bool BinaryCrystalFile;
	G4String BinaryCrystalFile_Name;
	
	bool BinaryClusterFile;
	G4String BinaryClusterFile_Name;
	
	bool AddPrimaryGammas;
	bool AddTimeOfGammas;
	int AutoSave;
	bool RootFile;
	bool ClusterMGatedHists;
	bool CrystalMGatedHists;
	bool LiHShell[2];
	double LiHDensity;
	double LiHInnerRadius;
	double LiHOuterRadius;
	bool BeamPipe[2];
	bool SupportingSphere[2];
	bool Holders[2];
	bool PMLense[2];
	int ShowCrystals[2];
	bool PVCShell;
		
	bool ExternalInput;
	G4String ExternalInputFile;

	bool MacroInput;
	G4String MacroInputFile;
	
	double SingleGateECrystal_dE[3];	
	double SingleGateECrystal_E[3];	
	double SingleGateECluster_dE[3];	
	double SingleGateECluster_E[3];	

	G4String RootFileName;
	G4String MasterInputFile;
		
	double E_threshold;
	bool BeamPipeVacuum;
	double CrystalDistance;
	bool EndOfExternalInput;

	int MultiplicityOfPrimaryEvent;
	
	int PrimaryTrackID[30];
	double PrimaryEnergy[30];
	int Prim_Mult;
	
	double Q_energy;
	double dQ_energy;
	
	ofstream *binary_crystal_file;


};

