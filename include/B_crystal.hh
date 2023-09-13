
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


class B_crystal{

public:

B_crystal();
~B_crystal();

// user defined

// here crystal B

G4double pDz_B1,pTheta_B1,pPhi_B1;
G4double pDy1_B1,pDx1_B1,pDx2_B1,pAlp1_B1; 
G4double BaF_B_ratio,pDy2_B1,pDx3_B1,pDx4_B1,pAlp2_B1;
G4double pDz_B2,pTheta_B2,pPhi_B2;
G4double pDy1_B2,pDx1_B2,pDx2_B2,pAlp1_B2; 
G4double pDy2_B2,pDx3_B2,pDx4_B2,pAlp2_B2;

G4double BaF_B_frontdist;
G4double BaF_B_length;

G4double trans,trans_loc;
// constructor

B_crystal(G4double a1,G4double a2,G4double a3,G4double a4,G4double a5,G4double a6,G4double a7){

// a1 - length
// a2-a6 
// a7 - front distance
	BaF_B_frontdist=a7;
	BaF_B_length=a1;

	pDz_B1=a1/2.;
	pTheta_B1=0.;
	pPhi_B1=90*deg;
 
	pDy1_B1=a2;
	pDx1_B1=a3;
	pDx2_B1=a4; 
	pAlp1_B1=0.; 

	BaF_B_ratio=(1.+a1/BaF_B_frontdist);

	pDy2_B1=BaF_B_ratio*pDy1_B1;
	pDx3_B1=BaF_B_ratio*pDx1_B1;
	pDx4_B1=BaF_B_ratio*pDx2_B1;
	pAlp2_B1=0.;

	pTheta_B1=-atan(pDy1_B1/BaF_B_frontdist)*rad;

	pDz_B2=pDz_B1;
	pTheta_B2=0.;
	pPhi_B2=90*deg;
 
	pDy1_B2=a5; //0.1140406/2.*RI;
	pDx1_B2=pDx2_B1;
	pDx2_B2=a6; 
	pAlp1_B2=0.; 

	pDy2_B2=BaF_B_ratio*pDy1_B2;
	pDx3_B2=BaF_B_ratio*pDx1_B2;
	pDx4_B2=BaF_B_ratio*pDx2_B2;
	pAlp2_B2=0.;

	pTheta_B2=atan(pDy1_B2/BaF_B_frontdist)*rad;

	trans_loc=0.5*(pDy1_B1+pDy2_B1+pDy1_B2+pDy2_B2);
	trans=-0.5*(pDy1_B1+pDy2_B1);
//	trans=-0.5*(pDy1_B2+pDy2_B2);
};

G4UnionSolid *Build_B_crystal(){

	G4Trap * BaF_B1_crystal=new G4Trap("BaF_B1_crystal",pDz_B1, pTheta_B1,pPhi_B1, pDy1_B1,pDx1_B1, pDx2_B1,pAlp1_B1,pDy2_B1,pDx3_B1, pDx4_B1,pAlp2_B1);
	G4Trap * BaF_B2_crystal=new G4Trap("BaF_B2_crystal",pDz_B2, pTheta_B2,pPhi_B2, pDy1_B2,pDx1_B2, pDx2_B2,pAlp1_B2,pDy2_B2,pDx3_B2, pDx4_B2,pAlp2_B2);

	G4UnionSolid *BaF_B_crystal=new G4UnionSolid("BaF_B_crystal",BaF_B1_crystal,BaF_B2_crystal,0,G4Vector3D(0,trans_loc,0));
	
	return BaF_B_crystal;
};

};
