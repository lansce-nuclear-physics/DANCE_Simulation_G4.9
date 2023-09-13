
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


class C_crystal{

public:
C_crystal(){};
~C_crystal();


// user defined 

G4double BaF_C_length;				//=1.0*BaF_length/2.; // length
G4double BaF_C_front_th1;			//=0.1749639/2.*RI;   // thickness of the middle part 
G4double BaF_C_front_th2;			//=0.0977372/2.*RI; // thickness of the left part 
G4double BaF_C_front_th3;			//=0.028076016/2.*RI; // thickness of the right part 


G4double BaF_C_front_height1;			//=0.19423*sin(59.7*deg)*RI;
G4double BaF_C_front_height2;			//0.16622*cos(9.7*deg)*RI;


G4double BaF_C_frontdist;			//=17.*cm; //BaF_frontdist;

G4double shift_C;				//=0.9*cm;

// derived

G4double BaF_C_ratio,BaF_C_back_th1,BaF_C_back_th2,BaF_C_back_th3,BaF_C_back_height1,BaF_C_back_height2;
G4double frontshift_C;
G4double mid_pDz,mid_pTheta,mid_pPhi,mid_pDy1,mid_pDx1,mid_pDx2,mid_pAlp1; 
G4double mid_pDy2,mid_pDx3,mid_pDx4,mid_pAlp2,mid_trans;

G4double left_pDz1,left_pPhi1,left_pDy1,left_pDx1,left_pDx2,left_pAlp1; 
G4double left_pDy2,left_pDx3,left_pDx4,left_pAlp2;
G4double left_pTheta1,left_trans;

G4double right_pDz1,right_pPhi1,right_pDy1,right_pDx1,right_pDx2,right_pAlp1; 
G4double right_pDy2,right_pDx3,right_pDx4,right_pAlp2;
G4double right_pTheta1,right_trans;

///// basic shapes at 0,0,0

C_crystal(G4double a1,G4double a2,G4double a3,G4double a4,G4double a5,G4double a6,G4double a7,G4double a8){

	BaF_C_length=a1;
	BaF_C_front_th1=a2;
	BaF_C_front_th2=a3;
	BaF_C_front_th3=a4;
	BaF_C_front_height1=a5;
	BaF_C_front_height2=a6;
	BaF_C_frontdist=a7;
	shift_C=a8;


// derived dimensions

	BaF_C_ratio=(1.+2.*BaF_C_length/BaF_C_frontdist);

	BaF_C_back_th1=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_th1;
	BaF_C_back_th2=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_th2;
	BaF_C_back_th3=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_th3;

	BaF_C_back_height1=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_height1;
	BaF_C_back_height2=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_height2;


// definition of parts of C irregular hexagon
// here we assume some basic properties of crystal
// prolonged edges will meet at 0,0,0 
// at BaF_C_frontdist-half of the length of crystal (this can be fixed)

// middle part of the crystal C dimensions
//    /----\
//   /      \
//  /________\


	frontshift_C=shift_C;

//shift_C*(1.+2.*BaF_C_length/BaF_C_frontdist);

	mid_pDz=BaF_C_length;
	mid_pTheta=0.;
	mid_pPhi=0.;
 
	mid_pDy1=BaF_C_front_th1;
	mid_pDx1=BaF_C_front_height1; 
	 mid_pDx2=BaF_C_front_height2;
	mid_pAlp1=0.; 

	mid_pDy2=BaF_C_back_th1;
	mid_pDx3=BaF_C_back_height1; 
	mid_pDx4=BaF_C_back_height2;
	mid_pAlp2=0.;


// shifting focus of the crystal
	mid_pPhi=90*deg;
	mid_pTheta=(atan(0.5*frontshift_C/(BaF_C_frontdist/2.)))*rad;


	mid_trans=-2*BaF_C_length*tan(mid_pTheta);
// left part of the crystal C dimensions

	left_pDz1=BaF_C_length;
	left_pPhi1=-90*deg;
 
	left_pDy1=BaF_C_front_th2;
	left_pDx1=0.01; 
	left_pDx2=BaF_C_front_height1;
	left_pAlp1=0.; 

	left_pDy2=BaF_C_back_th2;
	left_pDx3=(BaF_C_ratio*0.01); 
	left_pDx4=BaF_C_back_height1;
	left_pAlp2=0.;

	left_pTheta1=atan(-0.5*(left_pDy1-left_pDy2+mid_pDy1-mid_pDy2-mid_trans)/BaF_C_length)*rad;
	left_trans=-0.5*(left_pDy2+left_pDy1+mid_pDy2+mid_pDy1)-mid_trans;


///// ADD right piece

	right_pDz1=BaF_C_length;
	right_pPhi1=90.*deg;
 
	right_pDy1=BaF_C_front_th3;
	right_pDx1=BaF_C_front_height2; 
	right_pDx2=0.01;
	right_pAlp1=0.; 

	right_pDy2=BaF_C_back_th3;
	right_pDx3=BaF_C_back_height2; 
	right_pDx4=(BaF_C_ratio*0.01);
	right_pAlp2=0.;

	right_pTheta1=atan(0.5*(right_pDy2-right_pDy1+mid_pDy2-mid_pDy1-mid_trans)/BaF_C_length)*rad;
	right_trans=0.5*(right_pDy2+right_pDy1+mid_pDy2+mid_pDy1)-mid_trans;;
};


G4UnionSolid *Build_C_crystal(){

G4Trap * BaF_C1_crystal=	new G4Trap("BaF_C1_crystal",mid_pDz, mid_pTheta,mid_pPhi, mid_pDy1,mid_pDx1, mid_pDx2,mid_pAlp1,mid_pDy2,mid_pDx3, mid_pDx4,mid_pAlp2);
G4Trap * BaF_C2_crystal=	new G4Trap("BaF_C2_crystal",left_pDz1, left_pTheta1,left_pPhi1, left_pDy1,left_pDx1, left_pDx2,left_pAlp1,left_pDy2,left_pDx3, left_pDx4,left_pAlp2);
G4Trap * BaF_C3_crystal=	new G4Trap("BaF_C3_crystal",right_pDz1, right_pTheta1,right_pPhi1, right_pDy1,right_pDx1, right_pDx2,right_pAlp1,right_pDy2,right_pDx3, right_pDx4,right_pAlp2);

G4UnionSolid *BaF_C_crystal_1=new G4UnionSolid("BaF_C_crystal_1",BaF_C1_crystal,BaF_C2_crystal,0,G4Vector3D(0,0+left_trans+mid_trans,0));
G4UnionSolid *BaF_C_crystal=new G4UnionSolid("BaF_C_crystal",BaF_C_crystal_1,BaF_C3_crystal,0,G4Vector3D(0,0+right_trans+mid_trans,0));

return BaF_C_crystal;

};


};
