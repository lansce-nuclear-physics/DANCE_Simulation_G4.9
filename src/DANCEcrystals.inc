
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


//------------------------------ a BaF_C crystal  Detector

G4double RI=17.*cm;
//G4double RA=32*cm;

G4double BaF_length=15.*cm;
G4double BaF_frontdist=17.*cm;
G4double pi=acos(-1.);
G4double BaF_ratio=(1.+BaF_length/BaF_frontdist);

//////////////////////////////
// HERE CRYSTAL A           \\
//////////////////////////////

//G4int BaF_A_crystal_numRZ=5;
G4double BaF_A_crystal_R1=RI*0.5*(0.21436)/tan(36.*deg);
G4double BaF_A_crystal_R2=BaF_ratio*RI*0.5*(0.21436)/tan(36.*deg);

//G4double BaF_A_crystal_r[2]={BaF_A_crystal_R1,BaF_A_crystal_R2};
//G4double BaF_A_crystal_z[2]={0.,BaF_length};


G4double phiStart = 0;
G4double phiTotal= 2*acos(-1.);
G4int numSide = 5;
G4int numZPlanes = 2;
//G4double rInner[7]={0,0,0,0,0,0,0};
//G4double rOuter[7]= { 0, 15, 15, 4, 4, 10, 10 };
//G4double zPlane[7]={0,5,8,13,30,32,35};

G4double rInner[2]={0,0};
G4double rOuter[2]= {BaF_A_crystal_R1,BaF_A_crystal_R2};
G4double zPlane[2]={0,BaF_length};

G4Polyhedra* BaF_A_crystal=new G4Polyhedra("BaF_A_crystal",phiStart,phiTotal,numSide,numZPlanes,zPlane,rInner,rOuter);


// wrapping for crystal A
G4double thickness=1*mm;
G4double rInner_thicker[2]={0,0};
G4double rOuter_thicker[2]= {BaF_A_crystal_R1+thickness,BaF_A_crystal_R2+thickness};
G4double zPlane_thicker[2]={0,BaF_length+thickness};

G4Polyhedra* BaF_A_crystal_thicker=new G4Polyhedra("BaF_A_crystal_thicker",phiStart,phiTotal,numSide,numZPlanes,zPlane_thicker,rInner_thicker,rOuter_thicker);
G4SubtractionSolid *BaF_A_wrap=new G4SubtractionSolid("BaF_A_wrap",BaF_A_crystal_thicker,BaF_A_crystal,0,G4ThreeVector());

// holder for crystal A
G4double inch=2.54*cm;
G4double rInnerHol[2]={BaF_A_crystal_R2,2.375*inch};
//G4double rInnerHol[2]={0.5*1.5*inch,0.5*3.91*inch};
G4double rOuterHol[2]= {BaF_A_crystal_R2+1.*inch,1.*3.0*inch};
G4double zPlaneHol[2]={0,3.*inch};

G4Polyhedra* BaF_A_crystal_holder_part1=new G4Polyhedra("BaF_A_crystal_holder_part1",phiStart,phiTotal,numSide,numZPlanes,zPlaneHol,rInnerHol,rOuterHol);
G4Cons *BaF_A_crystal_holder_part2=new G4Cons("BaF_A_crystal_holder_part2",1.86*inch,2.375*inch,1.86*inch,2.375*inch,3*inch,0.*deg,360.*deg);
//G4Tubs *BaF_A_crystal_holder_part3=new G4Tubs("BaF_A_crystal_holder_part3",1.86*inch,2.375*inch,1.5*inch,0.*deg,360.*deg);
//G4Cons *BaF_A_crystal_holder_part3=new G4Cons("BaF_A_crystal_holder_part3",1.86*inch,2.375*inch,1.86*inch,2.375*inch,1.5*inch,0.*deg,360.*deg);

G4SubtractionSolid *BaF_A_crystal_holder=new G4SubtractionSolid("BaF_A_crystal_holder_1",BaF_A_crystal_holder_part2,BaF_A_crystal_holder_part1,0,G4ThreeVector(0,0,-3.*inch));
//G4UnionSolid *BaF_A_crystal_holder=new G4UnionSolid("BaF_A_crystal_holder",BaF_A_crystal_holder_1,BaF_A_crystal_holder_part3,0,G4ThreeVector(0,0,3*inch));


//////////////////////////////
// here CRYSTAL D
//////////////////////////////

// this heaxagon is strange - paper dimensions show 0.17742 side

//G4int BaF_D_crystal_numRZ=12;
G4double BaF_D_crystal_R1=(0.17442)*RI*cos(30.*deg);
G4double BaF_D_crystal_R2=BaF_ratio*(0.15)*RI;

//G4double BaF_D_crystal_r[2]={BaF_D_crystal_R1,BaF_D_crystal_R2};
//G4double BaF_D_crystal_z[2]={0.,BaF_length};

G4int numSide_D = 6;
G4int numZPlanes_D = 2;
G4double rInner_D[2]={0,0};
G4double rOuter_D[2]= {BaF_D_crystal_R1,BaF_D_crystal_R2};
G4double zPlane_D[2]={0,BaF_length};


G4Polyhedra* BaF_D_crystal=new G4Polyhedra("BaF_D_crystal",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlane_D,rInner_D,rOuter_D);
//BaF_D_crystal_log= new G4LogicalVolume(BaF_D_crystal,BaF,"BaF_A_crystal_log",0,0,0);
//G4Polyhedra* BaF_D_crystal_wrap=new G4Polyhedra("BaF_D_crystal",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlane_D,rInner_D,rOuter_D);

// wrapping
G4double rInner_D_thicker[2]={0,0};
G4double rOuter_D_thicker[2]= {BaF_D_crystal_R1+thickness,BaF_D_crystal_R2+thickness};
G4double zPlane_D_thicker[2]={0,BaF_length+thickness};

G4Polyhedra* BaF_D_crystal_thicker=new G4Polyhedra("BaF_D_crystal_thicker",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlane_D_thicker,rInner_D_thicker,rOuter_D_thicker);
G4SubtractionSolid *BaF_D_wrap=new G4SubtractionSolid("BaF_D_wrap",BaF_D_crystal_thicker,BaF_D_crystal,0,G4ThreeVector());

// holder for crystal D

//G4double rInnerDHol[2]={BaF_D_crystal_R2,2.375*inch};
G4double rInnerDHol[2]={BaF_D_crystal_R2,2.375*inch};
G4double rOuterDHol[2]= {BaF_D_crystal_R2+3*inch,BaF_D_crystal_R2+3*inch};
G4double zPlaneDHol[2]={-0.02*inch,3.0*inch};

G4Polyhedra* BaF_D_crystal_holder_part1=new G4Polyhedra("BaF_D_crystal_holder_part1",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlaneDHol,rInnerDHol,rOuterDHol);
//G4Polyhedra* BaF_D_crystal_holder=new G4Polyhedra("BaF_D_crystal_holder",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlaneDHol,rInnerDHol,rOuterDHol);

//G4Cons *BaF_D_crystal_holder_part2=new G4Cons("BaF_D_crystal_holder_part2",0*inch,1.9*inch,0*inch,2.375*inch,1.5*inch,0.*deg,360.*deg);
G4Cons *BaF_D_crystal_holder_part2=new G4Cons("BaF_D_crystal_holder_part2",1.75*inch,2.125*inch,2.00*inch,2.375*inch,3.0*inch,0.*deg,360.*deg);


//G4Tubs *BaF_D_crystal_holder_part3=new G4Tubs("BaF_D_crystal_holder_part3",2.00*inch,2.375*inch,1.5*inch,0.*deg,360.*deg);

G4SubtractionSolid *BaF_D_crystal_holder=new G4SubtractionSolid("BaF_D_crystal_holder_1",BaF_D_crystal_holder_part2,BaF_D_crystal_holder_part1,0,G4ThreeVector(0,0,-3.*inch));
//G4UnionSolid *BaF_D_crystal_holder=new G4UnionSolid("BaF_D_crystal_holder",BaF_D_crystal_holder_1,BaF_D_crystal_holder_part3,0,G4ThreeVector(0,0,3*inch));



/*

// here crystal B

G4double pDz_B1=BaF_length/2.;
G4double pTheta_B1=0.;
G4double pPhi_B1=90*deg;
 
G4double pDy1_B1=0.144641/2.*RI;
G4double pDx1_B1=0.21428*RI/2.;
G4double pDx2_B1=0.20783*RI; 
G4double pAlp1_B1=0.; 

G4double BaF_B_ratio=(1.+BaF_length/BaF_frontdist);

G4double pDy2_B1=BaF_B_ratio*pDy1_B1;
G4double pDx3_B1=BaF_B_ratio*pDx1_B1;
G4double pDx4_B1=BaF_B_ratio*pDx2_B1;
G4double pAlp2_B1=0.;

pTheta_B1=-atan(pDy1_B1/BaF_frontdist)*rad;

G4double pDz_B2=BaF_length/2.;
G4double pTheta_B2=0.;
G4double pPhi_B2=90*deg;
 
G4double pDy1_B2=0.5*RI*0.1663*sin(43.295*deg); //0.1140406/2.*RI;
G4double pDx1_B2=0.20783*RI;
G4double pDx2_B2=0.17358*RI/2.; 
G4double pAlp1_B2=0.; 

G4double pDy2_B2=BaF_B_ratio*pDy1_B2;
G4double pDx3_B2=BaF_B_ratio*pDx1_B2;
G4double pDx4_B2=BaF_B_ratio*pDx2_B2;
G4double pAlp2_B2=0.;

pTheta_B2=atan(pDy1_B2/BaF_frontdist)*rad;


G4Trap * BaF_B1_crystal=new G4Trap("BaF_B1_crystal",pDz_B1, pTheta_B1,pPhi_B1, pDy1_B1,pDx1_B1, pDx2_B1,pAlp1_B1,pDy2_B1,pDx3_B1, pDx4_B1,pAlp2_B1);
G4Trap * BaF_B2_crystal=new G4Trap("BaF_B2_crystal",pDz_B2, pTheta_B2,pPhi_B2, pDy1_B2,pDx1_B2, pDx2_B2,pAlp1_B2,pDy2_B2,pDx3_B2, pDx4_B2,pAlp2_B2);

//BaF_B1_log= new G4LogicalVolume(BaF_B1_crystal,BaF,"BaF_B1_log",0,0,0);
//BaF_B2_log= new G4LogicalVolume(BaF_B2_crystal,BaF,"BaF_B2_log",0,0,0);

*/

// HERE BEAMPIPE - Aluminium 2.1825 - 2.5 cm


G4double pRMin=0.5*1.75*inch; //2.1825*cm;
G4double pRMax=0.5*1.83*inch; //2.5*cm;
G4double pDz=50*cm;
G4double pSPhi=0.*deg;
G4double pDPhi=360.*deg;

G4Tubs *BeamPipe=new G4Tubs("BeamPipe",pRMin,pRMax,pDz,pSPhi,pDPhi);

G4Tubs *RTH=new G4Tubs("RTH",2.*cm,pRMin,0.94*inch,pSPhi,pDPhi);

// LiH6 sphere

//G4Sphere* LiHSphere=new G4Sphere("LiHSphere",10.5*cm,16.0*cm,0.*deg,360.*deg,0.*deg,180.*deg);

// new LiH sphere with holes
//G4Sphere *LiHSphere=new G4Sphere("LiHSphere",10.5*cm,16.5*cm,0.*deg,360.*deg,0.*deg,180.*deg);
G4Sphere *LiHSphere=new G4Sphere("LiHSphere",GPP->LiHInnerRadius*cm,GPP->LiHOuterRadius*cm,0.*deg,360.*deg,0.*deg,180.*deg);

G4Tubs *BP_subtract=new G4Tubs("Hole",0.,pRMax,17.5*cm,pSPhi,pDPhi);

G4SubtractionSolid *LiHSphere2= new G4SubtractionSolid("LiHHoles", LiHSphere, BP_subtract,0,G4ThreeVector(0,0,0));

// support ball

G4Sphere *SuppSphere=new G4Sphere("SuppSphere",19.5*inch,21.5*inch,0.*deg,360.*deg,0.*deg,180.*deg);
G4Tubs *tube_subtract=new G4Tubs("Hole",0.,0.5*3.91*inch,62.*inch,pSPhi,pDPhi);

G4double rInner_Dhole[2]={0,0};
G4double rOuter_Dhole[2]= {3.*inch,3.*inch};
G4double zPlane_Dhole[2]={0.,22.*inch};

G4Polyhedra* hexa_subtract=new G4Polyhedra("Hole2",phiStart,phiTotal,numSide_D,numZPlanes_D,zPlane_Dhole,rInner_Dhole,rOuter_Dhole);

G4RotationMatrix *BeamPipeRMhexa=new G4RotationMatrix();

BeamPipeRMhexa->rotateY(90.*deg);


G4SubtractionSolid *SuppSphereHoles[162];
SuppSphereHoles[0]=new G4SubtractionSolid("SuppHoles", SuppSphere, tube_subtract,BeamPipeRMhexa,G4ThreeVector(0,0,0));


/////////////////////////////////////////
///// USED PARTS
///// THROWN AWAY BUT RECYCLABLE
///// ///////////////////////////////////

/*
///////////////////////////////////////////////////////////////////////////////////////////
// Crystal C
//////////////////////////////
//RI=32*cm;

G4double BaF_C_length=1.0*BaF_length/2.; // length
G4double BaF_C_front_th1=0.1749639/2.*RI;   // thickness of the middle part 
G4double BaF_C_front_th2=0.0977372/2.*RI; // thickness of the left part 
G4double BaF_C_front_th3=0.028076016/2.*RI; // thickness of the right part 

// more precise
BaF_C_front_th1=0.087482477*RI;
BaF_C_front_th2=0.5*0.19423*cos(59.7*deg)*RI;
BaF_C_front_th3=0.5*0.16622*sin(9.7*deg)*RI;
//BaF_C_front_th3=

//G4double BaF_C_front_height2=0.163831699*RI;
//G4double BaF_C_front_height1=0.167847344*RI;

G4double BaF_C_front_height1=0.19423*sin(59.7*deg)*RI;
G4double BaF_C_front_height2=0.16622*cos(9.7*deg)*RI;


G4double BaF_C_frontdist=17.*cm; //BaF_frontdist;

G4double BaF_C_ratio=(1.+2.*BaF_C_length/BaF_C_frontdist);

G4double BaF_C_back_th1=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_th1;
G4double BaF_C_back_th2=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_th2;
G4double BaF_C_back_th3=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_th3;

G4double BaF_C_back_height1=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_height1;
G4double BaF_C_back_height2=(1.+2.*BaF_C_length/BaF_C_frontdist)*BaF_C_front_height2;


// definition of parts of C irregular hexagon
// here we assume some basic properties of crystal
// prolonged edges will meet at 0,0,0 
// at BaF_C_frontdist-half of the length of crystal (this can be fixed)

// middle part of the crystal C dimensions
//    /----\
//   /      \
//  /________\

G4double shift_C=0.9*cm;

G4double frontshift_C=shift_C;

//shift_C*(1.+2.*BaF_C_length/BaF_C_frontdist);

G4double mid_pDz=BaF_C_length;
G4double mid_pTheta=0.;
G4double mid_pPhi=0.;
 
G4double mid_pDy1=BaF_C_front_th1;
G4double mid_pDx1=BaF_C_front_height1; 
G4double mid_pDx2=BaF_C_front_height2;
G4double mid_pAlp1=0.; 

G4double mid_pDy2=BaF_C_back_th1;
G4double mid_pDx3=BaF_C_back_height1; 
G4double mid_pDx4=BaF_C_back_height2;
G4double mid_pAlp2=0.;


// shifting focus of the crystal
mid_pPhi=90*deg;
mid_pTheta=(atan(0.5*frontshift_C/(BaF_C_frontdist/2.)))*rad;
G4double mid_trans;

mid_trans=-2*BaF_C_length*tan(mid_pTheta);
// left part of the crystal C dimensions

printf("mid_trans= %lf \n",(mid_trans+1*cm)/(1*cm));
printf("BaF_C_length*tan(mid_pTheta)=%lf \n", BaF_C_length*tan(mid_pTheta));

//G4double backshift_C=(BaF_C_front_th1-shift_C)


/////////////////////////////////////////////////////////////////

G4double left_pDz1=BaF_C_length;
//G4double left_pTheta1=(atan (BaF_C_front_th1/2./BaF_C_length))*rad;
//G4double left_pTheta1=atan((mid_pDy2-mid_pDy1)/BaF_C_length)*rad;

G4double left_pPhi1=-90*deg;
 
G4double left_pDy1=BaF_C_front_th2;
G4double left_pDx1=0.01; 
G4double left_pDx2=BaF_C_front_height1;
G4double left_pAlp1=0.; 

G4double left_pDy2=BaF_C_back_th2;
G4double left_pDx3=(BaF_C_ratio*0.01); 
G4double left_pDx4=BaF_C_back_height1;
G4double left_pAlp2=0.;

// withiut shift
//G4double left_pTheta1=atan(0.5*(left_pDy2-left_pDy1+mid_pDy2-mid_pDy1)/BaF_C_length)*rad;
//G4double left_trans=-0.5*(left_pDy2+left_pDy1+mid_pDy2+mid_pDy1);

// with shift
G4double left_pTheta1=atan(-0.5*(left_pDy1-left_pDy2+mid_pDy1-mid_pDy2-mid_trans)/BaF_C_length)*rad;
G4double left_trans=-0.5*(left_pDy2+left_pDy1+mid_pDy2+mid_pDy1)-mid_trans;


///// ADD right piece

G4double right_pDz1=BaF_C_length;
G4double right_pPhi1=90.*deg;
 
G4double right_pDy1=BaF_C_front_th3;
G4double right_pDx1=BaF_C_front_height2; 
G4double right_pDx2=0.01;
G4double right_pAlp1=0.; 

G4double right_pDy2=BaF_C_back_th3;
G4double right_pDx3=BaF_C_back_height2; 
G4double right_pDx4=(BaF_C_ratio*0.01);
G4double right_pAlp2=0.;

//G4double right_pTheta1=atan(0.5*(right_pDy2-right_pDy1+mid_pDy2-mid_pDy1)/BaF_C_length)*rad;
//G4double right_trans=0.5*(right_pDy2+right_pDy1+mid_pDy2+mid_pDy1);

// with shift
G4double right_pTheta1=atan(0.5*(right_pDy2-right_pDy1+mid_pDy2-mid_pDy1-mid_trans)/BaF_C_length)*rad;
G4double right_trans=0.5*(right_pDy2+right_pDy1+mid_pDy2+mid_pDy1)-mid_trans;;

///// basic shapes at 0,0,0

G4Trap * BaF_C1_crystal=	new G4Trap("BaF_C1_crystal",mid_pDz, mid_pTheta,mid_pPhi, mid_pDy1,mid_pDx1, mid_pDx2,mid_pAlp1,mid_pDy2,mid_pDx3, mid_pDx4,mid_pAlp2);
//printf("BaF_tube=	new G4Trap( First Initialized \n");

G4Trap * BaF_C2_crystal=	new G4Trap("BaF_C2_crystal",left_pDz1, left_pTheta1,left_pPhi1, left_pDy1,left_pDx1, left_pDx2,left_pAlp1,left_pDy2,left_pDx3, left_pDx4,left_pAlp2);
//printf("BaF_tube2=	new G4Trap( First Initialized \n");

G4Trap * BaF_C3_crystal=	new G4Trap("BaF_C3_crystal",right_pDz1, right_pTheta1,right_pPhi1, right_pDy1,right_pDx1, right_pDx2,right_pAlp1,right_pDy2,right_pDx3, right_pDx4,right_pAlp2);
//printf("BaF_tube3=	new G4Trap( First Initialized \n");

G4UnionSolid *BaF_C_crystal_1=new G4UnionSolid("BaF_C_crystal_1",BaF_C1_crystal,BaF_C2_crystal,0,G4Vector3D(0,0+left_trans+mid_trans,0));

//G4UnionSolid *BaF_C_crystal=new G4UnionSolid("BaF_C_crystal",BaF_C_crystal_1,BaF_C3_crystal,0,G4Vector3D(0,0+right_trans+mid_trans,0));

//BaF_frontdist=17.*cm;

*/
