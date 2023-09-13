
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


#include "DANCEGlobal.hh"
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

DANCEGlobal::DANCEGlobal(){
}
DANCEGlobal::~DANCEGlobal(){
}

DANCEGlobal *DANCEGlobal::Instance()
{

	static DANCEGlobal *pinstance=new DANCEGlobal();
	return pinstance;

}
void DANCEGlobal::Initialize(){
// read in MasterInput.txt

ifstream *in;

in=new ifstream("MasterInput.txt");

if(!in->is_open()) {
	cout << "-------------------------------------------------" << endl;
	cout << "MasterInput.txt file not found, exiting ! " <<	endl;
	cout << "-------------------------------------------------" << endl;
	exit(0);
}

G4String gstring;
bool eof_flag;

while(!in->eof()){

	(*in) >> gstring;
	eof_flag=in->eof();


	if(!eof_flag){

		if(gstring=="BinaryCrystalFile"){
			BinaryCrystalFile=CheckForFlags(in);
			cout << "Setting: " << gstring << "	to: " << BinaryCrystalFile << endl;
			BinaryCrystalFile_Name=gstring;
			
			//binary_crystal_file=new ofstream(BinaryCrystalFile_Name.c_str(),ios::out|ios::binary);
			
		};
		if(gstring=="BinaryClusterFile"){
			BinaryClusterFile=CheckForFlags(in);
			cout << "Setting: " << gstring << "	to: " << BinaryClusterFile << endl;
		};
		if(gstring=="AddPrimaryGammas") {		
			AddPrimaryGammas=CheckForFlags(in);
			cout << "Setting: " << gstring << "	to: " << AddPrimaryGammas << endl;
		};
		if(gstring=="AddTimeOfGammas") {		
			AddTimeOfGammas=CheckForFlags(in);
			cout << "Setting: " << gstring << "	to: " << AddTimeOfGammas << endl;
		};
		if(gstring=="RootFile") {		
			RootFile=CheckForFlags(in);
			cout << "Setting: " << gstring << "	to: " << RootFile << endl;
		};

	}
}

exit(0);

}

void DANCEGlobal::Init(){
ifstream in2;
in2.open("config/SymmetryID.dat");
if(!in2) {
	cout << "-------------------------------------------------" << endl;
	cout << "config/SymmetryID.dat file not found, exiting ! " <<	endl;
	cout << "-------------------------------------------------" << endl;
	exit(0);
}

int temp;

for(int i=0;i<81;i++){

	in2 >> temp >> SymID[i];
//	cout << temp << "	" << SymID[i] << endl;
	if(temp!=i) {
		cout << " WRONG /config/SymmetryID.dat FILE !! " << endl;
		break;  
	}
}

in2.close();
// read in MasterInput.txt

ifstream in(MasterInputFile);
if(!in) {
	cout << "-------------------------------------------------" << endl;
	cout << "MasterInput.txt file not found, exiting ! " <<	endl;
	cout << "-------------------------------------------------" << endl;
	exit(0);
}

G4String gstring;
G4String flag;
G4String flag2;
bool eof_flag;

ShowCrystals[0]=0;
ShowCrystals[1]=162;

Q_energy=-10.;

EndOfExternalInput=false;

while(!in.eof()){
	in >> gstring;
	eof_flag=in.eof();
		
	if(!eof_flag){
//		cout << gstring << "	"  << endl;
		if(gstring=="ExternalInput") {		
			
			in >> flag;
			

			if(flag=="true") ExternalInput=true;
			if(flag=="false") ExternalInput=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool " << gstring << "	to " << flag << endl;
		}

		if(gstring=="PVCShell") {		
			
			in >> flag;
			

			if(flag=="true") PVCShell=true;
			if(flag=="false") PVCShell=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool " << gstring << "	to " << flag << endl;
		}

		if(gstring=="ExternalInputFile") {		
			
			in >> flag;			
			ExternalInputFile=flag;			
			cout << "Setting bool " << gstring << "	to " << ExternalInputFile << endl;
		}
		if(gstring=="MacroInput") {		
			
			in >> flag;
			

			if(flag=="true") MacroInput=true;
			if(flag=="false") MacroInput=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool " << gstring << "	to " << flag << endl;
		}
		if(gstring=="MacroInputFile") {		
			
			in >> flag;			
			MacroInputFile=flag;			
			cout << "Setting bool " << gstring << "	to " << flag << endl;
		}
		
		if(gstring=="BinaryCrystalFile") {		
			
			in >> flag >> BinaryCrystalFile_Name;
			
			if(flag=="true") {
				BinaryCrystalFile=true;
				//in >> flag;				
				//BinaryCrystalFile_Name=flag;
				
				binary_crystal_file=new ofstream(BinaryCrystalFile_Name.c_str(),ios::out|ios::binary);
				
			}
			if(flag=="false") BinaryCrystalFile=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool BinaryCrystalFile=" << flag << endl;
//			cout << BinaryCrystalFile << endl;		
		}
		if(gstring=="BinaryClusterFile") {		
			
			in >> flag;
			

			if(flag=="true") {
				BinaryClusterFile=true;
				in >> flag;				
				BinaryClusterFile_Name=flag;
			}
			if(flag=="false") BinaryClusterFile=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool BinaryClusterFile=" << flag << endl;
		}
		if(gstring=="AddPrimaryGammas") {		
			
			in >> flag;
			

			if(flag=="true") AddPrimaryGammas=true;
			if(flag=="false") AddPrimaryGammas=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool AddPrimaryGammas=" << flag << endl;
		}
		if(gstring=="AddTimeOfGammas") {		
			
			in >> flag;
			

			if(flag=="true") AddTimeOfGammas=true;
			if(flag=="false") AddTimeOfGammas=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool AddTimeOfGammas=" << flag << endl;
		}
		if(gstring=="RootFile") {		
			
			in >> flag;
			

			if(flag=="true") RootFile=true;
			if(flag=="false") RootFile=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool RootFile=" << flag << endl;
		}
		if(gstring=="RootFileName") {		
			
			in >> flag;
			

			RootFileName=flag;
			cout << "Setting bool RootFileName=" << flag << endl;
		}
		if(gstring=="RootAutoSave") {		
			
			in >> flag;
			

			if(AutoSave=atoi(flag)){
				cout << "Setting int RootAutoSave=" << AutoSave << endl;
			}
			else {
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << flag << "  does not seem to be integer. Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
	
			}
		}
		if(gstring=="ClusterMGatedHists") {		
			
			in >> flag;
			

			if(flag=="true") ClusterMGatedHists=true;
			if(flag=="false") ClusterMGatedHists=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool ClusterMGatedHists=" << flag << endl;
		}
		if(gstring=="CrystalMGatedHists") {		
			
			in >> flag;
			

			if(flag=="true") CrystalMGatedHists=true;
			if(flag=="false") CrystalMGatedHists=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool CrystalMGatedHists=" << flag << endl;
		}

		if(gstring=="SingleGateEcrystal_1" || gstring=="SingleGateEcrystal_2" ||gstring=="SingleGateEcrystal_3") {		
			
			double E,dE;
			double multiply=1.;
			int ID_Gate;
			
			if(gstring=="SingleGateEcrystal_1") ID_Gate=0;
			if(gstring=="SingleGateEcrystal_2") ID_Gate=1;
			if(gstring=="SingleGateEcrystal_3") ID_Gate=2;
			
			in >> E >> flag;
			
			if(E==0) {
				SingleGateECrystal_E[ID_Gate]=0.;
				
			}
			else{
			if(flag=="MeV") multiply=1.;
			if(flag=="keV") multiply=1e-3;
			
			if(flag!="MeV" && flag!="keV"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << E << " " << flag << endl;
				cout << "Unrecognized units. Please use \"keV\" or \"MeV\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);						
			}
			
			SingleGateECrystal_E[ID_Gate]=E*multiply;
			
			in >> dE >> flag2;

			if(flag2=="MeV") multiply=1.;
			if(flag2=="keV") multiply=1e-3;
			
			if(flag2!="MeV" && flag2!="keV"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set width of the gate" << gstring << "  to  " << dE << " " << flag2 << endl;
				cout << "Unrecognized units. Please use \"keV\" or \"MeV\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);						
			}
			SingleGateECrystal_dE[ID_Gate]=dE*multiply;
						
			cout << "Setting the CrystalEGate No." << ID_Gate << " E=" << SingleGateECrystal_E[ID_Gate] << " MeV +- " << SingleGateECrystal_dE[ID_Gate]
			<< " MeV" << endl;
			}
		}

		if(gstring=="SingleGateEcluster_1" || gstring=="SingleGateEcluster_2" ||gstring=="SingleGateEcluster_3") {		
			
			double E,dE;
			double multiply=1.;
			int ID_Gate;
			
			if(gstring=="SingleGateEcluster_1") ID_Gate=0;
			if(gstring=="SingleGateEcluster_2") ID_Gate=1;
			if(gstring=="SingleGateEcluster_3") ID_Gate=2;
			
			in >> E >> flag;
			
			if(E==0) {
				SingleGateECluster_E[ID_Gate]=0.;				
			}
			
			else{
			if(flag=="MeV") multiply=1.;
			if(flag=="keV") multiply=1e-3;
			
			if(flag!="MeV" && flag!="keV"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << E << " " << flag << endl;
				cout << "Unrecognized units. Please use \"keV\" or \"MeV\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);						
			}
			
			SingleGateECluster_E[ID_Gate]=E*multiply;
			
			in >> dE >> flag2;

			if(flag2=="MeV") multiply=1.;
			if(flag2=="keV") multiply=1e-3;
			
			if(flag2!="MeV" && flag2!="keV"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set width of the gate" << gstring << "  to  " << dE << " " << flag2 << endl;
				cout << "Unrecognized units. Please use \"keV\" or \"MeV\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);						
			}
			SingleGateECluster_dE[ID_Gate]=dE*multiply;
						
			cout << "Setting the ClusterEGate No." << ID_Gate << " E=" << SingleGateECluster_E[ID_Gate] << " MeV  +- " <<
			SingleGateECluster_dE[ID_Gate] << " MeV" << endl;
			}
		}
		
		if(gstring=="Comments") {
		
			break;
		
		}		

		if(gstring=="SupportingSphere") {		
			
			in >> flag >> flag2;
			

			if(flag=="true") SupportingSphere[0]=true;
			if(flag=="false") SupportingSphere[0]=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}

			if(flag2=="true") SupportingSphere[1]=true;
			if(flag2=="false") SupportingSphere[1]=false;
			if(flag2!="false" && flag2!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag2 << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << gstring << "=" << flag << "	vis= " << flag2 << endl;
		}

		if(gstring=="LiHShell") {		
			
			in >> flag >> flag2;
			

			if(flag=="true") LiHShell[0]=true;
			if(flag=="false") LiHShell[0]=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}

			if(flag2=="true") LiHShell[1]=true;
			if(flag2=="false") LiHShell[1]=false;
			if(flag2!="false" && flag2!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag2 << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "LiHShell=" << flag << "	vis= " << flag2 << endl;
		}
		if(gstring=="LiHDensity") {		
			
			in >> LiHDensity;
			

			cout << "LiH Density:  " << LiHDensity << endl;
		}
		if(gstring=="LiHInnerRadius") {		
			
			in >> LiHInnerRadius;
			

			cout << "LiH Inner Radius:  " << LiHInnerRadius << endl;
		}
		if(gstring=="LiHOuterRadius") {		
			
			in >> LiHOuterRadius;
			

			cout << "LiH Outer Radius:  " << LiHOuterRadius << endl;
		}
		
		if(gstring=="ShowCrystals") {		
			
			in >> ShowCrystals[0] >> ShowCrystals[1];
			

			cout << "Visible detectors for Crystal_ID:  " << ShowCrystals[0] << " - " << ShowCrystals[1] << endl;
		}
		if(gstring=="PMLense") {		
			
			in >> flag >> flag2;
			

			if(flag=="true") PMLense[0]=true;
			if(flag=="false") PMLense[0]=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}

			if(flag2=="true") PMLense[1]=true;
			if(flag2=="false") PMLense[1]=false;
			if(flag2!="false" && flag2!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag2 << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "PMLense=" << flag << "	vis= " << flag2 << endl;
			
		}
		if(gstring=="E_threshold") {		
			
			double multiply;
			
			in >> E_threshold >> flag2;
			
			if(flag2=="MeV") multiply=1.;
			if(flag2=="keV") multiply=1e-3;
			
			if(flag2!="MeV" && flag2!="keV"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set width of the gate" << gstring << "  to  " << E_threshold << " " << flag2 << endl;
				cout << "Unrecognized units. Please use \"keV\" or \"MeV\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);						
			}
			E_threshold=E_threshold*multiply;
			cout << "E_threshold:  " << E_threshold << " MeV" << endl;
		}
		if(gstring=="BeamPipeVacuum") {		
			
			in >> flag;
			
			if(flag=="true") BeamPipeVacuum=true;
			if(flag=="false") BeamPipeVacuum=false;
			if(flag!="false" && flag!="true"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set " << gstring << "  to  " << flag << endl;
				cout << "Unrecognized flag. Please use lower case \"true\" or \"false\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);			
			}
			cout << "Setting bool BeamPipeVacuum=" << flag << endl;
//			cout << BinaryCrystalFile << endl;		
		}
		if(gstring=="CrystalDistance") {		
			
			double multiply;
			
			in >> CrystalDistance >> flag2;
			
			if(flag2=="cm") multiply=1.;
			if(flag2=="m") multiply=10.;
			
			if(flag2!="cm" && flag2!="mm"){
				cout << "-----------------------------------------------------------------------" << endl;
				cout << "Trying to set width of the gate" << gstring << "  to  " << E_threshold << " " << flag2 << endl;
				cout << "Unrecognized units. Please use \"keV\" or \"MeV\". Exiting" << endl;
				cout << "-----------------------------------------------------------------------" << endl;
				exit(0);						
			}
			CrystalDistance=CrystalDistance*multiply;
			cout << "CrystalDistance=" << CrystalDistance << " cm" << endl;
		}
		if(gstring=="Q_energy") {		
			
			in >> Q_energy >> dQ_energy;
			

			cout << "Q of the reaction:  " << Q_energy << "	dQ gate: " << dQ_energy << endl;
		}

	}
}

in.close();

// read in symmetry ID file

//ifstream in("/home/jandel/geant4/MyGeant/DANCE-NEW/config/SymmetryID.dat");


}

bool DANCEGlobal::CheckForFlags(ifstream *in){

G4String flag1;

(*in) >> flag1;

if(flag1=="true") return true;
else return false;

}
