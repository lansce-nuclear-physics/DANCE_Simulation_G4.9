
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



#ifndef DANCEAnalysisManager_h
#define DANCEAnalysisManager_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

using namespace AIDA;

class AIDA::IAnalysisFactory;
class AIDA::ITree;
class AIDA::IHistogramFactory;
class AIDA::ITupleFactory;
class AIDA::IPlotter;

class G4Track;

class DANCEAnalysisManager {
public:

  virtual ~DANCEAnalysisManager();
  static DANCEAnalysisManager* getInstance();
  static void dispose();

  IHistogramFactory* getHistogramFactory();
  ITupleFactory* getTupleFactory();
  IPlotter* createPlotter();

private:

  DANCEAnalysisManager();
  static DANCEAnalysisManager* instance;

  IAnalysisFactory* analysisFactory;
  IHistogramFactory* hFactory;
  ITupleFactory* tFactory;
  ITree* tree;


};
#endif // G4ANALYSIS_USE

#endif
