
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



#ifdef G4ANALYSIS_USE

#include <fstream>

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include <AIDA/AIDA.h>

#include "DANCEAnalysisManager.hh"

DANCEAnalysisManager* DANCEAnalysisManager::instance = 0;

DANCEAnalysisManager::DANCEAnalysisManager()
:analysisFactory(0), hFactory(0), tFactory(0)
{
  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory)
  {
    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create("DANCE.aida","xml",false,true,"compress=yes");
    hFactory = analysisFactory->createHistogramFactory(*tree);
    tFactory = analysisFactory->createTupleFactory(*tree);
    delete treeFactory; // Will not delete the ITree.
  }
}

DANCEAnalysisManager::~DANCEAnalysisManager()
{
  if (analysisFactory)
  {
    if (!tree->commit()) G4cout << "Commit failed: no AIDA file produced!" << G4endl;
    delete tree;
    delete tFactory;
    delete hFactory;
    G4cout << "Warning: Geant4 will NOT exit unless you close the JAS-AIDA window." << G4endl;
    delete analysisFactory;
  }
}
IHistogramFactory* DANCEAnalysisManager::getHistogramFactory()
{
  return hFactory;
}
ITupleFactory* DANCEAnalysisManager::getTupleFactory()
{
  return tFactory;
}
IPlotter* DANCEAnalysisManager::createPlotter()
{
  if (analysisFactory)
  {
    IPlotterFactory* pf = analysisFactory->createPlotterFactory(0,0);
    if (pf) return pf->create("Plotter");
  }
  return 0;
}

DANCEAnalysisManager* DANCEAnalysisManager::getInstance()
{
  if (instance == 0) instance = new DANCEAnalysisManager();
  return instance;
}

void DANCEAnalysisManager::dispose()
{
  if (instance != 0)
  {
    delete instance;
    instance = 0;
  }
}

#endif // G4ANALYSIS_USE

