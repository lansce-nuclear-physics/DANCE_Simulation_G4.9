
#////////////////////////////////////////////////////////////////////////
#//                                                                    //
#//   Software Name: DANCE Data Acquisition and Analysis Package       //
#//     Subpackage: DANCE_Simulation_G4.9                              //
#//   Identifying Number: C18105                                       // 
#//                                                                    //
#////////////////////////////////////////////////////////////////////////
#//                                                                    //
#//                                                                    //
#// Copyright 2019.                                                    //
#// Triad National Security, LLC. All rights reserved.                 //
#//                                                                    //
#//                                                                    //
#//                                                                    //
#// This program was produced under U.S. Government contract           //
#// 89233218CNA000001 for Los Alamos National Laboratory               //
#// (LANL), which is operated by Triad National Security, LLC          //
#// for the U.S. Department of Energy/National Nuclear Security        //
#// Administration. All rights in the program are reserved by          //
#// Triad National Security, LLC, and the U.S. Department of           //
#// Energy/National Nuclear Security Administration. The Government    //
#// is granted for itself and others acting on its behalf a            //
#// nonexclusive, paid-up, irrevocable worldwide license in this       //
#// material to reproduce, prepare derivative works, distribute        //
#// copies to the public, perform publicly and display publicly,       //
#// and to permit others to do so.                                     //
#//                                                                    //
#// This is open source software; you can redistribute it and/or       //
#// modify it under the terms of the GPLv2 License. If software        //
#// is modified to produce derivative works, such modified             //
#// software should be clearly marked, so as not to confuse it         //
#// with the version available from LANL. Full text of the GPLv2       //
#// License can be found in the License file of the repository         //
#// (GPLv2.0_License.txt).                                             //
#//                                                                    //
#////////////////////////////////////////////////////////////////////////


# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := DICEBOX_DANCE
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

ifdef G4ANALYSIS_USE
   CPPFLAGS += `aida-config --include`
   LDLIBS += `aida-config --lib`
endif

########################### ROOT #################################

ifdef ROOTSYS
ifndef G4UI_USE_ROOT
  ROOTCPPFLAGS   = $(shell $(ROOTSYS)/bin/root-config --cflags)
  CPPFLAGS      += -DG4ANALYSIS_USE_ROOT $(ROOTCPPFLAGS)
  ROOTLIBS       = $(shell $(ROOTSYS)/bin/root-config --nonew --glibs)
  ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
  ROOTLIBS      += -lMinuit -lHtml
  LDLIBS        += $(ROOTLIBS)
endif
endif
