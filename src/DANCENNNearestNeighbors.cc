
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


//--------------------------------------------------------------------------------
// CLASS:	NNNearestNeighbors
// Purpose:	Use to map a specific electronics channel to a crystal.
// History: 2003-06-10 (jmw) - original.
//--------------------------------------------------------------------------------
#include <cstdio>
#include <iostream>

//#include "include/general/GLGlobalLogging.h"
#include "DANCENNNearestNeighbors.hh"

using namespace std;

//--------------------------------------------------------------------------------
/*!
 * \method	NNNearestNeighbors
 * \brief	Constructor for NNNearestNeighbors class.  Does nothing.
 * \note	
 * \history 2003-06-10 (jmw) original.
 */
//--------------------------------------------------------------------------------
NNNearestNeighbors::NNNearestNeighbors( void )
{
}

//--------------------------------------------------------------------------------
/*!
 * \method	~NNNearestNeighbors
 * \brief	Destructor - Does nothing.
 * \author	Jan M. Wouters
 * \history 2003-06-10 (jmw) original.
 */
//--------------------------------------------------------------------------------
NNNearestNeighbors::~NNNearestNeighbors( )
{
}

NNNearestNeighbors *NNNearestNeighbors::Instance()
{

	static  NNNearestNeighbors *pinstance=new NNNearestNeighbors();
	return pinstance;

}
//--------------------------------------------------------------------------------
/*!
 * \method	MPInitialize
 * \brief	Read in the Dance map relating electronic channel to crystal.  Then
 *			place data into crsytal map.
 * \param	asFilename		- Name of file listing nearest neighbors.
 * \note	
 * \history 2003-06-10 (jmw) original.
 */
//--------------------------------------------------------------------------------
bool NNNearestNeighbors::NNInitialize( const string &asFilename )
{
//    ostringstream   msg;
    int		nn1;
    int		nn2;
    int		nn3;
    int		nn4;
    int		nn5;
    int		nn6;
    int		crystal;
    bool	bRetVal = false;
	
//	asFilename;
	
// Open the file.
	FILE *filePtr = fopen( asFilename.c_str(), "r" );
	
// Error - could not open the file.
	if ( !filePtr )
	{
		bRetVal = false;
		cout << "Error: Could not read in nearest neighbors file: " << asFilename << endl;
	}
	
// Have input file so read in the data.
	else
	{
//		msg << "Nearest neighbor file opened: " << asFilename;
//		gStatusLogObj->LOLogMessage( 'I', msg.str(), "NNNearestNeighbors::NNInitialize" );

		for ( short i = 0; i < NNcMaxCrystals; i++ )
		{
			fscanf( filePtr, "%d %d %d %d %d %d %d", &crystal, &nn1, &nn2, &nn3, &nn4, &nn5, &nn6 );
			if ( crystal == i )
			{
				NNgNeighborMap[ i ][ 0 ] = crystal;
				NNgNeighborMap[ i ][ 1 ] = nn1;
				NNgNeighborMap[ i ][ 2 ] = nn2;
				NNgNeighborMap[ i ][ 3 ] = nn3;
				NNgNeighborMap[ i ][ 4 ] = nn4;
				NNgNeighborMap[ i ][ 5 ] = nn5;
				NNgNeighborMap[ i ][ 6 ] = nn6;
			}
			
// Handle bad record.
			else
			{
				bRetVal = false;
				cout << "Error: Bad record in neighbor file for crystal: " << i << endl;
				fclose( filePtr );
				return( bRetVal );
			}
				
			if ( feof( filePtr ) ) break;
		}
		fclose( filePtr);
        bRetVal = true;
	}	
    
// Debugging - print the input map.
//	if ( gStatusLogObj->LOIsMessageActive( 'F' ) )
//    	NNPrintMap();
    
    return( bRetVal );
}

//--------------------------------------------------------------------------------
/*!
 * \method	NNIsNeighbor
 * \brief	Determines if test crystal is neighbor to reference crystal.
 * \param	pRefCrystal			- Reference crystal.
 * \param	pTestCrystal		- Test crystal.  See if this crystal is
 *								  adjacent to reference crystal.
 * \return	NNcNoNeighbor		- no neighbor
 *			NNcIsNeighbor		- Is Neighbor
 *			NNcError			- Error, one or both crystal ids bad.
 * \note	Start search at index 1 since first index is crystal number.
 * \history 2003-06-10 (jmw) original.
 */
//--------------------------------------------------------------------------------
short NNNearestNeighbors::NNIsNeighbor( short aRefCrystal, short aTestCrystal )
{
        
	if ( aRefCrystal >= 0 && aRefCrystal < NNcMaxCrystals )
	{
		for ( short i = 1; i < NNcMaxNearestNeighbors; i++ )
		{
//            cout << " Compare Crystal: " << pRefCrystal << " ref ( " << i << " ): " 
//                 << NNgNeighborMap[ pRefCrystal ][ i ] << " Test: " << pTestCrystal << endl;
// Break out of test when pTestCrystal is less than neighbor in reference crystal since
// cannot ever be a neighbor after that point.
            if ( aTestCrystal < NNgNeighborMap[ aRefCrystal ][ i ] )
                break;
			if ( aTestCrystal == NNgNeighborMap[ aRefCrystal ][ i ] ) 
				return( NNcIsNeighbor );
		}
		
		return( NNcNoNeighbor );
	}
			
	return( NNcError );
}

//--------------------------------------------------------------------------------
/*!
 * \method	NNPrintMap
 * \brief	Print the nearest neighbor map.
 * \note	
 * \history 2003-12-1 (jmw) original.
 */
//--------------------------------------------------------------------------------
void NNNearestNeighbors::NNPrintMap( void )
{
//	ostringstream 	msg;
	
//    msg << endl;
//    msg << "**** Nearest Neighbor Map ****" << endl;
    
    for ( short i = 0; i < NNcMaxCrystals; i++ )
    {
//        msg << "Crystal: " << i << "  ";
        for ( short j = 0; j < NNcMaxNearestNeighbors; j++ )
        {
//            msg << NNgNeighborMap[ i ][ j ] << ", ";
        }
//        msg << endl;
    }
    
//    gStatusLogObj->LOLogMessage( 'F', msg.str(), "NNNearestNeighbors::NNPrintMap" );
}
