/*==============================================================

 waveclimate.cpp

Objects storing wave climate PDFs or data with functions
that can sample to get offshore waves to feed to the 
coastline object.

Developed by:
Martin D. Hurst
Andrew Barkwith
Michael A. Ellis
Christopher W. Thomas
A. Brad Murray

Martin D. Hurst, British Geological Survey, June 2013

Copyright (C) 2015, Martin Hurst

Developer can be contacted:
mhurst@bgs.ac.uk

Martin D. Hurst
British Geological Survey,
Environmental Science Centre,
Nicker Hill,
Keyworth,
Nottingham,
UK,
NG12 5GG

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=====================================================================
*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include "waveclimate.hpp"

using namespace std;

#ifndef waveclimate_CPP
#define wavecliamte_CPP


//=====================================================================
//
// Wave Object
//
//---------------------------------------------------------------------

	/*********************************************\
	| WaveClimate Object Initialisation Functions |
	\*********************************************/

	
void Wave::Initialise()
{
	//cout << "initialised an empty SingleWaveClimate object" << endl;
}
void Wave::Initialise( double WvPrd, double WvHgt, double WvDir )
{
	Dir = WvDir;
	Height = WvHgt;
	Period = WvPrd;
}

void Wave::AssignWaveDirection(double WvDir) { Dir = WvDir; }
void Wave::AssignWaveHeight(double WvHgt) { Height = WvHgt; }
void Wave::AssignWavePeriod(double WvPrd) { Period = WvPrd; }
	
double Wave::Get_WavePeriod() {return Period;}	//get Wave Period	
double Wave::Get_WaveDirection() {return Dir;}	//get Wave Direction
double Wave::Get_WaveHeight() {return Height;}	//get Wave Height
	
//=====================================================================
//
// Single Wave Climate Object
//
//---------------------------------------------------------------------

	/*********************************************\
	| WaveClimate Object Initialisation Functions |
	\*********************************************/

	
void SingleWaveClimate::Initialise()
{
	cout 	<< "SingleWaveClimate.Initialise: Error, "
			<< "initialised an empty SingleWaveClimate object" << endl;
	exit(EXIT_FAILURE);
}

void SingleWaveClimate::Initialise(double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight)
{
	Period = OffshoreWavePeriod;
	Dir = OffshoreWaveDirection;
	Height = OffshoreWaveHeight;	
}

double SingleWaveClimate::Get_WavePeriod() {return Period;}	//get Wave Period	
double SingleWaveClimate::Get_WaveDirection() {return Dir;}	//get Wave Direction
double SingleWaveClimate::Get_WaveHeight() {return Height;}	//get Wave Height
	
////=====================================================================
////
//// PDF Wave Climate Object
////
////---------------------------------------------------------------------

void GaussianWaveClimate::Initialise()
{
	cout 	<< "GaussianWaveClimate.Initialise: Warning, "
			<< "initialised an empty GaussianWaveClimate object,"
			<< "this is only OK if it is one of the two components of a bimodal wave climate!" << endl;
	//exit(EXIT_FAILURE);
}

void GaussianWaveClimate::Initialise(double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveDirection, double OffshoreStDWaveDirection,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight)
{
	Period_Mean = OffshoreMeanWavePeriod;
	Period_StD = OffshoreStDWavePeriod;
	Dir_Mean = OffshoreMeanWaveDirection;
	Dir_StD = OffshoreStDWaveDirection;
	Height_Mean = OffshoreMeanWaveHeight;
	Height_StD = OffshoreStDWaveHeight;
	//srand(time(NULL));
	srand(1);
}

Wave GaussianWaveClimate::Get_Wave()
{
	/* Get_Wave returns a wave object containing wave height, direction
	and period sampled from the RealWaveClimate data sequentially. If the
	end of the wave climate object is encountered then the sampling starts
	again from the begininning of the file.
	
	Get wave properties based on random sampling from a normal distribution
	described by a mean and standard deviation of offshore wave conditions
	
	Random normal distribution from Box-Muller transform 
	sqrt(-2*log((double)rand()/RAND_MAX))*sin(2*PI*((double)rand()/RAND_MAX)) */
		
 	//initialise a temprary wave object
	Wave TempWave = Wave();
	
	double WavePeriod, OffshoreWaveHeight, OffshoreWaveDirection, rand1, rand2;
		
	//Get two random numbers and generate wave data
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	WavePeriod = Period_Mean + Period_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	OffshoreWaveHeight = Height_Mean + Height_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	OffshoreWaveDirection = Dir_Mean + Dir_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
			
	//assign wave properties from wave climate
	TempWave.AssignWaveDirection(OffshoreWaveDirection);
	TempWave.AssignWaveHeight(OffshoreWaveHeight);
	TempWave.AssignWavePeriod(WavePeriod);
	
	//return the wave object
	return TempWave;
}

//=====================================================================

//=====================================================================
//
// Real Wave Climate Object
//
//---------------------------------------------------------------------

	/*********************************************\
	| WaveClimate Object Initialisation Functions |
	\*********************************************/

	
void RealWaveClimate::Initialise()
{
	cout 	<< "RealWaveClimate.Initialise: Error, "
			<< "initialised an empty RealWaveClimate object" << endl
			<< "RealWaveClimate requires an input *.txt containing wave data"<< endl;
	exit(EXIT_FAILURE);
}

void RealWaveClimate::Initialise(string WaveFileName)
{
	/*	Reads a wave data from a text file 
		File format is
								Headers
								Dir[0] | Period[0] | Height[0]
								Dir[1] | Period[1] | Height[1]	
								 ...   |    ...    |   ...		*/

	//declarations
	float T, D, H;
	string Temp, line;
	
	NoWaves = 0;
	WaveCounter = 0;
	
	ifstream ReadWaveFile;
	ReadWaveFile.open(WaveFileName.c_str());
	if (ReadWaveFile.is_open()) ReadWaveFile >> Temp >> Temp >> Temp;
	else
	{
		cout << "Coastline.ReadCoast: Error, the file " << WaveFileName << " has not been read correctly." << endl;
		exit(EXIT_FAILURE);
	}

	//Loop through file and read in wave data
	while (!ReadWaveFile.eof())
	{
		ReadWaveFile >> D >> T >> H; 
		Dir.push_back(D);
		Period.push_back(T);
		Height.push_back(H);
		NoWaves += 1;
	}
}

Wave RealWaveClimate::Get_Wave()
{
	/* Get_Wave returns a wave object containing wave height, direction
	and period sampled from the RealWaveClimate data sequentially. If the
	end of the wave climate object is encountered then the sampling starts
	again from the begininning of the file. */
 	
 	//initialise a temprary wave object
	Wave TempWave = Wave();
	
	//assign wave properties from wave climate
	TempWave.AssignWaveDirection(Dir[WaveCounter]);
	TempWave.AssignWaveHeight(Height[WaveCounter]);
	TempWave.AssignWavePeriod(Period[WaveCounter]);
	
	//update wave counter and reset if necessary
	WaveCounter += 1;
	if (WaveCounter >= NoWaves) WaveCounter = 0;
	
	//return the wave object
	return TempWave;
}

//=====================================================================

//=====================================================================
//
// U and A Wave Climate Object
//
//---------------------------------------------------------------------

	/*********************************************\
	| WaveClimate Object Initialisation Functions |
	\*********************************************/

	
void UAWaveClimate::Initialise()
{
	cout 	<< "UAWaveClimate.Initialise: Error, "
			<< "initialised an empty UAWaveClimate object" << endl
			<< "UAWaveClimate requires arguments U and A" << endl;
	exit(EXIT_FAILURE);
}

void UAWaveClimate::Initialise(	double input_U, double input_A, double Trend, 
											double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
											double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight)
{
	CoastTrend = Trend-90.;
	U = input_U;
	A = input_A;
	
	Period_Mean = OffshoreMeanWavePeriod;
	Period_StD = OffshoreStDWavePeriod;
	Height_Mean = OffshoreMeanWaveHeight;
	Height_StD = OffshoreStDWaveHeight;
	
	//srand(time(NULL));
	srand(1);
}

Wave UAWaveClimate::Get_Wave()
{
	/* Get_Wave returns a wave object containing wave height, direction
	and period sampled from the UAWaveClimate paramters.
	
	Get height and period wave properties based on random sampling from 
	a normal distribution described by a mean and standard deviation of 
	offshore wave conditions
	
	Random normal distribution from Box-Muller transform 
	sqrt(-2*log((double)rand()/RAND_MAX))*sin(2*PI*((double)rand()/RAND_MAX))
	
	Wave direction sampled following the U and A approach of Ashton et al
	2006. U is fraction of high angle waves and A is fraction of waves coming
	from the "left" relative to CoastTrend. */
	
 	//initialise a temprary wave object
	Wave TempWave = Wave();
	
	double WavePeriod, OffshoreWaveHeight, OffshoreWaveDirection, rand1, rand2;
		
	//Get two random numbers and generate wave data
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	WavePeriod = Period_Mean + Period_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	OffshoreWaveHeight = Height_Mean + Height_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	
	//Get wave Direction from U and A following Ashton & Murray (2006) CEM
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	OffshoreWaveDirection = rand1*45.;
	rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
	if (rand1 < U) OffshoreWaveDirection += 45.;
	if (rand2 < A) OffshoreWaveDirection *= -1;
	OffshoreWaveDirection += CoastTrend;
			
	//assign wave properties from wave climate
	TempWave.AssignWaveDirection(OffshoreWaveDirection);
	TempWave.AssignWaveHeight(OffshoreWaveHeight);
	TempWave.AssignWavePeriod(WavePeriod);
	
	//return the wave object
	return TempWave;
}

//=====================================================================



//=====================================================================
//
// Bimodal Wave Climate Object
//
//---------------------------------------------------------------------

	/*********************************************\
	| WaveClimate Object Initialisation Functions |
	\*********************************************/
	
void BimodalWaveClimate::Initialise()
{
	cout 	<< "BimodalWaveClimate.Initialise: Error, "
				<< "initialised an empty BimodalWaveClimate object" << endl
				<< "BimodalWaveClimate requires arguments..." << endl
				<< "\t- FractionWaveDirection1" << endl
				<< "\t- MeanWaveDirection1, StadardDeviationWaveDirection1" << endl
				<< "\t- MeanWavePeriod1, StadardDeviationWavePeriod1" << endl
				<< "\t- MeanWaveHeight1, StadardDeviationWaveHeight1" << endl
				<< "\t- MeanWaveDirection2, StadardDeviationWaveDirection2" << endl
				<< "\t- MeanWavePeriod2, StadardDeviationWavePeriod2" << endl
				<< "\t- MeanWaveHeight2, StadardDeviationWaveHeight2" << endl << endl;
	exit(EXIT_FAILURE);
}

void BimodalWaveClimate::Initialise(	double FractionWaveDirection1,	
											double OffshoreMeanWaveDirection1, double OffshoreStDWaveDirection1,
											double OffshoreMeanWavePeriod1, double OffshoreStDWavePeriod1,
											double OffshoreMeanWaveHeight1, double OffshoreStDWaveHeight1,
											double OffshoreMeanWaveDirection2, double OffshoreStDWaveDirection2,
											double OffshoreMeanWavePeriod2, double OffshoreStDWavePeriod2,
											double OffshoreMeanWaveHeight2, double OffshoreStDWaveHeight2)
{
	cout << "BimodalWaveClimate.Initialise: Initialising a BimodalWaveClimate object" << endl;
	
	WaveMode1 = GaussianWaveClimate(OffshoreMeanWaveDirection1, OffshoreStDWaveDirection1,
											OffshoreMeanWavePeriod1, OffshoreStDWavePeriod1,
											OffshoreMeanWaveHeight1, OffshoreStDWaveHeight1);
											
	WaveMode2 = GaussianWaveClimate(OffshoreMeanWaveDirection2, OffshoreStDWaveDirection2,
											OffshoreMeanWavePeriod2, OffshoreStDWavePeriod2,
											OffshoreMeanWaveHeight2, OffshoreStDWaveHeight2);
	
	//srand(time(NULL));
	srand(1);
}

Wave BimodalWaveClimate::Get_Wave()
{
	/* Get_Wave returns a wave object containing wave height, direction
	and period sampled from the UAWaveClimate paramters.
	
	Get height and period wave properties based on random sampling from 
	a normal distribution described by a mean and standard deviation of 
	offshore wave conditions
	
	Random normal distribution from Box-Muller transform 
	sqrt(-2*log((double)rand()/RAND_MAX))*sin(2*PI*((double)rand()/RAND_MAX))
	
	Wave direction sampled from two guassian distributions weighted according
	to Dir_Fraction1 */
	
 	//initialise a temprary wave object
	Wave TempWave = Wave();
	
	double WavePeriod, OffshoreWaveHeight, OffshoreWaveDirection, rand1, rand2;

	//Get a random number to determine which wave direction mode to use
	rand1 = (double)rand()/RAND_MAX;
	
	//for WaveMode #1
	if (rand1 < 0.5)
	{
		//get a wave direction
		rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
		OffshoreWaveDirection = WaveMode1.Dir_Mean + WaveMode1.Dir_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	
		//Get wave Period
		rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
		WavePeriod = WaveMode1.Period_Mean + WaveMode1.Period_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
		
		//get wave hieght
		rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
		OffshoreWaveHeight = WaveMode1.Height_Mean + WaveMode1.Height_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	}
	
	//for WaveMode #2
	else
	{
		//get a wave direction
		rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
		OffshoreWaveDirection = WaveMode2.Dir_Mean + WaveMode2.Dir_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	
		//Get wave Period
		rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
		WavePeriod = WaveMode2.Period_Mean + WaveMode2.Period_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
		
		//get wave hieght
		rand1 = (double)rand()/RAND_MAX; rand2 = (double)rand()/RAND_MAX;
		OffshoreWaveHeight = WaveMode2.Height_Mean + WaveMode2.Height_StD*sqrt(-2.*log(rand1))*cos(2.*M_PI*(rand2));
	}

	//assign wave properties from wave climate
	TempWave.AssignWaveDirection(OffshoreWaveDirection);
	TempWave.AssignWaveHeight(OffshoreWaveHeight);
	TempWave.AssignWavePeriod(WavePeriod);
	
	//return the wave object
	return TempWave;
}

#endif 
