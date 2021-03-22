/*==============================================================

spiral_bay_driver.hpp

A driver function to simulate the evolution of a crenulate or
spiral bay from a straight coastline, forced by a Gaussian
wave climate.

Developed by:
Martin D. Hurst
Andrew Barkwith
Michael A. Ellis
Christopher W. Thomas
A. Brad Murray

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
√MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

==============================================================*/

/** @file spiral_bay_driver.cpp
@author Martin D. Hurst, British Geological Survey
@brief driver file for evolving a spiral bay from a straight coastline
@details This file generates a coastline object as a straightline and
evolves it driven by a Gaussian distribution wave climate. The model
boundaries are fixed to represent headlands or seawalls. Sediment is
allowed to escape around these fixed ends, but none is provided into the
model domain.
@date 28/10/2015
*/

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include "../coastline.hpp"
#include "../waveclimate.hpp"

using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 9)
	{
		cout << "Program needs 6 input arguments:	 \n\t- MeanWavePeriod\n\t- StDWavePeriod"
												<<	"\n\t- MeanWaveHeight\n\t- StDWaveHeight"
												<<	"\n\t- MeanWaveDirection\n\t- StDWaveDirection"
												<<	"\n\t- TimeStep\n\t- EndTime" << endl;
		exit(EXIT_FAILURE);
	}

	//Declare parameter for wave conditions
	double 	OffshoreMeanWavePeriod, OffshoreStDWavePeriod, OffshoreMeanWaveDirection,
			OffshoreStDWaveDirection, OffshoreMeanWaveHeight, OffshoreStDWaveHeight;

	//create output filename based on input conditions
	string FileName = "Spiral";
	string Underscore = "_";
	string Extension = ".xy";
	string arg1 = argv[1];
	string arg2 = argv[2];
	string arg3 = argv[3];
	string arg4 = argv[4];
	string arg5 = argv[5];
	string arg6 = argv[6];
	string WriteCoastFile = FileName+Underscore+arg1+Underscore+arg2+Underscore+arg3+Underscore+arg4+Underscore+arg5+Underscore+arg6+Extension;

	//Read wave parameters from argvs
	OffshoreMeanWavePeriod = atof(argv[1]);
	OffshoreStDWavePeriod = atof(argv[2]);
	OffshoreMeanWaveHeight = atof(argv[3]);
	OffshoreStDWaveHeight = atof(argv[4]);
	OffshoreMeanWaveDirection = atof(argv[5]);
	OffshoreStDWaveDirection = atof(argv[6]);

	//declare time control paramters
	double EndTime, TimeStep, MaxTimeStep, WaveTimeDelta;
	EndTime = atof(argv[8]);
	TimeStep = atof(argv[7]);
	MaxTimeStep = atof(argv[7]);
	WaveTimeDelta = 2.*TimeStep;


	// int EndTime = 20.;						      // End time (years)
	double Time = 0.;							      // Start Time (years)
	double PrintTimeDelta = 36.5/365.;  // how often to print coastline (years)
	double PrintTime = PrintTimeDelta;	// Print time (years)
	double TempTime;

	// double WaveTimeDelta = 0.1;			// Frequency at which to sample new waves (days)
	double GetWaveTime = 0.0;				// Time to get a new wave (days)
	// double TimeStep = 0.05;					// Time step (days)
	// double MaxTimeStep = 0.05;			// Maximum timestep (days)

	//initialise coast as straight line with low amp noise
	int MeanNodeSpacing = 50;
	double CoastLength = 2000;
	double Trend = 140.;

	//boundary conditions are fixed
	int StartBoundary = 2;
	int EndBoundary = 2;

	//initialise wave climates
	GaussianWaveClimate WaveClimate = GaussianWaveClimate(	OffshoreMeanWavePeriod, OffshoreStDWavePeriod,
															OffshoreMeanWaveDirection, OffshoreStDWaveDirection,
															OffshoreMeanWaveHeight, OffshoreStDWaveHeight);

	Wave MyWave = Wave();
	MyWave = WaveClimate.Get_Wave();

	//initialise the coastline as a straight line
	Coastline CoastVector = Coastline(MeanNodeSpacing, CoastLength, Trend, StartBoundary, EndBoundary);

	//int a, b; //to hold random numbers
	int CERCFlag = 1;
	int RefDiffFlag = 1;
	double FluxFraction = 0.0;

	//initiate random seed
	//srand (time(NULL));
	srand(3);

	//loop through time and evolve the coast
	CoastVector.WriteCoast(WriteCoastFile, Time);
	double visT=0.;
	while (Time < EndTime)
	{
		if (Time >= visT )
		{
			cout << "Simulated time is " << Time << " years"<< endl;
			visT += 1.;
		}
		//Get a new wave?
		if (Time > GetWaveTime)
		{
			MyWave = WaveClimate.Get_Wave();
			GetWaveTime += WaveTimeDelta/365.;
		}

		//Evolve coast
		TempTime = MaxTimeStep;
		CoastVector.TransportSediment(TimeStep, MyWave, CERCFlag, RefDiffFlag, FluxFraction);
		Time += TimeStep/365.;
		if (TimeStep < MaxTimeStep)
		{
			TempTime = TimeStep;
			while (TempTime < MaxTimeStep)
			{
				//Evolve coast
				CoastVector.TransportSediment(TimeStep, MyWave, CERCFlag, RefDiffFlag, FluxFraction);
				Time += TimeStep/365.;
				TempTime += TimeStep;
			}
		}
		TimeStep = MaxTimeStep;

		//Write results to file
		//CoastVector.WriteCoast(WriteCoastFile, Time);
		if (Time >= PrintTime)
		{
			CoastVector.WriteCoast(WriteCoastFile, Time);
			PrintTime += PrintTimeDelta;
		}
	}

	cout 	<< endl << "Coastline: Results written to " << WriteCoastFile
			  << endl << "Coastline: Model run complete." << endl << endl;

	exit(EXIT_SUCCESS);
}
