/*==============================================================

waveclimate.hpp

Objects storing wave climate PDFs or data with functions
that can sample to get offshore waves to feed to the 
coastline object.

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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

==============================================================*/

/** @file waveclimate.hpp
@author Martin D. Hurst, British Geological Survey
@brief wave climate objects for offshore wave conditions
@details This file contains objects for representing offshore wave 
conditions for the COVE model. Wave climates are described using 
wave direction (azimuth), significant wave height and wave period. 
Wave climates can be intialised as a single wave conditions, a 
Gaussian distribution of direction, height and period, a bimodal 
distribution (i.e. two Gaussians) or using real wave climate data.
@date 26/10/2015
*/

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <cstring>

using namespace std;

#ifndef waveclimate_HPP
#define waveclimate_HPP

/// @brief Single wave object, containing a direction, height and period.
/// @author Martin D. Hurst
/// @date 26/10/2015
class Wave
{
	private:
	
	//data members
	double Dir;
	double Period;
	double Height;
	
	//Initialise functions
	void Initialise();
	void Initialise( double WvPrd,  double WvHgt, double WvDir );
		
	public:
		
	/// @brief Initialise function. This is an empty intialisation and thus throws an error
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	Wave() { Initialise(); }
	/// @brief Initialise function. Creates a wave object from wave parameters
	/// @return Wave
	/// @param WvPrd Double, the wave period (seconds).
	/// @param WvHgt Double, the wave height (metres).
	/// @param WvDir Double, the wave direction (azimuth).
	/// @author Martin D. Hurst	
	/// @date 26/10/2015
	Wave(double WvPrd, double WvHgt, double WvDir) { Initialise( WvPrd, WvHgt, WvDir );	}
	
	//Assign functions
	void AssignWaveDirection( double WvDir );
	void AssignWaveHeight( double WvHgt );
	void AssignWavePeriod( double WvPrd );
	
	// Get Functions
	/// @brief Gets wave period value.
	/// @return Double, wave period.	
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	double Get_WavePeriod();
	/// @brief Gets wave height value.
	/// @return Double, wave height.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	double Get_WaveDirection();
	/// @brief Gets wave direction value.
	/// @return Double, wave direction.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	double Get_WaveHeight();
};

/// @brief Single wave climate object, containing a direction, height and period.
/// @author Martin D. Hurst
/// @date 26/10/2015
class SingleWaveClimate
{
	private:

	// data members
	double Dir;
	double Period;
	double Height;
	
	// Initialise functions
	void Initialise();
	void Initialise(double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight);

	public:
	
 	/// @brief Initialise function. This is an empty intialisation and thus throws an error	
	/// @author Martin D. Hurst
	/// @date 26/10/2015
 	SingleWaveClimate()	{ Initialise(); }
	/// @brief Initialise function. Creates a single wave climate object from wave parameters
	/// @return SingleWaveClimate
	/// @param OffshoreWavePeriod Double, the wave period (seconds).
	/// @param OffshoreWaveHeight Double, the wave height (metres).
	/// @param OffshoreWaveDirection Double, the wave direction (azimuth).
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	SingleWaveClimate(double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight) 
	{	Initialise(OffshoreWavePeriod, OffshoreWaveDirection, OffshoreWaveHeight); }
	
	// Get Functions
	/// @brief Gets wave period value.
	/// @return Double, wave period.	
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	double Get_WavePeriod();
	/// @brief Gets wave height value.
	/// @return Double, wave height.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	double Get_WaveDirection();
	/// @brief Gets wave direction value.
	/// @return Double, wave direction.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	double Get_WaveHeight();
	
};

/// @brief Gaussian wave climate object, described by mean and standard deviation
/// wave direction, height and period.
/// @author Martin D. Hurst
/// @date 26/10/2015
class GaussianWaveClimate 
{
	
	//BimodalWaveClimate will need to be a friend since it uses two GaussianWaveClimates
	friend class BimodalWaveClimate;
	
	private:
	
	//data members
	double Period_Mean;
	double Period_StD;
	double Dir_Mean;
	double Dir_StD;
	double Height_Mean;
	double Height_StD;
	
	//Initialise functions
	void Initialise();
	void Initialise(	double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveDirection, double OffshoreStDWaveDirection,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight);											
	
	public:
	/// @brief Initialise function. This is an empty intialisation and thus throws an error.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	GaussianWaveClimate()	{ Initialise(); }
	/// @brief Initialise function. Creates a Gaussian wave climate object from wave parameters.
	/// @return GaussianWaveClimate
	/// @param OffshoreMeanWavePeriod Double, the mean wave period (seconds).
	/// @param OffshoreMeanWaveHeight Double, the mean wave height (metres).
	/// @param OffshoreMeanWaveDirection Double, the mean wave direction (azimuth).
	/// @param OffshoreStDWavePeriod Double, the standard deviation of wave period (seconds).
	/// @param OffshoreStDWaveHeight Double, the standard deviation of wave height (metres).
	/// @param OffshoreStDWaveDirection Double, the standard deviation of wave direction (azimuth).
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	GaussianWaveClimate(double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveDirection, double OffshoreStDWaveDirection,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight)
	{
		Initialise(	OffshoreMeanWavePeriod, OffshoreStDWavePeriod,
						OffshoreMeanWaveDirection, OffshoreStDWaveDirection,
						OffshoreMeanWaveHeight, OffshoreStDWaveHeight);
	}
	/// @brief Function to randomly sample a wave from the wave climate.
	/// @return Wave object
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	Wave Get_Wave();
};

/// @brief Real wave climate object, reading wave data from a text file
/// containing wave direction, height and period.
/// @author Martin D. Hurst
/// @date 26/10/2015
class RealWaveClimate 
{
	private:
	
	// data members
	vector<double> Period;
	vector<double> Dir;
	vector<double> Height;
	int NoWaves, WaveCounter;
	
	// initialise functions
	void Initialise();
	void Initialise(string WaveFileName);											
	
	public:
	
	/// @brief Initialise function. This is an empty intialisation and thus throws an error.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	RealWaveClimate()	{ Initialise(); }
	
	/// @brief Initialise function. Creates a real wave climate object from a specified text file.
	/// @details  Reads wave data from a text file. File format is three columns, the first is wave 
	/// directions (azimuth), the second is wave periods (seconds) and the third is wave height (metres).
	/// The function assumes that the file as 1 header line.
	/// @return RealWaveClimate
	/// @param WaveFileName String, name of the wave file
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	RealWaveClimate(string WaveFileName) { Initialise(WaveFileName); }
	
	/// @brief Function to sample a wave from the wave climate.
	/// @details Function updates an index counter so that waves are sampled from the RealWaveClimate
	/// object sequentially.
	/// @return Wave object
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	Wave Get_Wave();
};

/// @brief Four-bin PDF wave climate object following Ashton & Murray (2006).
/// @details Wave directions represented in a four-bin PDF of wave directions. Wave climate thus described
/// by two parameters: asymmetry (A) and highness (U). Assymetry is the proportion of waves approaching from
/// the upcoast direction (relative to the overall trend of the coastline). Highness is the proportion of waves
/// that impinge the coastline at high angle (>45^o). See Ashton, A.D., Murray, A.B., 2006. High-angle wave 
/// instability and emergent shoreline shapes: 1. Modeling of sand waves, flying spits, and capes. J. 
/// Geophys. Res., 111(F4), F04011, doi:<a href="https://dx.doi.org/10.1029/2005jf000422">10.1029/2005jf000422.</a>.
/// Wave height and period represented by Gaussian distributions.
/// @author Martin D. Hurst
/// @date 26/10/2015
class UAWaveClimate 
{
	private:
	
	//data members
	double Period_Mean;
	double Period_StD;
	double Height_Mean;
	double Height_StD;
	double U;
	double A;
	double CoastTrend;
	
	void Initialise();
	void Initialise( double input_U, double input_A, double Trend, 
							double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight);
	
	public:
	
	/// @brief Initialise function. This is an empty intialisation and thus throws an error.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	UAWaveClimate()	{ Initialise(); }
	
	/// @brief Initialise function. Creates a UA wave climate object from U, A and wave parameters.
	/// @return UAWaveClimate
	/// @param input_U Double, fraction of waves that are high angle (0-1).
	/// @param input_A Double, fraction of waves from the upcoast direction (0-1).
	/// @param Trend Double, the orientation of the coast (azimuth). Sea to the left as you look down the line.
	/// @param OffshoreMeanWaveHeight Double, the mean wave height (metres).
	/// @param OffshoreMeanWavePeriod Double, the mean wave period (seconds).
	/// @param OffshoreStDWaveHeight Double, the standard deviation of wave height (metres).
	/// @param OffshoreStDWavePeriod Double, the standard deviation of wave period (seconds).
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	UAWaveClimate(	double input_U, double input_A, double Trend, 
						double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
						double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight)
	{
		Initialise(	input_U, input_A,Trend, 
						OffshoreMeanWavePeriod, OffshoreStDWavePeriod, 
						OffshoreMeanWaveHeight, OffshoreStDWaveHeight);
	}
	/// @brief Function to sample a random wave from the wave climate.
	/// @return Wave object
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	Wave Get_Wave();											
};

/// @brief Bimodal wave climate object, described by two GaussianWaveClimate objects.
/// @author Martin D. Hurst
/// @date 26/10/2015
class BimodalWaveClimate 
{
	private:
	
	//data members
	GaussianWaveClimate WaveMode1;
	GaussianWaveClimate WaveMode2;
	
	void Initialise();
	void Initialise(double FractionWaveDirection1,
							double OffshoreMeanWavePeriod1, double OffshoreStDWavePeriod1,
							double OffshoreMeanWaveDirection1, double OffshoreStDWaveDirection1,
							double OffshoreMeanWaveHeight1, double OffshoreStDWaveHeight1,
							double OffshoreMeanWavePeriod2, double OffshoreStDWavePeriod2,
							double OffshoreMeanWaveDirection2, double OffshoreStDWaveDirection2,
							double OffshoreMeanWaveHeight2, double OffshoreStDWaveHeight2);											
	
	public:
	
	/// @brief Initialise function. This is an empty intialisation and thus throws an error.
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	BimodalWaveClimate() { Initialise(); }
	
	/// @brief Initialise function. Creates a Bimodal wave climate object from two GaussianWaveCliamte objects according to wave parameters.
	/// @return BimodalWaveClimate
	/// @param OffshoreMeanWavePeriod1 Double, the mean wave period for the first mode (seconds).
	/// @param OffshoreMeanWaveHeight1 Double, the mean wave height  for the first mode (metres).
	/// @param OffshoreMeanWaveDirection1 Double, the mean wave direction for the first mode (azimuth).
	/// @param OffshoreStDWavePeriod1 Double, the standard deviation of wave period for the first mode (seconds).
	/// @param OffshoreStDWaveHeight1 Double, the standard deviation of wave height for the first mode (metres).
	/// @param OffshoreStDWaveDirection1 Double, the standard deviation of wave direction for the first mode (azimuth).
	/// @param OffshoreMeanWavePeriod2 Double, the mean wave period for the second mode (seconds).
	/// @param OffshoreMeanWaveHeight2 Double, the mean wave height  for the second mode (metres).
	/// @param OffshoreMeanWaveDirection2 Double, the mean wave direction for the second mode (azimuth).
	/// @param OffshoreStDWavePeriod2 Double, the standard deviation of wave period for the second mode (seconds).
	/// @param OffshoreStDWaveHeight2 Double, the standard deviation of wave height for the second mode (metres).
	/// @param OffshoreStDWaveDirection2 Double, the standard deviation of wave direction for the second mode (azimuth).
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	
	BimodalWaveClimate(double FractionWaveDirection1,
							double OffshoreMeanWavePeriod1, double OffshoreStDWavePeriod1,
							double OffshoreMeanWaveDirection1, double OffshoreStDWaveDirection1,
							double OffshoreMeanWaveHeight1, double OffshoreStDWaveHeight1,
							double OffshoreMeanWavePeriod2, double OffshoreStDWavePeriod2,
							double OffshoreMeanWaveDirection2, double OffshoreStDWaveDirection2,
							double OffshoreMeanWaveHeight2, double OffshoreStDWaveHeight2)
	{
		Initialise(	FractionWaveDirection1,
						OffshoreMeanWavePeriod1, OffshoreStDWavePeriod1,
						OffshoreMeanWaveDirection1, OffshoreStDWaveDirection1,
						OffshoreMeanWaveHeight1, OffshoreStDWaveHeight1,
						OffshoreMeanWavePeriod2, OffshoreStDWavePeriod2,
						OffshoreMeanWaveDirection2, OffshoreStDWaveDirection2,
						OffshoreMeanWaveHeight2, OffshoreStDWaveHeight2);
	}
	
	/// @brief Function to sample a random wave from the wave climate.
	/// @return Wave object	
	/// @author Martin D. Hurst
	/// @date 26/10/2015
	Wave Get_Wave();
	
};
#endif
