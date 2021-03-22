/*==============================================================

coastline.hpp

The coastline object
coastline is a vector based object defining the location
and paramaters of a soft-sediment coast line.

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

/** @file coastline.hpp
@author Martin D. Hurst, British Geological Survey
@brief coastline object for a soft sediment coast
@details This file contains the coastline object for the COVE model,
evolving a vector representing a soft sediment coast. Alongshore
sediment transport is driven by the height and angle of breaking waves
following a flux law (e.g. the CERC equation).  
@date 28/10/2015
*/

/** @mainpage
This is the documentation for the COastal Vector Evolition (COVE) model.

@image html COVE_logo.jpg

It has been generated automatically using <a href="http://www.stack.nl/~dimitri/doxygen/">Doxygen</a>.
More instructions for using the model can be found the project website. <-- Add href.

@author Martin D. Hurst, British Geological Survey
*/

#include <cmath>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <queue>
#include "waveclimate.hpp"

using namespace std;

#ifndef coastline_HPP
#define coastline_HPP

/// @brief coastline class object for evolving a vector coastline.
/// @details The coastline class object is the main object used in 
/// evolving a one-line coast. It consists of paired X, Y vectors 
/// which describe the location of the coastline, a series of vectors
/// describing the coastline and coastal cell geometry, and a series 
/// of functions with which to transform waves impinging on the coast, 
/// calculate alongshore sediment transport and update the coastline 
/// position during a model timestep.
/// @author Martin D. Hurst
/// @date 28/10/2015
class Coastline
{
	private:
	//data members
	int NoNodes;										//Number of nodes along coastline
	vector<double> X;									//position in x (m) (private)
	vector<double> Y;									//position in y (m) (private)
	vector<double> Distance;						//Distance along the vector
	vector<double> Orientation;					//orientation/azimuth of shoreline across i-1, i+1 (private)
	vector<double> FluxOrientation;				//shoreline orientation bewteen i, i+1 (private)
	vector<double> LeftOrientation;			  //orientation for meshing left boundary
	vector<double> RightOrientation;			//orientation for meshing right boundary
	vector<double> CellWidth;						//Width of individual cells parallel to orientation
	vector<double> BreakingWaveHeight;			//Heights of breaking waves
	vector<double> BreakingWaveAngle;			//Angle of breaking waves
	vector<double> Shadows;							//Shadow Zone Vector
	vector<double> ShadowZoneWaveDirection;	//Wave direction modified for ref/diff in shadow zone
	vector<double> ShadowZoneWaveHeight;		//Wave height modified for ref/diff in shadow zone
	vector<double> LongshoreFlux;					//Volume of sediment transported alongshore in m3/day
	vector<double> VolumeChange;					//Volume change in each cell during a particular timestep (m3/day)
	vector<double> PositionChange;				//Magnitude of change in shoreline position
	vector<double> Volume;
	vector<double> Area;
	vector<double> e1;
	vector<double> e2;
	vector<double> XMidPoints;						//X position of downcoast midpoint to next node
	vector<double> YMidPoints;						//Y position of downcoast midpoint to next node
	vector<double> alpha_0;
	vector<double> X0;
	vector<double> Y0;
	vector<double> XL;
	vector<double> YL;
	vector<double> XR;
	vector<double> YR;
	vector<double> Dsf;
	vector<double> Alpha_0;
	
	vector<int> Fixed;						//Is coastline fixed?
	vector<int> NodesToRemove;				//For removing nodes that are getting squashed
	
	int StartBoundary;						//Boundary type periodic=1 /no flux=2 /continous flux=3 /sinks=4 (private)
	int EndBoundary;						//Boundary type periodic=1 /no flux=2 /continous flux=3 /sinks=4 (private)
	int MeanNodeSpacing;					//average spacing between nodes (m) (private)
	int DesiredNodeSpacing;					//for adding and removing nodes?
	int ShadowFlag;
	int TriangleFlag;
	int BuildCellsFlag;
	int NodeAddedFlag;
	
	double Trend;							//redundant?
	double TotalVolume;
	
	double OffshoreWavePeriod, OffshoreWaveHeight, OffshoreWaveDirection;
	double ClosureDepth;					// } This needs to moved at some point
	double ShorefaceSlope;					// }
	double BeachWidth;						// }
	
	/** @brief Calculates the morphological properties of the coastline
	@details Resets the coastline morphological vectors and recalculates a variety
	of geometric properties and metrics. Calculates Distance along the coast, Orientation 
	(3-cell) of the coastline, FluxOrientation (2-cell), the orientation of cell boundaries 
	(e1 and e2), the CellWidth, and the position of the cell edges for constructing 
	coastline cells (BuildCellGeometries).
	@author Martin D. Hurst
  @date 29/10/2015 */ 
	void CalculateMorphology();

  /* @brief Updates the morphological properties of the coastline
	@details Updates the coastline morphological vectors and recalculates a variety
	of geometric properties and metrics. This is called when wqe don't want to build
	new cells using BuildCellGeometries every singel timestep
	@author Martin D. Hurst
  @date 29/10/2015 */ 
	void UpdateMorphology();
	
	/* @brief Checks nodes haven't got too close together or too far apart.
	@details Checks to see whether there are nodes too close together 
	(<0.66*DesiredNodeSpacing) or too far apart (>1.5*DesiredNodeSpacing). 
	Nodes are added by interpolation or deleted as appropriate
	@author Martin D. Hurst
  @date 29/10/2015 */ 
	void CheckNodeSpacing();
	void CalculateMeanNodeSpacing();
	void Sinks();
	
	/* @brief Checks the coastline doesn't intersect itself anywhere.
	@author Martin D. Hurst
  @date 29/10/2015 */ 
	void IntersectionAnalysis();
		
	/* @brief Transforms wave from offshore to nearshore following Airy wave theory.
	@details Assumes shore parallel contours and calculates refraction and shoaling
	coefficients then updates wave hieght and angle. Wave breaking occurs when wave 
	height exceeds 0.8*WaterDepth. Breaking wave conditions (height and angle) are
	stored to be fed to the alongshore seidment tranpsort equation (see CalculateFlux).
	@author Martin D. Hurst
  @date 29/10/2015 */
  void TransformWaves();
  
  /* @brief Finds regions of the coast that are in shadow of offshore waves.
	@details Loops through coast from start to end and end to start to check 
	whether the coastline casts shadows given the offshore wave direction
	There are five codes to the Shadows vector:
		#0: No shadow
		#1: In Shadow
		#2: Casts a shadow forward along the coast
		#3: Casts a shadow backward along the coast
		#4: Special case for cell just downdrift of backward shadow (#3)
	These Codes are required for handling sediment transport and
	for the refraction/diffraction code 
	@author Martin D. Hurst
  @date 29/10/2015 */
	void GetShadows();
	
	/* @brief Computes wave angles and heights in the shadow zone.
	@details Wave angle modified in the shadow zone as 1.5*omega, where omega
	is the angle within the shadow zone (Kraus, 1984). Wave height is assumed 
	to be halved at the tip of the shadow zone due to wave spreading into the
	shadow zone. The distribution of wave height in and adjacent to the shadow 
	zone is there described as a sine function of omega, with wave heights also
	reduced outside the shadow zone to compensate.
	@author Martin D. Hurst
  @date 29/10/2015 */
	void RefractDiffractShadowZone();
	
	/* @brief Function to mesh the shoreface cells and calculate cell area.
	@details Function to build shoreline cells based on coastline vector 
	by starting in the most concave out parts of the shoreline and building 
  outward using a priority queue. Calculates the cell surface area to 
  facilitate calculation of changes in shoreline position following a 
  quadratic (SolveQuadratic) or cubic (SolveCubic) solution. 
  CellType flag  identifies the whether the cell has been processed or 
  what combination of adjacent cells have been processed in order to 
  determine how the cell is handled:
    1 = Unprocessed cell
    2 = Unprocessed cell where upcoast cell has been processed
    3 = Unprocessed cell where downcoast cell has been processed
    4 = Unprocessed cell where both upcoast and downcoast cells have 
        been processed
  This is the most cumbersome part of the model code. It may be possible 
  to make this more efficient, but at the moment I have just limited how 
  often it gets called, since it is not needed every single timestep.
  See TransportSediment for more details on this.
  @author Martin D. Hurst
  @date 29/10/2015 */
	void BuildCellGeometries();
		
	//Use priority queue to deal with triangle nodes in order of increasing shoreface depth?
	struct CoastNode
	{
		double ShorefaceDepth;
		int i;
	};

	class CompareNodes
	{
		public:
		bool operator()(CoastNode& Node1, CoastNode& Node2)
		{
			if (Node1.ShorefaceDepth > Node2.ShorefaceDepth) return true;
			else return false;
		}	
	};

	priority_queue<CoastNode, vector<CoastNode>, CompareNodes> CoastlineQueue;
	vector< vector<CoastNode> > Recievers;
	
	//Structure to handle single nodes
	struct Node
	{
		double X, Y;
	};
	vector< vector<Node> > Vertices;
	
	void Initialise();
	void Initialise(string xyfilename);
	void Initialise(string xyfilename, float StartTime);
	void Initialise(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary);
	
  void WriteXYFile();
  void WriteNodesFile();
  void WriteRecieversFile();
  void WriteShadowsFile();
  void PrintCoastlineQueue();

  public:
	
	/// @brief Initialise function. This is an empty intialisation and thus throws an error	
	/// @author Martin D. Hurst
	/// @date 28/10/2015
	Coastline()	{ Initialise(); }
	
	/// @brief Initialise function. Creates a Coastline object from a file containing x and y data
	/// @return Coastline
	/// @param xyfilename String, the file containing the coastline x and y coordinates
	/// @author Martin D. Hurst
	/// @date 28/10/2015
	Coastline(string xyfilename) { Initialise(xyfilename); }
	
	/// @brief Initialise function. Creates a Coastline object from a file containing x and y data at a specified time.
	/// @details File format is the same as that outputted by the WriteCoast method.
	/// @return Coastline
	/// @param xyfilename String, the file containing the coastline x and y coordinates
	/// @param StartTime Float, the time to get at which to get x and y data from file
	/// @author Martin D. Hurst
	/// @date 28/10/2015
	Coastline(string xyfilename, float StartTime) { Initialise(xyfilename, StartTime); }
	
	/// @brief Initialise function. Creates a Coastline object as a straight line following an azimuth.
	/// @return Coastline
	/// @param MeanNodeSpacing Int, the average distance between coastal nodes (metres)
	/// @param CoastLength Double, the length of the coastline to be created
	/// @param Trend Double, the azimuth direction in which to create the coast.
	/// @param StartBoundary, the boundary condition at the start of the line (1=fixed, 2=periodic)
	/// @param EndBoundary, the boundary condition at the end of the line (1=fixed, 2=periodic), must be same as StartBoundary (for the moment).
	/// @author Martin D. Hurst
	/// @date 28/10/2015
	///Initialise coastline as a straight line.
	Coastline(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary)
	{
		Initialise(MeanNodeSpacing, CoastLength, Trend, StartBoundary, EndBoundary);
	}
	
	/// @brief Function to read a coastline from a file at a given time.
	/// @details File format is the same as that outputted by the WriteCoast method.
	/// @param InputFileName String, the file containing the coastline x and y coordinates
	/// @param Time Double, the time to get at which to get x and y data from file
	/// @author Martin D. Hurst
	/// @date 28/10/2015	
	void ReadCoast(string InputFileName, double Time);
	
	/// @brief Function to write a coastline to a file at a given time.
	/// @param OutputFileName String, the file to write the coastline x and y coordinates
	/// @param Time Double, the model time written to file
	/// @author Martin D. Hurst
	/// @date 28/10/2015	
	void WriteCoast(string OutputFileName, double Time);
	
	/// @brief Create priority queue object for cell builder
	/// @details Called before BuildCellGeometries to initiate a priority queue
	/// Function to find shoreface triangles and create priority queue for building
	/// cells, starting with the most convex out parts of the coastline
	/// @author Martin D. Hurst
	/// @date 29/10/2015	
	void SetupQueue();

	/* @brief Perform sediment transport and evolve coast
	@details This is the main function call in the model, all else is called from inside 
	this function during a model run. The method begins by calculates the coastline morphology 
	(CalculateMorphology), generates the coastal cells (BuildCellGeometries), identifies 
	shadowed regions of the coast (GetShadows), diffracts waves into the shadowed regions
	(RefractDiffractShadowZone) and transforms offshore waves to wave breaking (TransformWaves).
	Subsequently the model calculates alongshore sediment transport following your favourite 
	alongshore sediment transport formula (CalculateFlux) to determine the volume change in 
	each coastal cell, and then calculates the position change of the coast using SolveCubic or 
	SolveQuadratic. Finally IntersectionAnalysis checks for intersections in the coast (e.g. 
	breaching/reattaching) and then the model checks that the spacing of coastal nodes hasn't
  gotten too big or small (CheckNodeSpacing).
  @param TimeDelta Double, the model timestep (days)
  @param TheWave Wave Object sampled from a waveclimate object.
  @param FluxType Int, 1 = CERC equation, others are yet to be implemented, default is 1.
  @param RefDiffFlag Int, 0 = Off, 1 = On, default is 1.
  @param FluxFraction Double, Proportion of sediment input at boundaries when fixed (untested, default is 0).
  @param LostFluxFraction Double, Proportion of sediment lost offshore during flux (default is 0).
  @author Martin D. Hurst
	@date 29/10/2015	
	*/
	void TransportSediment(double &TimeDelta, Wave TheWave, int FluxType=1, int RefDiffFlag=1, double FluxFraction=0, double LostFluxFraction=0);
		
	/*Flux Type 	
               0 = Simple Diffusion
					1 = CERC equation (default)
					2 = Kamphius equation
					3 = Bailard equation
					4 = Deigaard equation
	*/
  
  /// @brief Calculate volumetric sediment fluxes
  /// @details Evaluates the chosen sediment flux law at node i
  /// @param i Int, the node index
  /// @param FluxType Int, 1 = CERC equation, others are yet to be implemented, default is 1.
  /// @author Martin D. Hurst
  /// @date 29/10/2015	
	void CalculateFlux(int i, int FluxType);
		
	/// @brief Perform quadratic solution for position change of a cell
	/// @details Position change of the coast is found to be a quadratic
	/// function of volume change and the geometry of the cell when it 
	/// extends to the bottom of the shoreface (see Hurst et al. 2015?).
	/// @param i Int, the node index
	/// @author Martin D. Hurst
  /// @date 29/10/2015	
	double SolveQuadratic(int i);
	
	/// @brief Perform cubic solution for position change of a cell
	/// @details Position change of the coast is found to be a cubic
	/// function of volume change and the geometry of the cell when it 
	/// does not extend to the bottom of the shoreface (see Hurst et al. 2015?).
	/// @param i Int, the node index
	/// @author Martin D. Hurst
  /// @date 29/10/2015
	double SolveCubic(int i);
	
	/*****************************************\
	| Get Functions to return private members |
	\*****************************************/
	
	/// @brief Return number of nodes in coastline.
	/// @details Function to return the private member NoNodes value, 
	/// the number of nodes in the coastline object
	/// @returns Int, Number of nodes
	/// @author Martin D. Hurst
  /// @date 6/1/2014
	int get_NoNodes() const						{return NoNodes;}			//get Number of nodes along coastline
	
	/// @brief Return coastline X-values.
  /// @details Function to return the private member X values
	/// @returns vector<double>, X coordinate values of the coastline 
	/// @author Martin D. Hurst
  /// @date 7/1/2014
  vector<double> get_X() const				{return X;}					//get X vector position in x (m)
	
	/// @brief Return coastline Y-values.
  /// @details Function to return the private member Y values
	/// @returns vector<double>, Y coordinate values of the coastline 
	/// @author Martin D. Hurst
  /// @date 7/1/2014
	vector<double> get_Y() const				{return Y;}					//get Y vectir position in y (m)
	
	/// @brief Return coastline orientation.
	/// @details Function to return the private member Orientation, 
	/// the azimuth direction of the line connecting the two adjacent 
	/// nodes to the node of interest (3-cell orientation). The evolving 
	/// cosatal node advances or retreats perpendicular to this orientation.
	/// @returns vector<double>, Orientation values of the coastline 
	/// @author Martin D. Hurst
  /// @date 6/1/2014
	vector<double> get_Orientation() const 	 	{return Orientation;}		//get orientation of shoreline between i-1 and i+1
	
	/// @brief Return coastline flux orientation.
	/// @details Function to return the private member FluxOrientation, 
	/// the azimuth direction of the line connecting the node of interest 
	/// to the next node in the vector (2-cell orientation). This 
	/// FluxOrientation is used to transform waves and for determining 
	/// the amount of alongshore sediment transport.
	/// @returns vector<double>, FluxOrientation values of the coastline 
	/// @author Martin D. Hurst
  /// @date 6/1/2014
	vector<double> get_FluxOrientation() const 	{return FluxOrientation;}	//get shoreline orientation bewteen i and i+1

  /// @brief Return coastal shadows.
	/// @details Function to return the private member Shadows, an integer
	/// vector which describes whether or not a node is in shadow, or casts 
	/// a shadow, or is open coast.
	/// @returns vector<double>, Shadows values of the coastline 
	/// @author Martin D. Hurst
  /// @date 6/1/2014
	vector<double> get_Shadows() const			{return Shadows;}			//get Shadows	
	
	/// @brief Return start boundary condition.
	/// @Function to return the private member StartBoundary, the boundary 
	/// condition at the start end of the vector. Currently, StartBoundary 
	/// and EndBoundary must have the same value. 1 = Periodic boundary, 
	/// 2 = Fixed Boundary (with or without sediment flux).
	/// @returns Int, StartBoundary condition
	/// @author Martin D. Hurst
  /// @date 8/1/2014
	int get_StartBoundary() const				{return StartBoundary;}		//get Boundary type start
	
	/// @brief Return end boundary condition.
	/// @Function to return the private member EndBoundary, the boundary 
	/// condition at the end of the vector. Currently, StartBoundary 
	/// and EndBoundary must have the same value. 1 = Periodic boundary, 
	/// 2 = Fixed Boundary (with or without sediment flux).
	/// @returns Int, EndBoundary condition
	/// @author Martin D. Hurst
  /// @date 8/1/2014
	int get_EndBoundary() const					{return EndBoundary;}		//get Boundary type end

	/// @brief Return desired mean node spacing.
  /// @details Function to return the private member MeanNodeSpacing, 
  /// the desired mean node spacing for the model.
	/// @returns Int, MeanNodeSpacing
	/// @author Martin D. Hurst
  /// @date 8/1/2014
	int get_MeanNodeSpacing() const				{return MeanNodeSpacing;}	//get average spacing between nodes (m)

};

#endif
