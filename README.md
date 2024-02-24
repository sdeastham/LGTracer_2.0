# LGTracer: Mass Lagrangian point advection code

LGTracer is designed to enable point-oriented tracking of quantities in an atmospheric flow field.

Required packages:

* MathNet.Numerics

## Basic principle

There are three key classes:
* LGPoint: represents a single air parcel. Can carry a few simple quantities.
* LGPointManager: handles a collection of LGPoints, including their creation and destruction.
* DomainManager: stores domain properties including boundaries and meteorology.

## Running LGTracer

The code is not currently designed for user friendliness. It uses a single domain and a single set of points, and assumes that MERRA-2 meteorological data will be available in a hard-coded path. However, a few options can be easily modified by tweaking variables at the top of the main routine:

* Maximum number of points allowed (an exception is thrown if exceeded)
* Initial number of points (spread evenly over the domain)
* Domain size (x is longitude, y latitude, p pressure in Pa)
* Simulation length in days
* Simulation time step in seconds
* Output time step in seconds

## Output

All data will be stored in the file output.nc in a point-wise fashion. At each time point, a vector of each of the following is stored:

* Point X positions (longitudes)
* Point Y positions (latitudes)
* Point pressures (Pa)
* Point unique identifiers, or UIDs (these are unsigned ints)
* Any other properties (e.g. temperature)

To track a point from one time to the next, you need to find the correct UID in the output. UIDs cycle so it is possible - albeit very unlikely - for a point to seem to "teleport" from one output time to the next.