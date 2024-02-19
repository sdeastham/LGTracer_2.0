using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;

using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
//using MathNet.Numerics.Integration;
using MathNet.Numerics.Random;
//using MathNet.Numerics.Distributions;

using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

namespace LGTracer
{
    public class Program
    {
        private static void Main(string[] args)
        {
            /* LGTracer is a very simple test code designed to simulate movement of points
            through a simple 2D space under the influence of a constant wind field. */
            Console.WriteLine("Initiating simple LGTracer program");

            // Set up output file
            string fileName = "output.nc";
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.Create
            };

            // Number of Lagrangian points to track
            int nPoints = 100000;
            int nInitial = 1000; // Points to initially scatter randomly
            double pointRate = 50.0/600.0; // Number of new points to add in each second
            bool debug = false;

            // Specify the domains
            double[] lonLims = [-80.0,15.0];
            double[] latLims = [30.0,60.0];
            int readLevel = 30;
            int readTime = 0;
            string metFileNameA3 = "C:/Data/MERRA-2/2023/01/MERRA2.20230101.A3dyn.05x0625.nc4";
            string metFileNameI3 = "C:/Data/MERRA-2/2023/01/MERRA2.20230101.I3.05x0625.nc4";

            // Major simulation settings
            double nDays = 30.0; // Days to run
            double dt = 60.0 * 5.0; // Time step in seconds
            double dtStorage = 60.0*15.0; // // How often to save out data (seconds)

            // Point settings
            double pressureDelta = 10.0e2; // Layer pressure thickness in Pa. This helps to determine how many points to seed
            double kgPerPoint = 2.0e12; // Air mass represented by a single point (mass flows are huge - but 1e11 seems very high? Total atm mass 5.1e18!)


            // CODE STARTS HERE

            // Read in MERRA-2 data and use that to set domain
            (double[] lonEdge, double[] latEdge, int[] lonSet, int[] latSet ) = MERRA2.ReadLatLon( metFileNameA3, lonLims, latLims );
            // Now read in the U and V data
            (double[,]uWind, double[,]vWind ) = MERRA2.ReadA3( metFileNameA3, readTime, readLevel, lonSet, latSet );
            // Also read in PS, T, and QV
            (double[,] surfacePressure, double[,]griddedTemperature, double[,]griddedSpecificHumidity ) = MERRA2.ReadI3( metFileNameI3, readTime, readLevel, lonSet, latSet );

            // Set up the mesh (units of degrees)
            // The original lon/lat limits aren't important - need the mesh limits
            int xPosts = lonEdge.Length;
            int xCells = xPosts - 1;
            int yPosts = latEdge.Length;
            int yCells = yPosts - 1;

            double xMin  = lonEdge[0];
            double xMax  =  lonEdge[xPosts-1];
            double xSpan = xMax - xMin;
            double yMin  = latEdge[0];
            double yMax  = latEdge[yPosts-1];
            double ySpan = yMax - yMin;
            // Assume uniform spacing
            double dx = xSpan/(xPosts - 1);
            double dy = ySpan/(yPosts - 1);

            // For convenience
            double[] xLims = {xMin,xMax};
            double[] yLims = {yMin,yMax};

            // X edges (degrees)
            double[] xMesh = new double[xPosts];
            double[] xMid = new double[xCells];
            for (int i=0;i<xPosts;i++)
            {
                xMesh[i] = lonEdge[i];
                if (i > 0)
                {
                    xMid[i-1] = (xMesh[i] + xMesh[i-1])/2.0;
                }
            }

            // Y edges (degrees)
            double[] yMesh = new double[yPosts];
            double[] yMid  = new double[yCells];
            for (int i=0;i<yPosts;i++)
            {
                //yMesh[i] = yMin + (dy * i);
                yMesh[i] = latEdge[i];
                if (i > 0)
                {
                    yMid[i-1] = (yMesh[i] + yMesh[i-1])/2.0;
                }
            }

            // Boundary lengths (m)
            double[] boundaryLengths = new double[xCells*2 + yCells*2];
            double earthCircumference = 2.0 * Math.PI * LGConstants.EarthRadius;
            double edgeLength;

            // South boundary
            edgeLength = (dx/360.0) * earthCircumference * Math.Cos(LGConstants.Deg2Rad*latEdge[0]);
            for (int i=0;i<xCells;i++)
            {
                boundaryLengths[i] = edgeLength;
            }
            // North boundary
            edgeLength = (dx/360.0) * earthCircumference * Math.Cos(LGConstants.Deg2Rad*latEdge[yPosts-1]);
            for (int i=0;i<xCells;i++)
            {
                boundaryLengths[xCells + yCells + i] = edgeLength;
            }
            // All cells on the East and West boundaries have a constant length
            edgeLength = (dy/360.0) * earthCircumference;
            for (int i=0;i<yCells;i++)
            {
                boundaryLengths[xCells + i] = edgeLength;
                boundaryLengths[(xCells*2) + yCells + i] = edgeLength;
            }

            // Create the velocity array in m/s
            double[,] xSpeed = new double[yCells,xCells];
            double[,] ySpeed = new double[yCells,xCells];
            for (int i=0;i<xCells;i++)
            {
                for (int j=0;j<yCells;j++)
                {
                    //(xSpeed[j,i], ySpeed[j,i]) = VelocitySolidBody(xMid[i],yMid[j],omega);
                    xSpeed[j,i] = uWind[j,i];// * roughWindConversion;
                    ySpeed[j,i] = vWind[j,i];// * roughWindConversion;
                }
            }

            // Time handling
            double duration = 60.0 * 60.0 * 24.0 * nDays; // Simulation duration in seconds
            double tStart = 0.0;
            double tStop = tStart + duration;
            double tCurr = tStart;
            int iterMax = (int)Math.Ceiling((tStop - tStart)/dt);
            double tStorage = tStart; // Next time that we want storage to occur

            // Central RNG for random point seeding
            System.Random RNG = SystemRandomSource.Default;
            // Use this if debugging
            //int seed = 31567891;
            //System.Random RNG = new SystemRandomSource(seed);

            // The point manager holds all the actual point data and controls velocity calculations (in deg/s)
            Func<double, double, (double, double)> vCalc = (double x, double y) => VelocityFromFixedSpaceArray(x,y,xMin,xMax,dx,yMin,yMax,dy,xSpeed,ySpeed,false);
            PointManager pointManager = new PointManager(nPoints,vCalc);

            // Scatter N points randomly over the domain
            (double[] xInitial, double[] yInitial) = MapRandomToXY(xLims[0],xLims[1],yLims[0],yLims[1],nInitial,RNG);
            pointManager.CreatePointSet(xInitial,yInitial);

            // Define boundary edges, normals etc
            ( Vector2[] xyPosts, Vector2[] boundaryNormals) = createBoundary(xMesh,yMesh);

            // Estimate the boundary velocity, given the velocity array (now in m/s)
            Func<double, double, (double, double)> vCalcMPS = (double x, double y) => VelocityFromFixedSpaceArray(x,y,xMin,xMax,dx,yMin,yMax,dy,xSpeed,ySpeed,true);
            Vector2[] vBoundary = GetBoundaryVelocities(xyPosts, vCalcMPS);

            // Set up output
            // Add a 2D variable
            // Define dimensions (variables)
            List<double> time = [];

            // X of all points will be stored as a 2D array
            List<double[]> xHistory = [];
            List<double[]> yHistory = [];
            List<uint[]> UIDHistory = [];

            // Store initial conditions - keeping track of the largest number of points being
            // tracked at any given output time
            int maxActive = ArchiveConditions(time,xHistory,yHistory,UIDHistory,tCurr,pointManager);
            tStorage += dtStorage;
            int nStored = 1;

            int nNew, nAvailable;

            // We only add an integer number of points each time step
            // If the number of points to be added is non-integer, retain
            // the surplus and add it in at the next time step
            double nNewExact;
            double nSurplus = 0.0;
            double massSurplus = 0.0;
            
            // Set up timing
            int nSteps = 0;
            var watch = new Stopwatch();
            watch.Start();
            Console.WriteLine("Beginning trajectory calculation");
            for (int iter=0;iter<iterMax; iter++)
            {
                // How many new points will we add (allowing for variable dt)?
                nNewExact = (pointRate * dt) + nSurplus;
                nNew = (int) Math.Floor(nNewExact);
                nSurplus = nNewExact - (double)nNew;

                // We want to introduce nNew at the edge, but we can only go up to nInactive
                nAvailable = Math.Min(pointManager.MaxPoints - pointManager.NActive,nNew);

                // If we have enough points available, scatter them evenly over the edges of the domain
                // TODO: Make this only at locations where we have inbound flow?
                //(double[] xSet, double[] ySet) = SeedBoundaryUniform(nAvailable,xLims,yLims,RNG);
                (double[] xSet, double[] ySet, massSurplus) = SeedBoundary(kgPerPoint, boundaryLengths, pressureDelta, dt,
                    xyPosts, boundaryNormals, vBoundary, RNG, massSurplus);

                pointManager.CreatePointSet(xSet, ySet);

                // Do the actual work
                if (debug) {Console.WriteLine($"TIME: {tCurr,7:f2}");}
                pointManager.Advance(dt);

                Cull(xLims,yLims,pointManager);

                nSteps++;
                tCurr = (iter+1) * dt;

                // For diagnostics - must take place AFTER tCurr advances
                // Only store data every dtStorage seconds. Use a small offset
                // to compensate for imperfect float comparisons
                if (tCurr >= (tStorage - 1.0e-10))
                {
                    maxActive = Math.Max(maxActive,ArchiveConditions(time,xHistory,yHistory,UIDHistory,tCurr,pointManager));
                    tStorage += dtStorage;
                    nStored += 1;
                }
            }
            watch.Stop();
            long elapsedTimeLong = watch.ElapsedMilliseconds;
            double elapsedTime = (double)elapsedTimeLong;
            double msPerStep = elapsedTime/nSteps;
            Console.WriteLine($"{nSteps} steps completed in {elapsedTime/1000.0,6:f1} seconds ({msPerStep,6:f2} ms per step)");

            bool success = WriteToFile(dsUri,time,xHistory,yHistory,UIDHistory,maxActive);
            if (success)
            {
                Console.WriteLine($"Output data with {nStored} samples [max points stored: {maxActive}] successfully written to {fileName}");
            }
            else
            {
                Console.WriteLine($"Could not write output to {fileName}");
            }

        }
        private static bool WriteToFile(NetCDFUri dsUri, List<double> time, List<double[]> xHistory, List<double[]> yHistory, List<uint[]> UIDHistory, int maxActive)
        {
            bool success = true;

            // Get the output sizes
            int nPoints = Math.Min(maxActive,xHistory[0].Length);
            int nTimes = time.Count;

            int[] index = new int[nPoints];
            for (int i=0; i<nPoints; i++ )
            {
                index[i] = i;
            }
            
            // Convert the lists into conventional 2D arrays
            double[,] x2D = new double[nTimes, nPoints];
            double[,] y2D = new double[nTimes, nPoints];
            uint[,] UIDs = new uint[nTimes, nPoints];

            for (int i=0; i<nTimes; i++)
            {
                for (int j=0; j<nPoints; j++)
                {
                    x2D[i,j] = xHistory[i][j];
                    y2D[i,j] = yHistory[i][j];
                    UIDs[i,j] = UIDHistory[i][j];
                }
            }

            using (DataSet ds = DataSet.Open(dsUri))
            {
                ds.AddAxis("index","-",index);
                ds.AddAxis("time","seconds",time.ToArray());
                ds.AddVariable(typeof(double), "x", x2D, ["time","index"]);
                ds.AddVariable(typeof(double), "y", y2D, ["time","index"]);
                ds.AddVariable(typeof(uint), "UID", UIDs, ["time","index"]);
                ds.Commit();
            }
            
            return success;
        }
        private static int ArchiveConditions(List<double> time, List<double[]> xHistory, List<double[]> yHistory, List<uint[]> UIDHistory, double tCurr, PointManager pointManager)
        {
            int nPoints = pointManager.MaxPoints;
            double[] xPoints = new double[nPoints];
            double[] yPoints = new double[nPoints];
            uint[] UIDs = new uint[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                if (i<pointManager.NActive)
                {
                    LGPoint point = pointManager.ActivePoints[i];
                    xPoints[i] = point.X;
                    yPoints[i] = point.Y;
                    UIDs[i] = point.UID;
                }
                else
                {
                    xPoints[i] = double.NaN;
                    yPoints[i] = double.NaN;
                    UIDs[i] = 0;
                }
            }
            time.Add(tCurr);
            xHistory.Add(xPoints);
            yHistory.Add(yPoints);
            UIDHistory.Add(UIDs);
            return pointManager.NActive;
        }

        private static (double, double) VelocityConst( double x, double y, double xSpeed, double ySpeed)
        {
            // Return a fixed velocity
            return (xSpeed, ySpeed);
        }
        
        private static (double, double) VelocityFromFixedSpaceArray( double x, double y, double xMin, double xMax, double dx, double yMin, double yMax, double dy, double[,] xSpeedArray, double[,] ySpeedArray, bool noConvert=false)
        {
            // Extract the velocity vector from an array
            // Assumes constant X spacing and constant Y spacing
            int nX = xSpeedArray.GetLength(1);
            int nY = xSpeedArray.GetLength(0);
            double dxdt, dydt;

            // If we made it this far - we are within the domain
            int xIndex = Math.Min(Math.Max(0,(int)Math.Floor((x - xMin)/dx)),nX-1);
            int yIndex = Math.Min(Math.Max(0,(int)Math.Floor((y - yMin)/dy)),nY-1);

            // Values in m/s
            dxdt = xSpeedArray[yIndex,xIndex];
            dydt = ySpeedArray[yIndex,xIndex];
            if (noConvert)
            {
                return (dxdt,dydt);
            }

            // Convert from m/s to deg/s
            dxdt = LGConstants.Rad2Deg * dxdt / (LGConstants.EarthRadius * Math.Cos(LGConstants.Deg2Rad*y));
            dydt = LGConstants.Rad2Deg * dydt / LGConstants.EarthRadius;
            return (dxdt, dydt);
        }

        private static (double[], double[]) SeedBoundaryUniform(int nPoints, double[] xLims, double[] yLims, System.Random RNG)
        {
            // Seed the domain boundaries completely uniformly
            double xMin = xLims[0];
            double xMax = xLims[1];
            double xSpan = xMax - xMin;
            double yMin = yLims[0];
            double yMax = yLims[1];
            double ySpan = yMax - yMin;
            double xCurr, yCurr;
            double smallDelta = 1.0e-5;
            double randomVal;

            double[] xVals = new double[nPoints];
            double[] yVals = new double[nPoints];
            
            for (int i=0; i<nPoints; i++)
            {
                // Activate the point at a location randomly chosen from the domain edge
                // Algorithm below basically goes around the edges of the domain in order
                randomVal = RNG.NextDouble() * ((xSpan*2) + (ySpan*2));
                if (randomVal < xSpan)
                {
                    yCurr = yMin + smallDelta;
                    xCurr = xMin + randomVal;
                }
                else if (randomVal < (xSpan + ySpan))
                {
                    yCurr = yMin + (randomVal - xSpan);
                    xCurr = xMax - smallDelta;
                }
                else if (randomVal < (xSpan + ySpan + xSpan))
                {
                    yCurr = yMax - smallDelta;
                    xCurr = xMin + (randomVal - (xSpan + ySpan));
                }
                else
                {
                    yCurr = yMin + (randomVal - (xSpan + ySpan + xSpan));
                    xCurr = xMin + smallDelta;
                }
                xVals[i] = xCurr;
                yVals[i] = yCurr;
            }
            return (xVals, yVals);
        }

        private static (Vector2[], Vector2[]) createBoundary( double[] xMeshDouble, double[] yMeshDouble)
        {
            // Create two arrays of 2-element vectors representing the boundary edge locations and the boundary normal vectors

            // Number of individual edge cells + 1
            int xPosts = xMeshDouble.Length;
            int yPosts = yMeshDouble.Length;

            // Vector2 is float only - do some conversions
            float[] xMesh = new float[xPosts];
            float[] yMesh = new float[yPosts];
            for (int i=0; i<xPosts; i++)
            {
                xMesh[i] = (float)xMeshDouble[i];
            }
            for (int i=0; i<yPosts; i++)
            {
                yMesh[i] = (float)yMeshDouble[i];
            }

            // Duplicate first post for convenience
            int nPosts = xPosts + (yPosts-1) + (xPosts-1) + (yPosts-1);
            Vector2[] xyPosts = new Vector2[nPosts];
            int currPost = 0;
            // South boundary
            for (int i=0; i<xPosts; i++)
            {
                xyPosts[currPost] = new Vector2(xMesh[i],yMesh[0]);
                currPost++;
            }
            // East boundary
            for (int i=1; i<yPosts; i++)
            {
                xyPosts[currPost] = new Vector2(xMesh[xPosts-1],yMesh[i]);
                currPost++;
            }
            // North boundary
            for (int i=xPosts-2; i>=0; i--)
            {
                xyPosts[currPost] = new Vector2(xMesh[i],yMesh[yPosts-1]);
                currPost++;
            }
            // West boundary (assume closure - ie include a final post identical to the first)
            for (int i=yPosts-2; i>=0; i--)
            {
                xyPosts[currPost] = new Vector2(xMesh[0],yMesh[i]);
                currPost++;
            }
            
            // Now calculate the boundary normals
            // First and last posts are duplicates for convenience
            Vector2[] normals = new Vector2[nPosts - 1];
            Vector2 boundaryVector2;
            Vector3 boundaryVector3, vertVector3;
            // Find the normal of each boundary edge by taking the cross-product of the "up" vector
            // with that of the boundary - this will always face into the domain
            vertVector3 = new Vector3(0,0,1);
            Vector3 normVector3;
            for (int i=0; i<(nPosts-1); i++)
            {
                boundaryVector2 = xyPosts[i+1] - xyPosts[i];
                boundaryVector3 = new Vector3(boundaryVector2,0);
                normVector3 = Vector3.Normalize(Vector3.Cross(vertVector3,boundaryVector3));
                normals[i] = new Vector2(normVector3.X,normVector3.Y);
            }

            return (xyPosts,normals);
        }

        private static Vector2[] GetBoundaryVelocities(Vector2[] xyPosts, Func<double, double, (double, double)> vCalc)
        {
            int nCells = xyPosts.Length - 1;
            double u, v, x, y;
            Vector2 xyMid;
            Vector2[] vBoundary = new Vector2[nCells];
            for (int i=0; i<nCells; i++)
            {
                xyMid = xyPosts[i] + (0.5f * (xyPosts[i+1] - xyPosts[i]));
                x = (double)xyMid.X;
                y = (double)xyMid.Y;
                (u, v) = vCalc(x,y);
                vBoundary[i] = new Vector2((float)u,(float)v);
            }
            return vBoundary;
        }
        private static (double[], double[], double) SeedBoundary(double kgPerPoint, double[] boundaryLengths, double pressureDelta, double dt,
            Vector2[] boundaryPosts, Vector2[] boundaryNormals, Vector2[] vBoundary, System.Random RNG, double massSurplus = 0.0)
        {
            // Seed the domain boundaries proportional to mass flow rate
            // Position along boundary for each cell is random
            int nEdges = boundaryPosts.Length - 1; // First post is duplicated as last post
            double smallDelta = 1.0e-5;
            double massFlux, cellFrac, randomVal;
            double pointSurplus = 0.0;
            double nPoints;
            Vector2 pointLocation, boundaryVector;
            double vNorm;
            double[] massFluxes = new double[nEdges];
            double[] weighting;

            for (int i=0; i<nEdges; i++)
            {
                // Mass flow rate across the boundary in kg/s
                vNorm = (double)Vector2.Dot(vBoundary[i],boundaryNormals[i]);
                massFlux = dt * vNorm * boundaryLengths[i] * pressureDelta / LGConstants.gravConstantSurface;
                massFluxes[i] = Math.Max(0.0,massFlux);
            }

            // We now know the total mass flux - use to figure out weightings
            massFlux = massFluxes.Sum();
            weighting = new double[nEdges];
            for (int i=0; i<nEdges; i++)
            {
                // Cumulative
                weighting[i] = (massFluxes[i] / massFlux);
                if (i>0)
                {
                    weighting[i] += weighting[i-1];
                }
            }
            // Force to 1.0 to ensure we don't accidentally miss anything
            weighting[nEdges-1] = 1.0;

            // How many points will we add (add mass surplus)?
            nPoints = (massFlux + massSurplus)/kgPerPoint;
            int nPointsTotal = (int)Math.Floor(nPoints);

            // Unused mass
            massSurplus = (nPoints - (double)nPointsTotal) * kgPerPoint;

            //Console.WriteLine($"Expecting {nPointsTotal} new points");
            double[] xVals = new double[nPointsTotal];
            double[] yVals = new double[nPointsTotal];
            int iCell;
            for (int iPoint=0; iPoint < nPointsTotal; iPoint++)
            {
                randomVal = RNG.NextDouble();
                iCell = 0;
                while (weighting[iCell] < randomVal)
                {
                    iCell++;
                }
                boundaryVector = boundaryPosts[iCell+1] - boundaryPosts[iCell];
                // Remainder of randomVal used to figure out how far along the cell we go
                if (iCell > 0)
                {
                    cellFrac = (randomVal - weighting[iCell-1])/(weighting[iCell] - weighting[iCell-1]);
                }
                else
                { 
                    cellFrac = randomVal/weighting[0];
                }
                pointLocation = boundaryPosts[iCell] + ((float)cellFrac * boundaryVector) + ((float)smallDelta * boundaryNormals[iCell]);
                xVals[iPoint] = pointLocation.X;
                yVals[iPoint] = pointLocation.Y;
            }
            return (xVals, yVals, massSurplus);
        }

        private static (double[], double[]) MapRandomToXY( double xMin, double xMax, double yMin, double yMax, int nPoints, System.Random rng )
        {
            // Scatter randomly throughout domain
            double xSpan = xMax - xMin;
            double xStart = xMin;
            double ySpan = yMax - yMin;
            double yStart = yMin;

            double[] xInitial = new double[nPoints];
            double[] yInitial = new double[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xInitial[i] = rng.NextDouble()*xSpan + xStart;
                yInitial[i] = rng.NextDouble()*ySpan + yStart;
            }
            return (xInitial, yInitial);
        }

        private static void Cull(double[] xLims, double[] yLims, PointManager pointManager)
        {
            // Deactivate any points which are outside the domain
            // Could also do this with LINQ
            for (int i=pointManager.NActive-1; i>=0; i--)
            {
                LGPoint point = pointManager.ActivePoints[i];
                if (point.X < xLims[0] || point.X >= xLims[1] || point.Y < yLims[0] || point.Y >= yLims[1] )
                {
                    pointManager.DeactivatePoint(i);
                }
            }
        }
    }
}