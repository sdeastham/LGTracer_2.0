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
            Func<double, double, (double, double)> vCalc = (double x, double y) => DomainManager.VelocityFromFixedSpaceArray(x,y,xMin,dx,yMin,dy,xSpeed,ySpeed,false);
            PointManager pointManager = new PointManager(nPoints,vCalc);

            // Scatter N points randomly over the domain
            (double[] xInitial, double[] yInitial) = DomainManager.MapRandomToXY(xLims[0],xLims[1],yLims[0],yLims[1],nInitial,RNG);
            pointManager.CreatePointSet(xInitial,yInitial);

            foreach (LGPoint point in pointManager.ActivePoints)
            {
                point.SetTemperature(DomainManager.NearestNeighbor(point.X,point.Y,xMin,dx,yMin,dy,griddedTemperature));
                point.SetSpecificHumidity(DomainManager.NearestNeighbor(point.X,point.Y,xMin,dx,yMin,dy,griddedSpecificHumidity));
            }

            // Define boundary edges, normals etc
            ( Vector2[] xyPosts, Vector2[] boundaryNormals) = DomainManager.CreateBoundary(xMesh,yMesh);

            // Estimate the boundary velocity, given the velocity array (now in m/s)
            Func<double, double, (double, double)> vCalcMPS = (double x, double y) => DomainManager.VelocityFromFixedSpaceArray(x,y,xMin,dx,yMin,dy,xSpeed,ySpeed,true);
            Vector2[] vBoundary = DomainManager.GetBoundaryVelocities(xyPosts, vCalcMPS);

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
                (double[] xSet, double[] ySet, massSurplus) = DomainManager.SeedBoundary(kgPerPoint, boundaryLengths, pressureDelta, dt,
                    xyPosts, boundaryNormals, vBoundary, RNG, massSurplus);

                pointManager.CreatePointSet(xSet, ySet);

                // Do the actual work
                if (debug) {Console.WriteLine($"TIME: {tCurr,7:f2}");}
                pointManager.Advance(dt);

                DomainManager.Cull(xLims,yLims,pointManager);

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
    }
}