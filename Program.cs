using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.Integration;
using MathNet.Numerics.Random;
using MathNet.Numerics.Distributions;

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
            double pointRate = 50.0/3600.0; // Number of new points to add in each second
            bool debug = false;

            // Specify the domains
            double[] lonLims = [-60.0,60.0];
            double[] latLims = [-60.0,60.0];
            int readLevel = 30;
            int readTime = 0;
            string metFileName = "C:/Data/MERRA-2/2023/01/MERRA2.20230101.A3dyn.05x0625.nc4";

            // Major simulation settings
            double nDays = 30.0; // Days to run
            double dt = 60.0 * 5.0; // Time step in seconds
            double dtStorage = 60.0*15.0; // // How often to save out data (seconds)

            // Read in MERRA-2 data and use that to set domain
            (double[] lonEdge, double[] latEdge, double[,]uWind, double[,]vWind ) = ReadMERRA2( metFileName, readTime, readLevel, lonLims, latLims );

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
            Vector<double> xMesh = Vector<double>.Build.Dense(xPosts);
            Vector<double> xMid = Vector<double>.Build.Dense(xCells);
            for (int i=0;i<xPosts;i++)
            {
                //xMesh[i] = xMin + (dx * i);
                xMesh[i] = lonEdge[i];
                if (i > 0)
                {
                    xMid[i-1] = (xMesh[i] + xMesh[i-1])/2.0;
                }
            }

            // Y edges (degrees)
            Vector<double> yMesh = Vector<double>.Build.Dense(yPosts);
            Vector<double> yMid = Vector<double>.Build.Dense(yCells);
            for (int i=0;i<yPosts;i++)
            {
                //yMesh[i] = yMin + (dy * i);
                yMesh[i] = latEdge[i];
                if (i > 0)
                {
                    yMid[i-1] = (yMesh[i] + yMesh[i-1])/2.0;
                }
            }

            // Create the velocity array in m/s
            Matrix<double> xSpeed = Matrix<double>.Build.DenseOfArray(uWind);
            Matrix<double> ySpeed = Matrix<double>.Build.DenseOfArray(vWind);
            for (int i=0;i<xCells;i++)
            {
                for (int j=0;j<yCells;j++)
                {
                    //(xSpeed[j,i], ySpeed[j,i]) = VelocitySolidBody(xMid[i],yMid[j],omega);
                    xSpeed[j,i] = uWind[j,i];// * roughWindConversion;
                    ySpeed[j,i] = vWind[j,i];// * roughWindConversion;
                }
            }


            // CODE STARTS HERE
            double duration = 60.0 * 60.0 * 24.0 * nDays; // Simulation duration in seconds
            double tStart = 0.0;
            double tStop = tStart + duration;
            double tCurr = tStart;
            int iterMax = (int)Math.Ceiling((tStop - tStart)/dt);
            double tStorage = tStart; // Next time that we want storage to occur

            // The point manager holds all the actual point data and controls velocity calculations
            //static void vCalc(double x, double y) => VelocityFromFixedSpaceArray(x,y,xMin,xMax,dx,yMin,yMax,dy,xSpeed,ySpeed);
            //PointManager pointManager = new PointManager(nInitial,nPoints,vCalc);
            PointManager pointManager = new PointManager(nInitial,nPoints,xLims,yLims,
                (double x, double y) => VelocityFromFixedSpaceArray(x,y,xMin,xMax,dx,yMin,yMax,dy,xSpeed,ySpeed));
        
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
                //Console.WriteLine($"Active before seeding: {pointManager.NActive}");
                pointManager.SeedBoundary(nAvailable);
                //Console.WriteLine($"Active after seeding:  {pointManager.NActive}");

                // Do the actual work
                if (debug) {Console.WriteLine($"TIME: {tCurr,7:f2}");}
                pointManager.Advance(dt);

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
                Console.WriteLine($"Output data with {nStored} samples successfully written to {fileName}");
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

        private static (double, double) RThetaAnalytical( double xInitial, double yInitial, double omega, double tCurr )
        {
            ( double radiusInitial, double thetaInitial ) = RThetaFromYX(yInitial, xInitial);
            double theta = thetaInitial + (tCurr * omega);
            return (radiusInitial, theta);
        }

        private static (double, double, double, double) RThetaError( double xInitial, double yInitial, double omega, double tCurr, double x, double y)
        {
            (double radius, double theta) = RThetaFromYX(y,x);
            (double radiusAnalytical, double thetaAnalytical) = RThetaAnalytical(xInitial,yInitial,omega,tCurr);
            double rErrPcg,thetaErrDeg;
            if (radius < 1.0e-10)
            {
                rErrPcg = double.NaN;
            }
            else
            {
                rErrPcg = 100.0 * (radius/radiusAnalytical - 1.0);
            }
            thetaErrDeg = Math.Abs((180.0/Math.PI)*Atan2(Trig.Sin(thetaAnalytical-theta),Trig.Cos(thetaAnalytical-theta)));
            return (radius, rErrPcg, theta * 180.0/Math.PI, thetaErrDeg);
        }
        private static (double, double) VelocitySolidBody( double x, double y, double omega )
        {
            (double radius, double theta) = RThetaFromYX(y,x);
            double dxdt = radius * omega * Trig.Sin(theta) * -1.0;
            double dydt = radius * omega * Trig.Cos(theta);
            return (dxdt, dydt);
        }

        private static (double, double) VelocityConst( double x, double y, double xSpeed, double ySpeed)
        {
            // Return a fixed velocity
            return (xSpeed, ySpeed);
        }
        private static (double, double) VelocityFromFixedSpaceArray( double x, double y, double xMin, double xMax, double dx, double yMin, double yMax, double dy, Matrix<double> xSpeedArray, Matrix<double> ySpeedArray)
        {
            // Radius of Earth in meters
            const double rEarth = 6.371e6;
            const double deg2rad = Math.PI / 180.0;
            const double rad2deg = 180.0 / Math.PI;

            // Extract the velocity vector from an array
            // Assumes constant X spacing and constant Y spacing
            double dxdt, dydt;
            if ((x >= xMax) || (y >= yMax))
            {
                // Propel point out of the domain
                dxdt = 1.0;
                dydt = 1.0;
                return (dxdt, dydt);
            }
            if ((x <= xMin) || (y <= yMin))
            {
                dxdt = -1.0;
                dydt = -1.0;
                return (dxdt, dydt);
            }
            // If we made it this far - we are within the domain
            int xIndex = (int)Math.Floor((x - xMin)/dx);
            int yIndex = (int)Math.Floor((y - yMin)/dy);
            // Convert from m/s to deg/s
            dxdt = rad2deg * xSpeedArray[yIndex,xIndex] / (rEarth * Math.Cos(deg2rad*y));
            dydt = rad2deg * ySpeedArray[yIndex,xIndex] / rEarth;
            return (dxdt, dydt);
        }
        private static (double, double) RThetaFromYX(double y, double x)
        {
            return (Math.Sqrt(x * x + y * y), Atan2(y,x));
        }
        private static double Atan2(double y, double x)
        {
            if (x>0)
            {
                return Trig.Atan(y/x);
            }
            else if ((x < 0) && (y >= 0 ))
            {
                return Trig.Atan(y/x) + Constants.Pi;
            }
            else if ((x < 0) && (y < 0 ))
            {
                return Trig.Atan(y/x) - Constants.Pi;
            }
            else // x == 0
            {
                if (y > 0)
                {
                    return Constants.Pi / 2.0;
                }
                else if (y < 0)
                {
                    return -1.0 * Constants.Pi / 2.0;
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        private static (double [], double[], double[,], double[,] ) ReadMERRA2( string fileName, int time, int level, double[] lonLims, double[] latLims )
        {
            // Returns [lon_edge],[lat_edge],[u],[v]
            // Later extend to get T, QV
            // Open netCDF4 file
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.ReadOnly
            };

            //Func<double,double,double,int> findLower = (targetValue, lowerBound, spacing) => ((int)Math.Floor((targetValue - lowerBound)/spacing));
            double[] lonEdge,latEdge;
            float[] lonMids, latMids;
            double[,] u,v;
            int nLon, nLat, lonFirst, latFirst, lonLast, latLast;
            double dLon, dLat, lonBase, latBase;
            using (DataSet ds = DataSet.Open(dsUri))
            {
                // Get the full dimension vectors
                latMids = ds.GetData<float[]>("lat");
                lonMids = ds.GetData<float[]>("lon");

                // Figure out which cells we need to keep in order to get all the data we need
                // Assume a fixed cell spacing for now
                dLon = lonMids[1] - lonMids[0];
                lonBase = lonMids[0] - (dLon/2.0);
                // For latitude, be careful about half-polar grids
                dLat = latMids[3] - latMids[2];
                latBase = latMids[1] - (3.0*dLon/2.0);

                // These indices are for the first and last cell _inclusive_
                //latFirst = findLower(latLims[0],latBase,dLat);
                //latLast  = findLower(latLims[1],latBase,dLat);
                //lonFirst = findLower(lonLims[0],lonBase,dLon);
                //lonLast  = findLower(lonLims[1],lonBase,dLon);
                latFirst = (int)Math.Floor((latLims[0]-latBase)/dLat);
                latLast  = (int)Math.Floor((latLims[1]-latBase)/dLat);
                lonFirst = (int)Math.Floor((lonLims[0]-lonBase)/dLon);
                lonLast  = (int)Math.Floor((lonLims[1]-lonBase)/dLon);

                nLon = 1 + (lonLast - lonFirst);
                nLat = 1 + (latLast - latFirst);

                u = new double[nLat,nLon];
                v = new double[nLat,nLon];

                // Be lazy for the moment and access the whole array (not certain if this reads into memory or just makes it available?)
                float[,,,] uFull = ds.GetData<float[,,,]>("U");
                float[,,,] vFull = ds.GetData<float[,,,]>("V");

                for (int iLon=0;iLon<nLon;iLon++)
                {
                    for (int iLat=0;iLat<nLat;iLat++)
                    {
                        u[iLat,iLon] = (double)uFull[time,level,iLat + latFirst,iLon + lonFirst];
                        v[iLat,iLon] = (double)vFull[time,level,iLat + latFirst,iLon + lonFirst];
                    }
                }
            }
            // Create lon/lat edge vectors
            lonEdge = new double[nLon+1];
            latEdge = new double[nLat+1];
            lonEdge[0] = lonMids[lonFirst] - (dLon/2.0);
            latEdge[0] = latMids[latFirst] - (dLat/2.0);
            for (int i=0;i<nLon;i++)
            {
                lonEdge[i+1] = lonEdge[0] + (dLon * (i+1));
            }
            for (int i=0;i<nLat;i++)
            {
                latEdge[i+1] = latEdge[0] + (dLat * (i+1));
            }
            
            return (lonEdge, latEdge, u, v);
        }
    }
}