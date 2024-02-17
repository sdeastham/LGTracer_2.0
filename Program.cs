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
            int nPoints = 50000;
            int nInitial = 100; // Points to initially scatter randomly
            float pointRate = 50; // Number of new points to add in each second

            // Set up the mesh (units of meters)
            double xMin = -200.0;
            double xMax =  200.0;
            double xRange = xMax - xMin;
            int xPosts = 401;
            int xCells = xPosts - 1;
            double dx = xRange/(xPosts - 1);

            double yMin = -200.0;
            double yMax =  200.0;
            double yRange = yMax - yMin;
            int yPosts = 401;
            int yCells = yPosts - 1;
            double dy = yRange/(yPosts - 1);

            // Rate of solid body rotation
            double omega = 0.05; // rad/s

            // X edges (meters)
            Vector<double> xMesh = Vector<double>.Build.Dense(xPosts);
            Vector<double> xMid = Vector<double>.Build.Dense(xCells);
            for (int i=0;i<xPosts;i++)
            {
                xMesh[i] = xMin + (dx * i);
                if (i > 0)
                {
                    xMid[i-1] = (xMesh[i] + xMesh[i-1])/2.0;
                }
            }

            // Y edges (meters)
            Vector<double> yMesh = Vector<double>.Build.Dense(yPosts);
            Vector<double> yMid = Vector<double>.Build.Dense(yCells);
            for (int i=0;i<yPosts;i++)
            {
                yMesh[i] = yMin + (dy * i);
                if (i > 0)
                {
                    yMid[i-1] = (yMesh[i] + yMesh[i-1])/2.0;
                }
            }

            // "Wind speed" in the x direction: m/s
            Matrix<double> xSpeed = Matrix<double>.Build.Dense(yCells,xCells);

            // "Wind speed" in the y direction: m/s
            Matrix<double> ySpeed = Matrix<double>.Build.Dense(yCells,xCells);

            for (int i=0;i<xCells;i++)
            {
                for (int j=0;j<yCells;j++)
                {
                    (xSpeed[j,i], ySpeed[j,i]) = VelocitySolidBody(xMid[i],yMid[j],omega);
                }
            }

            double tStart = 0.0;
            double tStop = 100.0;
            double dt = 0.1; // Time step in seconds
            double tCurr = tStart;
            int iterMax = (int)Math.Ceiling((tStop - tStart)/dt);

            // How often to save out data?
            double dtStorage = 1.0; // Storage period (seconds)


            // CODE STARTS HERE

            double tStorage = 0.0; // Next time that we want storage to occur

            // This should all be put into a class for a single point and held in a list of points
            List<LGPoint> points = [];

            // Get random values between zero and 1
            System.Random rng = SystemRandomSource.Default;
            (double[] xInitial, double[] yInitial) = MapRandomToXY(xMin,xMax,yMin,yMax,rng,nPoints);
            uint index = 0;
            for (int i=0;i<nPoints;i++)
            {
                //points.Add(new LGPoint(xInitial[i],yInitial[i],(double x, double y) => VelocitySolidBody(x,y,omega),index));
                //points.Add(new LGPoint(xInitial[i],yInitial[i],(double x, double y) => VelocityConst(x,y,5.0,2.0),index));
                points.Add(new LGPoint(xInitial[i],yInitial[i],(double x, double y) => VelocityFromFixedSpaceArray(x,y,xMin,xMax,dx,yMin,yMax,dy,xSpeed,ySpeed),index));
                index += 1;
                if (i >= nInitial)
                {
                    points[i].Deactivate();
                }
            }
        
            // Set up output
            // Add a 2D variable
            // Define dimensions (variables)
            List<double> time = [];

            // X of all points will be stored as a 2D array
            // Allow for 
            List<double[]> xHistory = [];
            List<double[]> yHistory = [];
            List<uint[]> UIDHistory = [];

            // Store initial conditions
            ArchiveConditions(time,xHistory,yHistory,UIDHistory,tCurr,points);
            tStorage += dtStorage;

            bool loopOK = true;
            double xCurr, yCurr;
            double rErrPcg, thetaErrDeg, radius, thetaDeg;
            int iPoint, nInactive, nNew, nAvailable;

            // We only add an integer number of points each time step
            // If the number of points to be added is non-integer, retain
            // the surplus and add it in at the next time step
            double nNewExact;
            double nSurplus = 0.0;
            
            //Console.WriteLine($"Rotational speed: {dt * omega * 180.0/Math.PI,7:f2} deg/timestep");
            for (int iter=0;iter<iterMax; iter++)
            {
                // How many new points will we add (allowing for variable dt)?
                nNewExact = (pointRate * dt) + nSurplus;
                nNew = (int) Math.Floor(nNewExact);
                nSurplus = nNewExact - (double)nNew;

                // Counter of active points
                iPoint = 0; // Counter, why not

                // Audit points: what has left the domain?
                foreach (LGPoint point in points.Where(point => point.Active))
                {
                    // Get the current location of the point and check that it's within bounds
                    xCurr = point.X;
                    yCurr = point.Y;

                    if ((xCurr < xMin) || (xCurr >= xMax) || (yCurr < yMin) || (yCurr >= yMax))
                    {
                        // Retire the point and move on
                        point.Deactivate();
                    }

                    iPoint += 1;
                }
                //Console.WriteLine($"{iPoint,5:d} points of {nPoints,5:d} active");
                nInactive = nPoints - iPoint;

                // We want to introduce nNew at the edge, but we can only go up to nInactive
                nAvailable = Math.Min(nInactive,nNew);
                iPoint = 0;

                // If we have enough points available, scatter them evenly over the edges of the domain
                // TODO: Make this only at locations where we have inbound flow?
                foreach (LGPoint point in points.Where(point => !point.Active))
                {
                    // Stop if we are out of available points
                    if (iPoint >= nAvailable)
                    {
                        break;
                    }
                    // Activate the point at a location randomly chosen from the domain edge
                    // Algorithm below basically goes around the edges of the domain in order
                    xCurr = rng.NextDouble() * ((xRange*2) + (yRange*2));
                    if (xCurr < xRange)
                    {
                        yCurr = yMin + (dy/100.0);
                        xCurr += xMin;
                    }
                    else if (xCurr < (xRange + yRange))
                    {
                        yCurr = yMin + (xCurr - xRange);
                        xCurr = xMax - (dx/100.0);
                    }
                    else if (xCurr < (xRange + yRange + xRange))
                    {
                        yCurr = yMax - (dy/100.0);
                        xCurr = xMin + (xCurr - (xRange + yRange));
                    }
                    else
                    {
                        yCurr = yMin + (xCurr - (xRange + yRange + xRange));
                        xCurr = xMin + (dx/100.0);
                    }
                    point.Activate(xCurr,yCurr,index);
                    index += 1;
                    iPoint += 1;
                }

                // Do the actual work
                foreach (LGPoint point in points.Where(point => point.Active))
                {
                    point.Advance(dt);
                }

                tCurr = (iter+1) * dt;

                // For diagnostics - must take place AFTER tCurr advances
                // Only store data every dtStorage seconds. Use a small offset
                // to compensate for imperfect float comparisons
                if (tCurr >= (tStorage - 1.0e-10))
                {
                    ArchiveConditions(time,xHistory,yHistory,UIDHistory,tCurr,points);
                    tStorage += dtStorage;
                }
            }
            if (loopOK)
            {
                Console.WriteLine($"Stopped at time {tCurr}");
            }

            bool success = WriteToFile(dsUri,time,xHistory,yHistory,UIDHistory);
            if (success)
            {
                Console.WriteLine($"Output data successfully written to {fileName}");
            }
            else
            {
                Console.WriteLine($"Could not write output to {fileName}");
            }

        }
        private static bool WriteToFile(NetCDFUri dsUri, List<double> time, List<double[]> xHistory, List<double[]> yHistory, List<uint[]> UIDHistory)
        {
            bool success = true;

            // Get the output sizes
            int nPoints = xHistory[0].Length;
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
        private static void ArchiveConditions(List<double> time, List<double[]> xHistory, List<double[]> yHistory, List<uint[]> UIDHistory, double tCurr, List<LGPoint> points)
        {
            int nPoints = points.Count;
            double[] xPoints = new double[nPoints];
            double[] yPoints = new double[nPoints];
            uint[] UIDs = new uint[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xPoints[i] = points[i].X;
                yPoints[i] = points[i].Y;
                UIDs[i] = points[i].UID;
            }
            time.Add(tCurr);
            xHistory.Add(xPoints);
            yHistory.Add(yPoints);
            UIDHistory.Add(UIDs);
        }

        private static (double[], double[]) MapRandomToXY( double xMin, double xMax, double yMin, double yMax, System.Random rng, int nPoints )
        {
            // Scatter randomly throughout domain
            double xRange = xMax - xMin;
            double xStart = xMin;
            double xInner = xRange;
            double yRange = yMax - yMin;
            double yStart = yMin;
            double yInner = yRange;

            double[] xInitial = new double[nPoints];
            double[] yInitial = new double[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xInitial[i] = rng.NextDouble()*xInner + xStart;
                yInitial[i] = rng.NextDouble()*yInner + yStart;
            }
            return (xInitial, yInitial);
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
            dxdt = xSpeedArray[yIndex,xIndex];
            dydt = ySpeedArray[yIndex,xIndex];
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
    }
}