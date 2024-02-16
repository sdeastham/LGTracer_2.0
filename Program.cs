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

            // Set up the mesh (units of meters)
            double xMin = -100.0;
            double xMax =  100.0;
            double xRange = xMax - xMin;
            int xPosts = 201;
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
            double tStop = 1000.0;
            double dt = 1.0; // Time step in seconds
            double tCurr = tStart;
            int iterMax = (int)Math.Ceiling((tStop - tStart)/dt);

            // How often to save out data?
            double dtStorage = 10.0; // Storage period (seconds)
            double tStorage = 0.0; // Next time that we want storage to occur

            // Number of Lagrangian points to track
            int nPoints = 10;

            // This should all be put into a class for a single point and held in a list of points
            List<LGPoint> points = [];

            // Get random values between zero and 1
            System.Random rng = SystemRandomSource.Default;
            (double[] xInitial, double[] yInitial) = MapRandomToXY(xMin,xMax,yMin,yMax,rng,nPoints);
            for (int i=0;i<nPoints;i++)
            {
                points.Add(new LGPoint(xInitial[i],yInitial[i],(double x, double y) => VelocitySolidBody(x,y,omega)));
            }
        
            // Set up output
            // Add a 2D variable
            // Define dimensions (variables)
            List<double> time = [];

            // X of all points will be stored as a 2D array
            // Allow for 
            List<double[]> xHistory = [];
            List<double[]> yHistory = [];

            // Store initial conditions
            ArchiveConditions(time,xHistory,yHistory,tCurr,points);
            tStorage += dtStorage;

            bool loopOK = true;
            double xCurr, yCurr;
            double rErrPcg, thetaErrDeg, radius, thetaDeg;
            int iPoint;
            
            Console.WriteLine($"Rotational speed: {dt * omega * 180.0/Math.PI,7:f2} deg/timestep");
            for (int iter=0;iter<iterMax; iter++)
            {
                //for (int i=0;i<nPoints;i++)
                iPoint = 0; // Counter, why not
                foreach (LGPoint point in points)
                {
                    // Get the current location of the point and check that it's within bounds
                    xCurr = point.X;
                    yCurr = point.Y;

                    if ((xCurr < xMin) || (xCurr >= xMax) || (yCurr < yMin) || (yCurr >= yMax))
                    {
                        loopOK = false;
                        break;
                    }

                    // If so desired - evaluate current error value
                    if (iPoint == 0 && (iter%10) == 0)
                    {
                        (radius, rErrPcg, thetaDeg, thetaErrDeg) = RThetaError(point.InitialLocation[0],point.InitialLocation[1],omega,tCurr,point.X,point.Y);
                        Console.WriteLine($" Time {tCurr,5:f1} --> Point {iPoint,3:d} at X {point.X,6:f1}/ Y {point.Y,6:f1}/ R {radius,6:f1}/ THETA {thetaDeg,6:f1}/ R ERR {rErrPcg,6:f2}%/ THETA ERR {thetaErrDeg,8:f4}");
                    }

                    point.Advance(dt);
                    iPoint += 1;
                }
                if (!loopOK)
                {
                    Console.WriteLine($"Stopped early at time {tCurr} due to out-of-domain point");
                    break;
                }
                tCurr = (iter+1) * dt;

                // For diagnostics - must take place AFTER tCurr advances
                // Only store data every dtStorage seconds. Use a small offset
                // to compensate for imperfect float comparisons
                if (tCurr >= (tStorage - 1.0e-10))
                {
                    ArchiveConditions(time,xHistory,yHistory,tCurr,points);
                    tStorage += dtStorage;
                }
            }
            if (loopOK)
            {
                Console.WriteLine($"Stopped at time {tCurr}");
            }

            bool success = WriteToFile(dsUri,time,xHistory,yHistory);
            if (success)
            {
                Console.WriteLine($"Output data successfully written to {fileName}");
            }
            else
            {
                Console.WriteLine($"Could not write output to {fileName}");
            }

        }
        private static bool WriteToFile(NetCDFUri dsUri, List<double> time, List<double[]> xHistory, List<double[]> yHistory)
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

            for (int i=0; i<nTimes; i++)
            {
                for (int j=0; j<nPoints; j++)
                {
                    x2D[i,j] = xHistory[i][j];
                    y2D[i,j] = yHistory[i][j];
                }
            }

            using (DataSet ds = DataSet.Open(dsUri))
            {
                ds.AddAxis("index","-",index);
                ds.AddAxis("time","seconds",time.ToArray());
                ds.AddVariable(typeof(double), "x", x2D, ["time","index"]);
                ds.AddVariable(typeof(double), "y", y2D, ["time","index"]);
                ds.Commit();
            }
            
            return success;
        }
        private static void ArchiveConditions(List<double> time, List<double[]> xHistory, List<double[]> yHistory,double tCurr, List<LGPoint> points)
        {
            int nPoints = points.Count;
            double[] xPoints = new double[nPoints];
            double[] yPoints = new double[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xPoints[i] = points[i].X;
                yPoints[i] = points[i].Y;
            }
            time.Add(tCurr);
            xHistory.Add(xPoints);
            yHistory.Add(yPoints);
        }

        private static (double[], double[]) MapRandomToXY( double xMin, double xMax, double yMin, double yMax, System.Random rng, int nPoints )
        {
            // Just map to the inner 50% of each range
            double xRange = xMax - xMin;
            double xStart = xMin + (xRange/4.0);
            double xInner = xRange/2.0;
            double yRange = yMax - yMin;
            double yStart = yMin + (yRange/4.0);
            double yInner = yRange/2.0;

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