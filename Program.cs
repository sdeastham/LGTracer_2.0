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

namespace LGTracer
{
    public class Program
    {
        private static void Main(string[] args)
        {
            /* LGTracer is designed to simulate movement of points
            through space under the influence of a wind field. */
            Console.WriteLine("Initiating LGTracer program");

            // Number of Lagrangian points to track
            int nPoints = 100000;
            int nInitial = 1000; // Points to initially scatter randomly
            bool debug = false;

            // Specify the domains
            // Huge domain
            //double[] lonLims = [-80.0,15.0];
            //double[] latLims = [10.0,60.0];
            //double[] pLims = [85000.0, 20000.0];

            // Moderate domain
            //double[] lonLims = [-80.0,15.0];
            //double[] latLims = [10.0,60.0];
            //double[] pLims = [85000.0, 20000.0];

            // Tiny domain
            double[] lonLims = [-30.0,0.0];
            double[] latLims = [30.0,40.0];
            double[] pLims = [400.0 * 1.0e2, 200.0 * 1.0e2];

            int readTime = 0;
            string metFileNameA3 = "C:/Data/MERRA-2/2023/01/MERRA2.20230101.A3dyn.05x0625.nc4";
            string metFileNameI3 = "C:/Data/MERRA-2/2023/01/MERRA2.20230101.I3.05x0625.nc4";
            string outputFileName = "output.nc";

            // Major simulation settings
            double nDays = 30.0; // Days to run
            double dt = 60.0 * 5.0; // Time step in seconds
            double dtStorage = 60.0*15.0; // // How often to save out data (seconds)

            // Point settings
            double kgPerPoint = 2.0e12; // Air mass represented by a single point (mass flows are huge - but 1e11 seems very high? Total atm mass 5.1e18!)


            // CODE STARTS HERE

            // Read in MERRA-2 data and use that to set domain
            (double[] lonEdge, double[] latEdge, int[] lonSet, int[] latSet ) = MERRA2.ReadLatLon( metFileNameA3, lonLims, latLims );
            // Now read in the U and V data, as well as pressure velocities
            (double[,,]uWind, double[,,]vWind, double[,,] pressureVelocity ) = MERRA2.ReadA3( metFileNameA3, readTime, lonSet, latSet );
            // Also read in PS, T, and QV
            (double[,] surfacePressure, double[,,]griddedTemperature, double[,,]griddedSpecificHumidity ) = MERRA2.ReadI3( metFileNameI3, readTime, lonSet, latSet );
            
            // Set up the domain
            DomainManager domainManager = new DomainManager(lonEdge, latEdge, pLims, MERRA2.AP, MERRA2.BP);

            // Set up the domain meteorology
            domainManager.UpdateMeteorology(surfacePressure,uWind,vWind,pressureVelocity,griddedTemperature,griddedSpecificHumidity);

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
            PointManager pointManager = new PointManager(nPoints,domainManager);

            // Scatter N points randomly over the domain
            (double[] xInitial, double[] yInitial, double[] pInitial) = domainManager.MapRandomToXYP(nInitial,RNG);
            pointManager.CreatePointSet(xInitial,yInitial,pInitial);

            // Set up output
            // Store initial conditions
            pointManager.ArchiveConditions(tCurr);
            tStorage += dtStorage;
            int nStored = 1;

            // We only add an integer number of points each time step
            // If the number of points to be added is non-integer, retain
            // the surplus and add it in at the next time step
            double massSurplus = 0.0;
            
            // Set up timing
            int nSteps = 0;
            var watch = new Stopwatch();
            watch.Start();
            Console.WriteLine("Beginning trajectory calculation");
            for (int iter=0;iter<iterMax; iter++)
            {
                // If we have enough points available, scatter them evenly over the edges of the domain
                // WARNING: This does not yet handle a vertical domain which is not uniform!
                //(double[] xSet, double[] ySet, double[] pSet, massSurplus) = domainManager.SeedBoundary(kgPerPoint, dt, RNG, massSurplus);
                //pointManager.CreatePointSet(xSet, ySet, pSet);

                (double[] xSetV, double[] ySetV, double[] pSetV, massSurplus) = domainManager.SeedPressureBoundaries(kgPerPoint, dt, RNG, massSurplus);
                pointManager.CreatePointSet(xSetV, ySetV, pSetV);

                // Do the actual work
                if (debug) {Console.WriteLine($"TIME: {tCurr,7:f2}");}
                pointManager.Advance(dt);

                pointManager.Cull();

                nSteps++;
                tCurr = (iter+1) * dt;

                // For diagnostics - must take place AFTER tCurr advances
                // Only store data every dtStorage seconds. Use a small offset
                // to compensate for imperfect float comparisons
                if (tCurr >= (tStorage - 1.0e-10))
                {
                    pointManager.ArchiveConditions(tCurr);
                    tStorage += dtStorage;
                    nStored += 1;
                }
            }
            watch.Stop();
            long elapsedTimeLong = watch.ElapsedMilliseconds;
            double elapsedTime = (double)elapsedTimeLong;
            double msPerStep = elapsedTime/nSteps;
            Console.WriteLine($"{nSteps} steps completed in {elapsedTime/1000.0,6:f1} seconds ({msPerStep,6:f2} ms per step)");

            double xMin = domainManager.XMax;
            foreach (LGPoint point in pointManager.ActivePoints)
            {
                xMin = Math.Min(xMin,point.X);
            }

            bool success = pointManager.WriteToFile(outputFileName);
            if (success)
            {
                Console.WriteLine($"Output data with {nStored} samples [max points stored: {pointManager.MaxStoredPoints}] successfully written to {outputFileName}");
            }
            else
            {
                Console.WriteLine($"Could not write output to {outputFileName}");
            }
        }
    }
}