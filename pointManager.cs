//using MathNet.Numerics;
//using MathNet.Numerics.LinearAlgebra;

using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

namespace LGTracer
{
    public class PointManager
    {
        public List<LGPoint> ActivePoints
        { get; private set; }

        private List<LGPoint> InactivePoints
        { get; set; }

        private uint nextUID
        { get; set; }

        private uint NextUID
        { get {nextUID += 1; return nextUID - 1;}}

        public int NPoints
        { get => NActive + NInactive; }

        public int NActive
        { get => ActivePoints.Count; }

        public int NInactive
        { get => InactivePoints.Count; }

        public int MaxPoints
        { get; private set; }

        public Func<double, double, (double, double)> VelocityCalc
        { get; protected set; }

        private bool Debug
        { get; set; }

        public PointManager( int maxPoints, Func<double, double, (double, double)> vCalc, bool debug=false )
        {
            // UIDs start from 1 (0 reserved for inactive points)
            nextUID = 1;

            // Set the velocity calculation
            VelocityCalc = vCalc;

            // Limit on how many points can be managed
            MaxPoints = maxPoints;

            ActivePoints = [];
            InactivePoints = [];

            // Run in debug mode?
            Debug = debug;
        }

        private LGPoint AddPoint( double x, double y )
        {
            // Create a new point in the list
            // Start by creating an _inactive_ point
            LGPoint point = new LGPoint(VelocityCalc);
            InactivePoints.Add(point);
            // Activate a point (doesn't matter if it's the same one) and return it
            return ActivatePoint(x,y);
        }

        private LGPoint ActivatePoint( double x, double y )
        {
            // Reactivate the first available dormant point and assign it a new UID
            LGPoint point = InactivePoints[0];
            ActivePoints.Add(point);
            InactivePoints.RemoveAt(0);

            // Give the point its location and a new UID
            // Requesting the UID will automatically increment it
            point.Activate(x,y,NextUID);
            return point;
        }

        public LGPoint NextPoint( double x, double y )
        {
            // Function places a new point, taken from the inactive list if any available.
            // If no points are available, add one if possible; otherwise throw an exception

            // Are there any inactive points available?
            LGPoint point;
            if (NInactive > 0)
            {
                // Reactivate a dormant point
                point = ActivatePoint(x,y);
            }
            else if (NActive < MaxPoints)
            {
                // Add a new point
                point = AddPoint(x,y);
            }
            else
            {
                // No more points to return!
                //if (Debug) {Console.WriteLine("!!!");}
                throw new InvalidOperationException("Point maximum exceeded");
            }
            // Set point properties based on local values
            

            return point;
        }

        public void DeactivatePoint( int index )
        {
            // Deactivate point i of those present in ActivePoints
            LGPoint point = ActivePoints[index];
            InactivePoints.Add(point);
            point.Deactivate();
            ActivePoints.RemoveAt(index);
        }

        public void CreatePointSet( double[] x, double[] y )
        {
            // Create multiple points
            for (int i=0; i<x.Length; i++)
            {
                NextPoint(x[i],y[i]);
            }
        }

        public void Advance( double dt )
        {
            // Advances all active points one time step
            foreach (LGPoint point in ActivePoints)
            {
                //Console.WriteLine($"{point.UID,5:d}: {point.X,7:f2}/{point.Y,7:f2}");
                point.Advance(dt);
            }
        }

        public bool WriteToFile(string fileName, List<double> time, List<double[]> xHistory, List<double[]> yHistory, List<uint[]> UIDHistory, int maxActive)
        {
            bool success = true;

            // Set up output file
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.Create
            };

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

        public int ArchiveConditions(List<double> time, List<double[]> xHistory, List<double[]> yHistory, List<uint[]> UIDHistory, double tCurr)
        {
            int nPoints = this.MaxPoints;
            double[] xPoints = new double[nPoints];
            double[] yPoints = new double[nPoints];
            uint[] UIDs = new uint[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                if (i<this.NActive)
                {
                    LGPoint point = this.ActivePoints[i];
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
            return this.NActive;
        }
    }
}