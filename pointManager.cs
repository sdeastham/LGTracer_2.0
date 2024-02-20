//using MathNet.Numerics;
//using MathNet.Numerics.LinearAlgebra;

using System;
using System.Linq;
using System.Collections.Generic;

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

        // Calculate all values and return them as a struct
        public Func<double, double, InterpolatedProperties> ValueCalc
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
    }

    public struct InterpolatedProperties
    {
        public InterpolatedProperties(double temperature, double specificHumidity)
        {
            Temperature = temperature;
            SpecificHumidity = specificHumidity;
        }
        public double Temperature{ get; }
        public double SpecificHumidity{ get; }
    }
}