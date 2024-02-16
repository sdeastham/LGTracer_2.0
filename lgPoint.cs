//using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace LGTracer
{
    public class LGPoint
    {
        protected Vector<double> _location
        { get; set; }

        public Vector<double> InitialLocation
        { get; protected set; }

        public Func<double, double, (double, double)> VelocityCalc
        { get; protected set; }

        // Convenience properties
        public double X
        { 
            get => _location[0];
            set => _location[0] = value;
        }
        public double Y
        { 
            get => _location[1];
            set => _location[1] = value;
        }

        public LGPoint( double x, double y, Func<double, double, (double, double)> vCalc )
        {
            this._location = Vector<double>.Build.Dense(2);
            this.X = x;
            this.Y = y;
            this.VelocityCalc = vCalc;

            // Copy this data for later comparison
            this.InitialLocation = Vector<double>.Build.Dense(2);
            this.InitialLocation[0] = x;
            this.InitialLocation[1] = y;

        }

        public void Advance( double dt )
        {
            // Forward Euler
            /*
            // Get velocity vector at this location
            ( double xSpeed, double ySpeed ) = VelocityCalc( X, Y );
            // Simple forward Euler for now
            X += xSpeed * dt;
            Y += ySpeed * dt;
            */

            // RK4
            ( double dx1, double dy1 ) = VelocityCalc( X, Y );
            ( double dx2, double dy2 ) = VelocityCalc( X + dx1*dt/2.0, Y + dy1*dt/2.0);
            ( double dx3, double dy3 ) = VelocityCalc( X + dx2*dt/2.0, Y + dy2*dt/2.0);
            ( double dx4, double dy4 ) = VelocityCalc( X + dx3*dt, Y + dy3*dt);
            double xSpeed = (1.0/6.0) * (dx1 + 2.0*dx2 + 2.0*dx3 + dx4);
            double ySpeed = (1.0/6.0) * (dy1 + 2.0*dy2 + 2.0*dy3 + dy4);
            X += xSpeed * dt;
            Y += ySpeed * dt;
        }

    }
}