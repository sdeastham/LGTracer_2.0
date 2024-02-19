using System.Numerics;

namespace LGTracer
{
    public class LGPoint
    {
        protected Vector2 _location
        { get; set; }

        public Vector2 InitialLocation
        { get; protected set; }

        public Func<double, double, (double, double)> VelocityCalc
        { get; protected set; }

        public bool Active
        { get; protected set; }

        public uint UID
        { get; protected set; }

        // Convenience properties
        public double X
        { 
            get => _location.X;
            set => _location = new Vector2((float)value,_location.Y);
        }
        public double Y
        { 
            get => _location.Y;
            set => _location = new Vector2(_location.X,(float)value);
        }

        public LGPoint( Func<double, double, (double, double)> vCalc )
        {
            // Point starts inactive
            this._location = new Vector2(float.NaN,float.NaN);
            this.InitialLocation = new Vector2(float.NaN,float.NaN);
            this.VelocityCalc = vCalc;
            this.UID = 0; // Reserved for inactive points
            this.Active = false;
        }

        public void Activate( double x, double y, uint uniqueID )
        {
            // Change the particle from being inactive to active
            this.Active = true;
            this._location = new Vector2((float)x, (float)y);
            // Copy this data for later comparison
            this.InitialLocation = new Vector2((float)x, (float)y);
            this.UID = uniqueID;
        }

        public void Deactivate()
        {
            // Kill the particle and put it into storage
            this.Active = false;
            this.X = double.NaN;
            this.Y = double.NaN;
            this.UID = 0;
        }

        public void Advance( double dt )
        {
            if (!Active)
            {
                return;
            }

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