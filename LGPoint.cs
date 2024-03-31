using System.Numerics;

namespace LGTracer
{
    public class LGPoint
    {
        protected Vector3 _location
        { get; set; }

        public Vector3 InitialLocation
        { get; protected set; }

        public Func<double, double, double, (double, double, double)> VelocityCalc
        { get; protected set; }

        public bool Active
        { get; protected set; }

        public uint UID
        { get; protected set; }

        // This might need to be made more flexible
        public double Temperature
        { get; protected set; }

        public double SpecificHumidity // kg water vapor per kg air
        { get; protected set; }
        
        public double LiquidWaterContent // kg water per kg air
        { get; protected set; }
        
        public double IceWaterContent // kg ice per kg air
        { get; protected set; }

        public double TotalWaterContent => IceWaterContent + LiquidWaterContent + SpecificHumidity;

        public double RelativeHumidityLiquid => Physics.RelativeHumidityLiquid(Temperature, Pressure, SpecificHumidity);
        
        public double RelativeHumidityIce => Physics.RelativeHumidityIce(Temperature, Pressure, SpecificHumidity);

        public double Age
        { get; protected set; }

        // Convenience properties
        public double X
        { 
            get => _location.X;
            set => _location = new Vector3((float)value,_location.Y,_location.Z);
        }
        public double Y
        { 
            get => _location.Y;
            set => _location = new Vector3(_location.X,(float)value,_location.Z);
        }

        public double Pressure
        { 
            get => _location.Z;
            set => _location = new Vector3(_location.X,_location.Y,(float)value);
        }

        private bool IncludeCompression;

        public const double GammaRatio = 0.4/1.4;

        public LGPoint( Func<double, double, double, (double, double, double)> vCalc, bool includeCompression = false )
        {
            // Point starts inactive
            _location = new Vector3(float.NaN,float.NaN,float.NaN);
            InitialLocation = new Vector3(float.NaN,float.NaN,float.NaN);
            VelocityCalc = vCalc;
            // Set the rest of the properties by deactivating the point
            Deactivate();
            // Do we want to calculate the effect of adiabatic compression on temperature?
            IncludeCompression = includeCompression;
        }

        public virtual void Activate( double x, double y, double pressure, uint uniqueID )
        {
            // Change the particle from being inactive to active
            Active = true;
            _location = new Vector3((float)x, (float)y, (float)pressure);
            // Copy this data for later comparison
            InitialLocation = new Vector3((float)x, (float)y, (float)pressure);
            UID = uniqueID;
            Age = 0.0;
        }

        public virtual void Deactivate()
        {
            // Kill the particle and put it into storage
            Active = false;
            _location = new Vector3(float.NaN,float.NaN,float.NaN);
            UID = 0;
            // Properties
            Temperature = double.NaN;
            SpecificHumidity = double.NaN;
            Age = double.NaN;
        }

        public virtual double GetProperty(string property)
        {
            // Hard code for now - the more flexible this is made, the slower it will get
            switch (property.ToLower().Replace("_",""))
            {
                case "temperature":
                    return Temperature;
                case "specifichumidity":
                case "qv":
                    return SpecificHumidity;
                case "relativehumidityliquid":
                case "rhl":
                    return RelativeHumidityLiquid;
                case "relativehumidityice":
                case "rhi":
                    return RelativeHumidityIce;
                default:
                    throw new ArgumentException($"No property for LGPoint called {property}");
            }
        }

        public void SetTemperature( double value )
        {
            Temperature = value;
        }

        public void SetSpecificHumidity( double value )
        {
            SpecificHumidity = value;
        }

        public virtual void Advance( double dt )
        {
            if (!Active)
            {
                return;
            }

            // RK4
            ( double dx1, double dy1, double dp1 ) = VelocityCalc( X, Y, Pressure );
            ( double dx2, double dy2, double dp2 ) = VelocityCalc( X + dx1*dt/2.0, Y + dy1*dt/2.0, Pressure + dp1*dt/2.0);
            ( double dx3, double dy3, double dp3 ) = VelocityCalc( X + dx2*dt/2.0, Y + dy2*dt/2.0, Pressure + dp2*dt/2.0);
            ( double dx4, double dy4, double dp4 ) = VelocityCalc( X + dx3*dt, Y + dy3*dt, Pressure + dp3*dt);
            double xSpeed = (1.0/6.0) * (dx1 + 2.0*dx2 + 2.0*dx3 + dx4);
            double ySpeed = (1.0/6.0) * (dy1 + 2.0*dy2 + 2.0*dy3 + dy4);
            double pSpeed = (1.0/6.0) * (dp1 + 2.0*dp2 + 2.0*dp3 + dp4);
            X += xSpeed * dt;
            Y += ySpeed * dt;
            double oldPressure = Pressure;
            double newPressure = oldPressure + (pSpeed * dt);
            Pressure = newPressure;
            Age += dt;

            if (IncludeCompression)
            {
                Temperature = Temperature * Math.Pow(newPressure/oldPressure,GammaRatio);
            }
        }

    }
}