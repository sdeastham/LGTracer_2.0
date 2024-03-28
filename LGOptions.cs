namespace LGTracer
{
    public class LGOptions
    {
        public LGOptionsTiming Timing
        { get; private set; }
        
        public LGOptionsTimesteps Timesteps
        { get; private set; }
        
        public LGOptionsPoints Points
        { get; private set; }
        
        public LGOptionsDomain Domain
        { get; private set; }
        
        public LGOptionsIO InputOutput
        { get; private set; }

        public int? Seed 
        { get; private set; } = null;

        public bool TimeDependentMeteorology
        { get; private set; } = true;
        public bool Debug
        { get; private set; } = false;

        public LGOptions()
        {
            Points = new LGOptionsPoints();
            Domain = new LGOptionsDomain();
            Timesteps = new LGOptionsTimesteps();
            Timing = new LGOptionsTiming();
            InputOutput = new LGOptionsIO();
        }
    }

    public class LGOptionsTiming
    {
        public DateTime StartDate = new DateTime(2000, 1, 1);
        public DateTime EndDate = new DateTime(2001,1,1);
    }

    public class LGOptionsTimesteps
    {
        public double Simulation = 60.0;
        public double Reporting = 3600.0;
        public double Storage = 3600.0;
    }
    
    public class LGOptionsPoints
    {
        public long? Max = null;
        public long Initial = 0;
        public bool AdiabaticCompression = true;
        public double KgPerPoint = 1.0e16;
        public string OutputFilename = "default_output.nc";
    }

    public class LGOptionsDomain
    {
        public double EastLimit = double.NegativeInfinity;
        public double WestLimit = double.PositiveInfinity;
        public double SouthLimit = double.NegativeInfinity;
        public double NorthLimit = double.PositiveInfinity;
        public double PressureBase = 1000.0;
        public double PressureCeiling = 10.0;
        public double[] LonLimits => new double[] { EastLimit, WestLimit };
        public double[] LatLimits => new double[] { SouthLimit, NorthLimit };
    }

    public class LGOptionsIO
    {
        public string MetDirectory = "";
        public string OutputDirectory = "";
    }
}