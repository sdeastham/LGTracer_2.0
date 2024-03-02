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
        { get; private set; }
        
        public bool TimeDependentMeteorology
        { get; private set; }
        
        public bool Debug
        { get; private set; }
    }

    public struct LGOptionsTiming
    {
        public DateTime StartDate;
        public DateTime EndDate;
    }

    public struct LGOptionsTimesteps
    {
        public double Simulation;
        public double Reporting;
        public double Storage;
    }
    
    public struct LGOptionsPoints
    {
        public long? Max;
        public long Initial;
        public bool AdiabaticCompression;
        public double KgPerPoint;
        public string OutputFilename;
    }

    public class LGOptionsDomain
    {
        public double EastLimit;
        public double WestLimit;
        public double SouthLimit;
        public double NorthLimit;
        public double PressureBase;
        public double PressureCeiling;
        public double[] LonLimits => new double[] { EastLimit, WestLimit };
        public double[] LatLimits => new double[] { SouthLimit, NorthLimit };
    }

    public struct LGOptionsIO
    {
        public string MetDirectory;
        public string OutputDirectory;
    }
}