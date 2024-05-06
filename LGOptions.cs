namespace LGTracer;

public class LGOptions
{
    public LGOptionsTiming Timing
    { get; private set; }
        
    public LGOptionsTimesteps Timesteps
    { get; private set; }
        
    public LGOptionsPointsDense PointsDense
    { get; private set; }
        
    public LGOptionsPointsFlights PointsFlights
    { get; private set; }
        
    public LGOptionsDomain Domain
    { get; private set; }
        
    public LGOptionsIO InputOutput
    { get; private set; }

    public int? Seed 
    { get; private set; } = null;

    public bool TimeDependentMeteorology
    { get; private set; } = true;
    public bool Verbose
    { get; private set; } = false;

    public LGOptions()
    {
        PointsDense = new LGOptionsPointsDense();
        PointsFlights = new LGOptionsPointsFlights();
        Domain = new LGOptionsDomain();
        Timesteps = new LGOptionsTimesteps();
        Timing = new LGOptionsTiming();
        InputOutput = new LGOptionsIO();
    }

    public static double ParseHms(string hhmmss)
    {
        // Expects a string in the form hh:mm:ss
        TimeSpan ts = new TimeSpan(int.Parse(hhmmss.Split(':')[0]),  // hours
                                int.Parse(hhmmss.Split(':')[1]),   // minutes
                                int.Parse(hhmmss.Split(':')[2])); // seconds
        return ts.TotalSeconds;
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
    public string Output = "240000";
}

public abstract class LGOptionsPoints
{
    public long? Max = 0;
    public bool AdiabaticCompression = true;
    public bool WritePeriodic = true;
    public string OutputFilename = "default_output_{date}.nc";
    public string[] OutputVariables = [];
    public bool WriteTrajectories = false;
    public string TrajectoryFilename = "default_trajectory_{date}_{uid}.nc";
    public string[] TrajectoryVariables = [];

    public bool Active => (Max == null || Max > 0);
}
public class LGOptionsPointsDense : LGOptionsPoints
{
    public long Initial = 0;
    public double KgPerPoint = 1.0e16;
}
    
public class LGOptionsPointsFlights : LGOptionsPoints
{
    public string? ScheduleFilename = null;
    public string? AirportsFilename = null;
    public string? SegmentsFilename = null;
    public double PointSpacing = 60.0; // Seconds
    public string SegmentsOutputFilename = "segments_output.nc";
    public bool ContrailSimulation = false;
    public bool IncludeSettling = false;
    public bool ComplexContrails = false;
    public bool UseIcao = true;
    public double MinimumLifetime = 0.0;
    public bool SkipNewtonIterationForTlm = false;
    public bool UsePonaterTlc = false;
    public string TrajectoryIdentifierFilename = "flightmatch-{date}.parquet";
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
    public string InputDirectory = "";
    public bool SerialMetData = false;
    public string MetSource = "MERRA-2";
}