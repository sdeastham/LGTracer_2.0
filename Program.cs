using System.Diagnostics;
using System.Drawing;
using MathNet.Numerics.Random;

using YamlDotNet.Serialization;
using YamlDotNet.Serialization.NamingConventions;

namespace LGTracer;

public class Program
{
    private static void Main(string[] args)
    {
        /* LGTracer is designed to simulate movement of points
        through space under the influence of a wind field. */
        Console.WriteLine("Initiating LGTracer program");

        // Initialize timing
        var watch = new Stopwatch();
        watch.Start();
        Dictionary<string,Stopwatch> subwatches = [];
        foreach (string watchName in (string[])["Point seeding", "Point physics", "Point culling", "Met advance",
                     "Derived quantities", "Met interpolate", "Archiving", "File writing"])
        {
            subwatches.Add(watchName, new Stopwatch());   
        }
        
        // Read in first argument as configuration file (or use default)
        string configFile = "config.yaml";
        if (args.Length > 0)
        {
            configFile = args[0];
        }
        LGOptions configOptions = ReadConfig(configFile);
            
        // Extract and store relevant variables
        bool verbose = configOptions.Verbose;
        bool updateMeteorology = configOptions.TimeDependentMeteorology;

        // Specify the domain
        double[] lonLims = configOptions.Domain.LonLimits;
        double[] latLims = configOptions.Domain.LatLimits;
        double[] pLims   = [configOptions.Domain.PressureBase * 100.0,
            configOptions.Domain.PressureCeiling * 100.0];

        // Major simulation settings
        DateTime startDate = configOptions.Timing.StartDate;
        DateTime endDate = configOptions.Timing.EndDate;
        double dt = configOptions.Timesteps.Simulation; // Time step in seconds
        double dtStorage = 60.0 * configOptions.Timesteps.Storage; // How often to save out data (seconds)
        double dtReport = 60.0 * configOptions.Timesteps.Reporting; // How often to report to the user?
            
        DateTime currentDate = startDate; // DateTime is a value type so this creates a new copy
        
        // Check if the domain manager will need to calculate box heights (expensive)
        bool boxHeightsNeeded = configOptions.PointsFlights is { Active: true, ComplexContrails: true };

        // Set up the meteorology and domain
        MetManager meteorology = new MetManager(configOptions.InputOutput.MetDirectory, lonLims, latLims, startDate, 
            configOptions.InputOutput.SerialMetData, subwatches);
        (double[] lonEdge, double[] latEdge) = meteorology.GetXYMesh();
        DomainManager domainManager = new DomainManager(lonEdge, latEdge, pLims, MERRA2.AP, MERRA2.BP,
            meteorology, subwatches, boxHeightsNeeded);

        // Time handling
        double nDays = (endDate - startDate).TotalDays; // Days to run
        double duration = 60.0 * 60.0 * 24.0 * nDays; // Simulation duration in seconds
        double tStart = 0.0;
        double tStop = tStart + duration;
        double tCurr = tStart;
        int iterMax = (int)Math.Ceiling((tStop - tStart)/dt);
        double tStorage = tStart; // Next time that we want storage to occur
        double tReport = tStart; // Next time we want output to go to the user
            
        List<PointManager> pointManagers = [];

        // Use a master RNG to generate seeds predictably
        Random masterRNG;
        if (configOptions.Seed != null)
        {
            // Use this if debugging
            masterRNG = new SystemRandomSource((int)configOptions.Seed);
        }
        else
        {
            masterRNG = SystemRandomSource.Default;
        }
        List<int> seedsUsed = [];
            
        // Dense point managers need an RNG for random point seeding
        // This approach is designed to avoid two failure modes:
        // * The relationship between successive managers being consistent (avoided by using a master RNG)
        // * Seeds being reused (avoided by generating until you hit a new seed)
        // The generation-until-new-seed is in theory slow but that would only matter if we were generating
        // a large number (>>>10) of dense point managers, which is not expected to be the case
        if (configOptions.PointsDense.Active)
        {
            Random pmRNG = GetNextRNG(masterRNG, seedsUsed);

            double kgPerPoint = configOptions.PointsDense.KgPerPoint;

            // The point manager holds all the actual point data and controls velocity calculations (in deg/s)
            string outputFileName = Path.Join(configOptions.InputOutput.OutputDirectory,
                configOptions.PointsDense.OutputFilename);
            PointManager pointManager = new PointManagerDense(configOptions.PointsDense.Max, domainManager,
                outputFileName, includeCompression: configOptions.PointsDense.AdiabaticCompression,
                propertyNames: configOptions.PointsDense.OutputVariables, rng: pmRNG, kgPerPoint: kgPerPoint,
                verboseOutput: configOptions.Verbose);

            // Scatter N points randomly over the domain
            (double[] xInitial, double[] yInitial, double[] pInitial) =
                domainManager.MapRandomToXYP(configOptions.PointsDense.Initial, pmRNG);
            pointManager.CreatePointSet(xInitial, yInitial, pInitial);

            // Add to the list of _all_ point managers
            pointManagers.Add(pointManager);
        }

        // Now add plume point managers - contrail point managers, exhaust point managers...
        // Current proposed approach will be to do this via logical connections (i.e. one manager handles
        // all contrails) rather than e.g. one manager per flight
        // The point manager holds all the actual point data and controls velocity calculations (in deg/s)
        if (configOptions.PointsFlights.Active)
        {
            string outputFileName = Path.Join(configOptions.InputOutput.OutputDirectory,
                configOptions.PointsFlights.OutputFilename);
            // TODO: Just pass configOptions.PointsFlights to the point manager!
            double pointPeriod = configOptions.PointsFlights.PointSpacing;
            PointManagerFlight pointManager = new PointManagerFlight(configOptions.PointsFlights.Max, domainManager,
                outputFileName, startDate, pointPeriod, configOptions.PointsFlights.SegmentsOutputFilename,
                contrailSimulation: configOptions.PointsFlights.ContrailSimulation,
                includeSettling: configOptions.PointsFlights.IncludeSettling,
                includeCompression: configOptions.PointsFlights.AdiabaticCompression,
                propertyNames: configOptions.PointsFlights.OutputVariables,
                verboseOutput: configOptions.Verbose, useIcao: configOptions.PointsFlights.UseIcao);

            if (configOptions.PointsFlights.ScheduleFilename != null)
            {
                Debug.Assert(configOptions.PointsFlights.AirportsFilename != null,
                    "No airport file provided");
                string scheduleFileName = Path.Join(configOptions.InputOutput.InputDirectory,
                    configOptions.PointsFlights.ScheduleFilename);
                string airportFileName = Path.Join(configOptions.InputOutput.InputDirectory,
                    configOptions.PointsFlights.AirportsFilename);
                pointManager.ReadScheduleFile(scheduleFileName, airportFileName, startDate, endDate);
            }
                
            if (configOptions.PointsFlights.SegmentsFilename != null)
            {
                pointManager.ReadSegmentsFile(configOptions.PointsFlights.SegmentsFilename);
            }
                
            // Add to the list of _all_ point managers
            pointManagers.Add((PointManager)pointManager);
        }

        foreach (PointManager pointManager in pointManagers)
        {
            // Store initial conditions
            pointManager.ArchiveConditions(tCurr);
        }

        if (pointManagers.Count == 0)
        {
            throw new ArgumentException("No point managers enabled.");
        }

        tStorage += dtStorage;
        int nStored = 1;

        // Don't report at initialization
        tReport += dtReport;
        
        int nSteps = 0;
        int nReports = 0;
        double pointSum = 0;
        bool accuratePointAverage = false;
        
        Console.WriteLine("Beginning trajectory calculation");
        for (int iter=0;iter<iterMax; iter++)
        {
            // Update meteorological data
            if (updateMeteorology)
            {
                meteorology.AdvanceToTime(currentDate);
                // Calculate derived quantities
                domainManager.UpdateMeteorology();
            }
            
            foreach (PointManager pointManager in pointManagers)
            {
                // Seed new points
                subwatches["Point seeding"].Start();
                pointManager.Seed(dt);
                subwatches["Point seeding"].Stop();
                    
                // Do the actual work
                subwatches["Point physics"].Start();
                pointManager.Advance(dt);
                subwatches["Point physics"].Stop();

                // TODO: Allow for this to not happen every time step
                subwatches["Point culling"].Start();
                pointManager.Cull();
                subwatches["Point culling"].Stop();
            }

            nSteps++;
            tCurr = (iter+1) * dt;
            currentDate = currentDate.AddSeconds(dt);

            // For diagnostics - must take place AFTER tCurr advances
            // Only store data every dtStorage seconds. Use a small offset
            // to compensate for imperfect float comparisons
            if (tCurr >= (tStorage - 1.0e-10))
            {
                subwatches["Archiving"].Start();
                foreach (PointManager pointManager in pointManagers)
                {
                    pointManager.ArchiveConditions(tCurr);
                }
                subwatches["Archiving"].Stop();
                tStorage += dtStorage;
                nStored += 1;
            }
            if (tCurr >= (tReport - 1.0e-10))
            {
                long totalActive = pointManagers.Sum(pm => pm.NActive);
                Console.WriteLine($" --> Time at end of time step: {currentDate}. Point count across all managers: {totalActive,10:d}");
                tReport += dtReport;
                pointSum += totalActive;
                nReports++;
            }
            else if (accuratePointAverage)
            {
                long totalActive = pointManagers.Sum(pm => pm.NActive);
                pointSum += totalActive;
            }
        }
        subwatches["File writing"].Start();
        foreach (PointManager pointManager in pointManagers)
        {
            bool success = pointManager.WriteToFile();
            Console.WriteLine(
                success
                    ? $"Output data with {nStored} samples [max points stored: {pointManager.MaxStoredPoints}] successfully written to {pointManager.OutputFilename}"
                    : $"Could not write output to {pointManager.OutputFilename}");
        }
        subwatches["File writing"].Stop();
        
        watch.Stop();
        long elapsedTimeLong = watch.ElapsedMilliseconds;
        double elapsedTime = (double)elapsedTimeLong;
        double msPerStep = elapsedTime/nSteps;
        Console.WriteLine($"{nSteps} steps completed in {elapsedTime/1000.0,6:f1} seconds ({msPerStep,6:f2} ms per step)");
        int nPointSums = accuratePointAverage ? nSteps : nReports;
        if (nPointSums >= 1)
        {
            Console.WriteLine($"Simulation average point count (at reporting times): {pointSum / nPointSums}");
        }
        else
        {
            Console.WriteLine("No point average calculated (simulation too short).");
        }

        foreach (string watchName in subwatches.Keys)
        {
            string watchOutput;
            if (subwatches[watchName].IsRunning)
            {
                watchOutput = "error (watch still running)";
            }
            else
            {
                double subwatchTime = (double)subwatches[watchName].ElapsedMilliseconds;
                watchOutput =
                    $"{subwatchTime / 1000.0,12:f2} seconds ({100.0 * subwatchTime / elapsedTime,10:f2}% of total)";
            }
            Console.WriteLine($" --> {watchName,20:s}: {watchOutput}");
        }
    }

    private static LGOptions ReadConfig(string filename)
    {
        string yaml = File.ReadAllText(filename);
        var deserializer = new DeserializerBuilder()
            .WithNamingConvention(HyphenatedNamingConvention.Instance)
            .Build();

        return deserializer.Deserialize<LGOptions>(yaml);
    }
        
    private static Random GetNextRNG(Random masterRNG, List<int> seedsUsed)
    {
        int seed;
        do { seed = masterRNG.Next(); } while (seedsUsed.Contains(seed));
        Random pmRNG = new SystemRandomSource(seed);
        seedsUsed.Add(seed);
        return pmRNG;
    }
}