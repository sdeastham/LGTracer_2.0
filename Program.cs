using System.Diagnostics;
using MathNet.Numerics.Random;

using YamlDotNet.Serialization;
using YamlDotNet.Serialization.NamingConventions;

namespace LGTracer
{
    public class Program
    {
        private static void Main(string[] args)
        {
            /* LGTracer is designed to simulate movement of points
            through space under the influence of a wind field. */
            Console.WriteLine("Initiating LGTracer program");

            // Read in first argument as configuration file (or use default)
            string configFile = "config.yaml";
            if (args.Length > 0)
            {
                configFile = args[0];
            }
            LGOptions configOptions = ReadConfig(configFile);
            
            // Extract and store relevant variables
            bool debug = configOptions.Debug;
            bool updateMeteorology = configOptions.TimeDependentMeteorology;

            // Specify the domain
            double[] lonLims = configOptions.Domain.LonLimits;
            double[] latLims = configOptions.Domain.LatLimits;
            double[] pLims   = [configOptions.Domain.PressureBase * 100.0,
                                configOptions.Domain.PressureCeiling * 100.0];
            double kgPerPoint = configOptions.Points.KgPerPoint;

            // Major simulation settings
            DateTime startDate = configOptions.Timing.StartDate;
            DateTime endDate = configOptions.Timing.EndDate;
            double dt = configOptions.Timesteps.Simulation; // Time step in seconds
            double dtStorage = 60.0 * configOptions.Timesteps.Storage; // How often to save out data (seconds)
            double dtReport = 60.0 * configOptions.Timesteps.Reporting; // How often to report to the user?
            
            DateTime currentDate = startDate; // DateTime is a value type so this creates a new copy

            // Set up the meteorology and domain
            MetManager meteorology = new MetManager(configOptions.InputOutput.MetDirectory, lonLims, latLims, startDate);
            (double[] lonEdge, double[] latEdge) = meteorology.GetXYMesh();
            DomainManager domainManager = new DomainManager(lonEdge, latEdge, pLims, MERRA2.AP, MERRA2.BP, meteorology);

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

            // Which point properties will we output?
            string[] densePropertyNames = ["temperature", "relative_humidity_ice", "relative_humidity_liquid", "specific_humidity"];
            string[] flightPropertyNames = ["temperature", "relative_humidity_ice"];
            
            // Use a master RNG to generate seeds predictably
            System.Random masterRNG;
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
            
            for (int i = 0; i < 1; i++)
            {
                // Dense point managers need an RNG for random point seeding
                // This approach is designed to avoid two failure modes:
                // * The relationship between successive managers being consistent (avoided by using a master RNG)
                // * Seeds being reused (avoided by generating until you hit a new seed)
                // The generation-until-new-seed is in theory slow but that would only matter if we were generating
                // a large number (>>>10) of dense point managers, which is not expected to be the case
                int seed;
                do { seed = masterRNG.Next(); } while (seedsUsed.Contains(seed));
                System.Random pmRNG = new SystemRandomSource(seed);
                seedsUsed.Add(seed);
                
                // The point manager holds all the actual point data and controls velocity calculations (in deg/s)
                string outputFileName = Path.Join(configOptions.InputOutput.OutputDirectory,
                    configOptions.Points.OutputFilename);
                PointManager pointManager = new PointManagerDense(configOptions.Points.Max, domainManager,
                    outputFileName, includeCompression: configOptions.Points.AdiabaticCompression,
                    propertyNames: densePropertyNames, rng: pmRNG, kgPerPoint: kgPerPoint);

                // Scatter N points randomly over the domain
                (double[] xInitial, double[] yInitial, double[] pInitial) =
                    domainManager.MapRandomToXYP(configOptions.Points.Initial, pmRNG);
                pointManager.CreatePointSet(xInitial, yInitial, pInitial);
                
                // Store initial conditions
                pointManager.ArchiveConditions(tCurr);

                // Add to the list of _all_ point managers
                pointManagers.Add(pointManager);
            }
            
            // Now add plume point managers - contrail point managers, exhaust point managers...
            // Current proposed approach will be to do this via logical connections (i.e. one manager handles
            // all contrails) rather than e.g. one manager per flight
            for (int i = 0; i < 1; i++)
            {
                
                // The point manager holds all the actual point data and controls velocity calculations (in deg/s)
                string outputFileName = Path.Join(configOptions.InputOutput.OutputDirectory,
                    configOptions.Points.OutputFilename);
                PointManagerFlight pointManager = new PointManagerFlight(configOptions.Points.Max, domainManager,
                    outputFileName,startDate, 
                    includeCompression: configOptions.Points.AdiabaticCompression,
                    propertyNames: densePropertyNames, kgPerPoint: kgPerPoint);
                
                // Add some flights [TESTING]
                double machOneKPH = (3600.0/1000.0) * Math.Sqrt(1.4 * Physics.RGasUniversal * 200.0 / 28.97e-3);
                double lonBOS = -1.0*(71.0 +  0.0 / 60.0 + 23.0 / 3600.0 );
                double latBOS =       42.0 + 21.0 / 60.0 + 47.0 / 3600.0;
                double lonLHR = -1.0*( 0.0 + 27.0 / 60.0 + 41.0 / 3600.0 );
                double latLHR =       51.0 + 28.0 / 60.0 + 39.0 / 3600.0;
                double cruiseSpeedKPH = machOneKPH * 0.8;
                pointManager.SimulateFlight(lonBOS,latBOS,lonLHR,latLHR,startDate,
                    cruiseSpeedKPH,flightLabel: $"BOS_LHR_{startDate}_{endDate}");
                pointManager.SimulateFlight(lonLHR,latLHR,lonBOS,latBOS,startDate,
                    cruiseSpeedKPH,flightLabel: $"LHR_BOS_{startDate}_{endDate}");
                
                // Store initial conditions
                pointManager.ArchiveConditions(tCurr);

                // Add to the list of _all_ point managers
                pointManagers.Add((PointManager)pointManager);
            }

            tStorage += dtStorage;
            int nStored = 1;

            // Don't report at initialization
            tReport += dtReport;
            
            // Set up timing
            int nSteps = 0;
            var watch = new Stopwatch();
            watch.Start();
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
                    pointManager.Seed(dt);
                    
                    // Do the actual work
                    pointManager.Advance(dt);

                    // TODO: Allow for this to not happen every time step
                    pointManager.Cull();
                }

                nSteps++;
                tCurr = (iter+1) * dt;
                currentDate = currentDate.AddSeconds(dt);

                // For diagnostics - must take place AFTER tCurr advances
                // Only store data every dtStorage seconds. Use a small offset
                // to compensate for imperfect float comparisons
                if (tCurr >= (tStorage - 1.0e-10))
                {
                    foreach (PointManager pointManager in pointManagers)
                    {
                        pointManager.ArchiveConditions(tCurr);
                    }
                    tStorage += dtStorage;
                    nStored += 1;
                }
                if (tCurr >= (tReport - 1.0e-10))
                {
                    Console.WriteLine($" --> Time at end of time step: {currentDate}. Current point count in first manager: {pointManagers.First().NActive,10:d}");
                    tReport += dtReport;
                }
            }
            watch.Stop();
            long elapsedTimeLong = watch.ElapsedMilliseconds;
            double elapsedTime = (double)elapsedTimeLong;
            double msPerStep = elapsedTime/nSteps;
            Console.WriteLine($"{nSteps} steps completed in {elapsedTime/1000.0,6:f1} seconds ({msPerStep,6:f2} ms per step)");

            foreach (PointManager pointManager in pointManagers)
            {
                bool success = pointManager.WriteToFile();
                Console.WriteLine(
                    success
                        ? $"Output data with {nStored} samples [max points stored: {pointManager.MaxStoredPoints}] successfully written to {pointManager.OutputFilename}"
                        : $"Could not write output to {pointManager.OutputFilename}");
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
    }
}