using System.Diagnostics;
using MathNet.Numerics.Random;

using YamlDotNet.RepresentationModel;
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

            // Central RNG for random point seeding
            System.Random RNG;
            if (configOptions.Seed != null)
            {
                // Use this if debugging
                RNG = new SystemRandomSource((int)configOptions.Seed);
            }
            else
            {
                RNG = SystemRandomSource.Default;
            }
            
            List<PointManager> pointManagers = [];

            for (int i = 0; i < 1; i++)
            {
                // The point manager holds all the actual point data and controls velocity calculations (in deg/s)
                string outputFileName = Path.Join(configOptions.InputOutput.OutputDirectory,
                    configOptions.Points.OutputFilename);
                PointManager pointManager = new PointManager(configOptions.Points.Max, domainManager, outputFileName,
                    includeCompression: configOptions.Points.AdiabaticCompression);

                // Scatter N points randomly over the domain
                (double[] xInitial, double[] yInitial, double[] pInitial) =
                    domainManager.MapRandomToXYP(configOptions.Points.Initial, RNG);
                pointManager.CreatePointSet(xInitial, yInitial, pInitial);
                
                // Store initial conditions
                pointManager.ArchiveConditions(tCurr);

                // Add to the list of point managers
                pointManagers.Add(pointManager);
            }

            tStorage += dtStorage;
            int nStored = 1;

            // Don't report at initialization
            tReport += dtReport;

            // We only add an integer number of points each time step
            // If the number of points to be added is non-integer, retain
            // the surplus and add it in at the next time step
            double massSurplus = 0.0;
            
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

                // If we have enough points available, scatter them evenly over the edges of the domain
                // WARNING: In testing
                foreach (PointManager pointManager in pointManagers)
                {
                    (double[] xSet, double[] ySet, double[] pSet, massSurplus) =
                        domainManager.SeedBoundary(kgPerPoint, dt, RNG, massSurplus);
                    pointManager.CreatePointSet(xSet, ySet, pSet);

                    (double[] xSetV, double[] ySetV, double[] pSetV, massSurplus) =
                        domainManager.SeedPressureBoundaries(kgPerPoint, dt, RNG, massSurplus);
                    pointManager.CreatePointSet(xSetV, ySetV, pSetV);

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