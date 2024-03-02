using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

namespace LGTracer
{
    public class PointManager
    {
        public LinkedList<LGPoint> ActivePoints
        { get; private set; }

        private LinkedList<LGPoint> InactivePoints
        { get; set; }

        private uint nextUID
        { get; set; }

        private uint NextUID
        { get {nextUID += 1; return nextUID - 1;}}

        public long NPoints => NActive + NInactive;

        // Handling this explicitly is more efficient than constantly taking LongCounts
        public long NActive
        { get; private set; }

        public long NInactive
        { get; private set; }

        public long MaxPoints
        { get; private set; }

        public Func<double, double, double, (double, double, double)> VelocityCalc
        { get; protected set; }

        public DomainManager Domain
        { get; protected set; }

        private List<double[]> XHistory;

        private List<double[]> YHistory;

        private List<double[]> PressureHistory;

        private List<double[]> TemperatureHistory;

        private List<double[]> SpecificHumidityHistory;

        private List<uint[]> UIDHistory;

        private List<double> TimeHistory;

        public long MaxStoredPoints
        { get; private set; }

        private bool Debug;

        // Do we calculate the effect of compression on temperature?
        private bool IncludeCompression;

        public string OutputFilename
        { get; private set; }

        public double MassSurplus;

        public PointManager( long? maxPoints, DomainManager domain, string filename, bool debug=false, bool includeCompression=false )
        {
            // UIDs start from 1 (0 reserved for inactive points)
            nextUID = 1;

            VelocityCalc = VCalc;

            // No initial mass surplus
            MassSurplus = 0.0;
            
            // Are we calculating the effect of adiabatic compression?
            IncludeCompression = includeCompression;

            // Where to output data
            OutputFilename = filename;
            
            // Limit on how many points can be managed
            if (maxPoints == null)
            {
                // No limit (danger!)
                MaxPoints = long.MaxValue;
            }
            else
            {
                MaxPoints = (long)maxPoints;
            }

            // Max number of points stored out in any single sample (diagnostic only)
            MaxStoredPoints = 0;

            ActivePoints = [];
            InactivePoints = [];
            NActive = 0;
            NInactive = 0;

            // For output
            TimeHistory = [];
            XHistory = [];
            YHistory = [];
            PressureHistory = [];
            UIDHistory = [];

            TemperatureHistory = [];
            SpecificHumidityHistory = [];

            // Domain manager to use for culling etc
            Domain = domain;

            // Run in debug mode?
            Debug = debug;
            return;

            // Set the velocity calculation
            (double, double, double) VCalc(double x, double y, double pressure) => domain.VelocityFromFixedSpaceArray(x, y, pressure, false);
        }

        private LGPoint AddPoint( double x, double y, double pressure )
        {
            // Create a new point in the list
            // Start by creating an _inactive_ point
            LGPoint point = new LGPoint(VelocityCalc,IncludeCompression);
            InactivePoints.AddLast(point);
            NInactive++;
            // Activate a point (doesn't matter if it's the same one) and return it
            return ActivatePoint(x,y,pressure);
        }

        private LGPoint ActivatePoint( double x, double y, double pressure )
        {
            // Reactivate the first available dormant point and assign it a new UID
            System.Diagnostics.Debug.Assert(InactivePoints.First != null, "InactivePoints.First != null");
            LinkedListNode<LGPoint> node = InactivePoints.First;
            LGPoint point = node.Value;
            InactivePoints.Remove(node);
            ActivePoints.AddLast(node);

            // Give the point its location and a new UID
            // Requesting the UID will automatically increment it
            point.Activate(x,y,pressure,NextUID);
            NInactive--;
            NActive++;
            return point;
        }

        public LGPoint NextPoint( double x, double y, double pressure )
        {
            // Function places a new point, taken from the inactive list if any available.
            // If no points are available, add one if possible; otherwise throw an exception

            // Are there any inactive points available?
            LGPoint point;
            if (NInactive > 0)
            {
                // Reactivate a dormant point
                point = ActivatePoint(x,y,pressure);
            }
            else if (NActive < MaxPoints)
            {
                // Add a new point
                point = AddPoint(x,y,pressure);
            }
            else
            {
                // No more points to return!
                //if (Debug) {Console.WriteLine("!!!");}
                throw new InvalidOperationException("Point maximum exceeded");
            }
            // Set point properties based on local values#
            double temperature = Domain.NearestNeighbor3D(x,y,pressure,Domain.TemperatureXYP);
            double specificHumidity = Domain.NearestNeighbor3D(x,y,pressure,Domain.SpecificHumidityXYP);
            point.SetTemperature(temperature);
            point.SetSpecificHumidity(specificHumidity);

            return point;
        }

        public void DeactivatePoint( LinkedListNode<LGPoint> node )
        {
            // Deactivate point i of those present in ActivePoints
            LGPoint point = node.Value;
            point.Deactivate();
            ActivePoints.Remove(node);
            InactivePoints.AddLast(node);
            NInactive++;
            NActive--;
        }

        public void CreatePointSet( double[] x, double[] y, double[] pressure )
        {
            // Create multiple points
            for (long i=0; i<x.Length; i++)
            {
                NextPoint(x[i],y[i],pressure[i]);
            }
        }

        public void Advance( double dt )
        {
            // Advances all active points one time step
            foreach (LGPoint point in ActivePoints)
            {
                point.Advance(dt);
            }
        }

        public void Cull()
        {
            // Deactivate any points which are outside the domain
            LinkedListNode<LGPoint>? node = ActivePoints.First;
            // The structure below is necessary because we can't get the next node from a deactivated node
            while (node != null)
            {
                LinkedListNode<LGPoint>? nextNode = node.Next;
                LGPoint point = node.Value;
                if (point.X < Domain.XMin || point.X >= Domain.XMax || point.Y < Domain.YMin || point.Y >= Domain.YMax || point.Pressure > Domain.PBase || point.Pressure < Domain.PCeiling )
                {
                    DeactivatePoint(node);
                }
                node = nextNode;
            }
        }

        public bool WriteToFile()
        {
            return WriteToFile(OutputFilename);
        }
        public bool WriteToFile(string fileName)
        {
            bool success = true;

            // Set up output file
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.Create
            };

            // Get the output sizes
            long nPoints = MaxStoredPoints;//Math.Min(MaxStoredPoints,XHistory[0].Length);
            int nTimes = TimeHistory.Count;

            long[] index = new long[nPoints];
            for (long i=0; i<nPoints; i++ )
            {
                index[i] = i;
            }
            
            // Convert the lists into conventional 2D arrays
            double[,] x2D = new double[nTimes, nPoints];
            double[,] y2D = new double[nTimes, nPoints];
            double[,] p2D = new double[nTimes, nPoints];
            double[,] temperature2D = new double[nTimes, nPoints];
            double[,] specificHumidity2D = new double[nTimes, nPoints];
            uint[,] UIDs = new uint[nTimes, nPoints];
            int nCurrent;

            for (int i=0; i<nTimes; i++)
            {
                nCurrent = XHistory[i].Length;
                for (int j=0; j<nPoints; j++)
                {
                    if (j < nCurrent)
                    {
                        x2D[i,j] = XHistory[i][j];
                        y2D[i,j] = YHistory[i][j];
                        p2D[i,j] = PressureHistory[i][j];
                        temperature2D[i,j] = TemperatureHistory[i][j];
                        specificHumidity2D[i,j] = SpecificHumidityHistory[i][j];
                        UIDs[i,j] = UIDHistory[i][j];
                    }
                    else
                    {
                        x2D[i,j] = double.NaN;
                        y2D[i,j] = double.NaN;
                        p2D[i,j] = double.NaN;
                        temperature2D[i,j] = double.NaN;
                        specificHumidity2D[i,j] = double.NaN;
                        UIDs[i,j] = 0;
                    }
                }
            }

            using (DataSet ds = DataSet.Open(dsUri))
            {
                ds.AddAxis("index","-",index);
                ds.AddAxis("time","seconds",TimeHistory.ToArray());
                ds.AddVariable(typeof(double), "x", x2D, ["time","index"]);
                ds.AddVariable(typeof(double), "y", y2D, ["time","index"]);
                ds.AddVariable(typeof(double), "pressure", p2D, ["time","index"]);
                ds.AddVariable(typeof(double), "temperature", temperature2D, ["time","index"]);
                ds.AddVariable(typeof(double), "specific_humidity", specificHumidity2D, ["time","index"]);
                ds.AddVariable(typeof(uint), "UID", UIDs, ["time","index"]);
                ds.Commit();
            }
            
            return success;
        }

        public void ArchiveConditions(double tCurr)
        {
            MaxStoredPoints = Math.Max(MaxStoredPoints,NActive);
            long nPoints = NActive;
            double[] xPoints                = new double[nPoints];
            double[] yPoints                = new double[nPoints];
            double[] pressurePoints         = new double[nPoints];
            double[] temperaturePoints      = new double[nPoints];
            double[] specificHumidityPoints = new double[nPoints];
            uint[] UIDs = new uint[nPoints];
            long i=0;
            foreach (LGPoint point in ActivePoints)
            {
                xPoints[i] = point.X;
                yPoints[i] = point.Y;
                pressurePoints[i] = point.Pressure;
                temperaturePoints[i] = point.Temperature;
                specificHumidityPoints[i] = point.SpecificHumidity;
                UIDs[i] = point.UID;
                i++;
            }
            TimeHistory.Add(tCurr);
            XHistory.Add(xPoints);
            YHistory.Add(yPoints);
            PressureHistory.Add(pressurePoints);
            TemperatureHistory.Add(temperaturePoints);
            SpecificHumidityHistory.Add(specificHumidityPoints);
            UIDHistory.Add(UIDs);
        }
    }
}