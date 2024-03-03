using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

using System.Diagnostics.CodeAnalysis;

namespace LGTracer
{
    public class MetFile
    {
        // The MetFile class holds all MetData variables
        // read from a single file. It is responsible for
        // keeping track of the dataset handle and updating
        // the variables as time proceeds
        private DataSet DS;
        private NetCDFUri DSUri;
        private DateTime[] TimeVec;
        private string FileTemplate;
        private int SecondOffset; // Number of seconds to offset times which are read in
        public List<IMetData> DataVariables { get; private set; }
        public List<string> DataNames { get; private set; }
        private int[] XBounds, YBounds;
        private double[] XEdge, YEdge;
        private double[] XLim, YLim;
        private int NLevels;
        private int TimeIndex;
        private TimeSpan TimeDelta;

    public MetFile(string fileTemplate, DateTime firstTime, string[] dataFields2D, string[] dataFields3D, double[] xLim, double[] yLim, int secondOffset=0)
        {
            FileTemplate = fileTemplate;
            SecondOffset = secondOffset;
            // Set domain boundaries
            XLim = xLim;
            YLim = yLim;
            // First read will also identify where, in this specific file's data,
            // the X and Y bounds for reading are found
            ReadFile(firstTime,true); // << Seems that this is not establishing DS?
            int nTimes = TimeVec.Length;
            DataVariables = [];
            DataNames = [];
            double scaleValue = 1.0;
            double offsetValue = 0.0;
            TimeIndex = 0;
            // Assume uniform spacing
            TimeDelta = TimeVec[1] - TimeVec[0];
            foreach (string varName in dataFields2D)
            {
                IMetData metVar = new MetData2DFixed(varName, XBounds, YBounds, nTimes, scaleValue, offsetValue);
                // Need to update twice to fill the initial data array
                // and align to the first entry
                metVar.Update(DS);
                metVar.Update(DS);
                DataVariables.Add(metVar);
                DataNames.Add(varName);
            }
            foreach (string varName in dataFields3D)
            {
                IMetData metVar = new MetData3DFixed(varName, XBounds, YBounds, NLevels, nTimes, scaleValue, offsetValue);
                // Need to update twice to fill the initial data array
                // and align to the first entry
                metVar.Update(DS);
                metVar.Update(DS);
                DataVariables.Add(metVar);
                DataNames.Add(varName);
            }
            AdvanceToTime(firstTime);
        }
    
        public void AdvanceToTime(DateTime newTime)
        {
            // Scan through the current times
            // Track whether we need to update - this is important because in theory
            // we could end up with the same time index
            bool updateFiles = false;
            while (TimeVec[TimeIndex + 1] < newTime)
            {
                updateFiles = true;
                TimeIndex++;
                // If the time index now points to the final time in 
                // the vector, then we need to update the underlying
                // date structure
                if (TimeIndex >= (TimeVec.Length - 1))
                {
                    DateTime nextTime = TimeVec[TimeVec.Length - 1] + TimeDelta;
                    ReadFile(nextTime, false);
                    TimeIndex = 0;
                }

                // Update all the variables
                // An interface would be a good idea here...
                foreach (IMetData metVar in DataVariables)
                {
                    metVar.Update(DS);
                }
            }
            // To allow for interpolation. This could get quite expensive
            // so first verify that it's actually going to be necessary
            if (updateFiles)
            {
                double newTimeFraction = IntervalFraction(newTime);
                foreach (IMetData metVar in DataVariables)
                {
                    metVar.SetTimeFraction(newTimeFraction);
                }
            }
        }

        private double IntervalFraction(DateTime targetTime)
        {
            // How far through the current time interval is the proposed time?
            return (targetTime - TimeVec[TimeIndex]).TotalSeconds / TimeDelta.TotalSeconds;
        }
        private string FillTemplate(DateTime targetTime)
        {
            return string.Format(FileTemplate,targetTime.Year,targetTime.Month,targetTime.Day);
        }

        [MemberNotNull(nameof(DS),nameof(DSUri),nameof(TimeVec),nameof(XBounds),nameof(YBounds),nameof(XEdge),nameof(YEdge))]
        private void ReadFile(DateTime targetTime, bool firstRead)
        {
            string fileName = FillTemplate(targetTime);
            DSUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.ReadOnly
            };
            DS = DataSet.Open(DSUri);
            // Parse the time data
            int[] timeInts = DS.GetData<int[]>("time");
            string timeUnits = DS.GetAttr<string>("time","units");
            int nTimes = timeInts.Length;
            TimeVec = ParseFileTimes(timeUnits,timeInts,SecondOffset);
            // Philosophical question: who is handling all these boundaries? Feels like
            // this is the domain manager's job
            if (firstRead || XBounds == null || YBounds == null || XEdge == null || YEdge == null)
            {
                // Set up the domain too
                (XEdge, YEdge, XBounds, YBounds ) = ReadLatLon( DS, XLim, YLim );
                if (DS.Variables.Contains("lev"))
                {
                    float[] levels = DS.GetData<float[]>("lev");
                    NLevels = levels.Length;
                }
                else
                {
                    NLevels = 0;
                }
            }
        }

        private static DateTime[] ParseFileTimes(string units, int[] timeDeltas, int secondOffset=0)
        {
            // Reads a units string (e.g. "minutes since 2023-01-01 00:00:00.0")
            // and a series of integers, returns the corresponding vector of DateTimes
            int nTimes = timeDeltas.Length;
            int secondsMult;
            string[] substrings = units.Split(' ');
            string timeType = substrings[0].ToLower();
            switch (timeType)
            {
                case "seconds":
                    secondsMult = 1;
                    break;
                case "minutes":
                    secondsMult = 60;
                    break;
                case "hours":
                    secondsMult = 3600;
                    break;
                case "days":
                    secondsMult = 3600 * 24;
                    break;
                default:
                    throw new ArgumentException($"Invalid time units {timeType} in string {units}");
            }
            string ymd = substrings[2];
            string hms = substrings[3];
            string ymdhms = $"{ymd} {hms}";
            DateTime refTime = DateTime.Parse(ymdhms);
            DateTime[] timeVec = new DateTime[nTimes];
            for (int i=0; i<nTimes; i++)
            {
                timeVec[i] = refTime.AddSeconds(secondsMult * timeDeltas[i] + secondOffset);
            }
            return timeVec;
        }

        private static (double [], double[], int[], int[] ) ReadLatLon( DataSet ds, double[] lonLims, double[] latLims )
        {
            Func<double,double,double,int> findLower = (targetValue, lowerBound, spacing) => ((int)Math.Floor((targetValue - lowerBound)/spacing));
            double[] lonEdge,latEdge;
            float[] lonMids, latMids;
            int nLon, nLat, lonFirst, latFirst, lonLast, latLast;
            double dLon, dLat, lonBase, latBase;

            // Get the full dimension vectors
            latMids = ds.GetData<float[]>("lat");
            lonMids = ds.GetData<float[]>("lon");

            // Figure out which cells we need to keep in order to get all the data we need
            // Assume a fixed cell spacing for now
            dLon = lonMids[1] - lonMids[0];
            lonBase = lonMids[0] - (dLon/2.0);
            // For latitude, be careful about half-polar grids
            dLat = latMids[3] - latMids[2];
            latBase = latMids[1] - (3.0*dLon/2.0);

            // These indices are for the first and last cell _inclusive_
            latFirst = findLower(latLims[0],latBase,dLat);
            latLast  = findLower(latLims[1],latBase,dLat);
            lonFirst = findLower(lonLims[0],lonBase,dLon);
            lonLast  = findLower(lonLims[1],lonBase,dLon);

            nLon = 1 + (lonLast - lonFirst);
            nLat = 1 + (latLast - latFirst);
            
            // Create lon/lat edge vectors
            lonEdge = new double[nLon+1];
            latEdge = new double[nLat+1];
            lonEdge[0] = lonMids[lonFirst] - (dLon/2.0);
            latEdge[0] = latMids[latFirst] - (dLat/2.0);
            for (int i=0;i<nLon;i++)
            {
                lonEdge[i+1] = lonEdge[0] + (dLon * (i+1));
            }
            for (int i=0;i<nLat;i++)
            {
                latEdge[i+1] = latEdge[0] + (dLat * (i+1));
            }

            // To help in subsetting data, provide the lon/lat limits
            int[] lonSet = [lonFirst,lonLast+1];
            int[] latSet = [latFirst,latLast+1];
            return (lonEdge, latEdge, lonSet, latSet);
        }
        public (double[], double[]) GetXYMesh()
        {
            return (XEdge, YEdge);
        }
    }
}