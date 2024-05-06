using System.Diagnostics;
using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using SerializeNC;

namespace LGTracer;

public interface IMetFile
{
    public IMetData GetMetData(int i);
    public void AdvanceToTime(DateTime newTime);
    public (double[], double[]) GetXYMesh();
    public void Initialize(string[] dataFields2D, string[] dataFields3D, bool timeInterp);
    public void UpdateAllVars(bool readFile);
    public void RegisterStopwatches(Dictionary<string, Stopwatch> stopwatches);
}

public static class MetFileFactory
{
    public static IMetFile CreateMetFile(string fileTemplate, DateTime firstTime, string[] dataFields2D,
        string[] dataFields3D, double[] xLim, double[] yLim, Dictionary<string, Stopwatch> stopwatches,
        int secondOffset=0, bool timeInterp=true, bool useSerial=false)
    {
        IMetFile metFile;
        if (useSerial)
        {
            metFile = new MetFileSerial(fileTemplate, firstTime, xLim, yLim, secondOffset);
        }
        else
        {
            metFile = new MetFileNetCDF(fileTemplate, firstTime, xLim, yLim, secondOffset);
        }
        metFile.RegisterStopwatches(stopwatches);
        metFile.Initialize(dataFields2D, dataFields3D, timeInterp);
        return metFile;
    }
}

public abstract class MetFile : IMetFile
{
    // The MetFile class holds all MetData variables
    // read from a single file. It is responsible for
    // keeping track of the dataset handle and updating
    // the variables as time proceeds
    protected DateTime[] TimeVec;
    protected string FileTemplate;
    protected int SecondOffset; // Number of seconds to offset times which are read in
    public List<IMetData> DataVariables { get; private set; }
    public List<string> DataNames { get; private set; }
    protected int[] XBounds, YBounds;
    protected double[] XEdge, YEdge;
    protected double[] XLim, YLim;
    protected int NLevels;
    protected int TimeIndex;
    protected TimeSpan TimeDelta;
    protected DateTime FirstTime;
    protected double ScaleValue, OffsetValue;
    
    // Options to handle different data conventions
    protected bool SurfaceFirst; // Is the first index the surface or TOA?
    protected bool FlipLats; // Are the file latitudes reversed (i.e. running N->S)?

    protected DateTime LeftBracketTime => TimeVec[TimeIndex - 1];
    protected DateTime RightBracketTime => TimeVec[TimeIndex];

    private Dictionary<string, Stopwatch> Stopwatches;
    public MetFile(string fileTemplate, DateTime firstTime, double[] xLim, double[] yLim, int secondOffset=0)
    {
        FirstTime = firstTime;
        FileTemplate = fileTemplate;
        SecondOffset = secondOffset;
        // Set domain boundaries
        XLim = xLim;
        YLim = yLim;
        DataVariables = [];
        DataNames = [];
        ScaleValue = 1.0;
        OffsetValue = 0.0;
        // Corresponds to the _right bracket_ of the data (if interpolating)
        TimeIndex = 1;
    }

    public void RegisterStopwatches(Dictionary<string, Stopwatch> stopwatches)
    {
        Stopwatches = stopwatches;
    }
    
    public void Initialize(string[] dataFields2D, string[] dataFields3D, bool timeInterp=true)
    {
        // First read will also identify where, in this specific file's data,
        // the X and Y bounds for reading are found
        OpenFile(FirstTime,true); // << Seems that this is not establishing DS?
        int nTimes = TimeVec.Length - 1;
        // Assume uniform spacing
        TimeDelta = TimeVec[1] - TimeVec[0];
        // Use this to indicate to the met data that the latitudes need to be flipped
        int[] yBounds;
        if (FlipLats)
        {
            yBounds = [YBounds[1], YBounds[0]];
        }
        else
        {
            yBounds = YBounds;
        }
        foreach (string varName in dataFields2D)
        {
            IMetData metVar;
            if (timeInterp)
            {
                metVar = new MetData2DLinterp(varName, XBounds, yBounds, nTimes, ScaleValue, OffsetValue);
            }
            else
            {
                metVar = new MetData2DFixed(varName, XBounds, yBounds, nTimes, ScaleValue, OffsetValue);
            }
            DataVariables.Add(metVar);
            DataNames.Add(varName);
        }
        foreach (string varName in dataFields3D)
        {
            IMetData metVar;
            if (timeInterp)
            {
                metVar = new MetData3DLinterp(varName, XBounds, yBounds, NLevels, nTimes, ScaleValue,
                    OffsetValue, SurfaceFirst);
            }
            else
            {
                metVar = new MetData3DFixed(varName, XBounds, yBounds, NLevels, nTimes, ScaleValue,
                    OffsetValue, SurfaceFirst);
            }
            DataVariables.Add(metVar);
            DataNames.Add(varName);
        }
        // Set the time index to zero; this forces an update
        TimeIndex = 0;
        AdvanceToTime(FirstTime);
    }

    public void UpdateAllVars(bool readFile)
    {
        foreach (IMetData metVar in DataVariables)
        {
            try
            {
                UpdateVar(metVar, readFile);
            }
            catch
            {
                Console.WriteLine($"Error updating {metVar.GetName()}.");
                throw;
            }
        }
    }

    protected abstract void UpdateVar(IMetData metVar, bool readFile);

    public IMetData GetMetData(int i)
    {
        return DataVariables[i];
    }
    
    public IMetData GetMetData(string varName)
    {
        int i = GetVarIndex(varName);
        return DataVariables[i];
    }

    public int GetVarIndex(string varName)
    {
        return DataNames.FindIndex(a => a == varName);
    }

    public void AdvanceToTime(DateTime newTime)
    {
        // Scan through the current times
        // Track whether we need to update - this is important because in theory
        // we could end up with the same time index
        Stopwatches["Met advance"].Start();
        bool readFile = TimeIndex == 0;
        while (newTime >= TimeVec[TimeIndex]) // While the current time is AFTER the right bracket..
        {
            TimeIndex++;
            // If the time index now points past the final time in the vector, then we need to update the underlying
            // date structure
            if (TimeIndex >= TimeVec.Length)
            {
                DateTime nextTime = TimeVec[^1] + TimeDelta;
                OpenFile(nextTime, false);
                // Since zero corresponds to the last time from the previous file
                TimeIndex = 1;
                // Variables need to read in new data
                readFile = true;
            }

            // Update all the variables
            // An interface would be a good idea here...
            UpdateAllVars(readFile);
        }
        Stopwatches["Met advance"].Stop();
        // To allow for interpolation. This could get quite expensive, but only
        // variables which do interpolate will do anything
        Stopwatches["Met interpolate"].Start();
        double newTimeFraction = IntervalFraction(newTime);
        foreach (IMetData metVar in DataVariables)
        {
            metVar.SetTimeFraction(newTimeFraction);
        }
        Stopwatches["Met interpolate"].Stop();
    }

    protected double IntervalFraction(DateTime targetTime)
    {
        // How far through the current time interval is the proposed time?
        return (targetTime - TimeVec[TimeIndex-1]).TotalSeconds / TimeDelta.TotalSeconds;
    }
    protected string FillTemplate(DateTime targetTime)
    {
        return string.Format(FileTemplate,targetTime.Year,targetTime.Month,targetTime.Day);
    }

    public abstract void OpenFile(DateTime targetTime, bool firstRead);

    protected static DateTime[] ParseFileTimes(string units, int[] timeDeltas, int secondOffset=0)
    {
        // Reads a units string (e.g. "minutes since 2023-01-01 00:00:00.0")
        // and a series of integers, returns the corresponding vector of DateTimes
        int nTimes = timeDeltas.Length;
        // WARNING: This needs to be a long or you can get integer overflow when multiplying by timeOffsets
        long secondsMult;
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

    protected static (double [], double[], int[], int[] ) ParseLatLon( float[] lonMids, float[] latMids, double[] lonLims, double[] latLims )
    {
        Func<double,double,double,int> findLower = (targetValue, lowerBound, spacing) => ((int)Math.Floor((targetValue - lowerBound)/spacing));
        double[] lonEdge,latEdge;
        int nLon, nLat, lonFirst, latFirst, lonLast, latLast;
        double dLon, dLat, lonBase, latBase;
        // Longitude and latitude midpoints and limits must be monotonically increasing (W->E, S->N)

        // Figure out which cells we need to keep in order to get all the data we need
        // Assume a fixed cell spacing for now
        dLon = lonMids[1] - lonMids[0];
        lonBase = lonMids[0] - (dLon/2.0);
        // Need to deal with the possibility that the domain crosses the "edge". We handle this by ensuring that
        // data fetches are done module NX, at which point the only requirement is that all our boundaries are
        // monotonically increasing (i.e. longitude boundary 0 is less than boundary 1)
        while (lonBase > lonLims[0])
        {
            lonBase -= 360.0;
        }
        // These indices are for the first and last cell _inclusive_
        lonFirst = findLower(lonLims[0],lonBase,dLon);
        lonLast  = findLower(lonLims[1],lonBase,dLon);
        nLon = 1 + (lonLast - lonFirst);

        // For latitude, be careful about half-polar grids
        dLat = latMids[3] - latMids[2];
        latBase = latMids[1] - (3.0 * dLat / 2.0);
        
        latFirst = findLower(latLims[0],latBase,dLat);
        latLast  = findLower(latLims[1],latBase,dLat);
        
        nLat = 1 + (latLast - latFirst);
        
        // Create lon/lat edge vectors
        lonEdge = new double[nLon+1];
        latEdge = new double[nLat+1];
        latEdge[0] = latMids[latFirst] - (dLat/2.0);
        
        // Internally, we treat [-180,180) as the valid longitude range
        double lonEastBoundary = lonMids[lonFirst%lonMids.Length] - (dLon/2.0);
        while (lonEastBoundary < -180.0) { lonEastBoundary += 360.0; }
        while (lonEastBoundary >= 180.0) { lonEastBoundary -= 360.0; }
        lonEdge[0] = lonEastBoundary;
        
        // This ensures monotonicity
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

public class MetFileNetCDF : MetFile
{
    private DataSet DS;
    private NetCDFUri DSUri;
    private static readonly string[] pressureUnits = new[] {"mb","millibar","millibars","pa","pascal","pascals","hpa","hectopascal","hectopascals"};

    public MetFileNetCDF(string fileTemplate, DateTime firstTime, double[] xLim, double[] yLim, int secondOffset = 0) :
        base(fileTemplate, firstTime, xLim, yLim, secondOffset)
    {
        // No change
    }

    protected override void UpdateVar(IMetData metVar, bool readFile)
    {
        metVar.Update(DS,TimeIndex,readFile);
    }
    
    public override void OpenFile(DateTime targetTime, bool firstRead)
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
        TimeVec = new DateTime[nTimes + 1];
        DateTime[] timeVec = ParseFileTimes(timeUnits,timeInts,SecondOffset);
        for (int i = 0; i < nTimes; i++)
        {
            TimeVec[i + 1] = timeVec[i];
        }
        TimeVec[0] = TimeVec[1] - (timeVec[1] - timeVec[0]);
        // Philosophical question: who is handling all these boundaries? Feels like
        // this is the domain manager's job
        if (firstRead || XBounds == null || YBounds == null || XEdge == null || YEdge == null)
        {
            // Set up the domain too
            string? latVar, lonVar, levVar;
            lonVar = testVarPresence(["lon", "longitude"], DS);
            latVar = testVarPresence(["lat", "latitude"], DS);
            levVar = testVarPresence(["lev", "level"], DS);
            // If we are missing longitude or latitude, cannot proceed
            if (lonVar == null || latVar == null)
            {
                throw new ArgumentException($"No valid longitude or latitude variable found in file.");
            }
            // Missing level might be OK
            bool levFound = (levVar != null);
            float[] latMids = DS.GetData<float[]>(latVar);
            float[] lonMids = DS.GetData<float[]>(lonVar);
            
            // Some files have latitudes arranged N -> S; these will need to be flipped
            FlipLats = (latMids[1] < latMids[0]);
            if (FlipLats)
            {
                Array.Reverse(latMids);
            }
            
            (XEdge, YEdge, XBounds, YBounds ) = ParseLatLon( lonMids, latMids, XLim, YLim );
            if (levFound)
            {
                float[] levels;
                if (DS.Variables[levVar].TypeOfData == typeof(int))
                {
                    int[] intLevels = DS.GetData<int[]>(levVar);
                    levels = intLevels.Select(d => (float)d).ToArray();
                }
                else
                {
                    levels = DS.GetData<float[]>(levVar);
                }
                NLevels = levels.Length;
                string levUnits = Convert.ToString(DS.Variables[levVar].Metadata.ToDictionary()["units"]) ?? "none";
                if (pressureUnits.Contains(levUnits.ToLower()))
                {
                    // Pressure decreases away from the surface, so if the first value is greater than the second,
                    // we are closer to the surface with the first value
                    SurfaceFirst = levels[0] > levels[1];
                }
                else
                {
                    // Otherwise, no information - can only hope that lower numbers == closer to surface
                    SurfaceFirst = levels[1] > levels[0];
                }
            }
            else
            {
                NLevels = 0;
                SurfaceFirst = true; // No harm
            }
        }
    }

    private static string? testVarPresence(string[] varOpts, DataSet dataSet)
    {
        string? result = null;
        foreach (string testVar in varOpts)
        {
            if (dataSet.Variables.Contains(testVar))
            {
                result = testVar;
                break;
            }
        }
        return result;
    }
}

public class MetFileSerial : MetFile
{
    private string CurrentFileTemplate;
    public MetFileSerial(string fileTemplate, DateTime firstTime, double[] xLim, double[] yLim, int secondOffset = 0) :
        base(fileTemplate, firstTime, xLim, yLim, secondOffset)
    {
        // Nothing additional to do
    }

    protected override void UpdateVar(IMetData metVar, bool readFile)
    {
        metVar.Update(string.Format(CurrentFileTemplate,metVar.GetName()),TimeIndex,readFile);
    }

    private void FillCurrentTemplate(DateTime targetTime)
    {
        // Updates the file template so that the full path to the file for a given variable can be easily generated
        CurrentFileTemplate = string.Format(FileTemplate, targetTime.Year, targetTime.Month, targetTime.Day, "{0}");
    }

    private string VariableFilePath(string varName)
    {
        return string.Format(CurrentFileTemplate, varName);
    }
    
    public override void OpenFile(DateTime targetTime, bool firstRead)
    {
        FillCurrentTemplate(targetTime);
        string fileName = VariableFilePath("DIMENSIONS");
        (DateTime[]? timeVec, int? nLevels, float[]? latMids, float[]? lonMids, bool[] dimsFound) =
            NetcdfSerializer.DeserializeDimensions(fileName);
        int nTimes = timeVec?.Length ?? 1;
        TimeVec = new DateTime[nTimes + 1];
        for (int i = 0; i < nTimes; i++)
        {
            TimeVec[i+1] = timeVec[i] + TimeSpan.FromSeconds(SecondOffset);
        }
        TimeVec[0] = TimeVec[1] - (timeVec[1] - timeVec[0]);
        if (firstRead || XBounds == null || YBounds == null || XEdge == null || YEdge == null)
        {
            // Set up the domain too
            (XEdge, YEdge, XBounds, YBounds ) = ParseLatLon( lonMids, latMids, XLim, YLim );
            NLevels = nLevels ?? 1;
        }
    }
}