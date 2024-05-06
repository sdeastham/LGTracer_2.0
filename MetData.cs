using System.Numerics;
using System.Xml;
using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;
using SerializeNC;

namespace LGTracer;

public interface IMetData
{
    void Update<T>(T dataSource, int timeIndex, bool readFile);
    void SetTimeFraction(double timeFraction);
    string GetName();
}
public abstract class MetData<T> : IMetData
{
    public int NDimensions
    {get; protected set; }

    protected double TimeFraction;
    protected T[] FullData { get; set; }

    public abstract T CurrentData { get; protected set; }
    public T NextData => FullData[TimeIndex];
    public T PreviousData => FullData[TimeIndex-1];

    // TimeIndex is the RIGHT BRACKET
    protected int TimeIndex; // What was the last time index in FullData used to set NextData?
    protected int TimesPerFile; // How many times in each file?

    // Three other variables should be present but dimensionality-dependent:
    // --> PreviousData
    // --> NextData
    // --> CurrentData (which just accesses PreviousData for now)

    // The field name read from the file
    public string FieldName
    { get; protected set; }

    // The X and Y bounds to use when reading from the file
    protected int[] XBounds;
    protected int[] YBounds;
    public int NX
    { get; protected set; }
    public int NY
    { get; protected set; }

    protected abstract void ShuffleLastToFirst();
    protected abstract void ReadData(DataSet ds);
    protected abstract void ReadData(string fileTemplate);

    protected double ScaleValue;
    protected double OffsetValue;

    protected bool SerializedData;
    protected int Rank;

    private bool Initialized;

    protected bool FlipY;

    protected bool ConvertShort; // Does the data need to be converted from an array of short[...] to float[...]?
    protected float ShortScaling = 1.0f;
    protected float ShortOffset = 0.0f;

    protected MetData(string fieldName, int[] xBounds, int[] yBounds, int timesPerFile, double scaleValue=1.0, double offsetValue=0.0, bool serializedData=false)
    {
        // Bounds of domain to be read in
        XBounds = xBounds;
        NX = XBounds[1] - XBounds[0];

        FlipY = (yBounds[1] < yBounds[0]);
        YBounds = FlipY ? [yBounds[1],yBounds[0]] : yBounds;
        NY = YBounds[1] - YBounds[0];

        // The field to be read from the source file
        FieldName = fieldName;
        // Scaling to be applied to the variable
        ScaleValue = scaleValue;
        // Offset to be applied after scaling is applied
        OffsetValue = offsetValue;

        NDimensions = 2;
        TimesPerFile = timesPerFile;
        // At any given time (once initialized), an instance of this class will contain a jagged array
        // FullData which is size [N+1][,]. Entry 0 will be the *last entry from the previous read*;
        // entry 1 to N will be the entries from the current file read. Since TimeIndex refers to the
        // entry in FullData relevant to NextData, it gets reset to 1 each time a new file is read in.
        // If TimeIndex hits N+1, a new read operation is required.
        TimeIndex = 0; // Will be overwritten

        // Indicate that the first update needs to fill the data arrays
        Initialized = false;
        
        // Will we use netCDF data or serialized?
        SerializedData = serializedData;
        Rank = 0;
    }

    public string GetName()
    {
        return FieldName;
    }

    public void Update<T2>(T2 dataSource, int timeIndex, bool readFile)
    {
        // Structure changed to be top-down
        TimeIndex = timeIndex;
        if (Initialized)
        {
            ShuffleLastToFirst();
        }
        if (readFile || !Initialized)
        {
            switch (dataSource)
            {
                case DataSet:
                    if (!Initialized)
                    {
                        // ERA5 data is held as int16, i.e. short, rather than 32-bit floats, which are then converted
                        // SDSLite does not have elegant methods to handle this
                        Variable dsVariable = ((DataSet)(object)dataSource!).Variables[FieldName];
                        ConvertShort = dsVariable.TypeOfData == typeof(short);
                        if (ConvertShort)
                        {
                            var propertyDict = dsVariable.Metadata.AsDictionary();
                            ShortScaling = Convert.ToSingle(propertyDict["scale_factor"]);
                            ShortOffset = Convert.ToSingle(propertyDict["add_offset"]);
                        }
                    }
                    ReadData((DataSet)(object)dataSource!);
                    break;
                case string:
                    ConvertShort = false;
                    ReadData((string)(object)dataSource);
                    break;
                default:
                    throw new ArgumentException("Unrecognized data type for met update.");
            }
            //Console.WriteLine($"VARIABLE {GetName()} NEW READ-IN COMPLETE");
        }
        Initialized = true;
    }

    public virtual void SetTimeFraction(double timeFraction)
    {
        TimeFraction = timeFraction;
    }
}

public abstract class MetData2D : MetData<double[,]>
{
    public MetData2D(string fieldName, int[] xBounds, int[] yBounds, int timesPerFile, double scaleValue=1.0, double offsetValue=0.0) : base(fieldName,xBounds,yBounds,timesPerFile,scaleValue,offsetValue)
    {
        Rank = 2;
        FullData = new double[TimesPerFile+1][,];
        for (int i=0; i<(TimesPerFile+1); i++)
        {
            FullData[i] = new double[NY,NX];
        }
    }

    protected override void ShuffleLastToFirst()
    {
        // Cycle the last entry of FullData back to the first
        for (int i=0; i<NX; i++)
        {
            for (int j=0; j<NY; j++)
            {
                FullData[0][j,i] = FullData[TimesPerFile][j,i];
            }
        }
    }

    private void ScaleShiftData(float[,,] rawData)
    {
        FlipRawData(rawData);
        int nXRaw = rawData.GetLength(2);
        for (int t = 0; t < TimesPerFile; t++)
        {
            for (int i = 0; i < NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    FullData[t+1][j,i] = (rawData[t,j+YBounds[0],(i+XBounds[0])%nXRaw] * ScaleValue) + OffsetValue;
                }
            }
        }
    }
    
    protected override void ReadData(DataSet ds)
    {
        // ds.Variables[FieldName].GetData() WILL return an array, regardless of type
        float[,,] rawData;
        if (ConvertShort)
        {
            short[,,] tempShort = (short[,,])ds.Variables[FieldName].GetData();
            int nTRaw = tempShort.GetLength(0);
            int nYRaw = tempShort.GetLength(1);
            int nXRaw = tempShort.GetLength(2);
            rawData = new float[nTRaw, nYRaw, nXRaw];
            for (int t = 0; t < nTRaw; t++)
            {
                for (int j = 0; j < nYRaw; j++)
                {
                    for (int i = 0; i < nXRaw; i++)
                    {
                        rawData[t, j, i] = ShortScaling * tempShort[t, j, i] + ShortOffset;
                    }                
                }
            }
        }
        else
        {
            rawData = ds.GetData<float[,,]>(FieldName);
        }
        ScaleShiftData(rawData);
    }

    protected override void ReadData(string fileTemplate)
    {
        (_, float[,,] rawData) = NetcdfSerializer.Deserialize2D(string.Format(fileTemplate, FieldName));
        ScaleShiftData(rawData);
    }

    protected void FlipRawData(float[,,] rawData)
    {
        // We need to be able to handle data where the latitudes are reversed relative to what LGTracer uses
        // (South -> North)
        if (!FlipY) { return; }
        int nYRaw = rawData.GetLength(1);
        int nXRaw = rawData.GetLength(2);
        for (int t = 0; t < rawData.GetLength(0); t++)
        {
            for (int j = 0; j < nYRaw; j++)
            {
                int jFlip = nYRaw - (j + 1);
                for (int i = 0; i < nXRaw; i++)
                {
                    // Swaps the variables
                    (rawData[t, j, i], rawData[t, jFlip, i]) = (rawData[t, jFlip, i], rawData[t, j, i]);
                }
            }
        }
    }
}
    
public abstract class MetData3D : MetData<double[,,]>
{
    public int NZ
    { get; protected set; }

    private bool SurfaceFirst;

    public MetData3D(string fieldName, int[] xBounds, int[] yBounds, int nLevels, int timesPerFile, double scaleValue=1.0, double offsetValue=0.0,
        bool surfaceFirst=true) : base(fieldName,xBounds,yBounds,timesPerFile,scaleValue,offsetValue)
    {
        NZ = nLevels;
        SurfaceFirst = surfaceFirst;
        FullData = new double[TimesPerFile+1][,,];
        for (int i=0; i<(TimesPerFile+1); i++)
        {
            FullData[i] = new double[NZ,NY,NX];
        }
    }
    protected override void ShuffleLastToFirst()
    {
        // Cycle the last entry of FullData back to the first
        for (int i=0; i<NX; i++)
        {
            for (int j=0; j<NY; j++)
            {
                for (int k = 0; k < NZ; k++)
                {
                    FullData[0][k, j, i] = FullData[TimesPerFile][k, j, i];
                }
            }
        }
    }

    private void ScaleShiftData(float[,,,] rawData)
    {
        FlipRawData(rawData);
        for (int t = 0; t < TimesPerFile; t++)
        {
            for (int i = 0; i < NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    for (int k = 0; k < NZ; k++)
                    {
                        FullData[t + 1][k, j, i] =
                            (rawData[t, k, j + YBounds[0], (i + XBounds[0])%NX] * ScaleValue) + OffsetValue;
                    }
                }
            }
        }
    }
    
    protected void FlipRawData(float[,,,] rawData)
    {
        // We need to be able to handle data where the latitude and/or altitude are reversed relative to what
        // LGTracer uses (surface = 0, South -> North)
        if (!(FlipY || SurfaceFirst)) { return; }
        int nZRaw = rawData.GetLength(1);
        int nYRaw = rawData.GetLength(2);
        int nXRaw = rawData.GetLength(3);
        for (int t = 0; t < rawData.GetLength(0); t++)
        {
            for (int k = 9; k < nZRaw; k++)
            {
                int kFlip;
                if (SurfaceFirst)
                {
                    kFlip = k;
                }
                else
                {
                    kFlip = nZRaw - (k + 1);
                }
                for (int j = 0; j < nYRaw; j++)
                {
                    int jFlip;
                    if (FlipY)
                    {
                        jFlip = nYRaw - (j + 1);
                    }
                    else
                    {
                        jFlip = j;
                    }
                    for (int i = 0; i < nXRaw; i++)
                    {
                        // Swaps the variables
                        (rawData[t, k, j, i], rawData[t, kFlip, jFlip, i]) =
                            (rawData[t, kFlip, jFlip, i], rawData[t, k, j, i]);
                    }
                }
            }
        }
    }
    
    protected override void ReadData(DataSet ds)
    {
        float[,,,] rawData;
        if (ConvertShort)
        {
            short[,,,] tempShort = (short[,,,])ds.Variables[FieldName].GetData();
            int nTRaw = tempShort.GetLength(0);
            int nZRaw = tempShort.GetLength(1);
            int nYRaw = tempShort.GetLength(2);
            int nXRaw = tempShort.GetLength(3);
            rawData = new float[nTRaw, nZRaw, nYRaw, nXRaw];
            for (int t = 0; t < nTRaw; t++)
            {
                for (int i = 0; i < nXRaw; i++)
                {
                    for (int j = 0; j < nYRaw; j++)
                    {
                        for (int k = 0; k < nZRaw; k++)
                        {
                            rawData[t, k, j, i] = ShortScaling * tempShort[t, k, j, i] + ShortOffset;
                        }
                    }
                }
            }
        }
        else
        {
            rawData = ds.GetData<float[,,,]>(FieldName);
        }
        ScaleShiftData(rawData);
    }

    protected override void ReadData(string fileTemplate)
    {
        ScaleShiftData(NetcdfSerializer.Deserialize3D(string.Format(fileTemplate, FieldName)).Item2);
    }
}
    
// Fixed implementations
public class MetData3DFixed : MetData3D
{
    public override double[,,] CurrentData
    {
        get => PreviousData;
        protected set => throw new InvalidOperationException("CurrentData field should not be set directly for a fixed met data object");
    }

    public MetData3DFixed(string fieldName, int[] xBounds, int[] yBounds, int nLevels, int timesPerFile,
        double scaleValue = 1.0, double offsetValue = 0.0, bool surfaceFirst = true) : base(fieldName, xBounds, yBounds, nLevels, timesPerFile,
        scaleValue, offsetValue, surfaceFirst) {}
}

public class MetData2DFixed : MetData2D
{
    public override double[,] CurrentData
    {
        get => PreviousData;
        protected set => throw new InvalidOperationException("CurrentData field should not be set directly for a fixed met data object");
    }

    public MetData2DFixed(string fieldName, int[] xBounds, int[] yBounds, int timesPerFile,
        double scaleValue = 1.0, double offsetValue = 0.0) : base(fieldName, xBounds, yBounds, timesPerFile,
        scaleValue, offsetValue) {}
}
    
// The Linterp approach interpolates the field when the time is updated, and only then
public class MetData3DLinterp : MetData3D
{
    public override double[,,] CurrentData { get; protected set; }

    public override void SetTimeFraction(double timeFraction)
    {
        TimeFraction = timeFraction;
        double previousFraction = 1.0 - timeFraction;
        double nextFraction = timeFraction;
        CurrentData = new double[NZ, NY, NX];
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                for (int k = 0; k < NZ; k++)
                {
                    CurrentData[k, j, i] = (previousFraction * PreviousData[k, j, i]) +
                                           (nextFraction * NextData[k, j, i]);
                }
            }
        }
    }

    public MetData3DLinterp(string fieldName, int[] xBounds, int[] yBounds, int nLevels, int timesPerFile,
        double scaleValue = 1.0, double offsetValue = 0.0, bool surfaceFirst = true) : base(fieldName, xBounds, yBounds, nLevels, timesPerFile,
        scaleValue, offsetValue, surfaceFirst) {}
}
    
public class MetData2DLinterp : MetData2D
{
    public override double[,] CurrentData { get; protected set; }
        
    public override void SetTimeFraction(double timeFraction)
    {
        TimeFraction = timeFraction;
        double previousFraction = 1.0 - timeFraction;
        double nextFraction = timeFraction;
        CurrentData = new double[NY, NX];
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                CurrentData[j, i] = (previousFraction * PreviousData[j, i]) + 
                                    (nextFraction * NextData[j, i]);
            }
        }
    }
        
    public MetData2DLinterp(string fieldName, int[] xBounds, int[] yBounds, int timesPerFile,
        double scaleValue = 1.0, double offsetValue = 0.0) : base(fieldName, xBounds, yBounds, timesPerFile,
        scaleValue, offsetValue) {}
}