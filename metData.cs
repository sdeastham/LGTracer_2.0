using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

namespace LGTracer
{
    public abstract class MetData
    {
        public int NDimensions
        {get; protected set; }

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

        protected abstract void ReadData(DataSet ds);

        protected double ScaleValue;
        protected double OffsetValue;

        protected MetData(string fieldName, int[] xBounds, int[] yBounds, int timesPerFile, double scaleValue=1.0, double offsetValue=0.0)
        {
            // Bounds of domain to be read in
            XBounds = xBounds;
            YBounds = yBounds;
            NX = XBounds[1] - XBounds[0];
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
            TimeIndex = TimesPerFile; // Forces read on next update
        }

        public void Update(DataSet ds)
        {
            if (TimeIndex >= TimesPerFile)
            {
                TimeIndex = 1;
                ReadData(ds);
            }
            else
            {
                TimeIndex++;
            }
        }
    }

    public class MetData2D : MetData
    {
        protected double[][,] FullData;
        protected double[,] NextData
        { get { return FullData[TimeIndex%TimesPerFile]; } }
        protected double[,] PreviousData
        { get { return FullData[(TimeIndex-1)%TimesPerFile]; } }
        public double[,] CurrentData
        { get { return PreviousData; } }

        public MetData2D(string fieldName, int[] xBounds, int[] yBounds, int timesPerFile, double scaleValue=1.0, double offsetValue=0.0) : base(fieldName,xBounds,yBounds,timesPerFile,scaleValue,offsetValue)
        {
            FullData = new double[TimesPerFile+1][,];
            for (int i=0; i<=TimesPerFile; i++)
            {
                FullData[i] = new double[NY,NX];
            }
        }
        protected override void ReadData(DataSet ds)
        {
            // Reads in all time slices from the file
            int lonFirst = XBounds[0];
            int latFirst = YBounds[0];

            // Cycle the last entry of FullData back to the first
            for (int i=0; i<NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    FullData[0][j,i] = FullData[TimesPerFile][j,i];
                }
            }

            float[,,] rawData = ds.GetData<float[,,]>(FieldName);
            for (int t=0; t<TimesPerFile; t++)
            {
                for (int i=0; i<NX; i++)
                {
                    for (int j=0; j<NY; j++)
                    {
                        FullData[t+1][j,i] = (rawData[t,j+YBounds[0],i+XBounds[0]] * ScaleValue) + OffsetValue;
                    }
                }
            }
        }
    }

    public class MetData3D : MetData
    {
        protected double[][,,] FullData;
        protected double[,,] NextData
        { get { return FullData[TimeIndex%TimesPerFile]; } }
        protected double[,,] PreviousData
        { get { return FullData[(TimeIndex-1)%TimesPerFile]; } }
        public double[,,] CurrentData
        { get { return PreviousData; } }

        public int NZ
        { get; protected set; }

        public MetData3D(string fieldName, int[] xBounds, int[] yBounds, int nLevels, int timesPerFile, double scaleValue=1.0, double offsetValue=0.0) : base(fieldName,xBounds,yBounds,timesPerFile,scaleValue,offsetValue)
        {
            NZ = nLevels;
            FullData = new double[TimesPerFile+1][,,];
            for (int i=0; i<=TimesPerFile; i++)
            {
                FullData[i] = new double[NZ,NY,NX];
            }
        }
        protected override void ReadData(DataSet ds)
        {
            // Reads in all time slices from the file
            // Cycle the last entry of FullData back to the first
            for (int i=0; i<NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    for (int k=0; k<NZ; k++)
                    {
                        FullData[0][k,j,i] = FullData[TimesPerFile][k,j,i];
                    }
                }
            }

            float[,,,] rawData = ds.GetData<float[,,,]>(FieldName);
            for (int t=0; t<TimesPerFile; t++)
            {
                for (int i=0; i<NX; i++)
                {
                    for (int j=0; j<NY; j++)
                    {
                        for (int k=0; k<NZ; k++)
                        {
                            FullData[t+1][k,j,i] = (rawData[t,k,j+YBounds[0],i+XBounds[0]] * ScaleValue) + OffsetValue;
                        }
                    }
                }
            }
        }
    }
}