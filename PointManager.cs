using System.ComponentModel;
using System.Globalization;
using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;
using Parquet;
using Parquet.Schema;

namespace LGTracer;

public abstract class PointManager
{
    public LinkedList<IAdvected> ActivePoints { get; private set; }

    protected LinkedList<IAdvected> InactivePoints { get; set; }

    private uint nextUID { get; set; }

    private uint NextUID
    {
        get
        {
            nextUID += 1;
            return nextUID - 1;
        }
    }

    public long NPoints => NActive + NInactive;

    // Handling this explicitly is more efficient than constantly taking LongCounts
    public long NActive { get; private set; }

    public long NInactive { get; private set; }

    public long MaxPoints { get; private set; }

    public Func<double, double, double, (double, double, double)> VelocityCalc { get; protected set; }

    public DomainManager Domain { get; protected set; }

    private List<double[]> XHistory;

    private List<double[]> YHistory;

    private List<double[]> PressureHistory;

    private List<double[]> AgeHistory;

    private List<uint[]> UIDHistory;

    private List<double> TimeHistory;

    private Dictionary<string, List<double[]>> PropertyHistory;

    private List<string> PropertyNames => PropertyHistory.Keys.ToList();
    
    public long MaxStoredPoints
    { get; private set; }

    protected bool VerboseOutput;

    // Do we calculate the effect of compression on temperature?
    protected bool IncludeCompression;

    public bool WriteOutput
    { get; private set; }
    public string OutputFilename
    { get; private set; }
    public string OutputDirectory
    { get; private set; }
    
    public bool WriteTrajectories
    { get; private set; }
    public string TrajectoryFilename
    { get; private set; }

    private DateTime StorageStartTime;

    private DateTime CurrentTime;

    public PointManager(DomainManager domain, LGOptions configOptions, LGOptionsPoints configSubOptions)
    {
        // UIDs start from 1 (0 reserved for inactive points)
        nextUID = 1;

        VelocityCalc = VCalc;

        long? maxPoints = configSubOptions.Max;
            
        // Are we calculating the effect of adiabatic compression?
        IncludeCompression = configSubOptions.AdiabaticCompression;

        // Where to output data
        WriteOutput = configSubOptions.WritePeriodic; // Output data on all points to netCDF files periodically
        OutputDirectory = configOptions.InputOutput.OutputDirectory;
        OutputFilename = configSubOptions.OutputFilename; // Must contain {date} to prevent overwrites
        WriteTrajectories = configSubOptions.WriteTrajectories; // Each point writes trajectory to unique parquet file
        TrajectoryFilename = configSubOptions.TrajectoryFilename; // Must contain {date} and {uid} to prevent overwrites
        
        // The start of the current storage period
        StorageStartTime = configOptions.Timing.StartDate;
        CurrentTime = configOptions.Timing.StartDate;
            
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

        ActivePoints = [];
        InactivePoints = [];
        NActive = 0;
        NInactive = 0;
        InitializeHistory(configSubOptions.OutputVariables);

        // Domain manager to use for culling etc
        Domain = domain;

        // Run in debug mode?
        VerboseOutput = configOptions.Verbose;
        return;

        // Set the velocity calculation
        (double, double, double) VCalc(double x, double y, double pressure) => domain.VelocityFromFixedSpaceArray(x, y, pressure, false);
    }

    // Must be overridden!
    public abstract void Seed(double dt);

    // Use this to abstract the underlying type of the points and allow inheriting classes to change them
    protected virtual IAdvected CreatePoint()
    {
        return new LGPoint(VelocityCalc);
    }

    private void InitializeHistory(string[]? propertyNames)
    {
        PropertyHistory = [];
        if (propertyNames != null)
        {
            foreach (string property in propertyNames)
            {
                List<double[]> newProperty = [];
                PropertyHistory.Add(property,newProperty);
            }
        }
        InitializeHistory();
    }
    private void InitializeHistory()
    {
        // Max number of points stored out in any single sample (diagnostic only)
        MaxStoredPoints = 0;

        // For output
        TimeHistory = [];
        XHistory = [];
        YHistory = [];
        PressureHistory = [];
        AgeHistory = [];
        UIDHistory = [];

        List<string> propertyNames = PropertyHistory.Keys.ToList();
        PropertyHistory.Clear();
        foreach (string property in propertyNames)
        {
            List<double[]> newProperty = [];
            PropertyHistory[property] = newProperty;
        }
    }
    
    private IAdvected AddPoint( double x, double y, double pressure )
    {
        // Create a new point in the list
        // Start by creating an _inactive_ point
        IAdvected point = CreatePoint();
        InactivePoints.AddLast(point);
        NInactive++;
        if (WriteTrajectories)
        {
            // Tell the point what data it will need to archive
            point.SetupHistory(PropertyNames);
        }
        // Activate a point (doesn't matter if it's the same one) and return it
        return ActivatePoint(x,y,pressure);
    }

    private IAdvected ActivatePoint( double x, double y, double pressure )
    {
        // Reactivate the first available dormant point and assign it a new UID
        System.Diagnostics.Debug.Assert(InactivePoints.First != null, "InactivePoints.First != null");
        LinkedListNode<IAdvected> node = InactivePoints.First;
        IAdvected point = node.Value;
        InactivePoints.Remove(node);
        ActivePoints.AddLast(node);

        // Give the point its location and a new UID
        // Requesting the UID will automatically increment it
        point.Activate(x,y,pressure,NextUID,CurrentTime);
        NInactive--;
        NActive++;
        return point;
    }

    public virtual IAdvected NextPoint( double x, double y, double pressure )
    {
        // Function places a new point, taken from the inactive list if any available.
        // If no points are available, add one if possible; otherwise throw an exception

        // Are there any inactive points available?
        IAdvected point;
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
        return point;
    }

    public void DeactivatePoint( LinkedListNode<IAdvected> node )
    {
        // Deactivate point i of those present in ActivePoints
        IAdvected point = node.Value;
        
        // Write history to file if requested and the point made it past its first time step
        if (WriteTrajectories)
        {
            Dictionary<string, List<double>> pointHistory = point.GetHistory();
            if (pointHistory.Count > 0 && pointHistory["age"].Count > 1)
            {
                // Add the location at "death" 
                if (Math.Abs(pointHistory["age"].Last() - point.GetAge()) > 1.0e-3)
                {
                    point.ArchiveConditions();
                }
                // Add this point's data to the running record
                //AddTrajectoryToBuffer(point.GetUID(),pointHistory);
            }
        }
        point.Deactivate();
        ActivePoints.Remove(node);
        InactivePoints.AddLast(node);
        NInactive++;
        NActive--;
    }

    public void DeactivateAllPoints()
    {
        // Copied from Cull
        LinkedListNode<IAdvected>? node = ActivePoints.First;
        // The structure below is necessary because we can't get the next node from a deactivated node
        while (node != null)
        {
            LinkedListNode<IAdvected>? nextNode = node.Next;
            DeactivatePoint(node);
            node = nextNode;
        }
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
        foreach (IAdvected point in ActivePoints)
        {
            point.Advance(dt, Domain);
        }
        CurrentTime += TimeSpan.FromSeconds(dt);
    }

    public virtual void Cull()
    {
        // Deactivate any points which are outside the domain
        LinkedListNode<IAdvected>? node = ActivePoints.First;
        // The structure below is necessary because we can't get the next node from a deactivated node
        while (node != null)
        {
            LinkedListNode<IAdvected>? nextNode = node.Next;
            IAdvected point = node.Value;
            (double x, double y, double p) = point.GetLocation();
            // Multiple possible reasons for invalidity
            // CheckValid is designed to see if the point is invalid for physics reasons (e.g. it is too diffuse to
            // track) whereas the domain checks are universal
            if (!(point.CheckValid() && Domain.InDomainXYP(x,y,p)))
            {
                DeactivatePoint(node);
            }
            node = nextNode;
        }
    }

    public virtual string WriteToFile(DateTime currentTime, bool reset=true)
    {
        // Start by assuming success
        bool success = true;
        
        // Parse the file name
        string fileNameShort = OutputFilename.Replace("{date}", StorageStartTime.ToString("yyyyMMddTHHmmss"));
        string fileName = Path.Join(OutputDirectory, fileNameShort);
        
        // Advance the time to be used for storage next time around
        StorageStartTime = currentTime;

        // Set up output file
        var dsUri = new NetCDFUri
        {
            FileName = fileName,
            OpenMode = ResourceOpenMode.Create
        };

        // Get the output sizes
        long nPoints = MaxStoredPoints;
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
        double[,] age2D = new double[nTimes, nPoints];
        uint[,] UIDs = new uint[nTimes, nPoints];
        int nProperties = PropertyNames.Count();
        double[,,] properties2D = new double[nProperties, nTimes, nPoints];
            
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
                    age2D[i,j] = AgeHistory[i][j];
                    UIDs[i,j] = UIDHistory[i][j];
                    //for (int k = 0; k < nProperties; k++)
                    //{
                    //    properties2D[k,i,j] = PropertyHistory[k][i][j];
                    //}
                    int k = 0;
                    foreach (string property in PropertyNames)
                    {
                        properties2D[k,i,j] = PropertyHistory[property][i][j];
                        k++;
                    }
                }
                else
                {
                    x2D[i,j] = double.NaN;
                    y2D[i,j] = double.NaN;
                    p2D[i,j] = double.NaN;
                    age2D[i,j] = double.NaN;
                    for (int k = 0; k < nProperties; k++)
                    {
                        properties2D[k,i,j] = double.NaN;
                    }
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
            ds.AddVariable(typeof(double), "age", age2D, ["time","index"]);
            ds.AddVariable(typeof(uint), "UID", UIDs, ["time","index"]);
            for (int k = 0; k < nProperties; k++)
            {
                double[,] property2D = new double[nTimes, nPoints];
                for (int i = 0; i < nTimes; i++)
                {
                    for (int j = 0; j < nPoints; j++)
                    {
                        property2D[i, j] = properties2D[k, i, j];
                    }
                }
                ds.AddVariable(typeof(double), PropertyNames[k], property2D, ["time","index"]);
            }
                
            ds.Commit();
        }

        if (reset)
        {
            MaxStoredPoints = 0;
            InitializeHistory();
        }

        if (WriteTrajectories)
        {
            WriteTrajectoriesToFile();
        }
        
        return fileName;
    }

    public void ArchiveConditions(double tCurr)
    {
        MaxStoredPoints = Math.Max(MaxStoredPoints,NActive);
        long nPoints = NActive;
        double[] ages           = new double[nPoints];
        double[] xPoints        = new double[nPoints];
        double[] yPoints        = new double[nPoints];
        double[] pressurePoints = new double[nPoints];
        uint[] UIDs = new uint[nPoints];
        List<double[]> properties = [];
        int nProperties = PropertyNames.Count();
        for (int k = 0; k < nProperties; k++)
        {
            double[] vec = new double[nPoints];
            properties.Add(vec);
        }
        long i=0;
        foreach (IAdvected point in ActivePoints)
        {
            ages[i] = point.GetAge();
            (xPoints[i], yPoints[i], pressurePoints[i]) = point.GetLocation();
            UIDs[i] = point.GetUID();
            // Get any remaining properties from the point
            // Need to allow for different point classes
            // It would be better here to register the get methods at manager initialization, but
            // that may not be straightforward
            int k = 0;
            foreach (string property in PropertyNames)
            {
                properties[k][i] = point.GetProperty(property);
                k++;
            }

            if (WriteTrajectories)
            {
                point.ArchiveConditions();
            }
            i++;
        }
        TimeHistory.Add(tCurr);
        XHistory.Add(xPoints);
        YHistory.Add(yPoints);
        PressureHistory.Add(pressurePoints);
        UIDHistory.Add(UIDs);
        AgeHistory.Add(ages);
        int iProperty = 0;
        foreach (string property in PropertyNames)
        {
            PropertyHistory[property].Add(properties[iProperty]);
            iProperty++;
        }
    }
    
    private ParquetSchema? Schema = null;
    private async void WriteTrajectoriesToFile()
    {
        //Console.WriteLine($"Writing point information to {Filename}");
        // Start lazy..
        //ParquetSerializer.SerializeAsync(History, Filename);
        /*
        if (Schema == null)
        {
            // Only define this once
            Field[] dataFields = new Field[History.Count];
            int iField = 0;
            foreach (string property in History.Keys)
            {
                dataFields[iField] = new DataField<double>(property); 
                iField++;
            }
            Schema = new ParquetSchema(dataFields);
        }
        using (Stream fs = System.IO.File.OpenWrite(Filename))
        {
            using (ParquetWriter writer = await ParquetWriter.CreateAsync(Schema,fs))
            {
                using (ParquetRowGroupWriter groupWriter = writer.CreateRowGroup())
                {
                    int iColumn = 0;
                    foreach (string property in History.Keys)
                    {
                        var column = new DataColumn(Schema.DataFields[iColumn], History[property].ToArray());
                        await groupWriter.WriteColumnAsync(column);
                        iColumn++;
                    }
                }
            }
        }
        // Clear the trajectory history
        */
    }
}