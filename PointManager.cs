using System.ComponentModel;
using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

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

    private List<List<double[]>> PropertyHistory;

    private List<string> PropertyNames;

    public long MaxStoredPoints
    { get; private set; }

    private bool VerboseOutput;

    // Do we calculate the effect of compression on temperature?
    protected bool IncludeCompression;

    public string OutputFilename
    { get; private set; }

    public PointManager( long? maxPoints, DomainManager domain, string filename, bool verboseOutput=false, bool includeCompression=false, string[]? propertyNames=null )
    {
        // UIDs start from 1 (0 reserved for inactive points)
        nextUID = 1;

        VelocityCalc = VCalc;
            
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
        AgeHistory = [];
        UIDHistory = [];

        PropertyNames = [];
        PropertyHistory = [];
        if (propertyNames != null)
        {
            foreach (string property in propertyNames)
            {
                PropertyNames.Add(property);
                List<double[]> newProperty = [];
                PropertyHistory.Add(newProperty);
            }
        }

        // Domain manager to use for culling etc
        Domain = domain;

        // Run in debug mode?
        VerboseOutput = verboseOutput;
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
        
    private IAdvected AddPoint( double x, double y, double pressure )
    {
        // Create a new point in the list
        // Start by creating an _inactive_ point
        IAdvected point = CreatePoint();
        InactivePoints.AddLast(point);
        NInactive++;
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
        point.Activate(x,y,pressure,NextUID);
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
        foreach (IAdvected point in ActivePoints)
        {
            LGPoint lgp = (LGPoint)point;
            if (double.IsNaN(lgp.X))
            {
                bool what;
            }
            point.Advance(dt, Domain);
            if (double.IsNaN(lgp.X))
            {
                bool what;
            }
        }
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
            if ((!point.CheckValid()) || x < Domain.XMin || x >= Domain.XMax || y < Domain.YMin || y >= Domain.YMax || p > Domain.PBase || p < Domain.PCeiling )
            {
                DeactivatePoint(node);
            }
            node = nextNode;
        }
    }

    public virtual bool WriteToFile()
    {
        return WriteToFile(OutputFilename);
    }
    private bool WriteToFile(string fileName)
    {
        bool success = true;

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
                    for (int k = 0; k < nProperties; k++)
                    {
                        properties2D[k,i,j] = PropertyHistory[k][i][j];
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
            
        return success;
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
            for (int k = 0; k < nProperties; k++)
            {
                properties[k][i] = point.GetProperty(PropertyNames[k]);
            }
            i++;
        }
        TimeHistory.Add(tCurr);
        XHistory.Add(xPoints);
        YHistory.Add(yPoints);
        PressureHistory.Add(pressurePoints);
        UIDHistory.Add(UIDs);
        AgeHistory.Add(ages);
        for (int k = 0; k < nProperties; k++)
        {
            PropertyHistory[k].Add(properties[k]);
        }
    }

    /*
    public virtual double GetPromotedProperty(IAdvected point, string property)
    {
        return point.GetProperty(property);
    }
    */
}