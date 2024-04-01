using System.Data;
using System.Reflection.Metadata;
using AtmosTools;
using LGTracer;
using YamlDotNet.Serialization.ObjectGraphTraversalStrategies;

public class LGSegment
{
    public LGPointConnected Head { get; protected set; }
    public LGPointConnected Tail { get; protected set; }
    public double InitialSegmentLength { get; protected set; }
    private bool HeadUpdated;
    private bool TailUpdated;
    private string Source;
    private uint UniqueID;
    
    // Avoid unnecessary recalculations
    private double _SegmentLength;

    public double SegmentLength
    {
        get
        {
            UpdateProperties();
            return _SegmentLength;
        }
    }
    
    private double _XMidpoint;
    public double XMidpoint
    {
        get
        {
            UpdateProperties();
            return _XMidpoint;
        }
    }

    private double _YMidpoint;
    public double YMidpoint
    {
        get
        {
            UpdateProperties();
            return _YMidpoint;
        }
    }

    protected bool Stale { get; private set; }
    
    public LGSegment(LGPointConnected head, LGPointConnected tail, uint uniqueId, string source="UNKNOWN")
    {
        Head = head;
        Tail = tail;
        Stale = true;
        UpdateProperties();
        InitialSegmentLength = SegmentLength;
        HeadUpdated = false;
        TailUpdated = false;
        Stale = false;
        UniqueID = uniqueId;
        Source = source;
    }

    public void Advance(LGPointConnected source)
    {
        // Only perform segment calculations once both the head and tail have been updated to their new positions
        if (source == Head)
        {
            HeadUpdated = true;
        }
        else if (source == Tail)
        {
            TailUpdated = true;
        }
        if (!(HeadUpdated && TailUpdated)) { return; }

        // Properties will need to be updated
        Stale = true;
        
        // Prepare for the next time step
        HeadUpdated = false;
        TailUpdated = false;
    }

    private void UpdateProperties()
    {
        // If the properties are still up to date,
        // don't bother re-calculating (expensive!)
        if (!Stale) return;
        // Update the segment's length
        CalculateSegmentLength();
        // Update the midpoint
        CalculateSegmentMidpoint();
        Stale = false;
    }

    private void CalculateSegmentLength()
    {
        _SegmentLength = Geodesy.GreatCircleDistance(Head.X, Head.Y, Tail.X, Tail.Y);
    }

    private void CalculateSegmentMidpoint()
    {
        // Get the midpoint
        (double[] lons, double[] lats, double[] lengths) = Geodesy.GreatCircleWaypointsByCount(Tail.X, Tail.Y,
            Head.X, Head.Y, 3);
        _XMidpoint = lons[1];
        _YMidpoint = lats[1];
    }

    public double GetProperty(string property)
    {
        switch (property.ToLower().Replace("_",""))
        {
            case "stretch":
                return SegmentLength / InitialSegmentLength;
            case "length":
                return SegmentLength;
            case "tailx":
            case "taillon":
                return Tail.X;
            case "taily":
            case "taillat":
                return Tail.Y;
            case "tailp":
            case "tailpressure":
                return Tail.Pressure;
            case "headx":
            case "headlon":
                return Head.X;
            case "heady":
            case "headlat":
                return Head.Y;
            case "headp":
            case "headpressure":
                return Head.Pressure;
            default:
                throw new ArgumentException($"No property for LGSegment called {property}");
        }
    }
}