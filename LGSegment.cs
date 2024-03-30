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
    public double SegmentLength { get; protected set; }
    private bool HeadUpdated;
    private bool TailUpdated;
    public double XMidpoint { get; protected set; }
    public double YMidpoint { get; protected set; }
    
    public LGSegment(LGPointConnected head, LGPointConnected tail)
    {
        Head = head;
        Tail = tail;
        UpdateProperties();
        InitialSegmentLength = SegmentLength;
        HeadUpdated = false;
        TailUpdated = false;
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

        UpdateProperties();
        
        // Prepare for the next time step
        HeadUpdated = false;
        TailUpdated = false;
    }

    private void UpdateProperties()
    {
        // Update the segment's length
        CalculateSegmentLength();
        
        // Update the midpoint
        CalculateSegmentMidpoint();
    }

    private void CalculateSegmentLength()
    {
        SegmentLength = Geodesy.GreatCircleDistance(Head.X, Head.Y, Tail.X, Tail.Y);
    }

    public void CalculateSegmentMidpoint()
    {
        // Get the midpoint
        (double[] lons, double[] lats, double[] lengths) = Geodesy.GreatCircleWaypointsByCount(Tail.X, Tail.Y,
            Head.X, Head.Y, 3);
        XMidpoint = lons[1];
        YMidpoint = lats[1];
    }
}