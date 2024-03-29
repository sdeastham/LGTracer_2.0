using System.Diagnostics;

namespace LGTracer;

public class LGPointConnected(
    Func<double, double, double, (double, double, double)> vCalc,
    bool includeCompression = false)
    : LGPoint(vCalc, includeCompression)
{
    /* LGPointConnected objects are intended to represent vertices along a chain of points, e.g. representing a
     trajectory  
     These have two major differences from standard LGPoints:
     - Each is associated with a predecessor point (if one exists)
     - When an LGPoint has a predecessor, it will also have a "segment" associated [splitting?]
    */
    private LGPointConnected? Previous;
    private LGPointConnected? Next;
    //private LGSegment? Segment;

    public void Activate(double x, double y, double pressure, uint uniqueID, LGPointConnected? predecessor)
    {
        Previous = predecessor;
        if (Previous != null)
        {
            Previous.Next = this;
        }
        CreateSegment();
        base.Activate(x, y, pressure, uniqueID);
    }

    public override void Deactivate()
    {
        if (Previous != null)
        {
            // Deactivate the segment
            DestroySegment();
            // We don't want to just stitch these together - with the death of this point the chain is broken
            Previous.Next = null;
        }

        if (Next != null)
        {
            // Tell the next node that its segment is dead
            Next.DestroySegment();
            // We don't want to just stitch these together - with the death of this point the chain is broken
            Next.Previous = null;
        }
        base.Deactivate();
    }

    public void CreateSegment()
    {
        if (Previous == null)
        {
            return;
        }
        //Debug.Assert(Segment != null,"Segment already exists");
        // Create segment and associate it to this class
    }

    public void DestroySegment()
    {
        if (Previous == null)
        {
            return;
        }
        // Delete the segment
    }
}