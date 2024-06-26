﻿using System.Diagnostics;

namespace LGTracer;

public class LGPointConnected(
    Func<double, double, double, (double, double, double)> vCalc, double minimumPointLifetime=0.0)
    : LGPoint(vCalc)
{
    /*
     LGPointConnected objects are intended to represent vertices along a chain of points, e.g. representing a
     trajectory  
     These have two major differences from standard LGPoints:
     - Each is associated with a predecessor point (if one exists)
     - When an LGPoint has a predecessor, it will also have a "segment" associated [splitting?]
    */
    protected LGPointConnected? Previous;
    protected LGPointConnected? Next;
    protected LGSegment? Segment;
    // Cleanup allows callbacks - it will be called with this point as an argument when the point is deactivated
    protected Action<LGPointConnected>? Cleanup;
    // Indicates whether this was the most recently created point in its chain
    protected bool IsLeader = false;
    // How long before we let a point die even if it is no longer connecteD?
    protected double MinimumPointLifetime = minimumPointLifetime * 3600.0;

    public void Connect(LGPointConnected? predecessor, uint? segmentID = null, string segmentSource = "UNKNOWN", Action<LGPointConnected>? cleanup=null)
    {
        Previous = predecessor;
        // Action to call on death
        Cleanup = cleanup;
        // Indicate that this is the most recently-created point
        IsLeader = true;
        if (Previous == null) { return; }
        // Update the previous connected point
        Previous.Next = this;
        // Cannot create a segment if the previous point is inactive
        if (!Previous.Active) return;
        // Use the point ID if no segment ID is given
        CreateSegment(segmentID ?? this.UID, segmentSource);
        Previous.IsLeader = false;
    }

    public void MakeFollower()
    {
        IsLeader = false;
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
        Cleanup?.Invoke(this);
        base.Deactivate();
    }

    public void CreateSegment(uint uniqueID, string source)
    {
        if (Previous == null) { return; }
        Debug.Assert(Segment == null,"Segment already exists");
        // Create segment and associate it to this class
        Segment = new LGSegment(this,Previous,uniqueID,source);
    }

    public void DestroySegment()
    {
        // Delete the segment
        Segment = null;
    }

    public override void Advance(double dt, DomainManager domain)
    {
        if (!Active) return;
        base.Advance(dt, domain);
        // Segments advance once both the head and tail indicate that they are updated
        // Update the segment we own (are head of) if not null
        Segment?.Advance(this);
        // If we are tail of a segment, also update that
        Next?.Segment?.Advance(this);
    }

    public override double GetProperty(string property)
    {
        if (!property.StartsWith("segment", StringComparison.CurrentCultureIgnoreCase))
        {
            return base.GetProperty(property);
        }
        if (Segment == null)
        {
            return double.NaN;
        }
        return Segment.GetProperty(property.Substring(7));
    }

    public bool HasNeighbours()
    {
        return IsLeader || !(Next is not { Active: true } || Previous is not { Active: true });
    }

    public override bool CheckValid()
    {
        // Prevent point from being deactivated if it has not yet reached minimum lifetime
        return (Age < MinimumPointLifetime) || base.CheckValid();
    }
}