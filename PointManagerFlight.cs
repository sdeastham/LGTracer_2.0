using System.Collections;
using System.Diagnostics;
using System.Xml.Serialization;
using AtmosTools;

namespace LGTracer;

public class PointManagerFlight : PointManager
{
    protected LinkedList<FlightSegment> FlightSegments;
    protected DateTime LastSeedTime;
    protected Dictionary<string,LinkedList<FlightSegment>> FlightTable;

    public PointManagerFlight(long? maxPoints, DomainManager domain, string filename, DateTime initialSeedTime,
        bool debug = false, bool includeCompression = false, string[]? propertyNames = null,
        double kgPerPoint = 1.0e12) : base(maxPoints, domain, filename, debug, includeCompression, propertyNames)
    {
        FlightSegments = [];
        LastSeedTime = initialSeedTime;
        FlightTable = [];
    }

    public void SimulateFlight(double originLon, double originLat, double destinationLon, double destinationLat,
        DateTime takeoffTime, double cruiseSpeedKPH, string? flightLabel = null, double pointPeriod = 60.0 * 5.0)
    {
        // Crude flight simulation between two airports. Currently only handles cruise
        double cruiseAltitude = 10.0; // km
        double flightDistance = Geodesy.GreatCircleDistance(originLon, originLat, destinationLon, destinationLat);
        DateTime endTime = takeoffTime + TimeSpan.FromSeconds(3600.0 * flightDistance / cruiseSpeedKPH);
        AddFlight([originLon,destinationLon],[originLat,destinationLat],[cruiseAltitude,cruiseAltitude],
            [takeoffTime,endTime], flightLabel: flightLabel, pointPeriod: pointPeriod);
    }

    public void PrintFlights()
    {
        foreach (FlightSegment segment in FlightSegments)
        {
            Console.WriteLine($"Flight segment {segment.FlightID} with {segment.WaypointsRemaining} waypoints");
            foreach (FlightSegment.Waypoint waypoint in segment.Waypoints)
            {
                Console.WriteLine($" --> Waypoint: {waypoint.Lat,9:f2}N/{waypoint.Lon,9:f2} to deploy at {waypoint.DeploymentTime}");
            }
        }
    }

    public void AddFlight(double[] lons, double[] lats, double[] pressureAltitudes, DateTime[] dateTimes,
        string? flightLabel=null, double pointPeriod = 60.0 * 5.0)
    {
        /* Add a flight as a series of waypoints
         Vector inputs (one value per waypoint):
          - lons                Longitude of each waypoint in degrees East
          - lats                Latitudes in degrees North
          - pressureAltitudes   Pressure altitudes in km
          - dateTimes           Times
          Scalar inputs (optional)
          - flightLabel         A unique string to identify the flight*
          - pointPeriod         Seconds between seeds to be dropped along the track
          
          If a flight label is not given then on is generated based on the origin location, end location, start time,
          and end time. This should be reasonably reliable as there should be essentially no times that two aircraft
          start from the same point go to the same point at exactly the same times.
          
          A word of caution - this function does NOT simulate cruise or similar, so if you want to simulate flight 
          between two airports either do the simulation properly or give the pressures as being at cruise already. 
        */
        int nSegments = lons.Length - 1;
        string flightLabelSafe;
        if (flightLabel == null)
        {
            flightLabelSafe = $"{lons[0]}/{lats[0]}_{lons[nSegments]}/{lats[nSegments]}_{dateTimes[0]}";
        }
        else
        {
            flightLabelSafe = flightLabel;
        }
        for (int i = 0; i < nSegments; i++)
        {
            double gcd = Geodesy.GreatCircleDistance(lons[i], lats[i], lons[i + 1], lats[i + 1]);
            double segmentDuration = (dateTimes[i + 1] - dateTimes[i]).TotalSeconds;
            double flightSpeed = 3600.0 * gcd / segmentDuration;
            AddFlightSegmentOrdered(lats[i], lons[i], dateTimes[i],
                lats[i + 1], lons[i + 1], dateTimes[i + 1],
                0.5 * (pressureAltitudes[i] + pressureAltitudes[i + 1]), flightSpeed, pointPeriod, flightLabelSafe);
        }
    }

    public void AddFlightSegmentOrdered(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod, string flightLabel)
    {
        // Convenience function - as long as flights are being added in full and their segments are added in order,
        // this will ensure that they are linked together
        LinkedListNode<FlightSegment>? node = FlightSegments.Last;
        AddFlightSegment(startLatitude, startLongitude, startDateTime, endLatitude, endLongitude,
            endDateTime, cruisePressureAltitude, flightSpeed, pointPeriod, flightLabel);
        if (FlightSegments.Last == null)
        {
            return;
        }
        if ((node != null) && (node.Value.FlightID == flightLabel))
        {
            FlightSegments.Last!.Value.PreviousSegment = node.Value;
        }
    }

    public void AddFlightSegment(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod, string flightLabel,
        FlightSegment? previousSegment = null)
    {
        // Define a new flight segment, to take place at constant cruise pressure altitude and flight speed
        // Units:
        // * start/end latitude         Degrees North
        // * start/end longitude        Degrees East
        // * start/end date time        C# date time objects
        // * altitude                   Kilometers, explicitly representing a pressure altitude
        // * pointPeriod                Seconds between production of a new point
        // The first point will be produced in the first timestep which contains the start time
        FlightSegment seg = new FlightSegment(startLatitude, startLongitude, startDateTime, endLatitude, endLongitude,
            endDateTime, cruisePressureAltitude, flightSpeed, pointPeriod, flightLabel, previousSegment);
        // Remove any waypoints which would already have been deployed
        seg.CullByDatetime(LastSeedTime);
        FlightSegments.AddLast(seg);
        // Add to the table of flights
        if (!FlightTable.ContainsKey(flightLabel))
        {
            FlightTable.Add(flightLabel, []);
        }
        FlightTable[flightLabel].AddLast(seg);
    }

    public override void Seed(double dt)
    {
        DateTime endTime = LastSeedTime + TimeSpan.FromSeconds(dt);
        // For each flight segment, check if there are any new seeds which should be created
        foreach (FlightSegment flightSegment in FlightSegments)
        {
            FlightSegment.Waypoint[] seedVector = flightSegment.DeploySeeds(LastSeedTime, endTime);
            foreach (FlightSegment.Waypoint seed in seedVector)
            {
                LGPoint newPoint = NextPoint(seed.Lon, seed.Lat, seed.Pressure);
                // TODO: Put all of this in the override for NextPoint
                // TODO: Determine if the new point's predecessor from the same flight exists
                // TODO: Write a plume segment class which represents a segment length
                // TODO: Create a segment between the predecessor and this point
                // TODO: Override LGPoint with LGPlume, which can carry a segment
                // TODO: Assign the segment to this point
                // TODO: Allow the point to carry/retrieve segment properties (nested class?)
            }
        }
        LastSeedTime = endTime;
    }

    public override void Cull()
    {
        // Kill generated points as usual
        base.Cull();
        // Also kill flight segments which are completed
        LinkedListNode<FlightSegment>? node = FlightSegments.First;
        while (node != null)
        {
            LinkedListNode<FlightSegment>? nextNode = node.Next;
            FlightSegment flightSegment = node.Value;
            if (flightSegment.WaypointsRemaining == 0)
            {
                FlightSegments.Remove(node);
            }
            node = nextNode;
        }
        // TODO: Go through subsegments and see if any have lost their head or tail
    }
}

public class FlightSegment
{
    public double[] StartLatLon { get; protected set; }
    public double[] EndLatLon { get; protected set; }
    public DateTime StartDateTime { get; protected set; }
    public DateTime EndDateTime { get; protected set; }

    // Altitudes are "pressure altitudes", meaning they are the pressure the ISA says corresponds to a given altitude
    private double StartAltitudePa; 
    private double EndAltitudePa;
    public LinkedList<Waypoint> Waypoints { get; protected set; }
    public int WaypointsRemaining => Waypoints.Count;
    public FlightSegment? PreviousSegment;
    public string FlightID { get; private set; }

    public FlightSegment(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod,
        string flightLabel = "UNKNOWN",FlightSegment? previousSegment = null)
    {
        StartDateTime = startDateTime;
        EndDateTime = endDateTime;
        StartLatLon = [startLatitude, startLongitude];
        EndLatLon = [endLatitude, endLongitude];
        // Convert from pressure altitude in km to pressure in Pa using the International Standard Atmosphere
        double cruisePressurePa = ISAtmos.AltitudeToPressure(cruisePressureAltitude);
        StartAltitudePa = cruisePressurePa;
        EndAltitudePa = cruisePressurePa;
        PreviousSegment = previousSegment;
        // A (hopefully!) unique identifier
        FlightID = flightLabel;

        // How far do we want to space the waypoints (in kilometers)?
        double waypointSpacing = flightSpeed * pointPeriod / 3600.0;
        
        // The waypoints (or seeds) are tracked as a linked list, where we will delete entries as they are dropped.
        // Once the list is empty, the segment can be safely deleted
        
        // Calculate how many points we will need to seed, when they should be seeded, and where
        // Get the locations of the seeds
        (double[] waypointLons, double[] waypointLats, double[] initialSubsegmentLengths) = Geodesy.GreatCircleWaypoints(startLongitude,
            startLatitude, endLongitude, endLatitude, waypointSpacing);
        int nWaypoints = waypointLons.Length;
        // Something to track which seeds have been dropped
        // Now the times at which the seeds should be created
        DateTime[] waypointTimes = new DateTime[nWaypoints];
        // All but the last seed are spaced evenly
        for (int iWaypoint = 0; iWaypoint < nWaypoints; iWaypoint++)
        {
            waypointTimes[iWaypoint] = StartDateTime + TimeSpan.FromSeconds(iWaypoint * pointPeriod);
        }
        waypointTimes[nWaypoints - 1] = EndDateTime;

        Waypoints = [];
        for (int iWaypoint = 0; iWaypoint < nWaypoints; iWaypoint++)
        {
            double segmentDistance = 0.0;
            if (iWaypoint < (nWaypoints - 1))
            {
                segmentDistance = initialSubsegmentLengths[iWaypoint];
            }
            // The leader flag is used to help determine whether the predecessor waypoint is from this flight segment
            // or from a previous flight segment
            Waypoints.AddLast(new Waypoint(waypointLons[iWaypoint], waypointLats[iWaypoint], cruisePressurePa,
                waypointTimes[iWaypoint], segmentDistance, iWaypoint == 0));
        }
    }

    public void CullByDatetime(DateTime now)
    {
        // Delete seeds which should already have been deployed
        LinkedListNode<Waypoint>? node = Waypoints.First;
        while (node != null)
        {
            LinkedListNode<Waypoint>? nextNode = node.Next;
            Waypoint waypoint = node.Value;
            if (waypoint.Deploy(now))
            {
                Waypoints.Remove(node);
            }
            node = nextNode;
        }
    }

    public Waypoint[] DeploySeeds(DateTime startTime, DateTime endTime)
    {
        // Returns lists of waypoints which are eligible for deployment in a given time period
        // Also removes those waypoints from the flight segment's overarching list so they'll
        // not be redeployed in the future
        LinkedListNode<Waypoint>? node = Waypoints.First;
        List<Waypoint> seedList = [];
        while (node != null)
        {
            LinkedListNode<Waypoint>? nextNode = node.Next;
            if (node.Value.Deploy(startTime,endTime))
            {
                seedList.Add(node.Value);
                Waypoints.Remove(node);
            }

            node = nextNode;
        }
        return seedList.ToArray();
    }

    // Nested class because no-one else needs this
    public class Waypoint(double lon, double lat, double pressure, DateTime deploymentTime, double segmentDistance = 0.0, bool leader = false)
    {
        public double Lon = lon;
        public double Lat = lat;
        public double Pressure = pressure;
        public double SegmentDistance = segmentDistance;
        public bool Leader = leader;
        public DateTime DeploymentTime = deploymentTime;

        public bool Deploy(DateTime now)
        {
            return DeploymentTime < now;
        }
        public bool Deploy(DateTime startTime, DateTime endTime)
        {
            return (DeploymentTime >= startTime && DeploymentTime < endTime);
        }
    }
}