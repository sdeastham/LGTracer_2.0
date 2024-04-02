using System.Diagnostics;
using System.Globalization;
using Microsoft.VisualBasic.FileIO; // For parsing CSVs
using AtmosTools;

namespace LGTracer;

public class PointManagerFlight : PointManager
{
    protected DateTime LastSeedTime;
    protected Dictionary<string, Flight> FlightTable;
    protected double PointPeriod;
    public string SegmentsOutputFilename { get; protected set; }
    public bool ContrailSimulation { get; private set; }

    private bool IncludeSettling;

    public PointManagerFlight(long? maxPoints, DomainManager domain, string filename, DateTime initialSeedTime,
        double pointPeriod, string segmentsOutputFilename,
        bool verboseOutput = false, bool includeCompression = false, bool includeSettling = false, string[]? propertyNames = null,
        bool contrailSimulation = false) : base(maxPoints, domain, filename, verboseOutput, includeCompression,
        propertyNames)
    {
        LastSeedTime = initialSeedTime;
        FlightTable = [];
        PointPeriod = pointPeriod;
        SegmentsOutputFilename = segmentsOutputFilename;
        ContrailSimulation = contrailSimulation;
        IncludeSettling = includeSettling;
        // Run contrail test suite
        /*
        if (ContrailSimulation && !LGContrail.TestSAC(true))
        {
            throw new InvalidOperationException("Failed contrail test suite!");
        }
        */
    }

    protected override IAdvected CreatePoint()
    {
        if (ContrailSimulation)
        {
            return new LGContrail(VelocityCalc, IncludeCompression, IncludeSettling);
        }
        else
        {
            return new LGPointConnected(VelocityCalc);
        }
    }

    public void SimulateFlight(double originLon, double originLat, double destinationLon, double destinationLat,
        DateTime takeoffTime, double cruiseSpeedKPH, string? flightLabel = null, double pointPeriod = 60.0 * 5.0)
    {
        // Crude flight simulation between two airports. Currently only handles cruise
        double cruiseAltitude = 10.0; // km
        double flightDistance = Geodesy.GreatCircleDistance(originLon, originLat, destinationLon, destinationLat);
        DateTime endTime = takeoffTime + TimeSpan.FromSeconds(3600.0 * flightDistance / cruiseSpeedKPH);
        AddFlight([originLon, destinationLon], [originLat, destinationLat], [cruiseAltitude, cruiseAltitude],
            [takeoffTime, endTime], flightLabel: flightLabel, pointPeriod: pointPeriod);
    }

    public void PrintFlights()
    {
        foreach (string flightName in FlightTable.Keys)
        {
            Flight flight = FlightTable[flightName];
            Console.WriteLine($"Flight {flightName}");
            foreach (FlightSegment segment in flight.SegmentList)
            {
                Console.WriteLine(
                    $" --> Flight segment {segment.FlightID} with {segment.WaypointsRemaining} waypoints");
                foreach (FlightSegment.Waypoint waypoint in segment.Waypoints)
                {
                    Console.WriteLine(
                        $" ----> Waypoint: {waypoint.Lat,9:f2}N/{waypoint.Lon,9:f2} to deploy at {waypoint.DeploymentTime}");
                }
            }
        }
    }

    public void AddFlight(double[] lons, double[] lats, double[] pressureAltitudes, DateTime[] dateTimes,
        string? flightLabel = null, double pointPeriod = 60.0 * 5.0)
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

          If a flight label is not given then one is generated based on the origin location, end location, start time,
          and end time. This should be reasonably reliable as there should be essentially no times that two aircraft
          start from the same point go to the same point at exactly the same times.

          A word of caution - this function does NOT simulate cruise or similar, so if you want to simulate flight
          between two airports either do the simulation properly or give the pressures as being at cruise already.
        */
        int nSegments = lons.Length - 1;
        // Auto-generate a label if one is not provided
        var flightLabelSafe = flightLabel ?? $"{lons[0]}/{lats[0]}_{lons[nSegments]}/{lats[nSegments]}_{dateTimes[0]}";

        // Check that we don't already have this flight in the table
        if (FlightTable.ContainsKey(flightLabelSafe))
        {
            throw new ArgumentException($"Flight label {flightLabelSafe} duplicated");
        }
        
        // Assume that the flight date times are in order
        FlightTable.Add(flightLabelSafe, new Flight(dateTimes[0]));

        for (int i = 0; i < nSegments; i++)
        {
            double gcd = Geodesy.GreatCircleDistance(lons[i], lats[i], lons[i + 1], lats[i + 1]);
            double segmentDuration = (dateTimes[i + 1] - dateTimes[i]).TotalSeconds;
            double flightSpeed = 3600.0 * gcd / segmentDuration;
            AddFlightSegmentOrdered(lats[i], lons[i], dateTimes[i],
                lats[i + 1], lons[i + 1], dateTimes[i + 1],
                0.5 * (pressureAltitudes[i] + pressureAltitudes[i + 1]), flightSpeed, pointPeriod, flightLabelSafe);
        }
        
        // If you want to be very safe, you can run FlightTable[flightLabelSafe].SetTakeoff()
    }

    public void AddFlightSegmentOrdered(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod, string flightLabel)
    {
        // Convenience function - as long as flights are being added in full and their segments are added in order,
        // this will ensure that they are linked together
        Flight flight = FlightTable[flightLabel];
        LinkedListNode<FlightSegment>? node = flight.SegmentList.Last;
        AddFlightSegment(startLatitude, startLongitude, startDateTime, endLatitude, endLongitude,
            endDateTime, cruisePressureAltitude, flightSpeed, pointPeriod, flightLabel);
    }

    public void AddFlightSegment(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod, string flightLabel)
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
            endDateTime, cruisePressureAltitude, flightSpeed, pointPeriod, flightLabel);
        // Remove any waypoints which would already have been deployed
        seg.CullByDatetime(LastSeedTime);
        // Add to the table of flights
        if (!FlightTable.ContainsKey(flightLabel))
        {
            FlightTable.Add(flightLabel, new Flight());
        }

        FlightTable[flightLabel].SegmentList.AddLast(seg);
    }

    public override IAdvected NextPoint(double x, double y, double pressure)
    {
        LGContrail point = (LGContrail)base.NextPoint(x, y, pressure);
        point.Temperature = Domain.NearestNeighbor3D(x,y,pressure,Domain.TemperatureXYP);
        point.AmbientSpecificHumidity = Domain.NearestNeighbor3D(x,y,pressure,Domain.SpecificHumidityXYP); 
        point.InitiateContrail(0.35, 1.0e14, 1.0, 265.0, 0.7,
            waterVapourEmissionsIndex: 1.223, lowerHeatingValue: 43.2e6);
        return point;
    }

public override void Seed(double dt)
    {
        DateTime endTime = LastSeedTime + TimeSpan.FromSeconds(dt);
        foreach (string flightLabel in FlightTable.Keys)
        {
            Flight flight = FlightTable[flightLabel];
            if (flight.Takeoff > endTime) { continue; }
            int nWaypointsLeft = 0;
            // For each flight segment, check if there are any new seeds which should be created
            foreach (FlightSegment flightSegment in flight.SegmentList)
            {
                FlightSegment.Waypoint[] seedVector = flightSegment.DeploySeeds(LastSeedTime, endTime);
                foreach (FlightSegment.Waypoint seed in seedVector)
                {
                    LGPointConnected newPoint = (LGPointConnected)NextPoint(seed.Lon, seed.Lat, seed.Pressure);
                    newPoint.Connect(flight.LastPoint, null, flightLabel);
                    flight.UpdateLast(newPoint);
                }
                nWaypointsLeft += flightSegment.WaypointsRemaining;
            }
            // Don't need to keep this around if there are no waypoints left to seed
            if (nWaypointsLeft == 0) { FlightTable.Remove(flightLabel); }
        }

        LastSeedTime = endTime;
    }

    public override void Cull()
    {
        // Kill generated points as usual
        base.Cull();
        // Also kill flight segments which are completed
        foreach (string flightLabel in FlightTable.Keys)
        {
            Flight flight = FlightTable[flightLabel];
            LinkedListNode<FlightSegment>? node = flight.SegmentList.First;
            while (node != null)
            {
                LinkedListNode<FlightSegment>? nextNode = node.Next;
                FlightSegment flightSegment = node.Value;
                if (flightSegment.WaypointsRemaining == 0)
                {
                    flight.SegmentList.Remove(node);
                }
                node = nextNode;
            }
            // If no segments left, delete the flight from the table
            if (flight.SegmentList.First == null)
            {
                FlightTable.Remove(flightLabel);
            }
        }
    }

    public void ReadSegmentsFile(string segmentsFilePath)
    {
        // This will be tricky - need to find all the segments for one flight, sort them by generation time,
        // and then add new flights
        throw new NotImplementedException("Segment file reading not yet implemented");
    }

    public void ReadScheduleFile(string scheduleFilePath, string airportFilePath, DateTime? simulationStart=null, DateTime? simulationEnd=null )
    {
        bool cullByStart = (simulationStart != null);
        bool cullByEnd = (simulationEnd != null);
        // First read in the airport data file
        //List<Airport> airports = [];
        Dictionary<string, Airport> airports = [];
        using (TextFieldParser csvParser = new TextFieldParser(airportFilePath))
        {
            csvParser.CommentTokens = ["#"];
            csvParser.SetDelimiters([","]);
            csvParser.HasFieldsEnclosedInQuotes = true;

            // Get the columns
            string[] colNames = csvParser.ReadFields();
            int icaoIndex = Array.IndexOf(colNames, "ICAO name");
            int iataIndex = Array.IndexOf(colNames, "IATA name");
            int latIndex = Array.IndexOf(colNames, "Latitude");
            int lonIndex = Array.IndexOf(colNames, "Longitude");

            while (!csvParser.EndOfData)
            {
                string[] fields = csvParser.ReadFields();
                string airportName = fields[iataIndex]; // Use IATA for now
                airports.Add(airportName, new Airport(double.Parse(fields[lonIndex]),
                    double.Parse(fields[latIndex]), fields[iataIndex], fields[icaoIndex]));
            }
        }

        // Now the schedule file
        int nEntries = 0;
        int nFlights = 0;
        using (TextFieldParser csvParser = new TextFieldParser(scheduleFilePath))
        {
            csvParser.CommentTokens = ["#"];
            csvParser.SetDelimiters([","]);
            csvParser.HasFieldsEnclosedInQuotes = true;

            // Get the columns
            string[] colNames = csvParser.ReadFields();
            int originIndex = Array.IndexOf(colNames, "Origin airport");
            int destinationIndex = Array.IndexOf(colNames, "Destination airport");
            int startDateIndex = Array.IndexOf(colNames, "First date");
            int stopDateIndex = Array.IndexOf(colNames, "Last date");
            int weekdayIndex = Array.IndexOf(colNames, "Weekdays");
            int timeIndex = Array.IndexOf(colNames, "Takeoff time");
            // Also Airline and Equipment, but not currently used

            // Not entirely sure why this is needed
            CultureInfo cultureInfo = new CultureInfo("en-US");
            while (!csvParser.EndOfData)
            {
                string[] fields = csvParser.ReadFields();
                string origin = fields[originIndex];
                string destination = fields[destinationIndex];
                Airport originAirport = airports[origin];
                Airport destinationAirport = airports[destination];
                if (originAirport == destinationAirport)
                {
                    continue;
                }
                // Keep track of how many schedule entries we have read
                nEntries++;
                // Date of first takeoff
                string dtFull = $"{fields[startDateIndex]} {fields[timeIndex]}Z";
                DateTime startDate = DateTime.ParseExact(dtFull, "u", cultureInfo);
                // Date from which there will be no more takeoffs (file lists the last takeoff rather than the stop date
                // so we add one day to compensate
                dtFull = $"{fields[stopDateIndex]} {fields[timeIndex]}Z";
                DateTime endDate = DateTime.ParseExact(dtFull, "u", cultureInfo) + TimeSpan.FromDays(1);
                if ((cullByStart && endDate < simulationStart) || (cullByEnd && startDate >= simulationEnd))
                {
                    continue;
                }
                // Get the weekday indicators
                string weekdaysString = fields[weekdayIndex];
                bool[] weekdays = new bool[7];
                for (int i = 0; i < 7; i++)
                {
                    // NB: First entry means Sunday
                    weekdays[i] = weekdaysString[i] ==  'T';
                }
                // Add flight each valid day
                for (int iDay = 0; iDay < (endDate - startDate).TotalDays; iDay++)
                {
                    DateTime currentDate = startDate + TimeSpan.FromDays(iDay);
                    if (!weekdays[(int)currentDate.DayOfWeek]) { continue; }
                    if ((cullByStart && currentDate < simulationStart) || (cullByEnd && currentDate >= simulationEnd)) { continue; }
                    SimulateFlight(originAirport.Longitude, originAirport.Latitude,
                        destinationAirport.Longitude, destinationAirport.Latitude,
                        currentDate, 820.0, null, PointPeriod);
                    nFlights++;
                }
            }
        }

        if (VerboseOutput)
        {
            Console.WriteLine($"Schedule parsed. Found {nFlights} flights in {nEntries} schedule entries.");
        }
    }

    public override bool WriteToFile()
    {
        bool pointFileSuccess = base.WriteToFile();
        // Also write segments to file!
        bool segmentFileSuccess = true;
        return (segmentFileSuccess && pointFileSuccess);
    }

    private class Airport(double longitude, double latitude, string nameIATA, string nameICAO, double elevation = 0.0)
    {
        public double Longitude = longitude;
        public double Latitude = latitude;
        public string NameIATA = nameIATA;
        public string NameICAO = nameICAO;
        public double Elevation = elevation;
    }
    protected class Flight(DateTime? takeoff = null)
    {
        public LinkedList<FlightSegment> SegmentList = [];
        public LGPointConnected? LastPoint { get; protected set; } = null;
        public DateTime? Takeoff = takeoff;

        public void SetTakeoff()
        {
            foreach (FlightSegment segment in SegmentList)
            {
                DateTime segmentTime = segment.StartDateTime;
                if (Takeoff == null || segmentTime < Takeoff)
                {
                    Takeoff = (DateTime)segmentTime;
                }
            }
        }
        
        public void UpdateLast(LGPointConnected newPoint)
        {
            LastPoint = newPoint;
        }
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
    public string FlightID { get; private set; }
    
    public FlightSegment(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod,
        string flightLabel = "UNKNOWN")
    {
        StartDateTime = startDateTime;
        EndDateTime = endDateTime;
        StartLatLon = [startLatitude, startLongitude];
        EndLatLon = [endLatitude, endLongitude];
        // Convert from pressure altitude in km to pressure in Pa using the International Standard Atmosphere
        double cruisePressurePa = ISAtmos.AltitudeToPressure(cruisePressureAltitude);
        StartAltitudePa = cruisePressurePa;
        EndAltitudePa = cruisePressurePa;
        // A (hopefully!) unique identifier
        FlightID = flightLabel;

        // How far do we want to space the waypoints (in kilometers)?
        double waypointSpacing = flightSpeed * pointPeriod / 3600.0;
        
        // The waypoints (or seeds) are tracked as a linked list, where we will delete entries as they are dropped.
        // Once the list is empty, the segment can be safely deleted
        
        // Calculate how many points we will need to seed, when they should be seeded, and where
        // Get the locations of the seeds
        (double[] waypointLons, double[] waypointLats, _) = 
            Geodesy.GreatCircleWaypointsByLength(startLongitude, startLatitude, endLongitude, endLatitude, waypointSpacing);
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
            Waypoints.AddLast(new Waypoint(waypointLons[iWaypoint], waypointLats[iWaypoint], cruisePressurePa,
                waypointTimes[iWaypoint]));
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
    public class Waypoint(double lon, double lat, double pressure, DateTime deploymentTime)
    {
        public double Lon = lon;
        public double Lat = lat;
        public double Pressure = pressure;
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