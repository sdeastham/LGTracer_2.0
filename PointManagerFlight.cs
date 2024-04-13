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
    public bool ContrailSimulation { get; private set; }
    private bool IncludeSettling;
    private string AirportNameField;
    protected Random? RandomNumberGenerator;
    protected double MinimumPointLifetime;

    public PointManagerFlight( DomainManager domain, LGOptions configOptions, LGOptionsPointsFlights configSubOptions, Random rng ) : base(
        domain, configOptions, configSubOptions )
    {
        RandomNumberGenerator = rng;
        LastSeedTime = configOptions.Timing.StartDate;
        FlightTable = [];
        PointPeriod = configSubOptions.PointSpacing;
        ContrailSimulation = configSubOptions.ContrailSimulation;
        IncludeSettling = configSubOptions.IncludeSettling;
        MinimumPointLifetime = configSubOptions.MinimumLifetime;
        // Run contrail test suite
        /*
        if (ContrailSimulation && !LGContrail.TestSAC(true))
        {
            throw new InvalidOperationException("Failed contrail test suite!");
        }
        */
        AirportNameField = configSubOptions.UseIcao ? "ICAO" : "IATA";
    }

    protected override IAdvected CreatePoint()
    {
        return ContrailSimulation ? new LGContrail(VelocityCalc, IncludeCompression, IncludeSettling, MinimumPointLifetime) : new LGPointConnected(VelocityCalc, MinimumPointLifetime);
    }

    public bool SimulateFlight(double originLon, double originLat, double destinationLon, double destinationLat,
        DateTime takeoffTime, double cruiseSpeedKPH, string? flightLabel = null, double pointPeriod = 60.0 * 5.0,
        IAircraft? equipment = null)
    {
        // Crude flight simulation between two airports. Currently only handles cruise
        double cruiseAltitude = RandomNumberGenerator == null ? 10.0 : 9.0 + RandomNumberGenerator.NextDouble() * (12.0 - 9.0);
        double flightDistance = Geodesy.GreatCircleDistance(originLon, originLat, destinationLon, destinationLat);
        DateTime endTime = takeoffTime + TimeSpan.FromSeconds(3600.0 * flightDistance / cruiseSpeedKPH);
        return AddFlight([originLon, destinationLon], [originLat, destinationLat], [cruiseAltitude, cruiseAltitude],
               [takeoffTime, endTime], flightLabel: flightLabel, pointPeriod: pointPeriod, equipment);
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

    public bool AddFlight(double[] lons, double[] lats, double[] pressureAltitudes, DateTime[] dateTimes,
        string? flightLabel = null, double pointPeriod = 60.0 * 5.0, IAircraft? equipment = null)
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
            return false;
        }
        
        // Store the flight data
        Flight flight = new Flight(dateTimes,lons,lats,pressureAltitudes,pointPeriod,flightLabelSafe,LastSeedTime,equipment);
        
        // Assume that the flight date times are in order
        FlightTable.Add(flightLabelSafe, flight);
        return true;
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
            if (!flight.SegmentsReady)
            {
                flight.PrepSegments(LastSeedTime,endTime);
            }
            int nWaypointsLeft = 0;
            // For each flight segment, check if there are any new seeds which should be created
            foreach (FlightSegment flightSegment in flight.SegmentList)
            {
                FlightSegment.Waypoint[] seedVector = flightSegment.DeploySeeds(LastSeedTime, endTime);
                foreach (FlightSegment.Waypoint seed in seedVector)
                {
                    if (!Domain.InDomainXYP(seed.Lon, seed.Lat, seed.Pressure)) { continue; }
                    LGPointConnected newPoint = (LGPointConnected)NextPoint(seed.Lon, seed.Lat, seed.Pressure);
                    newPoint.Connect(flight.LastPoint, null, flightLabel, flight.ClearPoint);
                    flight.UpdateLast(newPoint);
                }
                nWaypointsLeft += flightSegment.WaypointsRemaining;
            }
            // Don't need to keep this around if there are no waypoints left to seed
            if (nWaypointsLeft == 0)
            {
                // Tell the final seeded point that no more will follow
                if (flight.LastPoint != null && flight.LastPoint.Active)
                {
                    flight.LastPoint.MakeFollower();
                }
                FlightTable.Remove(flightLabel);
            }
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
            if (flight is { SegmentsReady: true, SegmentList.First: null })
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
            int nameIndex = Array.IndexOf(colNames, $"{AirportNameField} name");

            while (!csvParser.EndOfData)
            {
                string[] fields = csvParser.ReadFields();
                string airportName = fields[nameIndex];
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
            int equipmentIndex = Array.IndexOf(colNames, "Equipment");
            // Also Airline, but not currently used

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
                // Add flight that takes off on each valid day
                // Start from the day before, as we need to carry over flights which took off the previous day
                DateTime? startCutoff = cullByStart ? simulationStart - TimeSpan.FromDays(1) : null;
                for (int iDay = -1; iDay < (endDate - startDate).TotalDays; iDay++)
                {
                    DateTime currentDate = startDate + TimeSpan.FromDays(iDay);
                    if (!weekdays[(int)currentDate.DayOfWeek]) { continue; }
                    if ((cullByStart && currentDate < startCutoff) || (cullByEnd && currentDate >= simulationEnd)) { continue; }
                    bool flightOK = SimulateFlight(originAirport.Longitude, originAirport.Latitude,
                        destinationAirport.Longitude, destinationAirport.Latitude,
                        currentDate, 820.0, null, PointPeriod);
                    if (flightOK)
                    {
                        nFlights++;
                    }
                }
            }
        }

        if (VerboseOutput)
        {
            Console.WriteLine($"Schedule parsed. Found {nFlights} flights in {nEntries} schedule entries.");
        }
    }

    private class Airport(double longitude, double latitude, string nameIATA, string nameICAO, double elevation = 0.0)
    {
        public double Longitude = longitude;
        public double Latitude = latitude;
        public string NameIATA = nameIATA;
        public string NameICAO = nameICAO;
        public double Elevation = elevation;
    }
    protected class Flight
    {
        public LinkedList<FlightSegment> SegmentList;
        public LGPointConnected? LastPoint { get; protected set; } = null;
        public DateTime Takeoff => WaypointTimes[0];
        public DateTime Landing => WaypointTimes[^1];
        public bool SegmentsReady;
        
        // Waypoint data
        private double[] WaypointLons;
        private double[] WaypointLats;
        private double[] WaypointAltitudes;
        private DateTime[] WaypointTimes;
        private string FlightLabel;
        private IAircraft? Equipment;
        private double PointPeriod;

        public Flight(DateTime[] waypointTimes, double[] lons, double[] lats, double[] pressureAltitudes,
            double pointPeriod, string flightLabel, DateTime cullTime, IAircraft? equipment)
        {
            WaypointTimes = waypointTimes;
            WaypointLons = lons;
            WaypointLats = lats;
            WaypointAltitudes = pressureAltitudes;
            Equipment = equipment;
            FlightLabel = flightLabel;
            PointPeriod = pointPeriod;
            SegmentsReady = false;
            SegmentList = [];
        }

        public void PrepSegments(DateTime stepStart, DateTime stepEnd)
        {
            if (stepEnd < Takeoff)
            {
                return;
            }
            int nSegments = WaypointLons.Length - 1;
            for (int i = 0; i < nSegments; i++)
            {
                if (stepStart > WaypointTimes[i+1] || stepEnd < WaypointTimes[i])
                {
                    continue;
                }
                double gcd = Geodesy.GreatCircleDistance(WaypointLons[i], WaypointLats[i], WaypointLons[i + 1], WaypointLats[i + 1]);
                double segmentDuration = (WaypointTimes[i + 1] - WaypointTimes[i]).TotalSeconds;
                double flightSpeed = 3600.0 * gcd / segmentDuration;
                FlightSegment seg = new FlightSegment(WaypointLats[i], WaypointLons[i], WaypointTimes[i],
                    WaypointLats[i + 1], WaypointLons[i + 1], WaypointTimes[i + 1],
                    0.5 * (WaypointAltitudes[i] + WaypointAltitudes[i + 1]), flightSpeed, PointPeriod,
                    FlightLabel, Equipment);
                // Cull waypoints in the segment that have no business existing
                seg.CullByDatetime(stepStart);
                SegmentList.AddLast(seg);
                SegmentsReady = true;
            }
        }
        
        public void UpdateLast(LGPointConnected newPoint)
        {
            LastPoint = newPoint;
        }

        public void ClearPoint(LGPointConnected caller)
        {
            if (caller == LastPoint)
            {
                LastPoint = null;
            }
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

    public IAircraft? Equipment;
    
    public FlightSegment(double startLatitude, double startLongitude, DateTime startDateTime,
        double endLatitude, double endLongitude, DateTime endDateTime,
        double cruisePressureAltitude, double flightSpeed, double pointPeriod,
        string flightLabel = "UNKNOWN", IAircraft? equipment = null)
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
        Equipment = equipment;

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