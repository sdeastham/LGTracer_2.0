using AtmosTools;
namespace LGTracer;

public interface IAircraft
{
    public void Advance(double dt);
    public double InitialWakeCirculation(double airDensity);
    public double FlightSpeed();
    public double WaterEmissionsPerMeter();
    public double FuelFlowRate();
    public double NonVolatileNumberEmissionsPerMeter(); // Particles per second
    public double PlumeArea();
    public double OverallEfficiency()
    {
        return 0.35;
    }
    public void InitializeFlight(double distance, double loadFactor)
    {
        // Default to doing nothing
    }
}

public static class AircraftFactory
{
    public static IAircraft CreateAircraftU2016(string aircraftName, bool simpleDefault=false)
    {
        double keroseneDensity = 0.82e3; // kg/m3
        double kgPerLiterKerosene = 1e-3 * keroseneDensity;
        switch (aircraftName.ToUpper())
        {
            case "CRJ":
                return new Aircraft(maximumTakeoffMass: 24000.0, wingspan: 21.2, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(1.77e-3),
                    maximumFuelMass: 6489.0, operatingEmptyMass: 13835.0,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 6124.0, flightPayloadFraction: 0.8,
                    cruiseCeiling: 12.496, cruiseMach: 0.79, range: 3000.0, waterVapourEmissionsIndex: 1.223);
            case "A320":
                // U2016 uses 34.4 m as wingspan
                // Data here from A320 specification access 2024-04-21
                // https://web.archive.org/web/20120124123133/http://www.airbus.com/aircraftfamilies/passengeraircraft/a320family/a320/specifications/
                // Cruise ceiling from EASA type certification data sheet for A319/320/321 access 2024-04-21
                // (https://www.easa.europa.eu/en/downloads/16507/en)
                return new Aircraft(maximumTakeoffMass: 78000.0, wingspan: 34.1, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(3.70e-3),
                    maximumFuelMass: 30191 * kgPerLiterKerosene, operatingEmptyMass: 42600.0,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 16600.0,
                    cruiseCeiling: 12.500, cruiseMach: 0.82, range: 6150.0, waterVapourEmissionsIndex: 1.223);
            case "B737":
            case "737":
            case "B737-800":
            case "737-800":
                // U2016 uses 34.4 m as wingspan
                // Data here is for 737-800 with winglets
                return new Aircraft(maximumTakeoffMass: 79016.0, wingspan: 35.79, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(3.70e-3),
                    maximumFuelMass: 26022 * kgPerLiterKerosene, operatingEmptyMass: 41413.0,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 12856.0,
                    cruiseCeiling: 12.497, cruiseMach: 0.789, range: 5436.0, waterVapourEmissionsIndex: 1.223);
            case "A300":
                // A300-600R
                return new Aircraft(maximumTakeoffMass: 171700.0, wingspan: 44.84, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(7.26e-3),
                    maximumFuelMass: 53505.0, operatingEmptyMass: 88626.0,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 41374.0,
                    cruiseCeiling: 12.192, cruiseMach: 0.78, range: 7500.0, waterVapourEmissionsIndex: 1.223);
            case "B767":
                // B767-300
                // At 12 km, speed of sound is 274.07 m/s; so if the given speed is 850 kph, that means Mach 0.86
                // Range is actually given as being for 269 pax and an OEW of 85.2 tonnes..
                return new Aircraft(maximumTakeoffMass: 158800.0, wingspan: 47.57, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(7.26e-3),
                    maximumFuelMass: 73400.0, operatingEmptyMass: 86.1e3,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 40.0e3,
                    cruiseCeiling: 13.1, cruiseMach: 0.86, range: 7200.0, waterVapourEmissionsIndex: 1.223);
            case "A350":
                // A350-900
                return new Aircraft(maximumTakeoffMass: 283.0e3, wingspan: 64.75, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(15.0e-3),
                    maximumFuelMass: 110.5e3, operatingEmptyMass: 142.4e3,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 53.3e3,
                    cruiseCeiling: 13.1, cruiseMach: 0.85, range: 15372.0, waterVapourEmissionsIndex: 1.223);
            case "B777":
                // B777-300ER
                // Data taken from Boeing for a B777-300ER with GE engines (performance summary)
                // https://www.boeing.com/content/dam/boeing/boeingdotcom/company/about_bca/startup/pdf/historical/777_passenger.pdf
                // Using GE90-115BL engines.
                // Estimated max payload mass as max zero fuel weight minus operating empty weight.
                return new Aircraft(maximumTakeoffMass: 351530.0, wingspan: 64.8, numberOfEngines: 2,
                    cruiseFuelFlowRate: DeriveFuelFlow(15.0e-3),
                    maximumFuelMass: 181280.0 * kgPerLiterKerosene, operatingEmptyMass: 168780.0,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 68.9e3,
                    cruiseCeiling: 13.1, cruiseMach: 0.84, range: 14685.0, waterVapourEmissionsIndex: 1.223);
            case "B747":
                // B747-400 (not the 400ER)
                // Range is specifically with max payload
                // Data from https://www.airlines-inform.com/commercial-aircraft/boeing-747-400.html
                return new Aircraft(maximumTakeoffMass: 396.9e3, wingspan: 64.4, numberOfEngines: 4,
                    cruiseFuelFlowRate: DeriveFuelFlow(13.8e-3),
                    maximumFuelMass: 216840.0 * kgPerLiterKerosene, operatingEmptyMass: 181.120e3,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 70.620e3,
                    cruiseCeiling: 13.75, cruiseMach: 0.855, range: 13430.0, waterVapourEmissionsIndex: 1.223);
            case "A380":
                return new Aircraft(maximumTakeoffMass: 575.0e3, wingspan: 79.75, numberOfEngines: 4,
                    cruiseFuelFlowRate: DeriveFuelFlow(20.0e-3),
                    maximumFuelMass: 253983.0, operatingEmptyMass: 277.0e3,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 84.0e3,
                    cruiseCeiling: 13.15, cruiseMach: 0.85, range: 14800.0, waterVapourEmissionsIndex: 1.223);
            case "B787":
                // B787-9
                // Fuel flow rate very roughly estimated based on 5500 kg of fuel burn per hour
                return new Aircraft(maximumTakeoffMass: 254.7e3, wingspan: 60.12, numberOfEngines: 2,
                    cruiseFuelFlowRate: 1.53,
                    maximumFuelMass: 101.456e3, operatingEmptyMass: 129.0e3,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 53.0e3,
                    cruiseCeiling: 13.1, cruiseMach: 0.85, range: 14140.0, waterVapourEmissionsIndex: 1.223);
            /*
            case "??":
                return new Aircraft(maximumTakeoffMass: 0.0, wingspan: 0.0, numberOfEngines: 0,
                    cruiseFuelFlowRate: DeriveFuelFlow(0.0e-3),
                    maximumFuelMass: 0.0, operatingEmptyMass: 0.0,
                    nonVolatileNumberEmissionsIndex: 1.0e14, maximumPayloadMass: 0.0,
                    cruiseCeiling: 0.0, cruiseMach: 0.0, range: 0.0, waterVapourEmissionsIndex: 1.223);
            */
            default:
                if (!simpleDefault) throw new ArgumentException($"Aircraft {aircraftName} not recognized.");
                return new SimpleAircraft(wingspan: 60.0, numberOfEngines: 2);;
        }
    }

    private static double DeriveFuelFlow(double waterVapourEmissionsRate, double cruiseSpeed=230.0, double eiH2O=1.25)
    {
        // Derive kg/s of fuel flow given the water vapour emissions rate and cruise speed from Unterstrasser 2016.
        // An EI of 1.25 is used by default as this was the value used in the original paper when deriving the water
        // vapour emissions rate. This function inverts equation A8 of said paper.
        return cruiseSpeed * waterVapourEmissionsRate / eiH2O;
    }

    private static double DeriveAircraftCruisingMass(double initialCirculation, double wingSpan, double cruiseSpeed=230.0)
    {
        // Equations A1 and A2 of Unterstrasser 2016
        double vortexSeparation = Math.PI * wingSpan / 4.0;
        return initialCirculation * 0.4 * vortexSeparation * cruiseSpeed / Physics.G0;
    }
}

public class Aircraft : IAircraft
{
    private double MaximumTakeoffMass; // kg
    private double MaximumFuelMass;
    private double MaximumPayloadMass;
    private readonly double Wingspan; // m
    private readonly double InitialVortexSeparation; // m
    private readonly double CruiseCeiling; // m
    private readonly double CruiseMach; // Mach number
    private readonly double Range; // m
    private readonly double TypicalCruiseSpeed; // m/s
    private readonly double CruiseFuelFlowRate; // kg/s, summed over all engines
    private int NumberOfEngines; // integer
    private double NonVolatileNumberEmissionsIndex; // Particles per kg at cruise
    private readonly double EmptyMass;
    private double CurrentMass => EmptyMass + PayloadMass + FuelMass;
    private double WaterVapourEmissionsIndex;
    // We want the following to be updateable externally so that new aircraft can be instantiated and then modified
    public double FuelMass; // kg
    public double PayloadMass; // kg
    private double PerformanceFactor; // Overall efficiency * L/D
    // SReserve requirement is 30 minutes of flight at 1,500 ft; only have cruise fuel flow rate so use that
    private double ReserveFuelMass => 1800.0 * CruiseFuelFlowRate;

    public Aircraft(double maximumTakeoffMass, double maximumFuelMass, double maximumPayloadMass, double operatingEmptyMass,
        double wingspan, int numberOfEngines, double cruiseFuelFlowRate,
        double cruiseCeiling, double cruiseMach, double range,
        double nonVolatileNumberEmissionsIndex = 1.0e14, double flightPayloadFraction=0.8,
        double flightInitialFuelFraction = 1.0, double waterVapourEmissionsIndex=1.223)
    {
        /*
         * Range is assumed to be the 
         * Mass calculations:
         *  -> maximumTakeoffMass           Maximum takeoff mass (kg)
         *  -> maximumFuelMass              Maximum fuel mass (kg)
         *  -> maximumPayloadMass           Maximum payload mass (kg)
         *  -> operatingEmptyMass           Mass without fuel or payload (kg)
         *  -> flightPayloadFraction        The fraction of the max payload mass used on this flight (-)
         *  -> flightInitialFuelFraction    The fraction of the maximum fuel present when flight is instantiated (-)
         * IMPORTANT: MTOM is LESS than OEM + max fuel + max payload (as max payload cannot be carried max range)
         * Range and cruise ceiling are to be supplied in km but will be stored in meters
         */
        // Aircraft characteristics
        Wingspan = wingspan;
        NumberOfEngines = numberOfEngines;
        MaximumTakeoffMass = maximumTakeoffMass;
        MaximumFuelMass = maximumFuelMass;
        MaximumPayloadMass = maximumPayloadMass;
        FuelMass = maximumFuelMass * flightInitialFuelFraction;
        PayloadMass = maximumPayloadMass * flightPayloadFraction;
        // EmptyMass is the operating empty weight
        EmptyMass = operatingEmptyMass;
        
        // Cruise characteristics
        CruiseFuelFlowRate = cruiseFuelFlowRate; // All engines combined, kg/s
        // Calculate cruise speed given typical conditions at 10 km altitude and the cruise Mach
        CruiseMach = cruiseMach;
        double cruiseTemperature = ISAtmos.AltitudeToTemperature(Math.Min(cruiseCeiling,10.0));
        TypicalCruiseSpeed = IdealGases.SpeedOfSound(cruiseTemperature) * cruiseMach;
        Range = range * 1.0e3; // m
        CruiseCeiling = cruiseCeiling * 1.0e3; // m
        
        // Engine properties
        NonVolatileNumberEmissionsIndex = nonVolatileNumberEmissionsIndex;
        // Equation A2 of Unterstrasser (2016)
        InitialVortexSeparation = Math.PI * Wingspan / 4.0;
        WaterVapourEmissionsIndex = waterVapourEmissionsIndex;
        
        // Derive the performance factor (efficiency * L/D) so that we can later estimate fuel mass requirements
        PerformanceFactor = EstimatePerformance();
    }

    public void Advance(double dt)
    {
        FuelMass = double.Max(0.0, FuelMass - (CruiseFuelFlowRate * dt));
    }

    public double InitialWakeCirculation(double airDensity)
    {
        // Equation A1 of Unterstrasser (2016)
        return PhysConstants.G0 * CurrentMass / (airDensity * InitialVortexSeparation * FlightSpeed());
    }

    public double FuelFlowRate()
    {
        return CruiseFuelFlowRate;
    }
    
    public static double VortexPairDescentDistance(double circulation, double vortexSeparation)
    {
        return circulation / (2.0 * Math.PI * vortexSeparation);
    }

    public double FlightSpeed()
    {
        return TypicalCruiseSpeed;
    }

    public double WaterEmissionsPerMeter()
    {
        return FuelFlowRate() * WaterVapourEmissionsIndex / FlightSpeed();
    }
    
    public double NonVolatileNumberEmissionsPerMeter()
    {
        return FuelFlowRate() * NonVolatileNumberEmissionsIndex / FlightSpeed();
    }
    
    public double PlumeArea()
    {
        // Equations A6 and A7 of Unterstrasser 2016
        double plumeRadius = 1.5 + 0.314 * Wingspan;
        return 4.0 * Math.PI * plumeRadius * plumeRadius;
    }

    private double EstimatePerformance()
    {
        // Assume that the maximum range given is the range at maximum payload
        // Invert the Breguet range equation to estimate the product of overall efficiency and L/D
        const double fuelLowerHeatingValue = 46.2e6; // J/kg
        // Assume range is based on max fuel load and max payload
        // Assume that all aircraft land with 30 minutes of fuel remaining
        double weightRatio = MaximumTakeoffMass / (MaximumTakeoffMass - (MaximumPayloadMass + ReserveFuelMass));
        return PhysConstants.G0 * Range / (fuelLowerHeatingValue * Math.Log(weightRatio));
    }

    public double EstimateFuelRequirement(double distance, double fuelLowerHeatingValue = 46.2e6)
    {
        // Given a distance in kilometers, how much fuel is needed?
        // IMPORTANT: Payload fraction should be set before using this!
        // Assume again that aircraft always land with 30 minutes of fuel
        double finalWeight = EmptyMass + ReserveFuelMass + PayloadMass;
        // Based on inverting the Breguet range equation again
        double fractionalChange =
            (Math.Exp(PhysConstants.G0 * distance * 1.0e3 / (PerformanceFactor * fuelLowerHeatingValue)) - 1.0);
        return finalWeight * fractionalChange;
    }

    private void InitializeByDistance(double distance, double loadFactor, double fuelMargin)
    {
        // Allow MTOW to be exceeded by this factor
        const double maxExcess = 0.05;
        PayloadMass = loadFactor * MaximumPayloadMass;
        FuelMass = ReserveFuelMass + (1.0 + fuelMargin) * EstimateFuelRequirement(distance);
        if (CurrentMass > (1.0 + maxExcess) * MaximumTakeoffMass)
        {
            throw new ArgumentOutOfRangeException(nameof(distance),$"Takeoff mass ({CurrentMass*1.0e-3:f2} tonnes) for a {distance:f2} km flight exceeds maximum take-off mass ({MaximumTakeoffMass*1.0e-3:f2} tonnes).");
        }
    }
    
    private void InitializeByFuelMass(double fuelMass, double loadFactor)
    {
        PayloadMass = loadFactor * MaximumPayloadMass;
        FuelMass = fuelMass;
        if (CurrentMass > MaximumTakeoffMass)
        {
            throw new ArgumentOutOfRangeException(nameof(fuelMass),"Takeoff mass exceeds maximum take off mass.");
        }
    }
    
    // Use "initialize by distance" by default
    public void InitializeFlight(double distance, double loadFactor=0.8)
    {
        InitializeByDistance(distance, loadFactor, fuelMargin: 0.05);
    }
}

public class SimpleAircraft : IAircraft
{
    private readonly double Wingspan;
    private readonly double WingspanSquared;
    private readonly double NonVolatileNumberEmissionsIndex;
    public SimpleAircraft(double wingspan, int numberOfEngines, double typicalCruiseSpeed = 230.0,
        double nonVolatileNumberEmissionsIndex = 1.0e14)
    {
        Wingspan = wingspan;
        WingspanSquared = wingspan * wingspan;
        NonVolatileNumberEmissionsIndex = nonVolatileNumberEmissionsIndex;
    }
    public double InitialWakeCirculation(double _)
    {
        // Equation A5 of Unterstrasser (2016)
        return (10.0 * Wingspan - 70.0);
    }

    public void Advance(double dt)
    {
        // Null - nothing to advance
    }

    public double FuelFlowRate()
    {
        return WaterEmissionsPerMeter() / 1.233;
    }

    public double FlightSpeed()
    {
        return 230.0;
    }
    
    public double WaterEmissionsPerMeter()
    {
        // Equation A8 of Unterstrasser (2016)
        return 0.020 * WingspanSquared / 6400.0;
    }
    
    public double NonVolatileNumberEmissionsPerMeter()
    {
        return FuelFlowRate() * NonVolatileNumberEmissionsIndex;
    }
    
    public double PlumeArea()
    {
        // Equations A6 and A7 of Unterstrasser 2016
        double plumeRadius = 1.5 + 0.314 * Wingspan;
        return 4.0 * Math.PI * plumeRadius * plumeRadius;
    }
}