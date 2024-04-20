using AtmosTools;
namespace LGTracer;

public interface IAircraft
{
    public void Advance(double dt);
    public double InitialWakeCirculation(double airDensity);
    public double FlightSpeed();
    public double WaterEmissionsPerMeter();
    public double NonVolatileNumberEmissionsPerMeter(); // Particles per second
    public double PlumeArea();
}

public class Aircraft : IAircraft
{
    private double MaximumTakeoffMass; // kg
    private readonly double Wingspan; // m
    private readonly double InitialVortexSeparation; // m
    private readonly double TypicalCruiseSpeed; // m/s
    private readonly double CruiseFuelFlowRate; // kg/s, summed over all engines
    private int NumberOfEngines; // integer
    private double FuelMass; // kg
    private double NonVolatileNumberEmissionsIndex; // Particles per kg at cruise
    private readonly double PayloadMass;
    private readonly double EmptyMass;
    private double CurrentMass => EmptyMass + PayloadMass + FuelMass;
    private double WaterVapourEmissionsIndex;

    public Aircraft(double maximumTakeoffMass, double wingspan, int numberOfEngines, double cruiseFuelFlowRate, 
        double maximumFuelFraction, double typicalCruiseSpeed = 230.0, double initialFuelFraction = 1.0,
        double nonVolatileNumberEmissionsIndex = 1.0e14, double maximumPayloadMass=0.0, double flightPayloadFraction=0.0,
        double waterVapourEmissionsIndex=1.223)
    {
        MaximumTakeoffMass = maximumTakeoffMass;
        Wingspan = wingspan;
        double maximumFuelMass = maximumTakeoffMass * maximumFuelFraction;
        FuelMass = maximumFuelMass * initialFuelFraction;
        CruiseFuelFlowRate = cruiseFuelFlowRate;
        TypicalCruiseSpeed = typicalCruiseSpeed;
        NumberOfEngines = numberOfEngines;
        NonVolatileNumberEmissionsIndex = nonVolatileNumberEmissionsIndex;
        PayloadMass = maximumPayloadMass * flightPayloadFraction;
        EmptyMass = MaximumTakeoffMass - (maximumPayloadMass + maximumFuelMass);
        // Equation A2 of Unterstrasser (2016)
        InitialVortexSeparation = Math.PI * Wingspan / 4.0;
        WaterVapourEmissionsIndex = waterVapourEmissionsIndex;
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

    private double FuelFlowRate()
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

    private double FuelFlowRate()
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