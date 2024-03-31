namespace LGTracer;

public class LGContrail(
    Func<double, double, double, (double, double, double)> vCalc)
    : LGPointConnected(vCalc)
{
    // Trapezoidal model for the contrail?
    private double CrystalRadius;
    private double CrystalCount;
    private double CrossSectionArea;
    private double Depth;
    public double SpecificHumidity => WaterVapourMass / AirMass;
    private double AirMass; // Should be derived
    private double WaterVapourMass;
    private double IceMass;
    
    // Do we need a new pointmanager?
    // Need to make an LGPointAirMass (make LGPoint lighter weight)
    // Reverse of base?
}