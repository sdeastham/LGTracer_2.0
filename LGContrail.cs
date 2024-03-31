namespace LGTracer;

public class LGContrail : LGPointConnected
{
    // Trapezoidal model for the contrail?
    private readonly bool IncludeCompression;
    private double CrystalRadius;
    private double CrystalCount;
    private double CrossSectionArea;
    private double Depth;
    public double SpecificHumidity => WaterVapourMass / AirMass;
    private double AirMass; // Should be derived
    private double WaterVapourMass;
    private double IceMass;
    private double Temperature;

    private const double GammaRatio = 0.4 / 1.4;

    public LGContrail(Func<double, double, double, (double, double, double)> vCalc, bool includeCompression) :
        base(vCalc)
    {
        IncludeCompression = includeCompression;
        IceMass = 0.0;
        WaterVapourMass = 0.0;
        AirMass = 0.0;
        Depth = 0.0;
        CrossSectionArea = 0.0;
        CrystalCount = 1.0;
        CrystalRadius = 1.0;
    }
    
    public override void Advance(double dt, DomainManager domain)
    {
        if (!Active) return;
        double oldPressure = Pressure;
        base.Advance(dt, domain);
        // Update temperature based on adiabatic compression
        if (IncludeCompression)
        {
            Temperature = Temperature * Math.Pow(Pressure/oldPressure,GammaRatio);
        }
    }

    public override bool CheckValid()
    {
        return CrystalCount > 1.0e-10;
    }
}