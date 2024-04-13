namespace LGTracer;

public class PointManagerDense : PointManager
{
    // A "dense" point manager is designed to represent ALL air within the domain
    protected double MassSurplus;
    protected readonly Random Rng;
    protected readonly double KgPerPoint;

    public PointManagerDense( DomainManager domain, LGOptions configOptions, LGOptionsPointsDense configSubOptions,
        Random rng) : base(domain, configOptions, configSubOptions)
    {
        // No initial mass surplus
        MassSurplus = 0.0;
        KgPerPoint = configSubOptions.KgPerPoint;
        Rng = rng;
    }
        
    protected override IAdvected CreatePoint()
    {
        return new LGAirMass(VelocityCalc,IncludeCompression);
    }
        
    public override void Seed(double dt)
    {
        (double[] xSet, double[] ySet, double[] pSet, MassSurplus) =
            Domain.SeedBoundary(KgPerPoint, dt, Rng, MassSurplus);
        CreatePointSet(xSet, ySet, pSet);

        (double[] xSetV, double[] ySetV, double[] pSetV, MassSurplus) =
            Domain.SeedPressureBoundaries(KgPerPoint, dt, Rng, MassSurplus);
        CreatePointSet(xSetV, ySetV, pSetV);
    }

    public override IAdvected NextPoint(double x, double y, double pressure)
    {
        LGAirMass point = (LGAirMass)base.NextPoint(x, y, pressure);
        // Set point properties based on local values#
        double temperature = Domain.NearestNeighbor3D(x,y,pressure,Domain.TemperatureXYP);
        double specificHumidity = Domain.NearestNeighbor3D(x,y,pressure,Domain.SpecificHumidityXYP);
        point.Temperature = temperature;
        point.SpecificHumidity = specificHumidity;
        return point;
    }
}