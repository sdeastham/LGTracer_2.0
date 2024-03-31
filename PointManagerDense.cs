namespace LGTracer;

public class PointManagerDense : PointManager
{
    // A "dense" point manager is designed to represent ALL air within the domain
    public double MassSurplus;
    private Random Rng;
    private double KgPerPoint;

    public PointManagerDense( long? maxPoints, DomainManager domain, string filename, bool verboseOutput=false, 
        bool includeCompression=false, string[]? propertyNames=null, Random? rng=null,
        double kgPerPoint=1.0e12 ) : base(maxPoints,domain,filename,verboseOutput,includeCompression,propertyNames)
    {
        // No initial mass surplus
        MassSurplus = 0.0;
        KgPerPoint = kgPerPoint;
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

    public override LGPoint NextPoint(double x, double y, double pressure)
    {
        LGAirMass point = (LGAirMass)base.NextPoint(x, y, pressure);
        // Set point properties based on local values#
        double temperature = Domain.NearestNeighbor3D(x,y,pressure,Domain.TemperatureXYP);
        double specificHumidity = Domain.NearestNeighbor3D(x,y,pressure,Domain.SpecificHumidityXYP);
        point.Temperature = temperature;
        point.SpecificHumidity = specificHumidity;
        return point;
    }
        
    public override double GetPromotedProperty(IAdvected point, string property)
    {
        return ((LGAirMass)point).GetProperty(property);
    }
}