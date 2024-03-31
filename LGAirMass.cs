namespace LGTracer;

public class LGAirMass : LGPoint
{
    // This might need to be made more flexible
    public double Temperature;
    public double SpecificHumidity;
        
    public double LiquidWaterContent // kg water per kg air
    { get; protected set; }
        
    public double IceWaterContent // kg ice per kg air
    { get; protected set; }

    public double TotalWaterContent => IceWaterContent + LiquidWaterContent + SpecificHumidity;

    public double RelativeHumidityLiquid => Physics.RelativeHumidityLiquid(Temperature, Pressure, SpecificHumidity);
        
    public double RelativeHumidityIce => Physics.RelativeHumidityIce(Temperature, Pressure, SpecificHumidity);
    
    private bool IncludeCompression;

    public const double GammaRatio = 0.4/1.4;

    public LGAirMass(Func<double, double, double, (double, double, double)> vCalc, bool includeCompression) : base(vCalc)
    {
        IncludeCompression = includeCompression;
    }
    
    public override void Deactivate()
    {
        base.Deactivate();
        // Properties
        Temperature = double.NaN;
        SpecificHumidity = double.NaN;
    }

    public override void Activate( double x, double y, double pressure, uint uniqueID )
    {
        base.Activate( x,y,pressure,uniqueID );
        Temperature = double.NaN;
        SpecificHumidity = double.NaN;
    }

    public override double GetProperty(string property)
    {
        switch (property.ToLower().Replace("_",""))
        {
            case "temperature":
                return Temperature;
            case "specifichumidity":
            case "qv":
                return SpecificHumidity;
            case "relativehumidityliquid":
            case "rhl":
                return RelativeHumidityLiquid;
            case "relativehumidityice":
            case "rhi":
                return RelativeHumidityIce;
            default:
                throw new ArgumentException($"No property for LGAirMass called {property}");
        }
    }

    public override void Advance(double dt, DomainManager domain)
    {
        if (!Active) return;
        double oldPressure = Pressure;
        base.Advance(dt, domain);
        // Update temperature based on adiabatic compression
        if (IncludeCompression)
        {
            Temperature *= Math.Pow(Pressure/oldPressure,GammaRatio);
        }
    }
}