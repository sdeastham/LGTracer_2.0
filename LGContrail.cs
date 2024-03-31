namespace LGTracer;

public class LGContrail : LGPointConnected
{
    // Trapezoidal model for the contrail?
    private readonly bool IncludeCompression;
    private double CrystalRadius; // Meters
    private double CrystalCount; // Crystals per meter
    private double CrossSectionArea; // Meters squared
    private double Depth; // Meters
    public double SpecificHumidity => WaterVapourMass / AirMass;
    public double RelativeHumidityLiquid => Physics.RelativeHumidityLiquid(Temperature, Pressure, SpecificHumidity);
    public double RelativeHumidityIce => Physics.RelativeHumidityIce(Temperature, Pressure, SpecificHumidity);
    private double AirMass => AirDensity * CrossSectionArea; // Kilograms per meter
    private double WaterVapourMass; // Kilograms per meter
    private double IceMass => CrystalCount * CrystalRadius * CrystalRadius * CrystalRadius * Math.PI * (4.0/3.0);
    private double Temperature;
    private double AirDensity => Pressure * 28.97e-3 / (Physics.RGasUniversal * Temperature); // kg/m3

    private const double GammaRatio = 0.4 / 1.4;

    public LGContrail(Func<double, double, double, (double, double, double)> vCalc, bool includeCompression) :
        base(vCalc)
    {
        IncludeCompression = includeCompression;
        WaterVapourMass = 0.0;
        Depth = 0.0;
        CrossSectionArea = 0.0;
        CrystalCount = 1.0;
        CrystalRadius = 1.0;
    }

    public void InitiateContrail(double efficiency, double numberEmissionsIndex, double activationFraction,
        double airSpeed, double fuelFlowRate, double waterVapourEmissionsIndex=1.223, double lowerHeatingValue = 43.2e3)
    {
        if (!CheckSAC(efficiency))
        {
            return;
        }
        // Fuel burned per meter travelled
        double fuelPerMeter = fuelFlowRate / airSpeed;
        // Particles emitted per meter travelled
        double particleCount = numberEmissionsIndex * fuelPerMeter;
        // For now, just assume some activation and survival fraction
        CrystalCount = particleCount * activationFraction;
        CrystalRadius = 1.0e-9; // Just assume 1 nm to begin with. Assume soot is negligible
        CrossSectionArea = 100.0; // Assume crystals are spread over 100 m2 to begin with
        WaterVapourMass = fuelPerMeter * waterVapourEmissionsIndex;
    }

    private bool CheckSAC(double efficiency)
    {
        return SchmidtApplemanCriterion(Pressure, Temperature, efficiency);
    }

    private static bool SchmidtApplemanCriterion(double pressure, double temperature, double efficiency)
    {
        throw new NotImplementedException("SAC calculation not yet implemented");
    }

    private static double MixingLineGradient(double pressure, double efficiency, double waterVapourEmissionsIndex=1.223,
        double lowerHeatingValue = 43.2e3)
    {
        /*
         * pressure         Pressure in Pa
         * efficiency       Between 0 and 1 (typically 0.3 - 0.4)
         * waterVapour..    kg H2O per kg fuel (1.223 for kerosene)
         * lowerHeating...  J per kg (43.2e3 for kerosene)
         Value returned is in units of partial pressure of H2O (Pa) per K
         */
        const double cpExhaust = 1004; // J kg-1 K-1
        const double molarMassRatio = 0.622; // Ratio between molar masses of water vapor (~18e-3) and air (~29e-3)
        return waterVapourEmissionsIndex * pressure * cpExhaust /
               (molarMassRatio * lowerHeatingValue * (1.0 - efficiency));
    }
    
    public override void Advance(double dt, DomainManager domain)
    {
        if (!Active) return;
        double oldPressure = Pressure;
        bool contrailActive = (Segment != null) && CheckValid();
        double oldLength = 1.0;
        double oldTemperature = 1.0;
        if (contrailActive)
        {
            oldTemperature = Temperature;
            oldLength = Segment.SegmentLength;
        }
        base.Advance(dt, domain);
        // Update temperature based on adiabatic compression
        if (IncludeCompression)
        {
            Temperature *= Math.Pow(Pressure/oldPressure,GammaRatio);
        }
        if (!contrailActive) return;
        // Old air density (in mol/m3) divided by new air density
        double densityRatio = (oldPressure * Temperature)/(Pressure * oldTemperature);
        // For now, assume that the change in air density manifests only as a change in area
        CrossSectionArea *= densityRatio;
        // Apply the effect of segment stretching
        CrossSectionArea *= oldLength / Segment.SegmentLength;
        // TODO: Incorporate simple diffusion and mixing
        // TODO: Incorporate ice crystal microphysics
        // TODO: Incorporate settling
    }

    public override bool CheckValid()
    {
        return CrystalCount > 1.0e-10;
    }

    public override double GetProperty(string property)
    {
        if (property.StartsWith("Contrail", StringComparison.CurrentCultureIgnoreCase))
        {
            switch (property.Substring(8).ToLower().Replace("_", ""))
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
                case "crosssectionarea":
                case "xsa":
                    return CrossSectionArea;
                case "crystalcount":
                    return CrystalCount;
                case "crystalnumberdensity":
                    return CrystalCount / CrossSectionArea;
                case "crystalradius":
                    return CrystalRadius;
                case "icemass":
                    return IceMass;
                default:
                    throw new ArgumentException($"LGContrail does not contain property {property}");
            }
        }
        else
        {
            return base.GetProperty(property);
        }
    }
}