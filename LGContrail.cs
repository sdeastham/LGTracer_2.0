using System.Reflection.PortableExecutable;
using System.Runtime.CompilerServices;
using MathNet.Numerics;

namespace LGTracer;

public class LGContrail : LGPointConnected
{
    // Explicitly-controlled valuables
    private readonly bool IncludeCompression;
    private double CrystalRadius; // Meters
    private double CrystalCount; // Crystals per meter
    private double CrossSectionArea; // Meters squared
    private double UpperWidth; // Meters
    private double Depth; // Meters, vertical
    private double LowerWidth; // Meters
    private double SkewOffset; // Meters
    public double AmbientSpecificHumidity; // Fraction
    public double Temperature;
    private double WaterVapourMass; // Kilograms per meter
    // Derived properties for the plume itself
    public double SpecificHumidity => WaterVapourMass / AirMass;
    public double RelativeHumidityLiquid => Physics.RelativeHumidityLiquid(Temperature, Pressure, SpecificHumidity);
    public double RelativeHumidityIce => Physics.RelativeHumidityIce(Temperature, Pressure, SpecificHumidity);
    private double AirDensity => Pressure * 28.97e-3 / (Physics.RGasUniversal * Temperature); // kg/m3
    private double AirMass => AirDensity * CrossSectionArea; // Kilograms per meter
    private double IceMass => CrystalCount * CrystalRadius * CrystalRadius * CrystalRadius * Math.PI * (4.0/3.0);
    // Derived plume quantities relating to the total amount of water available
    private double TotalWaterMass => IceMass + WaterVapourMass;
    public double TotalRelativeHumidityIce => Physics.RelativeHumidityIce(Temperature, Pressure, TotalWaterMass/AirMass);
    // Derived ambient air properties
    public double AmbientRelativeHumidityLiquid => Physics.RelativeHumidityLiquid(Temperature, Pressure, AmbientSpecificHumidity);
    public double AmbientRelativeHumidityIce => Physics.RelativeHumidityIce(Temperature, Pressure, AmbientSpecificHumidity);
    // Derived geometric quantities
    private double Circumference => CalculateCircumference();
    // Diagnostic quantities only
    private double LastSettlingVelocity;
    // Constants
    private const double GammaRatio = 0.4 / 1.4;
    private double WaterVapourEmissionsIndex = 1.223; // kg H2O per kg fuel
    private double LowerHeatingValue = 43.2e6; // J/kg
    private bool SkipNewtonIterationForTlm;
    private bool UsePonaterTlc;
    
    public Func<double, double, double, (double, double, double)> VelocityCalcNoSettling { get; protected set; }

    public LGContrail(Func<double, double, double, (double, double, double)> vCalc, bool includeCompression,
        bool includeSettling, double minimumPointLifetime=0.0, bool skipNewtonIterationForTlm=false,
        bool usePonaterTlc=false) : base(vCalc, minimumPointLifetime)
    {
        // Override vCalc to allow for inclusion of a settling speed
        VelocityCalcNoSettling = vCalc;
        if (includeSettling)
        {
            VelocityCalc = VelocityCalcWithSettling;
        }
        else
        {
            VelocityCalc = VelocityCalcNoSettling;
        }
        // Do we use newton iteration to evaluate threshold temperature?
        SkipNewtonIterationForTlm = skipNewtonIterationForTlm;
        // Do we use the Ponater approach for SAC (skipping critical temperature calculation)?
        UsePonaterTlc = usePonaterTlc;
        IncludeCompression = includeCompression;
        ZeroContrail();
        // Contrail points are NOT valid by default
        DefaultValidity = false;
    }

    private (double, double, double) VelocityCalcWithSettling(double x, double y, double pressure)
    {
        (double dxdt, double dydt, double dpdt) = VelocityCalcNoSettling(x, y, pressure);
        if (!ContrailLive()) { return (dxdt,dydt,dpdt); }
        // Calculate settling velocity in m/s
        double dzdtSettling = CalculateSettlingVelocity(CrystalRadius, CalculateDynamicViscosity(Temperature));
        // Store for diagnostics
        LastSettlingVelocity = dzdtSettling;
        // Convert to Pa/s - use hydrostatic assumption and assume changes are small during the given time period (!)
        // dp/dz = -rho * g
        // dp/dt = dz/dt * (-rho * g)
        double dpdtSettling = dzdtSettling * Physics.G0 * -1.0 * AirDensity;
        return (dxdt, dydt, dpdt + dpdtSettling);
    }
    
    public void InitiateContrail(double efficiency, double numberEmissionsIndex, double activationFraction,
        double airSpeed, double fuelFlowRate, double waterVapourEmissionsIndex=1.223, double lowerHeatingValue = 43.2e6)
    {
        // Set up the contrail key values
        LowerHeatingValue = lowerHeatingValue;
        WaterVapourEmissionsIndex = waterVapourEmissionsIndex;
        // Do we pass the Schmidt-Appleman Criterion?
        bool contrailFormed;
        if (UsePonaterTlc)
        {
            // Positive means we are above (i.e. humid enough)
            double criticalRelativeHumidityDelta = CriticalRelativeHumidityDelta(efficiency);
            contrailFormed = criticalRelativeHumidityDelta >= 0.0;
        }
        else
        {
            // Negative means that we are below (i.e. cold enough)
            double criticalTemperatureDelta = CriticalTemperatureDelta(efficiency);
            contrailFormed = criticalTemperatureDelta < 0.0;
        }
        if (!contrailFormed)
        {
            // Local temperature is above critical
            ZeroContrail();
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
        double fuelWaterVapourMass = fuelPerMeter * waterVapourEmissionsIndex;
        double ambientWaterVapourMass = AmbientSpecificHumidity * AirMass;
        WaterVapourMass = fuelWaterVapourMass + ambientWaterVapourMass;
        // Subtract ice already on the crystals
        WaterVapourMass -= IceMass;
        LastSettlingVelocity = 0.0;
    }

    private double CalculateCircumference()
    {
        double depthSquared = Depth * Depth;
        double rightBase = LowerWidth + SkewOffset - UpperWidth;
        return LowerWidth + UpperWidth + Math.Sqrt(depthSquared + SkewOffset * SkewOffset) +
               Math.Sqrt(depthSquared + rightBase * rightBase);
    }
    private static double CalculateDynamicViscosity(double temperature)
    {
        // Sutherland (1893), "The viscosity of gases and molecular force"
        // Units are Pa.s
        return 1.458e-6 * Math.Pow(temperature, 1.5) / (temperature + 110.4);
    }
    
    private static double CalculateSettlingVelocity(double particleRadius, double dynamicViscosity)
    {
        // Stokes law for terminal velocity.
        // Valid for r < 0.03 mm (Lohmann et al 2016)
        // particleRadius in m
        // dynamic viscosity in Pa.s
        // Value returned is in meters per second
        const double particleDensity = 1000.0; // kg/m3
        const double preFactor = ( 2.0 / 9.0 ) * particleDensity * Physics.G0;
        return preFactor * particleRadius * particleRadius / dynamicViscosity;
    }

    private void ZeroContrail()
    {
        // Resets the contrail properties and prevents contrail calculations from occuring until correctly initialized
        WaterVapourMass = double.NaN;
        Depth = double.NaN;
        CrossSectionArea = double.NaN;
        CrystalCount = double.NaN;
        CrystalRadius = double.NaN;
        LastSettlingVelocity = double.NaN;
    }

    private double CriticalTemperatureDelta(double efficiency)
    {
        return Temperature - CalculateCriticalTemperature(Pressure, efficiency, AmbientRelativeHumidityLiquid,
            WaterVapourEmissionsIndex,LowerHeatingValue,SkipNewtonIterationForTlm);
    }

    private double CriticalRelativeHumidityDelta(double efficiency)
    {
        return AmbientRelativeHumidityLiquid - CalculateCriticalRelativeHumidityLiquid(Pressure, Temperature, efficiency,
            WaterVapourEmissionsIndex,LowerHeatingValue,SkipNewtonIterationForTlm);
    }

    private static double CalculateCriticalTemperature(double pressure, double efficiency, double relativeHumidity,
        double waterVapourEmissionsIndex, double lowerHeatingValue, bool approximateTlm=false)
    {
        double mixingLineGradient = MixingLineGradient(pressure, efficiency, waterVapourEmissionsIndex, lowerHeatingValue);
        double thresholdTemperature = approximateTlm ? EstimateLiquidThresholdTemperature(mixingLineGradient) : NewtonIterTlm(mixingLineGradient);
        double criticalTemperature = NewtonIterTlc(mixingLineGradient, thresholdTemperature, relativeHumidity);
        return criticalTemperature;
    }
    
    private static double CalculateCriticalRelativeHumidityLiquid(double pressure, double temperature, double efficiency,
        double waterVapourEmissionsIndex, double lowerHeatingValue, bool approximateTlm=false)
    {
        // G
        double mixingLineGradient = MixingLineGradient(pressure, efficiency, waterVapourEmissionsIndex, lowerHeatingValue);
        // Tlm
        double thresholdTemperature = approximateTlm ? EstimateLiquidThresholdTemperature(mixingLineGradient) : NewtonIterTlm(mixingLineGradient);
        double pSatTlm = Physics.SaturationPressureLiquid(thresholdTemperature);
        double pSat = Physics.SaturationPressureLiquid(temperature);
        // Equation 12 from Schumann (2012)
        return (mixingLineGradient * (temperature - thresholdTemperature) + pSatTlm) / pSat;
    }

    private static double MixingLineGradient(double pressure, double efficiency, double waterVapourEmissionsIndex=1.223,
        double lowerHeatingValue = 43.2e6)
    {
        /*
         * pressure         Pressure in Pa
         * efficiency       Between 0 and 1 (typically 0.3 - 0.4)
         * waterVapour..    kg H2O per kg fuel (1.25 for kerosene)
         * lowerHeating...  J per kg (43.0e3 for kerosene)
         Value returned is in units of partial pressure of H2O (Pa) per K. Default values for kerosene above taken from
         Schumann (1996) Table 1. Defaults for cp and the molar mass ratio likewise.
         */
        const double cpExhaust = Physics.CpAir; // J kg-1 K-1
        const double molarMassRatio = 1.0/Physics.WaterMolarConversion; // Ratio between molar masses of water vapor (~18e-3) and air (~29e-3) (0.622)
        return waterVapourEmissionsIndex * pressure * cpExhaust /
               (molarMassRatio * lowerHeatingValue * (1.0 - efficiency));
    }

    private static double EstimateLiquidThresholdTemperature(double mixingLineGradient)
    {
        // Mixing line gradient in Pa/K
        // Value returned in K
        // Approximations below from Gierens (2021), "Theory of Contrail Formation for Fuel Cells"
        if (mixingLineGradient <= 2.0)
        {
            double logG = Math.Log(mixingLineGradient - 0.053);
            return 226.69 + 9.43 * logG + 0.720 * logG * logG;
        }
        else
        {
            double logG = Math.Log(mixingLineGradient);
            return 226.031 + 10.2249 * logG + 0.335372 * logG * logG + 0.0642105 * logG * logG * logG;
        }
    }
    
    private static double EstimateCriticalTemperature(double thresholdTemperature, double mixingLineGradient,
        double relativeHumidity)
    {
        // Notation here is from Schumann (1996). Implements the estimation approach in Appendix 2 to provide a first
        // guess for subsequent Newton-Raphson iteration. Updated to follow Schumann (2012).
        // * thresholdTemperature     T_LM in Schumann (1996), T_Max in Gierens (2021) [K]
        // * mixingLineGradient       G, in Pa/K
        // * relativeHumidity         RH with respect to liquid as a fraction (0-1)
        double dT = 1.0e-4;
        double pSat = Physics.SaturationPressureLiquid(thresholdTemperature);
        if (relativeHumidity <= 0.0)
        {
            return thresholdTemperature - (pSat / mixingLineGradient);
        }
        double pSatFraction = pSat * relativeHumidity;
        double d2PSatByD2T = (Physics.SaturationPressureLiquid(thresholdTemperature + dT) - 2.0 * pSat +
                              Physics.SaturationPressureLiquid(thresholdTemperature - dT)) / (dT * dT);
        // This is a 2nd-order Taylor expansion around T_LM to find T_LC
        // Updated based on Schumann (2012), which states that U**2 in equation 34 (for commonFactor) should be U
        double commonFactor = 1.0 / (relativeHumidity * d2PSatByD2T);
        double aFactor = (1.0 - relativeHumidity) * mixingLineGradient * commonFactor;
        double bFactor = (pSat - pSatFraction) * commonFactor;
        // offsetTemperature is x in Schumann (1996); just the difference between T_LC and T_LM
        double offsetTemperature = (-1.0 * aFactor) + Math.Sqrt(aFactor * aFactor + 2.0 * bFactor);
        return thresholdTemperature - offsetTemperature;
    }

    private static double EstimateDPSatIceByDT(double temperature, double dT = 2.0e-10)
    {
        return (Physics.SaturationPressureIce(temperature + dT / 2.0) -
                Physics.SaturationPressureIce(temperature - dT / 2.0)) / dT;
    }
    
    private static double EstimateDPSatByDT(double temperature, double dT = 2.0e-10)
    {
        return (Physics.SaturationPressureLiquid(temperature + dT / 2.0) -
                Physics.SaturationPressureLiquid(temperature - dT / 2.0)) / dT;
    }

    private static double EstimateD2PSatByDTDTlm(double temperature, double dT = 2.0e-2)
    {
        return (EstimateDPSatByDT(temperature + dT / 2.0) - EstimateDPSatByDT(temperature - dT / 2.0)) / dT;
    }

    private static double NewtonIterTlm(double mixingLineGradient, double? initialGuess = null)
    {
        // Use Newton-Raphson iteration to refine the estimated temperature
        // Maximum allowable difference between G and the saturation pressure
        const double errorThreshold = 1.0e-5;
        double temperatureThreshold;
        if (initialGuess == null)
        {
            temperatureThreshold = EstimateLiquidThresholdTemperature(mixingLineGradient);
        }
        else
        {
            temperatureThreshold = (double)initialGuess;
        }
        double dPSatByDt = EstimateDPSatByDT(temperatureThreshold);
        double signedError = dPSatByDt - mixingLineGradient;
        int nAttempts = 0;
        const int nMax = 500;
        while (Math.Abs(signedError) > errorThreshold)
        {
            if (nAttempts >= nMax)
            {
                throw new NonConvergenceException("Failed to converge on value for threshold temperature");
            }
            // Finding the roots of the function f = (dp/dt - G) = f(T)
            double error = dPSatByDt - mixingLineGradient;
            double errorDerivative = EstimateD2PSatByDTDTlm(temperatureThreshold);
            // Refine the estimate of T
            temperatureThreshold -= error / errorDerivative;
            // Update the error estimate
            dPSatByDt = EstimateDPSatByDT(temperatureThreshold);
            signedError = dPSatByDt - mixingLineGradient;
            nAttempts++;
        }
        return temperatureThreshold;
    }

    private static double NewtonIterTlc(double mixingLineGradient, double thresholdTemperature, double relativeHumidity, double? initialGuess=null)
    {
        const double errorThreshold = 1.0e-5;
        double criticalTemperature;
        if (initialGuess == null)
        {
            criticalTemperature =
                EstimateCriticalTemperature(thresholdTemperature, mixingLineGradient, relativeHumidity);
        }
        else
        {
            criticalTemperature = (double)initialGuess;
        }

        if (relativeHumidity <= 0.0 || relativeHumidity >= 1.0)
        {
            return criticalTemperature;
        }
        double pSat = Physics.SaturationPressureLiquid(thresholdTemperature);
        double updatedTemperature = thresholdTemperature -
                                    (pSat - relativeHumidity * Physics.SaturationPressureLiquid(criticalTemperature)) /
                                    mixingLineGradient;
        int nAttempts = 0;
        int nMax = 150;
        while (Math.Abs(criticalTemperature - updatedTemperature) > errorThreshold)
        {
            if (nAttempts > nMax)
            {
                throw new NonConvergenceException("Failed to converge on value for critical temperature");
            }
            double f = criticalTemperature - updatedTemperature;
            double f_prime = 1.0 + (relativeHumidity / mixingLineGradient) * EstimateDPSatByDT(criticalTemperature);
            criticalTemperature -= f / f_prime;
            updatedTemperature = thresholdTemperature -
                                 (pSat - relativeHumidity * Physics.SaturationPressureLiquid(criticalTemperature)) /
                                 mixingLineGradient;
            nAttempts++;
        }
        return criticalTemperature;
    }

    public static bool TestSAC(bool iterate=true)
    {
        // Verify that the SAC routines are doing what they should be comparison to Schumann (1996) known results.
        // NB: Passing iterate=true means that N-R iteration will be used to refine guesses. This seems to
        double pressure = 220.0;
        double temperature = -59.0;
        double rh = 42.0;
        double efficiency = 0.308;
        double maxError = 0.01; // Generous because only 3 sig figs in Schumann (1996)
        
        bool success = true;
        success = success && Evaluate("Kerosene", 43.0, 1.25, 1.49, -50.6, -42.9,
            pressure, temperature, rh, efficiency,maxError,iterate);
        success = success && Evaluate("Methane", 50.0, 2.24, 2.31, -46.3, -38.3,
            pressure, temperature, rh, efficiency,maxError,iterate);
        success = success && Evaluate("Hydrogen", 120.0, 8.94, 3.82, -41.2, -32.7,
            pressure, temperature, rh, efficiency,maxError,iterate);

        // Verify that both approaches to checking SAC give the same result for a sweep which spans from
        // SAC = True to SAC = False
        double tMin = -52.0;
        double tMax = -49.0;
        double dt = 0.02;
        int nChecks = (int)Math.Ceiling((tMax - tMin) / dt);
        double tCurr = tMin;
        bool agreedTlcUlc = true;
        while (tCurr <= tMax)
        {
            agreedTlcUlc = agreedTlcUlc && CompareTlcUlc(43.0, 1.25, pressure, tCurr, rh, 
                efficiency, iterate, tCurr <= tMin + dt/2.0 || tCurr >= tMax - dt);
            tCurr += dt;
        }
        Console.WriteLine($"Agreement between SAC calculation approaches: {agreedTlcUlc}");
        return success && agreedTlcUlc;
    }
    
    public static bool Evaluate(string name, double lhv, double eiH2O, double refG, double refTLC, double refTLM,
        double ambientPressure, double ambientTemperature, double ambientRH, double efficiency,
        double maxError = 0.1, bool iterate=true)
    {
        // Input units follow Schumann (1996), i.e. LHV in MJ/kg, temperatures in C
        double temperature = ambientTemperature + 273.15; // Convert C to K
        double pressure = ambientPressure * 100.0; // Convert hPa to Pa 
        double rh = ambientRH * 0.01; // Convert % to fraction
        double mixingLineGradient = MixingLineGradient(pressure, efficiency, eiH2O, lhv*1.0e6);
        double pSatAmbient = Physics.SaturationPressureLiquid(temperature);
        double thresholdTemperature = EstimateLiquidThresholdTemperature(mixingLineGradient);
        if (iterate)
        {
            thresholdTemperature = NewtonIterTlm(mixingLineGradient, thresholdTemperature);
        }
        double criticalTemperature = EstimateCriticalTemperature(thresholdTemperature, mixingLineGradient, rh);
        if (iterate)
        {
            criticalTemperature =
                NewtonIterTlc(mixingLineGradient, thresholdTemperature, rh, criticalTemperature);
        }
        double criticalRelativeHumidity = CalculateCriticalRelativeHumidityLiquid(pressure, temperature, efficiency,
            eiH2O, lhv*1.0e6, approximateTlm: iterate);
        double mixingLineError = mixingLineGradient / refG - 1.0;
        double tlmError = thresholdTemperature / (273.15+refTLM) - 1.0;
        double tlcError = criticalTemperature / (273.15+refTLC) - 1.0;
        Console.WriteLine($"Evaluation for {name} (iteration: {iterate}) [T={ambientTemperature,9:f3} C, RHl={ambientRH,9:f3}%]:");
        Console.WriteLine($" --> Mixing line gradient G     : {mixingLineGradient,9:f3} vs {refG,9:f3} Pa/K, {100.0*mixingLineError,9:f4}% error");
        Console.WriteLine($" --> Threshold temperature T_LM : {thresholdTemperature-273.15,9:f3} vs {refTLM,9:f3} C,    {100.0*tlmError,9:f4}% error");
        Console.WriteLine($" --> Critical temperature T_LC  : {criticalTemperature-273.15,9:f3} vs {refTLC,9:f3} C,    {100.0*tlcError,9:f4}% error");
        Console.WriteLine($" --> Critical RH U_LC (Ponater) : {criticalRelativeHumidity*100.0,9:f3}%");
        Console.WriteLine($" -----> Formation (Schumann/Ponater): {temperature<=criticalTemperature}/{rh>=criticalRelativeHumidity}");
        return (mixingLineError < maxError && tlmError < maxError && tlcError < maxError);
    }
    
    public static bool CompareTlcUlc(double lhv, double eiH2O,
        double ambientPressure, double ambientTemperature, double ambientRH, double efficiency, bool iterate=true, bool verbose=true)
    {
        // Input units follow Schumann (1996), i.e. LHV in MJ/kg, temperatures in C
        double temperature = ambientTemperature + 273.15; // Convert C to K
        double pressure = ambientPressure * 100.0; // Convert hPa to Pa 
        double rh = ambientRH * 0.01; // Convert % to fraction
        double mixingLineGradient = MixingLineGradient(pressure, efficiency, eiH2O, lhv*1.0e6);
        double pSatAmbient = Physics.SaturationPressureLiquid(temperature);
        double thresholdTemperature = EstimateLiquidThresholdTemperature(mixingLineGradient);
        if (iterate)
        {
            thresholdTemperature = NewtonIterTlm(mixingLineGradient, thresholdTemperature);
        }
        double criticalTemperature = EstimateCriticalTemperature(thresholdTemperature, mixingLineGradient, rh);
        if (iterate)
        {
            criticalTemperature =
                NewtonIterTlc(mixingLineGradient, thresholdTemperature, rh, criticalTemperature);
        }
        double criticalRelativeHumidity = CalculateCriticalRelativeHumidityLiquid(pressure, temperature, efficiency,
            eiH2O, lhv*1.0e6, approximateTlm: iterate);
        bool tlcCheck = temperature <= criticalTemperature;
        bool rhCheck = rh >= criticalRelativeHumidity;
        if (!verbose) { return tlcCheck == rhCheck; }
        Console.WriteLine($"P/S comparison (iteration: {iterate}) [T={ambientTemperature,9:f3} C, RHl={ambientRH,9:f3}%]:");
        Console.WriteLine($" --> Mixing line gradient G     : {mixingLineGradient,9:f3} Pa/K");
        Console.WriteLine($" --> Threshold temperature T_LM : {thresholdTemperature-273.15,9:f3} C");
        Console.WriteLine($" --> Critical temperature T_LC  : {criticalTemperature-273.15,9:f3} C");
        Console.WriteLine($" --> Critical RH U_LC (Ponater) : {criticalRelativeHumidity*100.0,9:f3}%");
        Console.WriteLine($" -----> Formation (T_LC/U_LC): {tlcCheck}/{rhCheck}");
        return tlcCheck == rhCheck;
    }

    public override void Deactivate()
    {
        IsLeader = false;
        ZeroContrail();
        base.Deactivate();
    }

    public override void Advance(double dt, DomainManager domain)
    {
        if (!Active) return;
        double oldPressure = Pressure;
        bool contrailActive = (Segment != null) && ContrailLive();
        double oldLength = 1.0;
        double oldTemperature = 1.0;
        double oldIceMass = IceMass;
        double oldAmbientSpecificHumidity = AmbientSpecificHumidity;
        if (contrailActive)
        {
            oldTemperature = Temperature;
            oldLength = Segment.SegmentLength;
        }
        else if (!double.IsNaN(CrystalCount))
        {
            // Only called if contrail is not active, but hasn't been zeroed
            // Really we want this to happen when the segment is nullified - need a callback from the segment on death
            ZeroContrail();
        }
        
        // Perform the standard advection calculation
        base.Advance(dt, domain);
        
        // Update temperature based on adiabatic compression
        if (IncludeCompression)
        {
            Temperature *= Math.Pow(Pressure/oldPressure,GammaRatio);
        }
        if (!contrailActive) return;
        
        // Apply expansion and stretching to the contrail
        // Old air density (in mol/m3) divided by new air density
        double densityRatio = (oldPressure * Temperature)/(Pressure * oldTemperature);
        // Need to conserve mass prior to diffusion calculation
        double previousVolume = CrossSectionArea * oldLength;
        double newVolume = previousVolume * densityRatio;
        double oldArea = CrossSectionArea;
        CrossSectionArea = newVolume / Segment.SegmentLength;

        SimpleContrailDiffusionModel(dt, oldAmbientSpecificHumidity);
    }

    private readonly double SimpleGrowthConstant = 3600.0 / Math.Log(1.0 + 0.10); 
    private void SimpleContrailDiffusionModel(double dt, double oldAmbientSpecificHumidity)
    {
        // Mix in new air at the ambient humidity of the target location
        // Very, very simple mixing; assume a 10% increase in air mass per hour, and that the additional air is ambient
        double growthFactor = Math.Exp(dt / SimpleGrowthConstant); // Factor increase in air mass due to "mixing"
        // How much additional water will be brought in?
        double meanAmbientSpecificHumidity = 0.5 * (AmbientSpecificHumidity + oldAmbientSpecificHumidity);
        double newWater = meanAmbientSpecificHumidity * (growthFactor - 1.0) * AirMass;
        WaterVapourMass += newWater;
        CrossSectionArea *= growthFactor;
            
        // Calculate what the relative humidity with respect to ice would be if all water mass was vapour. If this is
        // less than 1.0 then the contrail cannot be sustained
        double availableIce = (1.0 - (1.0 / TotalRelativeHumidityIce)) * TotalWaterMass;
        if (availableIce < 0.0)
        {
            ZeroContrail();
            return;
        }
        // Grow/shrink crystals accordingly - assumes monodisperse
        CrystalRadius = Math.Cbrt(0.75 * availableIce / (Math.PI * CrystalCount));
    }

    public bool ContrailLive()
    {
        return CrystalCount > 1.0e-10 || !double.IsNaN(CrossSectionArea);
    }

    public override bool CheckValid()
    {
        // Am I base valid? If not, cull!
        // Am I the leader such that I might end up as the tail for a contrail?
        // Do I have a contrail?
        // If not; am I the tail of another contrail?
        return base.CheckValid() || (IsLeader || ContrailLive() || (Next != null && ((LGContrail)Next).ContrailLive()));
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
                case "settlingspeed":
                case "settlingvelocity":
                    return LastSettlingVelocity;
                default:
                    throw new ArgumentException($"LGContrail does not contain property {property}");
            }
        }
        else
        {
            return base.GetProperty(property);
        }
    }
    
    // Functions for Unterstrasser 2016 (https://acp.copernicus.org/articles/16/2059/2016/)
    private static double VortexDisplacementLength(double circulation, double bruntVaisalaFrequency)
    {
        // Returns z_desc
        return Math.Sqrt(circulation * 8.0 / (Math.PI * bruntVaisalaFrequency));
    }
    
    private static double AtmosphericSaturationLength(double temperature, double relativeHumidityIce, double dryLapseRate=9.8e-3)
    {
        // Temperature in K, relativeHumidityIce in fraction
        // Dry lapse rate is in K/m
        // Returns z_atm
        double errorThreshold = 1.0e-6;
        double leftHandSide = relativeHumidityIce * Physics.SaturationPressureIce(temperature) / temperature;
        double zAtm = 1000.0 * Double.Max(0.0, relativeHumidityIce - 1.0);
        double temperatureOffset = temperature + dryLapseRate * zAtm;
        double pSatOffset = Physics.SaturationPressureIce(temperatureOffset);
        double f = leftHandSide - (pSatOffset / temperatureOffset);
        int nAttempts = 0;
        int nMax = 100;
        while (Math.Abs(f) > errorThreshold)
        {
            if (nAttempts > nMax)
            {
                throw new NonConvergenceException("Failed to converge on value for atmospheric length scale in contrail formation.");
            }
            double fPrime = (dryLapseRate / (temperatureOffset * temperatureOffset)) *
                             (pSatOffset - temperatureOffset * EstimateDPSatIceByDT(temperatureOffset));
            zAtm -= f / fPrime;
            // Update derived quantities
            temperatureOffset = temperature + dryLapseRate * zAtm;
            pSatOffset = Physics.SaturationPressureIce(temperatureOffset);
            f = leftHandSide - (pSatOffset / temperatureOffset);
            nAttempts++;
        }
        return zAtm;
    }

    private const double WaterGasConstant = 461.0; // J/kg/k for water vapour
    private static double PlumeSaturationLength(double temperature, double plumeArea, double waterEmissionRate, double dryLapseRate=9.8e-3)
    {
        // Temperature in K, wing span in m, water emission rate in kg/m
        // Dry lapse rate is in K/m
        // Returns z_emit
        double errorThreshold = 1.0e-8;
        double waterPerUnitArea = waterEmissionRate / plumeArea;
        double leftHandSide = (Physics.SaturationPressureIce(temperature) / (temperature * WaterGasConstant)) +
                              waterPerUnitArea;
        double zEmit = 100.0;
        double temperatureOffset = temperature + dryLapseRate * zEmit;
        double pSatOffset = Physics.SaturationPressureIce(temperatureOffset);
        double f = leftHandSide - (pSatOffset / (WaterGasConstant * temperatureOffset));
        int nAttempts = 0;
        int nMax = 100;
        while (Math.Abs(f) > errorThreshold)
        {
            if (nAttempts > nMax)
            {
                throw new NonConvergenceException("Failed to converge on value for plume water length scale in contrail formation.");
            }
            double fPrime = (dryLapseRate / (WaterGasConstant * temperatureOffset * temperatureOffset)) *
                            (pSatOffset - temperatureOffset * EstimateDPSatIceByDT(temperatureOffset));
            zEmit -= f / fPrime;
            // Update derived quantities
            temperatureOffset = temperature + dryLapseRate * zEmit;
            pSatOffset = Physics.SaturationPressureIce(temperatureOffset);
            f = leftHandSide - (pSatOffset / (WaterGasConstant * temperatureOffset));
            nAttempts++;
        }
        return zEmit;
    }
}