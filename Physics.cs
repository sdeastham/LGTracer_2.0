using AtmosTools;

namespace LGTracer;

public class Physics
{
    public const double EarthRadius = PhysConstants.EarthRadius * 1000.0; // m
    public const double Deg2Rad = Math.PI / 180.0;
    public const double Rad2Deg = 180.0 / Math.PI;
    public const double WaterMolarMass = 18.01528 * 1.0e-3; // kg/mol
    public const double AirMolarMass = 28.964 * 1.0e-3; // kg/mol
    public const double WaterMolarConversion = AirMolarMass / WaterMolarMass;
    public const double RGasUniversal = PhysConstants.RGasUniversal; // J/K/mol
    public const double G0 = PhysConstants.G0; // m/s2
    public const double Avogadro = PhysConstants.Avogadro; // molec/mol
    public const double CpAir = 1004.0; // J/kg/K
        
    public static double SaturationPressureLiquid(double temperature)
    {
       /*
        Code adapted from APCEMM (github.com/MIT-LAE/APCEMM)
        DESCRIPTION:
          Returns water liquid saturation pressure in Pascal.
          Source: Sonntag (1994)
        
        INPUT PARAMETERS:
          - double T :: temperature expressed in K

        OUTPUT PARAMETERS:
          - double :: H2O liquid saturation pressure in Pascal
        */

        return (100.0 * Math.Exp(-6096.9385 / temperature
                                 + 16.635794 - 0.02711193 * temperature
                                 + 1.673952e-5 * temperature * temperature
                                 + 2.433502 * Math.Log(temperature)));
    }

    public static double SaturationPressureIce(double temperature)
    {
       /*
        Code adapted from APCEMM (github.com/MIT-LAE/APCEMM)
        DESCRIPTION:
          Returns water solid saturation pressure in Pascal.
          Source: Sonntag (1990)

        INPUT PARAMETERS:
          - double T :: temperature expressed in K

          OUTPUT PARAMETERS:
          - double :: H2O solid saturation pressure in Pascal
        */
       
        return 100.0 * Math.Exp(- 6024.5282 / temperature 
                                + 24.7219 + 0.010613868 * temperature
                                - 1.3198825E-5 * temperature * temperature
                                - 0.49382577 * Math.Log( temperature ) );
    }
        
    public static double RelativeHumidityLiquid(double temperature, double pressure, double specificHumidity)
    {
        // Returns relative humidity with respect to liquid water (fractional)
        // Temperature in K, pressure in Pa, specific humidity in kg water/kg air
        return (pressure * specificHumidity * WaterMolarConversion / SaturationPressureLiquid(temperature));
    }
        
    public static double RelativeHumidityIce(double temperature, double pressure, double specificHumidity)
    {
        // Returns relative humidity with respect to ice (fractional)
        // Temperature in K, pressure in Pa, specific humidity in kg water/kg air
        return (pressure * specificHumidity * WaterMolarConversion / SaturationPressureIce(temperature));
    }
}