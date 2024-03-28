namespace LGTracer
{
    public class Physics
    {
        public static double SaturationPressureLiquid(double temperature)
        {
            /*
            Code adapted from APCEMM (github.com/MIT-LAE/APCEMM)
            DESCRIPTION:
             * Returns water liquid saturation pressure in Pascal.
             * Source: Sonntag (1994)

            INPUT PARAMETERS:
             * - double T :: temperature expressed in K
             *
             * OUTPUT PARAMETERS:
             * - double :: H2O liquid saturation pressure in Pascal
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
             * Returns water solid saturation pressure in Pascal.
             * Source: Sonntag (1990)

            INPUT PARAMETERS:
             * - double T :: temperature expressed in K
             *
             * OUTPUT PARAMETERS:
             * - double :: H2O solid saturation pressure in Pascal
            */
            return 100.0 * Math.Exp(- 6024.5282 / temperature 
                                    + 24.7219 + 0.010613868 * temperature
                                    - 1.3198825E-5 * temperature * temperature
                                    - 0.49382577 * Math.Log( temperature ) );
        }

        private const double WaterMolarMass = 18.01528 * 1.0e-3; // kg/mol
        private const double AirMolarMass = 28.964 * 1.0e-3; // kg/mol
        private const double WaterMolarConversion = AirMolarMass / WaterMolarMass;
        
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
}