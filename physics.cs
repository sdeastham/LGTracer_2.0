namespace LGTracer
{
    public class Physics
    {
        public static double SaturationPressureLiquid(double temperature)
        {
            /* DESCRIPTION:
             * Returns water liquid saturation pressure in Pascal.
             * Source: Sonntag (1994) */

            /* INPUT PARAMETERS:
             * - double T :: temperature expressed in K
             *
             * OUTPUT PARAMETERS:
             * - double :: H2O liquid saturation pressure in Pascal */

            return (100.0 * Math.Exp(-6096.9385 / temperature
                                     + 16.635794 - 0.02711193 * temperature
                                     + 1.673952e-5 * temperature * temperature
                                     + 2.433502 * Math.Log(temperature)));

        }

        private const double WaterMolarMass = 18.01528 * 1.0e-3; // kg/mol
        private const double AirMolarMass = 28.964 * 1.0e-3; // kg/mol
        private const double WaterMolarConversion = AirMolarMass / WaterMolarMass;
        
        public static double RelativeHumidityLiquid(double temperature, double pressure, double specificHumidity)
        {
            // Returns relative humidity (fractional)
            // Temperature in K, pressure in Pa, specific humidity in kg water/kg air
            return (pressure * specificHumidity * WaterMolarConversion / SaturationPressureLiquid(temperature));

        }
    }
}