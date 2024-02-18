using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.Integration;
using MathNet.Numerics.Random;
using MathNet.Numerics.Distributions;

namespace LGTracer
{
    public static class SolidRotation
    {
        private static (double, double) VelocitySolidBody( double x, double y, double omega )
        {
            (double radius, double theta) = RThetaFromXY(x,y);
            double dxdt = radius * omega * Trig.Sin(theta) * -1.0;
            double dydt = radius * omega * Trig.Cos(theta);
            return (dxdt, dydt);
        }
        private static (double, double) RThetaAnalytical( double xInitial, double yInitial, double omega, double tCurr )
        {
            // Analytical solution to solid body rotation
            ( double radiusInitial, double thetaInitial ) = RThetaFromXY(xInitial, yInitial);
            double theta = thetaInitial + (tCurr * omega);
            return (radiusInitial, theta);
        }

        private static (double, double, double, double) RThetaError( double xInitial, double yInitial, double omega, double tCurr, double x, double y)
        {
            // Calculate the error in radius and angle for a calculated point versus the analytical solution
            (double radius, double theta) = RThetaFromXY(x,y);
            (double radiusAnalytical, double thetaAnalytical) = RThetaAnalytical(xInitial,yInitial,omega,tCurr);
            double rErrPcg,thetaErrDeg;
            const double radiusThreshold = 1.0e-10;
            // Avoid divide-by-zero errors
            if (radius < radiusThreshold)
            {
                rErrPcg = double.NaN;
            }
            else
            {
                rErrPcg = 100.0 * (radius/radiusAnalytical - 1.0);
            }
            thetaErrDeg = Math.Abs((180.0/Math.PI)*Atan2(Trig.Sin(thetaAnalytical-theta),Trig.Cos(thetaAnalytical-theta)));
            return (radius, rErrPcg, theta * 180.0/Math.PI, thetaErrDeg);
        }

        private static (double, double) RThetaFromXY(double x, double y)
        {
            return (Math.Sqrt(x * x + y * y), Atan2(y,x));
        }

        private static double Atan2(double y, double x)
        {
            if (x>0)
            {
                return Trig.Atan(y/x);
            }
            else if ((x < 0) && (y >= 0 ))
            {
                return Trig.Atan(y/x) + Constants.Pi;
            }
            else if ((x < 0) && (y < 0 ))
            {
                return Trig.Atan(y/x) - Constants.Pi;
            }
            else // x == 0
            {
                if (y > 0)
                {
                    return Constants.Pi / 2.0;
                }
                else if (y < 0)
                {
                    return -1.0 * Constants.Pi / 2.0;
                }
                else
                {
                    return double.NaN;
                }
            }
        }
    }
}