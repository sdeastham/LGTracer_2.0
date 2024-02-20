using System;
using System.Linq;
using System.Collections.Generic;

using System.Numerics;

namespace LGTracer
{
    public class DomainManager
    {
        public double TestValue
        { get; set; }

        public DomainManager()
        {
            TestValue = double.NaN;
        }

        public static (double, double) VelocityFromFixedSpaceArray( double x, double y, double xMin, double dx, double yMin, double dy, double[,] xSpeedArray, double[,] ySpeedArray, bool noConvert=false)
        {
            // Extract the velocity vector from an array
            // Values in m/s
            // Inefficient as we repeat the neighbor calculation
            double dxdt = NearestNeighbor(x,y,xMin,dx,yMin,dy,xSpeedArray);
            double dydt = NearestNeighbor(x,y,xMin,dx,yMin,dy,ySpeedArray);
            if (noConvert)
            {
                return (dxdt,dydt);
            }

            // Convert from m/s to deg/s
            dxdt = LGConstants.Rad2Deg * dxdt / (LGConstants.EarthRadius * Math.Cos(LGConstants.Deg2Rad*y));
            dydt = LGConstants.Rad2Deg * dydt / LGConstants.EarthRadius;
            return (dxdt, dydt);
        }

        public static double NearestNeighbor( double x, double y, double xMin, double dx, double yMin, double dy, double[,] valueArray)
        {
            // Extract the velocity vector from an array
            // Assumes constant X spacing and constant Y spacing
            int nX = valueArray.GetLength(1);
            int nY = valueArray.GetLength(0);

            // If we made it this far - we are within the domain
            int xIndex = Math.Min(Math.Max(0,(int)Math.Floor((x - xMin)/dx)),nX-1);
            int yIndex = Math.Min(Math.Max(0,(int)Math.Floor((y - yMin)/dy)),nY-1);

            // Values in m/s
            return valueArray[yIndex,xIndex];
        }

        public static (double[], double[]) SeedBoundaryUniform(int nPoints, double[] xLims, double[] yLims, System.Random RNG)
        {
            // Seed the domain boundaries completely uniformly
            double xMin = xLims[0];
            double xMax = xLims[1];
            double xSpan = xMax - xMin;
            double yMin = yLims[0];
            double yMax = yLims[1];
            double ySpan = yMax - yMin;
            double xCurr, yCurr;
            double smallDelta = 1.0e-5;
            double randomVal;

            double[] xVals = new double[nPoints];
            double[] yVals = new double[nPoints];
            
            for (int i=0; i<nPoints; i++)
            {
                // Activate the point at a location randomly chosen from the domain edge
                // Algorithm below basically goes around the edges of the domain in order
                randomVal = RNG.NextDouble() * ((xSpan*2) + (ySpan*2));
                if (randomVal < xSpan)
                {
                    yCurr = yMin + smallDelta;
                    xCurr = xMin + randomVal;
                }
                else if (randomVal < (xSpan + ySpan))
                {
                    yCurr = yMin + (randomVal - xSpan);
                    xCurr = xMax - smallDelta;
                }
                else if (randomVal < (xSpan + ySpan + xSpan))
                {
                    yCurr = yMax - smallDelta;
                    xCurr = xMin + (randomVal - (xSpan + ySpan));
                }
                else
                {
                    yCurr = yMin + (randomVal - (xSpan + ySpan + xSpan));
                    xCurr = xMin + smallDelta;
                }
                xVals[i] = xCurr;
                yVals[i] = yCurr;
            }
            return (xVals, yVals);
        }

        public static (Vector2[], Vector2[]) CreateBoundary( double[] xMeshDouble, double[] yMeshDouble)
        {
            // Create two arrays of 2-element vectors representing the boundary edge locations and the boundary normal vectors

            // Number of individual edge cells + 1
            int xPosts = xMeshDouble.Length;
            int yPosts = yMeshDouble.Length;

            // Vector2 is float only - do some conversions
            float[] xMesh = new float[xPosts];
            float[] yMesh = new float[yPosts];
            for (int i=0; i<xPosts; i++)
            {
                xMesh[i] = (float)xMeshDouble[i];
            }
            for (int i=0; i<yPosts; i++)
            {
                yMesh[i] = (float)yMeshDouble[i];
            }

            // Duplicate first post for convenience
            int nPosts = xPosts + (yPosts-1) + (xPosts-1) + (yPosts-1);
            Vector2[] xyPosts = new Vector2[nPosts];
            int currPost = 0;
            // South boundary
            for (int i=0; i<xPosts; i++)
            {
                xyPosts[currPost] = new Vector2(xMesh[i],yMesh[0]);
                currPost++;
            }
            // East boundary
            for (int i=1; i<yPosts; i++)
            {
                xyPosts[currPost] = new Vector2(xMesh[xPosts-1],yMesh[i]);
                currPost++;
            }
            // North boundary
            for (int i=xPosts-2; i>=0; i--)
            {
                xyPosts[currPost] = new Vector2(xMesh[i],yMesh[yPosts-1]);
                currPost++;
            }
            // West boundary (assume closure - ie include a final post identical to the first)
            for (int i=yPosts-2; i>=0; i--)
            {
                xyPosts[currPost] = new Vector2(xMesh[0],yMesh[i]);
                currPost++;
            }
            
            // Now calculate the boundary normals
            // First and last posts are duplicates for convenience
            Vector2[] normals = new Vector2[nPosts - 1];
            Vector2 boundaryVector2;
            Vector3 boundaryVector3, vertVector3;
            // Find the normal of each boundary edge by taking the cross-product of the "up" vector
            // with that of the boundary - this will always face into the domain
            vertVector3 = new Vector3(0,0,1);
            Vector3 normVector3;
            for (int i=0; i<(nPosts-1); i++)
            {
                boundaryVector2 = xyPosts[i+1] - xyPosts[i];
                boundaryVector3 = new Vector3(boundaryVector2,0);
                normVector3 = Vector3.Normalize(Vector3.Cross(vertVector3,boundaryVector3));
                normals[i] = new Vector2(normVector3.X,normVector3.Y);
            }

            return (xyPosts,normals);
        }

        public static Vector2[] GetBoundaryVelocities(Vector2[] xyPosts, Func<double, double, (double, double)> vCalc)
        {
            int nCells = xyPosts.Length - 1;
            double u, v, x, y;
            Vector2 xyMid;
            Vector2[] vBoundary = new Vector2[nCells];
            for (int i=0; i<nCells; i++)
            {
                xyMid = xyPosts[i] + (0.5f * (xyPosts[i+1] - xyPosts[i]));
                x = (double)xyMid.X;
                y = (double)xyMid.Y;
                (u, v) = vCalc(x,y);
                vBoundary[i] = new Vector2((float)u,(float)v);
            }
            return vBoundary;
        }
        public static (double[], double[], double) SeedBoundary(double kgPerPoint, double[] boundaryLengths, double pressureDelta, double dt,
            Vector2[] boundaryPosts, Vector2[] boundaryNormals, Vector2[] vBoundary, System.Random RNG, double massSurplus = 0.0)
        {
            // Seed the domain boundaries proportional to mass flow rate
            // Position along boundary for each cell is random
            int nEdges = boundaryPosts.Length - 1; // First post is duplicated as last post
            double smallDelta = 1.0e-5;
            double massFlux, cellFrac, randomVal;
            double pointSurplus = 0.0;
            double nPoints;
            Vector2 pointLocation, boundaryVector;
            double vNorm;
            double[] massFluxes = new double[nEdges];
            double[] weighting;

            for (int i=0; i<nEdges; i++)
            {
                // Mass flow rate across the boundary in kg/s
                vNorm = (double)Vector2.Dot(vBoundary[i],boundaryNormals[i]);
                massFlux = dt * vNorm * boundaryLengths[i] * pressureDelta / LGConstants.gravConstantSurface;
                massFluxes[i] = Math.Max(0.0,massFlux);
            }

            // We now know the total mass flux - use to figure out weightings
            massFlux = massFluxes.Sum();
            weighting = new double[nEdges];
            for (int i=0; i<nEdges; i++)
            {
                // Cumulative
                weighting[i] = (massFluxes[i] / massFlux);
                if (i>0)
                {
                    weighting[i] += weighting[i-1];
                }
            }
            // Force to 1.0 to ensure we don't accidentally miss anything
            weighting[nEdges-1] = 1.0;

            // How many points will we add (add mass surplus)?
            nPoints = (massFlux + massSurplus)/kgPerPoint;
            int nPointsTotal = (int)Math.Floor(nPoints);

            // Unused mass
            massSurplus = (nPoints - (double)nPointsTotal) * kgPerPoint;

            //Console.WriteLine($"Expecting {nPointsTotal} new points");
            double[] xVals = new double[nPointsTotal];
            double[] yVals = new double[nPointsTotal];
            int iCell;
            for (int iPoint=0; iPoint < nPointsTotal; iPoint++)
            {
                randomVal = RNG.NextDouble();
                iCell = 0;
                while (weighting[iCell] < randomVal)
                {
                    iCell++;
                }
                boundaryVector = boundaryPosts[iCell+1] - boundaryPosts[iCell];
                // Remainder of randomVal used to figure out how far along the cell we go
                if (iCell > 0)
                {
                    cellFrac = (randomVal - weighting[iCell-1])/(weighting[iCell] - weighting[iCell-1]);
                }
                else
                { 
                    cellFrac = randomVal/weighting[0];
                }
                pointLocation = boundaryPosts[iCell] + ((float)cellFrac * boundaryVector) + ((float)smallDelta * boundaryNormals[iCell]);
                xVals[iPoint] = pointLocation.X;
                yVals[iPoint] = pointLocation.Y;
            }
            return (xVals, yVals, massSurplus);
        }

        public static (double[], double[]) MapRandomToXY( double xMin, double xMax, double yMin, double yMax, int nPoints, System.Random rng )
        {
            // Scatter randomly throughout domain
            double xSpan = xMax - xMin;
            double xStart = xMin;
            double ySpan = yMax - yMin;
            double yStart = yMin;

            double[] xInitial = new double[nPoints];
            double[] yInitial = new double[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xInitial[i] = rng.NextDouble()*xSpan + xStart;
                yInitial[i] = rng.NextDouble()*ySpan + yStart;
            }
            return (xInitial, yInitial);
        }

        public static void Cull(double[] xLims, double[] yLims, PointManager pointManager)
        {
            // Deactivate any points which are outside the domain
            // Could also do this with LINQ
            for (int i=pointManager.NActive-1; i>=0; i--)
            {
                LGPoint point = pointManager.ActivePoints[i];
                if (point.X < xLims[0] || point.X >= xLims[1] || point.Y < yLims[0] || point.Y >= yLims[1] )
                {
                    pointManager.DeactivatePoint(i);
                }
            }
        }
    }
}