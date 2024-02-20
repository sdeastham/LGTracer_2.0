using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;

using System.Numerics;

namespace LGTracer
{
    public class DomainManager
    {

        public double XMin
        { get; protected set; }
        
        public double YMin
        { get; protected set; }
        
        public double XMax
        { get; protected set; }

        public double YMax
        { get; protected set; }

        public int NX
        { get; protected set; }

        public int NY
        { get; protected set; }

        public double DX
        { get; protected set; }

        public double DY
        { get; protected set; }

        public double[] XLims
        { get; protected set; }

        public double[] YLims
        { get; protected set; }
        
        public double[] XMesh
        { get; protected set; }
        
        public double[] XMids
        { get; protected set; }

        public double[] YMesh
        { get; protected set; }

        public double[] YMids
        { get; protected set; }

        public double[] BoundaryLengths
        { get; protected set; }

        public Vector2[] BoundaryNormals
        { get; protected set; }

        public Vector2[] BoundaryPosts
        { get; protected set; }

        public DomainManager(double[] lonEdge, double[] latEdge)
        {
            // Set up the mesh (units of degrees)
            // The original lon/lat limits aren't important - need the mesh limits
            int xPosts = lonEdge.Length;
            NX = xPosts - 1;
            int yPosts = latEdge.Length;
            NY = yPosts - 1;

            XMin  = lonEdge[0];
            XMax  =  lonEdge[xPosts-1];
            double xSpan = XMax - XMin;

            YMin  = latEdge[0];
            YMax  = latEdge[yPosts-1];
            double ySpan = YMax - YMin;

            // Assume uniform spacing
            DX = xSpan/(xPosts - 1);
            DY = ySpan/(yPosts - 1);

            // For convenience
            XLims = new double[] {XMin,XMax};
            YLims = new double[] {YMin,YMax};

            // X edges (degrees)
            XMesh = new double[xPosts];
            XMids = new double[NX];
            for (int i=0;i<xPosts;i++)
            {
                XMesh[i] = lonEdge[i];
                if (i > 0)
                {
                    XMids[i-1] = (XMesh[i] + XMesh[i-1])/2.0;
                }
            }

            // Y edges (degrees)
            YMesh = new double[yPosts];
            YMids  = new double[NY];
            for (int i=0;i<yPosts;i++)
            {
                //yMesh[i] = yMin + (dy * i);
                YMesh[i] = latEdge[i];
                if (i > 0)
                {
                    YMids[i-1] = (YMesh[i] + YMesh[i-1])/2.0;
                }
            }

            // Boundary lengths (m)
            BoundaryLengths = new double[NX*2 + NY*2];
            double earthCircumference = 2.0 * Math.PI * LGConstants.EarthRadius;
            double edgeLength;

            // South boundary
            edgeLength = (DX/360.0) * earthCircumference * Math.Cos(LGConstants.Deg2Rad*latEdge[0]);
            for (int i=0;i<NX;i++)
            {
                BoundaryLengths[i] = edgeLength;
            }
            // North boundary
            edgeLength = (DX/360.0) * earthCircumference * Math.Cos(LGConstants.Deg2Rad*latEdge[yPosts-1]);
            for (int i=0;i<NX;i++)
            {
                BoundaryLengths[NX + NY + i] = edgeLength;
            }
            // All cells on the East and West boundaries have a constant length
            edgeLength = (DY/360.0) * earthCircumference;
            for (int i=0;i<NY;i++)
            {
                BoundaryLengths[NX + i] = edgeLength;
                BoundaryLengths[(NX*2) + NY + i] = edgeLength;
            }

            // Set up BoundaryNormals and BoundaryPosts
            CreateBoundary();
        }

        public (double, double) VelocityFromFixedSpaceArray( double x, double y, double[,] xSpeedArray, double[,] ySpeedArray, bool noConvert=false)
        {
            // Extract the velocity vector from an array
            // Values in m/s
            // Inefficient as we repeat the neighbor calculation
            double dxdt = NearestNeighbor(x,y,xSpeedArray);
            double dydt = NearestNeighbor(x,y,ySpeedArray);
            if (noConvert)
            {
                return (dxdt,dydt);
            }

            // Convert from m/s to deg/s
            dxdt = LGConstants.Rad2Deg * dxdt / (LGConstants.EarthRadius * Math.Cos(LGConstants.Deg2Rad*y));
            dydt = LGConstants.Rad2Deg * dydt / LGConstants.EarthRadius;
            return (dxdt, dydt);
        }

        public double NearestNeighbor( double x, double y, double[,] valueArray)
        {
            // Extract the velocity vector from an array
            // Assumes constant X spacing and constant Y spacing
            int xIndex = Math.Min(Math.Max(0,(int)Math.Floor((x - XMin)/DY)),NX-1);
            int yIndex = Math.Min(Math.Max(0,(int)Math.Floor((y - YMin)/DY)),NY-1);

            // Values in m/s
            return valueArray[yIndex,xIndex];
        }

        public (double[], double[]) SeedBoundaryUniform(int nPoints, System.Random RNG)
        {
            // Seed the domain boundaries completely uniformly
            double xSpan = XMax - XMin;
            double ySpan = YMax - YMin;
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
                    yCurr = YMin + smallDelta;
                    xCurr = XMin + randomVal;
                }
                else if (randomVal < (xSpan + ySpan))
                {
                    yCurr = YMin + (randomVal - xSpan);
                    xCurr = XMax - smallDelta;
                }
                else if (randomVal < (xSpan + ySpan + xSpan))
                {
                    yCurr = YMax - smallDelta;
                    xCurr = XMin + (randomVal - (xSpan + ySpan));
                }
                else
                {
                    yCurr = YMin + (randomVal - (xSpan + ySpan + xSpan));
                    xCurr = XMin + smallDelta;
                }
                xVals[i] = xCurr;
                yVals[i] = yCurr;
            }
            return (xVals, yVals);
        }

        [MemberNotNull(nameof(BoundaryPosts),nameof(BoundaryNormals))]
        private void CreateBoundary()
        {
            // Create two arrays of 2-element vectors representing the boundary edge locations and the boundary normal vectors
            // Vector2 is float only - do some conversions
            int xPosts = NX + 1;
            int yPosts = NY + 1;
            float[] xMeshFloat = new float[xPosts];
            float[] yMeshFloat = new float[yPosts];
            for (int i=0; i<xPosts; i++)
            {
                xMeshFloat[i] = (float)XMesh[i];
            }
            for (int i=0; i<yPosts; i++)
            {
                yMeshFloat[i] = (float)YMesh[i];
            }

            // Duplicate first post for convenience
            int nPosts = xPosts + (yPosts-1) + (xPosts-1) + (yPosts-1);
            BoundaryPosts = new Vector2[nPosts];
            int currPost = 0;
            // South boundary
            for (int i=0; i<xPosts; i++)
            {
                BoundaryPosts[currPost] = new Vector2(xMeshFloat[i],yMeshFloat[0]);
                currPost++;
            }
            // East boundary
            for (int i=1; i<yPosts; i++)
            {
                BoundaryPosts[currPost] = new Vector2(xMeshFloat[xPosts-1],yMeshFloat[i]);
                currPost++;
            }
            // North boundary
            for (int i=xPosts-2; i>=0; i--)
            {
                BoundaryPosts[currPost] = new Vector2(xMeshFloat[i],yMeshFloat[yPosts-1]);
                currPost++;
            }
            // West boundary (assume closure - ie include a final post identical to the first)
            for (int i=yPosts-2; i>=0; i--)
            {
                BoundaryPosts[currPost] = new Vector2(xMeshFloat[0],yMeshFloat[i]);
                currPost++;
            }
            
            // Now calculate the boundary normals
            // First and last posts are duplicates for convenience
            BoundaryNormals = new Vector2[nPosts - 1];
            Vector2 boundaryVector2;
            Vector3 boundaryVector3, vertVector3;
            // Find the normal of each boundary edge by taking the cross-product of the "up" vector
            // with that of the boundary - this will always face into the domain
            vertVector3 = new Vector3(0,0,1);
            Vector3 normVector3;
            for (int i=0; i<(nPosts-1); i++)
            {
                boundaryVector2 = BoundaryPosts[i+1] - BoundaryPosts[i];
                boundaryVector3 = new Vector3(boundaryVector2,0);
                normVector3 = Vector3.Normalize(Vector3.Cross(vertVector3,boundaryVector3));
                BoundaryNormals[i] = new Vector2(normVector3.X,normVector3.Y);
            }
        }

        public Vector2[] GetBoundaryVelocities(Func<double, double, (double, double)> vCalc)
        {
            int nCells = BoundaryPosts.Length - 1;
            double u, v, x, y;
            Vector2 xyMid;
            Vector2[] vBoundary = new Vector2[nCells];
            for (int i=0; i<nCells; i++)
            {
                xyMid = BoundaryPosts[i] + (0.5f * (BoundaryPosts[i+1] - BoundaryPosts[i]));
                x = (double)xyMid.X;
                y = (double)xyMid.Y;
                (u, v) = vCalc(x,y);
                vBoundary[i] = new Vector2((float)u,(float)v);
            }
            return vBoundary;
        }

        public (double[], double[], double) SeedBoundary(double kgPerPoint, double pressureDelta, double dt,
            Vector2[] vBoundary, System.Random RNG, double massSurplus = 0.0)
        {
            // Seed the domain boundaries proportional to mass flow rate
            // Position along boundary for each cell is random
            int nEdges = BoundaryPosts.Length - 1; // First post is duplicated as last post
            double smallDelta = 1.0e-5;
            double massFlux, cellFrac, randomVal;
            double nPoints;
            Vector2 pointLocation, boundaryVector;
            double vNorm;
            double[] massFluxes = new double[nEdges];
            double[] weighting;

            for (int i=0; i<nEdges; i++)
            {
                // Mass flow rate across the boundary in kg/s
                vNorm = (double)Vector2.Dot(vBoundary[i],BoundaryNormals[i]);
                //Console.WriteLine($"{i} -> {BoundaryLengths[i]}");
                massFlux = dt * vNorm * BoundaryLengths[i] * pressureDelta / LGConstants.gravConstantSurface;
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
                boundaryVector = BoundaryPosts[iCell+1] - BoundaryPosts[iCell];
                // Remainder of randomVal used to figure out how far along the cell we go
                if (iCell > 0)
                {
                    cellFrac = (randomVal - weighting[iCell-1])/(weighting[iCell] - weighting[iCell-1]);
                }
                else
                { 
                    cellFrac = randomVal/weighting[0];
                }
                pointLocation = BoundaryPosts[iCell] + ((float)cellFrac * boundaryVector) + ((float)smallDelta * BoundaryNormals[iCell]);
                xVals[iPoint] = pointLocation.X;
                yVals[iPoint] = pointLocation.Y;
            }
            return (xVals, yVals, massSurplus);
        }

        public (double[], double[]) MapRandomToXY( int nPoints, System.Random rng )
        {
            // Scatter randomly throughout domain
            double xSpan = XMax - XMin;
            double xStart = XMin;
            double ySpan = YMax - YMin;
            double yStart = YMin;

            double[] xInitial = new double[nPoints];
            double[] yInitial = new double[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xInitial[i] = rng.NextDouble()*xSpan + xStart;
                yInitial[i] = rng.NextDouble()*ySpan + yStart;
            }
            return (xInitial, yInitial);
        }
    }
}