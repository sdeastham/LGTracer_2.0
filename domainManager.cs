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

        public int NLevels
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

        public double[,,] XSpeedXYP
        { get; protected set; }

        public double[,,] YSpeedXYP
        { get; protected set; }
        
        public double[,,] PressureVelocityXYP
        { get; protected set; }

        public double[,] PressureVelocityXYBase
        { get; protected set; }

        public double[,] PressureVelocityXYCeiling
        { get; protected set; }

        public double[,] SurfacePressureXY
        { get; protected set; }

        public double[,,] PressureEdgeXYPe
        { get; protected set; }

        public double[,,] TemperatureXYP
        { get; protected set; }

        public double[,,] SpecificHumidityXYP
        { get; protected set; }

        public double PBase
        { get; protected set; }

        public double PCeiling
        { get; protected set; }

        protected double[] AP
        { get; set; }

        protected double[] BP
        { get; set; }

        public double[,] CellArea
        { get; protected set; }

        public DomainManager(double[] lonEdge, double[] latEdge, double[] pLimits, double[] pOffsets, double[] pFactors)
        {
            // Set up the vertical coordinates
            // NB: PBase and PCeiling indicate where we cull, not the vertical
            // extent of the meteorological data
            PBase = pLimits[0];
            PCeiling = pLimits[1];

            // Set up the factors needed to calculate pressure edges
            AP = pOffsets;
            BP = pFactors;
            NLevels = AP.Length - 1;

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

            // Cell areas
            SetCellAreas();

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

            // Set up the meteorological data
            XSpeedXYP = new double[NLevels,NY,NX];
            YSpeedXYP = new double[NLevels,NY,NX];
            SurfacePressureXY = new double[NY,NX];
            PressureEdgeXYPe = new double[NLevels+1,NY,NX];
            PressureVelocityXYP = new double[NLevels,NY,NX];
            TemperatureXYP = new double[NLevels,NY,NX];
            SpecificHumidityXYP = new double[NLevels,NY,NX];
            PressureVelocityXYBase = new double[NY,NX];
            PressureVelocityXYCeiling = new double[NY,NX];
        }

        [MemberNotNull(nameof(CellArea))]
        private void SetCellAreas()
        {
            double[] yRad = new double[NY + 1];
            for (int i=0; i<(NY+1); i++)
            {
                yRad[i] = YMesh[i] * LGConstants.Deg2Rad;
            }

            // Total surface area in each meridional band (assuming regular spacing)
            double bandArea = LGConstants.EarthRadius * LGConstants.EarthRadius * DX * LGConstants.Deg2Rad;
            CellArea = new double[NY,NX];
            double sinDiff, singleCellArea;
            for (int iLat=0; iLat<NY; iLat++)
            {
                // Fraction of a single band taken up by one cell
                sinDiff = Math.Sin(yRad[iLat+1]) - Math.Sin(yRad[iLat]);
                singleCellArea = sinDiff * bandArea;
                for (int iLon=0; iLon<NX; iLon++)
                {
                    CellArea[iLat,iLon] = singleCellArea;
                }
            }
        }

        public void UpdateMeteorology( double[,] surfacePressure, double[,,] xSpeed, double[,,] ySpeed, double[,,] omega, double[,,] temperature, double[,,] specificHumidity )
        {
            // Update stored meteorology and derived quantities
            int kBase, kCeiling;
            for (int i=0; i<NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    SurfacePressureXY[j,i] = surfacePressure[j,i];
                    for (int k=0; k<NLevels; k++)
                    {
                        XSpeedXYP[k,j,i] = xSpeed[k,j,i];
                        YSpeedXYP[k,j,i] = ySpeed[k,j,i];
                        PressureVelocityXYP[k,j,i] = omega[k,j,i];
                        TemperatureXYP[k,j,i] = temperature[k,j,i];
                        SpecificHumidityXYP[k,j,i] = specificHumidity[k,j,i];
                    }
                }
            }

            // Get the pressure edges too
            PressureEdgeXYPe = CalculatePressures(SurfacePressureXY);

            // Find the level which corresponds to the base and ceiling of the domain
            for (int i=0; i<NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    kBase = -1;
                    kCeiling = -1;
                    for (int k=0; k<NLevels; k++)
                    {
                        // If next edge would be above (at lower pressure) than
                        // the base pressure, we have found the lower bound
                        if (kBase < 0 && PressureEdgeXYPe[k+1,j,i] < PBase)
                        {
                            kBase = k;
                            // Set base quantities
                            PressureVelocityXYBase[j,i] = PressureVelocityXYP[kBase,j,i];
                        }
                        if (kCeiling < 0 && PressureEdgeXYPe[k+1,j,i] < PCeiling)
                        {
                            kCeiling = k;
                            // Set ceiling quantities
                            PressureVelocityXYCeiling[j,i] = PressureVelocityXYP[kCeiling,j,i];
                            // Since kCeiling > kBase, we don't need to keep going
                            break;
                        }
                    }
                }
            }
        }

        public (double, double, double) VelocityFromFixedSpaceArray( double x, double y, double pressure, bool noConvert=false)
        {
            // Extract the velocity vector from an array
            // Values in m/s
            // Inefficient as we repeat the neighbor calculation
            double dxdt = NearestNeighbor3D(x,y,pressure,XSpeedXYP);
            double dydt = NearestNeighbor3D(x,y,pressure,YSpeedXYP);
            double dpdt = NearestNeighbor3D(x,y,pressure,PressureVelocityXYP);
            //double dpdt = 0.0;
            if (noConvert)
            {
                return (dxdt,dydt,dpdt);
            }

            // Convert from m/s to deg/s (keep omega as Pa/s)
            dxdt = LGConstants.Rad2Deg * dxdt / (LGConstants.EarthRadius * Math.Cos(LGConstants.Deg2Rad*y));
            dydt = LGConstants.Rad2Deg * dydt / LGConstants.EarthRadius;
            return (dxdt, dydt, dpdt);
        }

        public double NearestNeighbor2D( double x, double y, double[,] valueArray)
        {
            // Ignore pressure for now
            // Extract the velocity vector from an array
            // Assumes constant X spacing and constant Y spacing
            int xIndex = Math.Min(Math.Max(0,(int)Math.Floor((x - XMin)/DX)),NX-1);
            int yIndex = Math.Min(Math.Max(0,(int)Math.Floor((y - YMin)/DY)),NY-1);

            // Values in m/s
            return valueArray[yIndex,xIndex];
        }

        public double NearestNeighbor3D( double x, double y, double pressure, double[,,] valueArray)
        {
            // Ignore pressure for now
            // Extract the velocity vector from an array
            // Assumes constant X spacing and constant Y spacing
            int xIndex = Math.Min(Math.Max(0,(int)Math.Floor((x - XMin)/DX)),NX-1);
            int yIndex = Math.Min(Math.Max(0,(int)Math.Floor((y - YMin)/DY)),NY-1);
            // Need to search through the pressure edges at this location
            double ps = SurfacePressureXY[yIndex,xIndex];
            int pIndex = 0;
            while (pressure < PressureEdgeXYPe[pIndex+1,yIndex,xIndex])
            {
                pIndex++;
            }
            return valueArray[pIndex,yIndex,xIndex];
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

        public Vector2[] GetBoundaryVelocities()
        {
            // u,v velocities at the cell edges
            int nCells = BoundaryPosts.Length - 1;
            double u, v, w, x, y, pressure;
            Vector2 xyMid;
            Vector2[] vBoundary = new Vector2[nCells];
            for (int i=0; i<nCells; i++)
            {
                xyMid = BoundaryPosts[i] + (0.5f * (BoundaryPosts[i+1] - BoundaryPosts[i]));
                x = (double)xyMid.X;
                y = (double)xyMid.Y;
                pressure = (PBase + PCeiling)/2.0;
                (u, v, w) = VelocityFromFixedSpaceArray(x,y,pressure,true);
                vBoundary[i] = new Vector2((float)u,(float)v);
            }
            return vBoundary;
        }

        public (double[], double[], double[], double) SeedBoundary(double kgPerPoint, double dt, System.Random RNG, double massSurplus = 0.0)
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
            double pressureDelta = PBase - PCeiling;

            Vector2[] vBoundary = GetBoundaryVelocities();

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
            double[] pressureVals = new double[nPointsTotal];
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
                // Initialize pressure randomly - unrelated to other concerns
                pressureVals[iPoint] = RNG.NextDouble() * (PBase - PCeiling) + PCeiling;
            }
            return (xVals, yVals, pressureVals, massSurplus);
        }

        public (double[], double[], double[], double) SeedPressureBoundaries(double kgPerPoint, double dt, System.Random RNG, double massSurplus = 0.0)
        {
            // Seed the upper and lower domain boundaries proportional to vertical mass flow rate
            // Position along boundary for each cell is random
            int nFaces = NX * NY;
            double smallDelta = 1.0e-5; // This is now in pressure terms
            double massFlux, cellFrac, randomVal;
            double nPoints;
            Vector2 pointLocation;
            double vNorm;
            // Double because we will deal with upper and lower boundary together
            double[] massFluxes = new double[nFaces*2];
            double[] weighting;
            double pressureDelta = PBase - PCeiling;

            int iFace = 0;
            for (int i=0; i<NX; i++)
            {
                for (int j=0; j<NY; j++)
                {
                    // Mass flow rate across the boundary in kg/s
                    // omega * area / g = kg/s
                    // (vNorm * omega) = dot product of velocity with boundary normal
                    // Lower boundary first: negative omega -> mass flow into domain
                    vNorm = -1.0;
                    massFlux =  dt * (vNorm * PressureVelocityXYBase[j,i]) * CellArea[j,i] / LGConstants.gravConstantSurface;
                    massFluxes[iFace] = Math.Max(0.0,massFlux);

                    // Now the upper boundary: positive omega -> mass flow into domain
                    vNorm = +1.0;
                    massFlux = dt * (vNorm * PressureVelocityXYCeiling[j,i]) * CellArea[j,i] / LGConstants.gravConstantSurface;
                    massFluxes[iFace + nFaces] = Math.Max(0.0,massFlux);
                    iFace++;
                }
            }

            // We now know the total mass flux - use to figure out weightings
            massFlux = massFluxes.Sum();
            weighting = new double[nFaces*2];
            for (int i=0; i<(nFaces*2); i++)
            {
                // Cumulative
                weighting[i] = (massFluxes[i] / massFlux);
                if (i>0)
                {
                    weighting[i] += weighting[i-1];
                }
            }
            // Force to 1.0 to ensure we don't accidentally miss anything
            weighting[(nFaces*2)-1] = 1.0;

            // How many points will we add (add mass surplus)?
            nPoints = (massFlux + massSurplus)/kgPerPoint;
            int nPointsTotal = (int)Math.Floor(nPoints);

            // Unused mass
            massSurplus = (nPoints - (double)nPointsTotal) * kgPerPoint;

            //Console.WriteLine($"Expecting {nPointsTotal} new points");
            double[] xVals = new double[nPointsTotal];
            double[] yVals = new double[nPointsTotal];
            double[] pressureVals = new double[nPointsTotal];
            int iCell, iCellOneFace, iX, iY;
            double pVal, xVal, yVal;
            for (int iPoint=0; iPoint < nPointsTotal; iPoint++)
            {
                xVal = RNG.NextDouble();
                yVal = RNG.NextDouble();
                randomVal = RNG.NextDouble();
                iCell = 0;
                while (weighting[iCell] < randomVal)
                {
                    iCell++;
                }
                // First N are lower boundary; second N are upper boundary
                if (iCell < nFaces)
                {
                    pVal = PBase - smallDelta;
                }
                else
                {
                    pVal = PCeiling + smallDelta;
                }
                iCellOneFace = iCell % nFaces;

                // Iterate over Y internally, X externally in weight-setting loop
                iY = iCellOneFace%NY;
                iX = (int)Math.Floor((double)(iCellOneFace - iY)/(double)NY);
                //Console.WriteLine($"{iCellOneFace,5:d} --> {iX,4:d}/{NX,4:d}, {iY,4:d}/{NY,4:d}");
                xVals[iPoint] = (XMesh[iX]) + xVal * DX;
                yVals[iPoint] = (YMesh[iY]) + yVal * DY;
                pressureVals[iPoint] = pVal;
            }
            return (xVals, yVals, pressureVals, massSurplus);
        }

        public (double[], double[], double[]) MapRandomToXYP( int nPoints, System.Random rng )
        {
            // Scatter randomly throughout domain
            double xSpan = XMax - XMin;
            double xStart = XMin;
            double ySpan = YMax - YMin;
            double yStart = YMin;
            double pStart = PCeiling;
            double pSpan = PBase - PCeiling;

            double[] xInitial = new double[nPoints];
            double[] yInitial = new double[nPoints];
            double[] pressureInitial = new double[nPoints];
            for (int i=0; i<nPoints; i++)
            {
                xInitial[i] = rng.NextDouble()*xSpan + xStart;
                yInitial[i] = rng.NextDouble()*ySpan + yStart;
                pressureInitial[i] = rng.NextDouble()*pSpan + pStart;
            }
            return (xInitial, yInitial, pressureInitial);
        }

        public double[,,] CalculatePressures(double[,] surfacePressure)
        {
            // Calculate pressures at grid cell edges given surface pressures
            // AP and surfacePressure must be in the same units
            int nY = surfacePressure.GetLength(0);
            int nX = surfacePressure.GetLength(1);
            int nLevelEdges = AP.Length;
            double[,,] pressureEdges = new double[nLevelEdges,nY,nX];
            double localPS;
            for (int i=0; i<nX; i++)
            {
                for (int j=0; j<nY; j++)
                {
                    // Surface pressure
                    localPS = surfacePressure[j,i];
                    // Hybrid eta calculation
                    for (int level=0; level<nLevelEdges; level++)
                    {
                        pressureEdges[level,j,i] = (localPS * BP[level]) + AP[level];
                    }
                }
            }            
            return pressureEdges;
        }
    }
}