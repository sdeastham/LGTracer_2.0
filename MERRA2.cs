using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

namespace LGTracer
{
    public static class MERRA2
    {
        public static (double [], double[], int[], int[] ) ReadLatLon( string fileName, double[] lonLims, double[] latLims )
        {
            // Returns [lon_edge],[lat_edge],[lon_set],[lat_set]
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.ReadOnly
            };

            Func<double,double,double,int> findLower = (targetValue, lowerBound, spacing) => ((int)Math.Floor((targetValue - lowerBound)/spacing));
            double[] lonEdge,latEdge;
            float[] lonMids, latMids;
            int nLon, nLat, lonFirst, latFirst, lonLast, latLast;
            double dLon, dLat, lonBase, latBase;
            using (DataSet ds = DataSet.Open(dsUri))
            {
                // Get the full dimension vectors
                latMids = ds.GetData<float[]>("lat");
                lonMids = ds.GetData<float[]>("lon");

                // Figure out which cells we need to keep in order to get all the data we need
                // Assume a fixed cell spacing for now
                dLon = lonMids[1] - lonMids[0];
                lonBase = lonMids[0] - (dLon/2.0);
                // For latitude, be careful about half-polar grids
                dLat = latMids[3] - latMids[2];
                latBase = latMids[1] - (3.0*dLon/2.0);

                // These indices are for the first and last cell _inclusive_
                latFirst = findLower(latLims[0],latBase,dLat);
                latLast  = findLower(latLims[1],latBase,dLat);
                lonFirst = findLower(lonLims[0],lonBase,dLon);
                lonLast  = findLower(lonLims[1],lonBase,dLon);

                nLon = 1 + (lonLast - lonFirst);
                nLat = 1 + (latLast - latFirst);
            }
            // Create lon/lat edge vectors
            lonEdge = new double[nLon+1];
            latEdge = new double[nLat+1];
            lonEdge[0] = lonMids[lonFirst] - (dLon/2.0);
            latEdge[0] = latMids[latFirst] - (dLat/2.0);
            for (int i=0;i<nLon;i++)
            {
                lonEdge[i+1] = lonEdge[0] + (dLon * (i+1));
            }
            for (int i=0;i<nLat;i++)
            {
                latEdge[i+1] = latEdge[0] + (dLat * (i+1));
            }

            // To help in subsetting data, provide the lon/lat limits
            int[] lonSet = [lonFirst,lonLast];
            int[] latSet = [latFirst,latLast];
            return (lonEdge, latEdge, lonSet, latSet);
        }

        public static (double[,], double[,] ) ReadA3( string fileName, int time, int level, int[] lonSet, int[] latSet )
        {
            // Returns [u],[v]
            // Later extend to get T, QV
            // Open netCDF4 file
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.ReadOnly
            };

            double[,] u,v;
            int lonFirst = lonSet[0];
            int lonLast  = lonSet[1];
            int latFirst = latSet[0];
            int latLast  = latSet[1];
            int nLon, nLat;

            nLon = 1 + (lonLast - lonFirst);
            nLat = 1 + (latLast - latFirst);

            u = new double[nLat,nLon];
            v = new double[nLat,nLon];
            
            using (DataSet ds = DataSet.Open(dsUri))
            {
                // Be lazy for the moment and access the whole array (not certain if this reads into memory or just makes it available?)
                float[,,,] uFull = ds.GetData<float[,,,]>("U");
                float[,,,] vFull = ds.GetData<float[,,,]>("V");

                for (int iLon=0;iLon<nLon;iLon++)
                {
                    for (int iLat=0;iLat<nLat;iLat++)
                    {
                        u[iLat,iLon] = (double)uFull[time,level,iLat + latFirst,iLon + lonFirst];
                        v[iLat,iLon] = (double)vFull[time,level,iLat + latFirst,iLon + lonFirst];
                    }
                }
            }
            
            return (u, v);
        }

        public static (double[,], double[,], double[,] ) ReadI3( string fileName, int time, int level, int[] lonSet, int[] latSet )
        {
            // Returns PS, T, QV
            // Open netCDF4 file
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.ReadOnly
            };

            int lonFirst = lonSet[0];
            int lonLast  = lonSet[1];
            int latFirst = latSet[0];
            int latLast  = latSet[1];
            int nLon, nLat;

            nLon = 1 + (lonLast - lonFirst);
            nLat = 1 + (latLast - latFirst);

            double[,] ps = new double[nLat,nLon];
            double[,] temperature = new double[nLat,nLon];
            double[,] qv = new double[nLat,nLon];
            
            using (DataSet ds = DataSet.Open(dsUri))
            {
                // Be lazy for the moment and access the whole array (not certain if this reads into memory or just makes it available?)
                float[,,] psFull           = ds.GetData<float[,,]>("PS");
                float[,,,] temperatureFull = ds.GetData<float[,,,]>("T");
                float[,,,] qvFull          = ds.GetData<float[,,,]>("QV");

                for (int iLon=0;iLon<nLon;iLon++)
                {
                    for (int iLat=0;iLat<nLat;iLat++)
                    {
                        ps[iLat,iLon]          = (double)psFull[time,iLat + latFirst,iLon + lonFirst];
                        temperature[iLat,iLon] = (double)temperatureFull[time,level,iLat + latFirst,iLon + lonFirst];
                        qv[iLat,iLon]          = (double)qvFull[time,level,iLat + latFirst,iLon + lonFirst];
                    }
                }
            }
            
            return (ps, temperature, qv);
        }

        public static double[,,] CalculatePressures(double[,] surfacePressure,double[] AP, double[] BP)
        {
            // Calculate pressures at grid cell edges given surface pressures
            // AP and surfacePressure must be in the same units
            int nY = surfacePressure.GetLength(0);
            int nX = surfacePressure.GetLength(1);
            int nLevels = AP.Length;
            double[,,] pressureEdges = new double[nLevels,nY,nX];
            double localPS;
            for (int i=0; i<nX; i++)
            {
                for (int j=0; j<nY; j++)
                {
                    // Surface pressure
                    localPS = surfacePressure[j,i];
                    pressureEdges[0,j,i] = localPS;
                    // Hybrid eta calculation
                    for (int level=1; level<nLevels; level++)
                    {
                        pressureEdges[level,j,i] = (localPS * BP[level]) + AP[level];
                    }
                }
            }            
            return pressureEdges;
        }
    }
}