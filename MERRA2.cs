using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;
using Microsoft.Research.Science.Data.NetCDF4;

namespace LGTracer
{
    public static class MERRA2
    {
        // 72 levels, so 73 level edges
        // Values are in Pa
        public static readonly double[] AP = [ 0.000000e+00, 4.804826e+00, 6.593752e+02, 1.313480e+03, 1.961311e+03, 2.609201e+03,
                              3.257081e+03, 3.898201e+03, 4.533901e+03, 5.169611e+03, 5.805321e+03, 6.436264e+03,
                              7.062198e+03, 7.883422e+03, 8.909992e+03, 9.936521e+03, 1.091817e+04, 1.189586e+04,
                              1.286959e+04, 1.429100e+04, 1.562600e+04, 1.696090e+04, 1.816190e+04, 1.930970e+04,
                              2.032590e+04, 2.121500e+04, 2.187760e+04, 2.238980e+04, 2.243630e+04, 2.168650e+04,
                              2.011920e+04, 1.769300e+04, 1.503930e+04, 1.278370e+04, 1.086630e+04, 9.236572e+03,
                              7.851231e+03, 6.660341e+03, 5.638791e+03, 4.764391e+03, 4.017541e+03, 3.381001e+03,
                              2.836781e+03, 2.373041e+03, 1.979160e+03, 1.645710e+03, 1.364340e+03, 1.127690e+03,
                              9.292942e+02, 7.619842e+02, 6.216801e+02, 5.046801e+02, 4.076571e+02, 3.276431e+02,
                              2.620211e+02, 2.084970e+02, 1.650790e+02, 1.300510e+02, 1.019440e+02, 7.951341e+01,
                              6.167791e+01, 4.758061e+01, 3.650411e+01, 2.785261e+01, 2.113490e+01, 1.594950e+01,
                              1.197030e+01, 8.934502e+00, 6.600001e+00, 4.758501e+00, 3.270000e+00, 2.000000e+00,
                              1.000000e+00 ];

        public static readonly double[] BP = [ 1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
                              8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
                              7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
                              5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
                              2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
                              6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                              0.000000e+00 ];

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

        public static (double[,], double[,], double[,] ) ReadA3( string fileName, int time, int level, int[] lonSet, int[] latSet )
        {
            // Returns [u],[v],[omega]
            // Open netCDF4 file
            var dsUri = new NetCDFUri
            {
                FileName = fileName,
                OpenMode = ResourceOpenMode.ReadOnly
            };

            double[,] u,v,w;
            int lonFirst = lonSet[0];
            int lonLast  = lonSet[1];
            int latFirst = latSet[0];
            int latLast  = latSet[1];
            int nLon, nLat;

            nLon = 1 + (lonLast - lonFirst);
            nLat = 1 + (latLast - latFirst);

            u = new double[nLat,nLon];
            v = new double[nLat,nLon];
            w = new double[nLat,nLon];
            
            using (DataSet ds = DataSet.Open(dsUri))
            {
                // Be lazy for the moment and access the whole array (not certain if this reads into memory or just makes it available?)
                float[,,,] uFull = ds.GetData<float[,,,]>("U");
                float[,,,] vFull = ds.GetData<float[,,,]>("V");
                float[,,,] wFull = ds.GetData<float[,,,]>("OMEGA");

                for (int iLon=0;iLon<nLon;iLon++)
                {
                    for (int iLat=0;iLat<nLat;iLat++)
                    {
                        u[iLat,iLon] = (double)uFull[time,level,iLat + latFirst,iLon + lonFirst];
                        v[iLat,iLon] = (double)vFull[time,level,iLat + latFirst,iLon + lonFirst];
                        w[iLat,iLon] = (double)wFull[time,level,iLat + latFirst,iLon + lonFirst];
                    }
                }
            }
            
            return (u, v, w);
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
                float[,,] psFull           = ds.GetData<float[,,]>("PS"); // Pa
                float[,,,] temperatureFull = ds.GetData<float[,,,]>("T"); // K
                float[,,,] qvFull          = ds.GetData<float[,,,]>("QV"); // kg/kg

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
    }
}