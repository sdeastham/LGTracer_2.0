namespace LGTracer
{
    public class PointManagerDense : PointManager
    {
        // A "dense" point manager is designed to represent ALL air within the domain
        public double MassSurplus;
        private System.Random Rng;
        private double KgPerPoint;

        public PointManagerDense( long? maxPoints, DomainManager domain, string filename, bool debug=false, 
            bool includeCompression=false, string[]? propertyNames=null, System.Random? rng=null,
            double kgPerPoint=1.0e12 ) : base(maxPoints,domain,filename,debug,includeCompression,propertyNames)
        {
            // No initial mass surplus
            MassSurplus = 0.0;
            KgPerPoint = kgPerPoint;
            Rng = rng;
        }

        public override void Seed(double dt)
        {
            (double[] xSet, double[] ySet, double[] pSet, MassSurplus) =
                Domain.SeedBoundary(KgPerPoint, dt, Rng, MassSurplus);
            CreatePointSet(xSet, ySet, pSet);

            (double[] xSetV, double[] ySetV, double[] pSetV, MassSurplus) =
                Domain.SeedPressureBoundaries(KgPerPoint, dt, Rng, MassSurplus);
            CreatePointSet(xSetV, ySetV, pSetV);
        }
    }
}