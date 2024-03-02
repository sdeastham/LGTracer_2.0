namespace LGTracer
{
    public class MetManagerFixed(string metDir, double[] lonLims, double[] latLims, DateTime startDate)
        : MetManager(metDir,
            lonLims, latLims, startDate)
    {
        // MetManagerFixed supplies data 
        public override double[,,] UWindXYP => MetFiles[A3DynIndex].DataVariables3D[UIndex].CurrentData;

        public override double[,,] VWindXYP => MetFiles[A3DynIndex].DataVariables3D[VIndex].CurrentData;

        public override double[,,] PressureVelocityXYP => MetFiles[A3DynIndex].DataVariables3D[OmegaIndex].CurrentData;

        public override double[,] SurfacePressureXY => MetFiles[I3Index].DataVariables2D[PSIndex].CurrentData;

        public override double[,,] SpecificHumidityXYP => MetFiles[I3Index].DataVariables3D[QVIndex].CurrentData; // kg water vapour per kg air

        public override double[,,] CloudIceXYP => MetFiles[A3CldIndex].DataVariables3D[QIIndex].CurrentData; // kg ice water per kg air
        public override double[,,] CloudWaterXYP => MetFiles[A3CldIndex].DataVariables3D[QLIndex].CurrentData; // kg liquid water per kg air
        public override double[,,] TemperatureXYP => MetFiles[I3Index].DataVariables3D[TIndex].CurrentData;
    }
}