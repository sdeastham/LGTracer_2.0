namespace LGTracer
{
    public class MetManager
    {
        // A MetManager handles all met data and exposes the met variables from all files

        protected List<MetFile> MetFiles;

        public double[,,] UWindXYP
        { get { return MetFiles[A3DynIndex].DataVariables3D[UIndex].CurrentData; } }
        public double[,,] VWindXYP
        { get { return MetFiles[A3DynIndex].DataVariables3D[VIndex].CurrentData; } }
        public double[,,] PressureVelocityXYP
        { get { return MetFiles[A3DynIndex].DataVariables3D[OmegaIndex].CurrentData; } }

        public double[,] SurfacePressureXY
        { get { return MetFiles[I3Index].DataVariables2D[PSIndex].CurrentData; } }
        public double[,,] SpecificHumidityXYP
        { get { return MetFiles[I3Index].DataVariables3D[QVIndex].CurrentData; } }
        public double[,,] TemperatureXYP
        { get { return MetFiles[I3Index].DataVariables3D[TIndex].CurrentData; } }
        private int I3Index, A3DynIndex;
        private int PSIndex, TIndex, QVIndex;
        private int UIndex, VIndex, OmegaIndex;

        public MetManager(string metDir, double[] lonLims, double[] latLims, DateTime startDate)
        {
            MetFiles = [];

            string[] A3DynVarList2D = {};
            string[] A3DynVarList3D = {"U","V","OMEGA"};
            int A3DynOffset = -30 * 60;

            string[] I3VarList2D = {"PS"};
            string[] I3VarList3D = {"QV","T"};
            int I3Offset = 0;
            MetFile currentFile;
            string currentTemplate;

            // 3-hour instantaneous
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.I3.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,I3VarList2D,I3VarList3D,lonLims,latLims,I3Offset);
            MetFiles.Add(currentFile);

            // Set up connections
            I3Index = MetFiles.Count() - 1;
            PSIndex = Array.FindIndex(I3VarList2D,element => element == "PS");
            QVIndex = Array.FindIndex(I3VarList3D,element => element == "QV");
            TIndex  = Array.FindIndex(I3VarList3D,element => element == "T");

            // 3-hour averaged, dynamics
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,A3DynVarList2D,A3DynVarList3D,lonLims,latLims,A3DynOffset);
            MetFiles.Add(currentFile);

            A3DynIndex = MetFiles.Count() - 1;
            UIndex      = Array.FindIndex(A3DynVarList3D,element => element == "U");
            VIndex      = Array.FindIndex(A3DynVarList3D,element => element == "V");
            OmegaIndex  = Array.FindIndex(A3DynVarList3D,element => element == "OMEGA");
        }

        public void AdvanceToTime(DateTime targetTime)
        {
            foreach (MetFile metFile in MetFiles)
            {
                metFile.AdvanceToTime(targetTime);
            }
        }

        public (double[], double[]) GetXYMesh()
        {
            // Return the X and Y edge vectors from the first file in our possession
            return MetFiles[0].GetXYMesh();
        }
    }
}