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
        public double[,,] SpecificHumidityXYP // kg water vapour per kg air
        { get { return MetFiles[I3Index].DataVariables3D[QVIndex].CurrentData; } }
        public double[,,] CloudIceXYP // kg ice water per kg air
        { get { return MetFiles[A3CldIndex].DataVariables3D[QIIndex].CurrentData; } }
        public double[,,] CloudWaterXYP // kg liquid water per kg air
        { get { return MetFiles[A3CldIndex].DataVariables3D[QLIndex].CurrentData; } }
        public double[,,] TemperatureXYP
        { get { return MetFiles[I3Index].DataVariables3D[TIndex].CurrentData; } }
        private int I3Index, A3DynIndex, A3CldIndex;
        private int PSIndex, TIndex, QVIndex, QIIndex, QLIndex;
        private int UIndex, VIndex, OmegaIndex;

        public MetManager(string metDir, double[] lonLims, double[] latLims, DateTime startDate)
        {
            MetFiles = [];
            
            // Read time offsets in seconds - A3 files are timestamped for
            // mid-way through the 3 hour averaging period, but we want to read
            // them at the start of their averaging period
            int A3Offset = -30 * 60;
            int I3Offset = 0;

            string[] A3DynVarList2D = {};
            string[] A3DynVarList3D = {"U","V","OMEGA"};

            string[] I3VarList2D = {"PS"};
            string[] I3VarList3D = {"QV","T"};
            
            string[] A3CldVarList2D = {};
            string[] A3CldVarList3D = {"CLOUD","QI","QL"};
            
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
            currentFile = new MetFile(currentTemplate,startDate,A3DynVarList2D,A3DynVarList3D,lonLims,latLims,A3Offset);
            MetFiles.Add(currentFile);

            A3DynIndex = MetFiles.Count() - 1;
            UIndex      = Array.FindIndex(A3DynVarList3D,element => element == "U");
            VIndex      = Array.FindIndex(A3DynVarList3D,element => element == "V");
            OmegaIndex  = Array.FindIndex(A3DynVarList3D,element => element == "OMEGA");
            
            // 3-hour averaged, cloud
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3cld.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,A3CldVarList2D,A3CldVarList3D,lonLims,latLims,A3Offset);
            MetFiles.Add(currentFile);

            A3CldIndex = MetFiles.Count() - 1;
            QIIndex = Array.FindIndex(A3CldVarList3D,element => element == "QI");
            QLIndex = Array.FindIndex(A3CldVarList3D,element => element == "QL");
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