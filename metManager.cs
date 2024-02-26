namespace LGTracer
{
    public class MetManager
    {
        // A MetManager handles all met data and exposes the met variables from all files

        protected List<MetFile> MetFiles;

        public double[,] SurfacePressureXY
        { get { return MetFiles[I3Index].DataVariables2D[PSIndex].CurrentData; } }
        //double[,,] = UWindXYP;
        //double[,,] = VWindXYP;
        private int I3Index;
        private int PSIndex;

        public MetManager(string metDir, double[] lonLims, double[] latLims, DateTime startDate)
        {
            MetFiles = [];

            string[] A3DynVarList2D = {};
            string[] A3DynVarList3D = {"U","V"};
            int A3DynOffset = -30 * 60;

            string[] I3VarList2D = {"PS"};
            string[] I3VarList3D = {"QV"};
            int I3Offset = 0;
            MetFile currentFile;
            string currentTemplate;
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.I3.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,I3VarList2D,I3VarList3D,lonLims,latLims,I3Offset);
            MetFiles.Add(currentFile);

            // Set up connections
            I3Index = MetFiles.Count() - 1;
            PSIndex = Array.FindIndex(I3VarList2D,element => element == "PS");

            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,A3DynVarList2D,A3DynVarList3D,lonLims,latLims,A3DynOffset);
            MetFiles.Add(currentFile);
        }

        public void AdvanceToTime(DateTime targetTime)
        {
            foreach (MetFile metFile in MetFiles)
            {
                metFile.AdvanceToTime(targetTime);
            }
        }
    }
}