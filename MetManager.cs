namespace LGTracer;

public class MetManager
{
    protected List<MetFile> MetFiles;
    public double[,,] UWindXYP => ((MetData3D)MetFiles[A3DynIndex].DataVariables[UIndex]).CurrentData;

    public double[,,] VWindXYP => ((MetData3D)MetFiles[A3DynIndex].DataVariables[VIndex]).CurrentData;

    public double[,,] PressureVelocityXYP => ((MetData3D)MetFiles[A3DynIndex].DataVariables[OmegaIndex]).CurrentData;

    public double[,] SurfacePressureXY => ((MetData2D)MetFiles[I3Index].DataVariables[PSIndex]).CurrentData;

    public double[,,] SpecificHumidityXYP => ((MetData3D)MetFiles[I3Index].DataVariables[QVIndex]).CurrentData; // kg water vapour per kg air

    public double[,,] CloudIceXYP => ((MetData3D)MetFiles[A3CldIndex].DataVariables[QIIndex]).CurrentData; // kg ice water per kg air
    public double[,,] CloudWaterXYP => ((MetData3D)MetFiles[A3CldIndex].DataVariables[QLIndex]).CurrentData; // kg liquid water per kg air
    public double[,,] TemperatureXYP => ((MetData3D)MetFiles[I3Index].DataVariables[TIndex]).CurrentData;

    // Locations in arrays
    protected int I3Index, A3DynIndex, A3CldIndex;
    protected int PSIndex, TIndex, QVIndex, QIIndex, QLIndex;
    protected int UIndex, VIndex, OmegaIndex;
        
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
            currentFile = new MetFile(currentTemplate,startDate,I3VarList2D,I3VarList3D,lonLims,latLims,I3Offset,timeInterp:true);
            MetFiles.Add(currentFile);

            // Set up connections
            I3Index = MetFiles.Count() - 1;
            PSIndex = currentFile.DataNames.FindIndex(element => element == "PS");
            QVIndex = currentFile.DataNames.FindIndex(element => element == "QV");
            TIndex  = currentFile.DataNames.FindIndex(element => element == "T");

            // 3-hour averaged, dynamics
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,A3DynVarList2D,A3DynVarList3D,lonLims,latLims,A3Offset,timeInterp: false);
            MetFiles.Add(currentFile);

            A3DynIndex = MetFiles.Count() - 1;
            UIndex      = currentFile.DataNames.FindIndex(element => element == "U");
            VIndex      = currentFile.DataNames.FindIndex(element => element == "V");
            OmegaIndex  = currentFile.DataNames.FindIndex(element => element == "OMEGA");
            
            // 3-hour averaged, cloud
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3cld.05x0625.nc4");
            currentFile = new MetFile(currentTemplate,startDate,A3CldVarList2D,A3CldVarList3D,lonLims,latLims,A3Offset,timeInterp: false);
            MetFiles.Add(currentFile);

            A3CldIndex = MetFiles.Count() - 1;
            QIIndex = currentFile.DataNames.FindIndex(element => element == "QI");
            QLIndex = currentFile.DataNames.FindIndex(element => element == "QL");
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