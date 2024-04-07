using System.Diagnostics;

namespace LGTracer;

public class MetManager
{
    protected List<MetFile> MetFiles;
    public double[,,] UWindXYP => ((MetData3D)MetFiles[A3DynIndex].GetMetData(UIndex)).CurrentData;

    public double[,,] VWindXYP => ((MetData3D)MetFiles[A3DynIndex].GetMetData(VIndex)).CurrentData;

    public double[,,] PressureVelocityXYP => ((MetData3D)MetFiles[A3DynIndex].GetMetData(OmegaIndex)).CurrentData;

    public double[,] SurfacePressureXY => ((MetData2D)MetFiles[I3Index].GetMetData(PSIndex)).CurrentData;

    public double[,,] SpecificHumidityXYP => ((MetData3D)MetFiles[I3Index].GetMetData(QVIndex)).CurrentData; // kg water vapour per kg air

    public double[,,] CloudIceXYP => ((MetData3D)MetFiles[A3CldIndex].GetMetData(QIIndex)).CurrentData; // kg ice water per kg air
    public double[,,] CloudWaterXYP => ((MetData3D)MetFiles[A3CldIndex].GetMetData(QLIndex)).CurrentData; // kg liquid water per kg air
    public double[,,] TemperatureXYP => ((MetData3D)MetFiles[I3Index].GetMetData(TIndex)).CurrentData;

    // Locations in arrays
    protected int I3Index, A3DynIndex, A3CldIndex;
    protected int PSIndex, TIndex, QVIndex, QIIndex, QLIndex;
    protected int UIndex, VIndex, OmegaIndex;

    private Dictionary<string, Stopwatch> Stopwatches;
        
    public MetManager(string metDir, double[] lonLims, double[] latLims, DateTime startDate, Dictionary<string, Stopwatch> stopwatches)
    {
        MetFiles = [];
        Stopwatches = stopwatches;
        
        // Read time offsets in seconds - A3 files are timestamped for
        // mid-way through the 3 hour averaging period, but we want to read
        // them at the start of their averaging period
        int A3Offset = -90 * 60;
        int I3Offset = 0;

        string[] A3DynVarList2D = {};
        string[] A3DynVarList3D = {"U","V","OMEGA"};

        string[] I3VarList2D = {"PS"};
        string[] I3VarList3D = {"QV","T"};
        
        string[] A3CldVarList2D = {};
        string[] A3CldVarList3D = {"CLOUD","QI","QL"};
        
        MetFile currentFile;
        string currentTemplate;
        bool useSerial = true;

        // 3-hour instantaneous
        if (useSerial)
        {
            currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.I3.{3}.05x0625.nc4");
        }
        else
        {
            currentTemplate = Path.Combine(metDir,"{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.I3.05x0625.nc4");
        }
        currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate,startDate,I3VarList2D,I3VarList3D,
            lonLims,latLims,Stopwatches,I3Offset,timeInterp: true,useSerial: useSerial);
        MetFiles.Add(currentFile);

        // Set up connections
        I3Index = MetFiles.Count() - 1;
        PSIndex = currentFile.DataNames.FindIndex(element => element == "PS");
        QVIndex = currentFile.DataNames.FindIndex(element => element == "QV");
        TIndex  = currentFile.DataNames.FindIndex(element => element == "T");

        // 3-hour averaged, dynamics
        if (useSerial)
        {
            currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.{3}.05x0625.nc4");
        }
        else
        {
            currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.05x0625.nc4");
        }
        currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate,startDate,A3DynVarList2D,A3DynVarList3D,
            lonLims,latLims,Stopwatches,A3Offset,timeInterp: false,useSerial: useSerial);
        MetFiles.Add(currentFile);

        A3DynIndex = MetFiles.Count() - 1;
        UIndex      = currentFile.DataNames.FindIndex(element => element == "U");
        VIndex      = currentFile.DataNames.FindIndex(element => element == "V");
        OmegaIndex  = currentFile.DataNames.FindIndex(element => element == "OMEGA");
        
        // 3-hour averaged, cloud
        if (useSerial)
        {
            currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3cld.{3}.05x0625.nc4");
        }
        else
        {
            currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3cld.05x0625.nc4");
        }
        currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate,startDate,A3CldVarList2D,A3CldVarList3D,
            lonLims,latLims,Stopwatches,A3Offset,timeInterp: false,useSerial: useSerial);
        MetFiles.Add(currentFile);

        A3CldIndex = MetFiles.Count() - 1;
        QIIndex = currentFile.GetVarIndex("QI");//DataNames.FindIndex(element => element == "QI"));
        QLIndex = currentFile.GetVarIndex("QL");//DataNames.FindIndex(element => element == "QL"));
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