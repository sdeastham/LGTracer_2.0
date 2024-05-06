using System.Diagnostics;

namespace LGTracer;

public class MetManager
{
    protected List<MetFile> MetFiles;
    public double[,,] UWindXYP => ((MetData3D)MetFiles[UFileIndex].GetMetData(UIndex)).CurrentData;

    public double[,,] VWindXYP => ((MetData3D)MetFiles[VFileIndex].GetMetData(VIndex)).CurrentData;

    public double[,,] PressureVelocityXYP => ((MetData3D)MetFiles[OmegaFileIndex].GetMetData(OmegaIndex)).CurrentData;

    public double[,] SurfacePressureXY => ((MetData2D)MetFiles[PSFileIndex].GetMetData(PSIndex)).CurrentData;

    public double[,,] SpecificHumidityXYP => ((MetData3D)MetFiles[QVFileIndex].GetMetData(QVIndex)).CurrentData; // kg water vapour per kg air

    public double[,,] CloudIceXYP => ((MetData3D)MetFiles[QIFileIndex].GetMetData(QIIndex)).CurrentData; // kg ice water per kg air
    public double[,,] CloudWaterXYP => ((MetData3D)MetFiles[QLFileIndex].GetMetData(QLIndex)).CurrentData; // kg liquid water per kg air
    public double[,,] TemperatureXYP => ((MetData3D)MetFiles[TFileIndex].GetMetData(TIndex)).CurrentData;

    // Locations in arrays
    //protected int I3Index, A3DynIndex, A3CldIndex;
    protected int PSIndex, TIndex, QVIndex, QIIndex, QLIndex;
    protected int UIndex, VIndex, OmegaIndex;
    protected int PSFileIndex, TFileIndex, QVFileIndex, QIFileIndex, QLFileIndex;
    protected int UFileIndex, VFileIndex, OmegaFileIndex;

    private Dictionary<string, Stopwatch> Stopwatches;
        
    public MetManager(string metDir, double[] lonLims, double[] latLims, DateTime startDate, bool useSerial, Dictionary<string, Stopwatch> stopwatches, string dataSource)
    {
        MetFiles = [];
        Stopwatches = stopwatches;

        if (dataSource == "MERRA-2")
        {
            // Read time offsets in seconds - A3 files are timestamped for
            // mid-way through the 3 hour averaging period, but we want to read
            // them at the start of their averaging period
            int A3Offset = -90 * 60;
            int I3Offset = 0;

            string[] A3DynVarList2D = { };
            string[] A3DynVarList3D = { "U", "V", "OMEGA" };

            string[] I3VarList2D = { "PS" };
            string[] I3VarList3D = { "QV", "T" };

            string[] A3CldVarList2D = { };
            string[] A3CldVarList3D = { "CLOUD", "QI", "QL" };

            MetFile currentFile;
            string currentTemplate;

            // 3-hour instantaneous
            if (useSerial)
            {
                currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.I3.{3}.05x0625.serial");
            }
            else
            {
                currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.I3.05x0625.nc4");
            }

            currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate, startDate, I3VarList2D, I3VarList3D,
                lonLims, latLims, Stopwatches, I3Offset, timeInterp: true, useSerial: useSerial);
            MetFiles.Add(currentFile);

            // Set up connections
            int I3Index = MetFiles.Count() - 1;
            PSIndex = currentFile.DataNames.FindIndex(element => element == "PS");
            QVIndex = currentFile.DataNames.FindIndex(element => element == "QV");
            TIndex = currentFile.DataNames.FindIndex(element => element == "T");
            PSFileIndex = I3Index;
            QVFileIndex = I3Index;
            TFileIndex = I3Index;

            // 3-hour averaged, dynamics
            if (useSerial)
            {
                currentTemplate = Path.Combine(metDir,
                    "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.{3}.05x0625.serial");
            }
            else
            {
                currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3dyn.05x0625.nc4");
            }

            currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate, startDate, A3DynVarList2D,
                A3DynVarList3D,
                lonLims, latLims, Stopwatches, A3Offset, timeInterp: false, useSerial: useSerial);
            MetFiles.Add(currentFile);

            int A3DynIndex = MetFiles.Count() - 1;
            UIndex = currentFile.DataNames.FindIndex(element => element == "U");
            VIndex = currentFile.DataNames.FindIndex(element => element == "V");
            OmegaIndex = currentFile.DataNames.FindIndex(element => element == "OMEGA");
            UFileIndex = A3DynIndex;
            VFileIndex = A3DynIndex;
            OmegaFileIndex = A3DynIndex;

            // 3-hour averaged, cloud
            if (useSerial)
            {
                currentTemplate = Path.Combine(metDir,
                    "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3cld.{3}.05x0625.serial");
            }
            else
            {
                currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/MERRA2.{0}{1,2:d2}{2,2:d2}.A3cld.05x0625.nc4");
            }

            currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate, startDate, A3CldVarList2D,
                A3CldVarList3D,
                lonLims, latLims, Stopwatches, A3Offset, timeInterp: false, useSerial: useSerial);
            MetFiles.Add(currentFile);

            int A3CldIndex = MetFiles.Count() - 1;
            QIIndex = currentFile.GetVarIndex("QI"); //DataNames.FindIndex(element => element == "QI"));
            QLIndex = currentFile.GetVarIndex("QL"); //DataNames.FindIndex(element => element == "QL"));
            QIFileIndex = A3CldIndex;
            QLFileIndex = A3CldIndex;
        }
        else if (dataSource == "ERA5")
        {
            // ERA5 on fixed pressure levels
            int timeOffset = 0;

            // Only two files
            string[] varList2D = { "sp" };
            string[] varList3D = { "q", "t", "u", "v", "w", "ciwc", "clwc" };

            MetFile currentFile;
            string currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/ERA5_surface_{0}{1,2:d2}{2,2:d2}.nc");

            currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate, startDate, varList2D, [],
                lonLims, latLims, Stopwatches, timeOffset, timeInterp: false, useSerial: false);
            MetFiles.Add(currentFile);

            int fileIndex = MetFiles.Count() - 1;
            PSIndex = currentFile.DataNames.FindIndex(element => element == "sp");
            PSFileIndex = fileIndex;
            
            currentTemplate = Path.Combine(metDir, "{0}/{1,2:d2}/ERA5_plevs_{0}{1,2:d2}{2,2:d2}.nc");
            currentFile = (MetFile)MetFileFactory.CreateMetFile(currentTemplate, startDate, [], varList3D,
                lonLims, latLims, Stopwatches, timeOffset, timeInterp: false, useSerial: false);
            MetFiles.Add(currentFile);
            
            // Set up connections
            fileIndex = MetFiles.Count() - 1;
            QVIndex = currentFile.DataNames.FindIndex(element => element == "q");
            TIndex = currentFile.DataNames.FindIndex(element => element == "t");
            UIndex = currentFile.DataNames.FindIndex(element => element == "u");
            VIndex = currentFile.DataNames.FindIndex(element => element == "v");
            OmegaIndex = currentFile.DataNames.FindIndex(element => element == "w");
            QIIndex = currentFile.GetVarIndex("ciwc");
            QLIndex = currentFile.GetVarIndex("clwc");
            QVFileIndex = fileIndex;
            TFileIndex = fileIndex;
            UFileIndex = fileIndex;
            VFileIndex = fileIndex;
            OmegaFileIndex = fileIndex;
            QIFileIndex = fileIndex;
            QLFileIndex = fileIndex;
        }
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
