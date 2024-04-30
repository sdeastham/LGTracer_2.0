namespace LGTracer;

public interface IAdvected
{
    public double GetProperty(string property);
    public void Activate( double x, double y, double pressure, uint uniqueID, DateTime initiationDate );
    public void Deactivate();
    public void Advance(double dt, DomainManager domain);

    public (double,double,double) GetLocation();
    public uint GetUID();
    public double GetAge();
    public bool CheckValid();
    public void ArchiveConditions(); // Append current information to history
    public void SetupHistory(IEnumerable<string> propertyNames);
    public Dictionary<string, List<double>> GetHistory();
}