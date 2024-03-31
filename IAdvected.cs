namespace LGTracer;

public interface IAdvected
{
    public double GetProperty(string property);
    public void Activate( double x, double y, double pressure, uint uniqueID );
    public void Deactivate();
    public void Advance(double dt);

    public (double,double,double) GetLocation();
    public uint GetUID();
    public double GetAge();
}