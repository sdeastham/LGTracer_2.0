using System.Numerics;
using System.Runtime.Intrinsics.X86;
using Apache.Arrow;
using Parquet;
using Parquet.Data;
using Parquet.Schema;
using Parquet.Serialization;
using Field = Parquet.Schema.Field;

namespace LGTracer;

public class LGPoint : IAdvected
{
    protected Vector3 _location
    { get; set; }

    public Vector3 InitialLocation
    { get; protected set; }

    public Func<double, double, double, (double, double, double)> VelocityCalc
    { get; protected set; }

    public bool Active
    { get; protected set; }

    public uint UID
    { get; protected set; }
        
    public double Age
    { get; protected set; }

    public double GetAge()
    {
        return Age;
    }

    private static readonly DateTime ReferenceTime = new DateTime(1970, 1, 1, 0, 0, 0);

    // Convenience properties
    public double X
    { 
        get => _location.X;
        set => _location = new Vector3((float)value,_location.Y,_location.Z);
    }
    public double Y
    { 
        get => _location.Y;
        set => _location = new Vector3(_location.X,(float)value,_location.Z);
    }

    public double Pressure
    { 
        get => _location.Z;
        set => _location = new Vector3(_location.X,_location.Y,(float)value);
    }
    
    public (double, double, double) GetLocation()
    {
        return( X, Y, Pressure );
    }

    public uint GetUID()
    {
        return UID;
    }

    // Facilitates some slightly screwy logic for points which can become invalid (classes inheriting from this one)
    protected bool DefaultValidity;

    public Dictionary<string, List<double>> History
    { get; protected set; }

    public DateTime InitiationDate;

    public LGPoint( Func<double, double, double, (double, double, double)> vCalc)
    {
        // Point starts inactive
        _location = new Vector3(float.NaN,float.NaN,float.NaN);
        InitialLocation = new Vector3(float.NaN,float.NaN,float.NaN);
        VelocityCalc = vCalc;
        DefaultValidity = true;
        History = [];
        // Set the rest of the properties by deactivating the point
        Deactivate();
    }

    public virtual void Activate( double x, double y, double pressure, uint uniqueID, DateTime initiationDate)
    {
        // Change the particle from being inactive to active
        Active = true;
        _location = new Vector3((float)x, (float)y, (float)pressure);
        // Copy this data for later comparison
        InitialLocation = new Vector3((float)x, (float)y, (float)pressure);
        UID = uniqueID;
        Age = 0.0;
        InitiationDate = initiationDate;
        ArchiveConditions(true);
    }

    public virtual void Deactivate()
    {
        // Kill the particle and put it into storage
        Active = false;
        _location = new Vector3(float.NaN,float.NaN,float.NaN);
        UID = 0;
        Age = double.NaN;
        // Clear history
        foreach (string property in History.Keys)
        {
            History[property].Clear();
        }
    }

    public virtual double GetProperty(string property)
    {
        switch (property.ToLower().Replace("_",""))
        {
            case "longitude":
            case "lon":
            case "x":
                return X;
            case "latitude":
            case "lat":
            case "y":
                return Y;
            case "pressure":
            case "p":
                return Pressure;
            case "age":
                return Age;
            case "time":
                return (InitiationDate + TimeSpan.FromSeconds(Age) - ReferenceTime).TotalSeconds;
            case "uid":
                return UID; // Casting from ulong to double..
            default:
                throw new ArgumentException($"No property for LGPoint called {property}");
        }
    }

    public virtual void Advance( double dt, DomainManager domain )
    {
        if (!Active)
        {
            return;
        }

        // RK4
        ( double dx1, double dy1, double dp1 ) = VelocityCalc( X, Y, Pressure );
        ( double dx2, double dy2, double dp2 ) = VelocityCalc( X + dx1*dt/2.0, Y + dy1*dt/2.0, Pressure + dp1*dt/2.0);
        ( double dx3, double dy3, double dp3 ) = VelocityCalc( X + dx2*dt/2.0, Y + dy2*dt/2.0, Pressure + dp2*dt/2.0);
        ( double dx4, double dy4, double dp4 ) = VelocityCalc( X + dx3*dt, Y + dy3*dt, Pressure + dp3*dt);
        double xSpeed = (1.0/6.0) * (dx1 + 2.0*dx2 + 2.0*dx3 + dx4);
        double ySpeed = (1.0/6.0) * (dy1 + 2.0*dy2 + 2.0*dy3 + dy4);
        double pSpeed = (1.0/6.0) * (dp1 + 2.0*dp2 + 2.0*dp3 + dp4);
        X += xSpeed * dt;
        Y += ySpeed * dt;
        Pressure += pSpeed * dt;
        Age += dt;
    }

    public virtual bool CheckValid()
    {
        return DefaultValidity;
    }

    public void SetupHistory(IEnumerable<string> propertyNames)
    {
        // First add the default properties
        foreach (string property in (string[])["longitude","latitude","pressure","age"])
        {
            History.Add(property, []);
        }
        // Now the non-default ones
        foreach (string property in propertyNames.Where(property => !History.ContainsKey(property)))
        {
            History.Add(property, []);
        }
    }

    public void ArchiveConditions(bool clear=false)
    {
        if (clear)
        {
            foreach (string property in History.Keys)
            {
                History[property].Clear();
            }
        }

        // Store information
        foreach (string property in History.Keys)
        {
            History[property].Add(GetProperty(property));
        }
    }

    public Dictionary<string, List<double>> GetHistory()
    {
        return History;
    }
}