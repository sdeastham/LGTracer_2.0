using System.Numerics;
using System.Runtime.Intrinsics.X86;

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

    protected bool DefaultValidity;

    public LGPoint( Func<double, double, double, (double, double, double)> vCalc)
    {
        // Point starts inactive
        _location = new Vector3(float.NaN,float.NaN,float.NaN);
        InitialLocation = new Vector3(float.NaN,float.NaN,float.NaN);
        VelocityCalc = vCalc;
        DefaultValidity = true;
        // Set the rest of the properties by deactivating the point
        Deactivate();
    }

    public virtual void Activate( double x, double y, double pressure, uint uniqueID )
    {
        // Change the particle from being inactive to active
        Active = true;
        _location = new Vector3((float)x, (float)y, (float)pressure);
        // Copy this data for later comparison
        InitialLocation = new Vector3((float)x, (float)y, (float)pressure);
        UID = uniqueID;
        Age = 0.0;
    }

    public virtual void Deactivate()
    {
        // Kill the particle and put it into storage
        Active = false;
        _location = new Vector3(float.NaN,float.NaN,float.NaN);
        UID = 0;
        Age = double.NaN;
    }

    public virtual double GetProperty(string property)
    {
        throw new ArgumentException($"Property {property} requested but base LGPoints have no properties");
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
}