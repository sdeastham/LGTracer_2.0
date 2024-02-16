namespace LGTracer
{
    public abstract class Integrator
    {
        public double State
        { get; protected set; }

        public double TimeStep
        { get; protected set; }

        public double CurrentTime
        { get; protected set; }

        public int StepCount
        { get; protected set; }

        public Func<double, double> RateCalc
        { get; private set; }

        public string IntMethod
        { get; protected set; }

        public Integrator(double state, double timeStep, Func<double, double> rateCalc)
        {
            this.State = state;
            this.TimeStep = timeStep;
            this.CurrentTime = 0.0;
            this.StepCount = 0;
            this.RateCalc = rateCalc;
            this.IntMethod = "Undefined";
        }

        public abstract void Advance();

        public override string ToString()
        {
            return $"{IntMethod} integrator";
        }
    }

    public class FwdEuler : Integrator
    {
        public FwdEuler(double state, double timeStep, Func<double, double> rateCalc) : base(state, timeStep, rateCalc)
        {
            this.IntMethod = "Forward Euler";
        }

        public override void Advance()
        {
            CurrentTime += TimeStep;
            State += RateCalc(State) * TimeStep;
        }
    }

    public class BwdEuler : Integrator
    {
        public BwdEuler(double state, double timeStep, Func<double, double> rateCalc) : base(state, timeStep, rateCalc)
        {
            this.IntMethod = "Backward Euler";
        }

        public override void Advance()
        {
            CurrentTime += TimeStep;
            State += RateCalc(State) * TimeStep;
        }
    }

    public class RK4 : Integrator
    {
        public RK4(double state, double timeStep, Func<double, double> rateCalc) : base(state, timeStep, rateCalc)
        {
            this.IntMethod = "Runge-Kutta 4th-Order";
        }

        public override void Advance()
        {
            CurrentTime += TimeStep;
            double k1 = RateCalc(State);
            double k2 = RateCalc(State + (TimeStep*k1/2.0));
            double k3 = RateCalc(State + (TimeStep*k2/2.0));
            double k4 = RateCalc(State + (TimeStep*k3));
            State += (TimeStep/6.0)*(k1+2.0*k2+2.0*k3+k4);
        }
    }
}