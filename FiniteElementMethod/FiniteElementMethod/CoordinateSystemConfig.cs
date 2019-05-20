namespace FiniteElementMethod
{
    public class CoordinateSystemConfig
    {
        // FIELDS
        private double a, b;
        private double c, d;

        private int n, m;

        // PEOPERTIES
        public double A => a;
        public double B => b;
        public double C => c;
        public double D => d;

        public int N => n;
        public int M => m;

        // CONSTRUCTORS
        public CoordinateSystemConfig(double a, double b, double c, double d, int n, int m)
        {
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;

            this.n = n;
            this.m = m;
        }
    }
}
