using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethod
{
    public class CylindricalShellConfig
    {
        // FIELDS
        private double e;// модуль Юнга, пружність
        private double v;// коефіцієнт Пуасона, міра деформації
        private double h;// товщина оболонки в метрах
        private double r;// радіус оболонки в метрах

        // PROPERTIES
        public double E => e;
        public double V => v;
        public double H => h;
        public double R => r;

        // CONSTRUCTORS
        public CylindricalShellConfig(double e, double v, double h, double r)
        {
            this.e = e;
            this.v = v;
            this.h = h;
            this.r = r;
        }
    }
}
