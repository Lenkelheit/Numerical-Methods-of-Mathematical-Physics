using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using FiniteElementMethod;
using FiniteElementMethod.Matrices;
using FiniteElementMethod.Additions;
using System.Globalization;

namespace MainProgram
{
    class Program
    {
        static void Main(string[] args)
        {
            double a = 0, b = 0.4, c = 0, d = Math.PI / 25;

            int nPartOnX = 4, mPartOnY = 1;

            CoordinateSystemConfig coordinateSystemConfig = new CoordinateSystemConfig(a, b, c, d, nPartOnX, mPartOnY);

            double R = 0.1, E = 80.1e10, v = 0.3, h = 0.005;

            CylindricalShellConfig cylindricalShellConfig = new CylindricalShellConfig(E, v, h, R);

            CoordinatePlane coordinatePlane = new CoordinatePlane(coordinateSystemConfig, cylindricalShellConfig);

            coordinatePlane.NodesMatrix.Show();

            coordinatePlane.ConnectivityMatrix.Show();

            SpecialData.ShowArrayWithIndices(coordinatePlane.BoundaryConditionMatrix);

            Console.WriteLine();

            //coordinatePlane.CreateGlobalMatrix(0).Show();


            Console.ReadLine();
        }
    }
}
