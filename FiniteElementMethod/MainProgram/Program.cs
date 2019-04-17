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

            //coordinatePlane.CreateDifferentialOperatorMatrixCl(0, 1, 0).Show();

            Matrix matrix = coordinatePlane.CreateGlobalMatrix(0);//.Show();

            Matrix fromFile = new Matrix(48, 48);

            string[] data = null;
            using (StreamReader streamReader = new StreamReader(@"..\..\..\K48.txt")) 
            {
                data = streamReader.ReadToEnd().Split(new char[] { ' ','\n','\r' }, StringSplitOptions.RemoveEmptyEntries);
            }

            int k = 0;
            for (int i = 0; i < 48; i++)
            {
                for (int j = 0; j < 48; j++)
                {
                    fromFile[i, j] = double.Parse(data[k++], System.Globalization.NumberStyles.Float, CultureInfo.InvariantCulture);
                }
            }

            bool isEqual = true;
            for (int i = 0; i < 48; i++)
            {
                for (int j = 0; j < 48; j++)
                {
                    if (Math.Round(matrix[i, j]) != fromFile[i, j]) 
                    {
                        isEqual = false;
                        Console.WriteLine($"{i} {j}: {matrix[i, j]} {fromFile[i, j]}");
                    }
                }
            }

            Console.WriteLine(isEqual);

            Console.ReadLine();
        }
    }
}
