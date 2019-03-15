using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FiniteElementMethod;
using FiniteElementMethod.Matrices;
using FiniteElementMethod.Additions;

namespace MainProgram
{
    class Program
    {
        static void Main(string[] args)
        {
            Coordinate2D pointOnX = new Coordinate2D { X = 2, Y = 0 };
            Coordinate2D pointOnY = new Coordinate2D { X = 0, Y = 3 };

            int nPartOnX = 1, mPartOnY = 2;
            CoordinatePlane coordinatePlane = new CoordinatePlane(pointOnX, pointOnY);

            coordinatePlane.SplitCoordinatePlane(nPartOnX, mPartOnY);

            coordinatePlane.NodesMatrix.Show();

            coordinatePlane.ConnectivityMatrix.Show();

            SpecialData.ShowArrayWithIndices(coordinatePlane.BoundaryConditionMatrix);

            Console.ReadLine();
        }
    }
}
