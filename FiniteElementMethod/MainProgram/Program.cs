using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FiniteElementMethod;
using FiniteElementMethod.Additions;

namespace MainProgram
{
    class Program
    {
        static void Main(string[] args)
        {
            Coordinate2D basePoint = new Coordinate2D { X = 0, Y = 0 };
            Coordinate2D pointOnX = new Coordinate2D { X = 3, Y = 0 };
            Coordinate2D pointOnY = new Coordinate2D { X = 0, Y = 4 };

            int nPartOnX = 3, mPartOnY = 4;
            CoordinatePlane coordinatePlane = new CoordinatePlane(basePoint, pointOnX, pointOnY);

            coordinatePlane.SplitCoordinatePlane(nPartOnX, mPartOnY);

            coordinatePlane.NodesMatrix.Show();

            Console.ReadLine();
        }
    }
}
