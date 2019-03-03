using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FiniteElementMethod.Additions;

namespace FiniteElementMethod
{
    public class CoordinatePlane
    {
        // PROPERTIES
        public Coordinate2D BasePoint { get; }
        public Coordinate2D PointOnX { get; }
        public Coordinate2D PointOnY { get; }
        public NodesMatrix NodesMatrix { get; private set; }

        // CONSTRUCTORS
        public CoordinatePlane(Coordinate2D basePoint, Coordinate2D pointOnX, Coordinate2D pointOnY)
        {
            this.BasePoint = basePoint;
            this.PointOnX = pointOnX;
            this.PointOnY = pointOnY;
        }

        // METHODS
        public void SplitCoordinatePlane(int nPartOnX, int mPartOnY)
        {
            // finds steps on X and Y
            double hOnX = (PointOnX.X - BasePoint.X) / nPartOnX;
            double hOnY = (PointOnY.Y - BasePoint.Y) / mPartOnY;

            double stepCoef = 0.5;
            hOnY *= stepCoef;
            hOnX *= stepCoef;

            NodesMatrix = new NodesMatrix(nPartOnX, mPartOnY);

            int mAmount = 0;
            // on M plane
            for (double i = BasePoint.Y; i <= PointOnY.Y; i += hOnY)  
            {
                int nAmount = 0;
                // on N plane
                for (double j = BasePoint.X; j <= PointOnX.X; j += hOnX) 
                {
                    NodesMatrix[mAmount, nAmount] = new Coordinate2D { X = j, Y = i };
                    ++nAmount;
                }
                stepCoef = Math.Pow(stepCoef, -1);
                hOnX *= stepCoef;
                ++mAmount;
            }
        }
    }
}
