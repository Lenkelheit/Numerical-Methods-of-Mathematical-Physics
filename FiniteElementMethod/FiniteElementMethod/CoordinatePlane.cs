using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FiniteElementMethod.Matrices;
using FiniteElementMethod.Additions;

namespace FiniteElementMethod
{
    public class CoordinatePlane
    {
        // PROPERTIES
        #region POINTS
        public Coordinate2D PointOnX { get; }
        public Coordinate2D PointOnY { get; }
        #endregion

        #region MATRICES
        public NodesMatrix NodesMatrix { get; private set; }
        public Matrix ConnectivityMatrix { get; private set; }
        public int[] BoundaryConditionMatrix { get; private set; }
        #endregion

        // CONSTRUCTORS
        public CoordinatePlane(Coordinate2D pointOnX, Coordinate2D pointOnY)
        {
            this.PointOnX = pointOnX;
            this.PointOnY = pointOnY;
        }

        // METHODS
        public void SplitCoordinatePlane(int nPartOnX, int mPartOnY)
        {
            CreateNodesMatrix(nPartOnX, mPartOnY);

            CreateConnectivityMatrix(nPartOnX, mPartOnY);

            CreateBoundaryConditionMatrix(nPartOnX, mPartOnY);
        }

        private void CreateNodesMatrix(int nPartOnX, int mPartOnY)
        {
            double stepCoef = 0.5;

            // finds steps on X and Y
            double hOnX = (PointOnX.X - PointOnY.X) / nPartOnX * stepCoef;
            double hOnY = (PointOnY.Y - PointOnX.Y) / mPartOnY * stepCoef;

            NodesMatrix = nPartOnX < mPartOnY ? new NodesMatrix(mPartOnY, nPartOnX) : new NodesMatrix(nPartOnX, mPartOnY);

            int mAmount = 0;
            // fills matrix of nodes
            // on M plane
            for (double i = nPartOnX < mPartOnY ? PointOnX.Y : PointOnY.X; i <= (nPartOnX < mPartOnY ? PointOnY.Y : PointOnX.X); i += nPartOnX < mPartOnY ? hOnY : hOnX)
            {
                int nAmount = 0;
                // on N plane
                for (double j = nPartOnX < mPartOnY ? PointOnY.X : PointOnX.Y; j <= (nPartOnX < mPartOnY ? PointOnX.X : PointOnY.Y); j += nPartOnX < mPartOnY ? hOnX : hOnY)
                {
                    NodesMatrix[mAmount, nAmount] = nPartOnX < mPartOnY ? new Coordinate2D { X = j, Y = i } : new Coordinate2D { X = i, Y = j };

                    ++nAmount;
                }
                // turns over stepCoef from 1/2 to 2 or backward
                stepCoef = Math.Pow(stepCoef, -1);
                if (nPartOnX < mPartOnY)
                {
                    hOnX *= stepCoef;
                }
                else
                {
                    hOnY *= stepCoef;
                }
                ++mAmount;
            }
        }

        private void CreateConnectivityMatrix(int nPartOnX, int mPartOnY)
        {
            int finiteElementsAmount = nPartOnX * mPartOnY;
            ConnectivityMatrix = new Matrix(finiteElementsAmount, SpecialData.NODES_NUMBER_IN_FINITE_ELEMENT);

            // onMFirst - first row index, onMSecond - second row index, onMThird - third row index
            int onMFirst = 0, onMSecond = onMFirst + 1, onMThird = onMSecond + 1;

            // onNFirst - column index on first row, onNSecond - column index on second row, onNThird - column index on third row
            int onNFirst = 0, onNSecond = 0, onNThird = 0;


            // fills connectivity matrix
            for (int i = 0; i < ConnectivityMatrix.RowsAmount; ++i)
            {
                int currentColumnInConnectMatrix = nPartOnX < mPartOnY ? 0 : SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT - 1;

                // goes to other finite element
                // skips 1 column
                onNThird += SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT - 1;
                // goes to the next column
                ++onNSecond;

                // executes when onNFirst is the last column index and then goes to the row with next finite element
                if (onNFirst >= NodesMatrix.EvenNodesOnMAmount - 1)
                {
                    onMFirst += SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT - 1;
                    onMSecond = onMFirst + 1;
                    onMThird = onMSecond + 1;

                    onNFirst = 0;

                    // skips 1 column
                    onNThird = SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT - 1;

                    // goes to the next column
                    onNSecond = 1;
                }
                // fills row in ConnectivityMatrix
                for (int j = 0; j < SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT; ++j)
                {
                    // fills from M first row in such way: --> if nPartOnX < mPartOnY,
                    // else fills from M third row in such way: -->
                    ConnectivityMatrix[i, currentColumnInConnectMatrix] = NodesMatrix.CreateGlobalIndex(nPartOnX < mPartOnY ? onMFirst : onMThird, onNFirst + j);

                    if (nPartOnX < mPartOnY)
                    {
                        // fills from M third row in such way: <--
                        ConnectivityMatrix[i, currentColumnInConnectMatrix + SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT + 1] = NodesMatrix.CreateGlobalIndex(onMThird, onNThird - j);
                    }
                    else
                    {
                        // fills from M first row in such way: <-- , where the last index j must be set in 0
                        ConnectivityMatrix[i, j == SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT - 1 ?
                            currentColumnInConnectMatrix = 0 : currentColumnInConnectMatrix + SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT + 1] = NodesMatrix.CreateGlobalIndex(onMFirst, onNThird - j);
                    }

                    ++currentColumnInConnectMatrix;
                }
                // fills from M second row
                // It is middle column in ConnectivityMatrix if nPartOnX < mPartOnY,
                // else it is second column in ConnectivityMatrix
                ConnectivityMatrix[i, currentColumnInConnectMatrix] = NodesMatrix.CreateGlobalIndex(onMSecond, nPartOnX < mPartOnY ? onNSecond : onNSecond - 1);

                // It is last column in ConnectivityMatrix if nPartOnX < mPartOnY,
                // else it is third from back column in ConnectivityMatrix
                ConnectivityMatrix[i, currentColumnInConnectMatrix + SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT + 1] = NodesMatrix.CreateGlobalIndex(onMSecond, nPartOnX < mPartOnY ? onNSecond - 1 : onNSecond);

                // goes to other finite element
                // skips 1 column which has been already visited
                onNFirst += SpecialData.NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT - 1;
            }
        }

        private void CreateBoundaryConditionMatrix(int nPartOnX, int mPartOnY)
        {
            int nodesAmount = (nPartOnX + 1) * (2 * mPartOnY + 1) + nPartOnX * (mPartOnY + 1);
            BoundaryConditionMatrix = new int[nodesAmount];

            for (int i = 0; i < nodesAmount; ++i) 
            {
                // left side = 1
                // right side = 4
                // top side = 2
                // bottom side = 7
                // left top corner = left + top = 1 + 2 = 3
                // left bottom corner = left + bottom = 1 + 7 = 8
                // right top corner = right + top = 4 + 2 = 6
                // right bottom corner = right + bottom = 4 + 7 = 11

                if (NodesMatrix[i].X == PointOnY.X) ++BoundaryConditionMatrix[i];

                if (NodesMatrix[i].X == PointOnX.X) BoundaryConditionMatrix[i] += 4;

                if (NodesMatrix[i].Y == PointOnY.Y) BoundaryConditionMatrix[i] += 2;

                if (NodesMatrix[i].Y == PointOnX.Y) BoundaryConditionMatrix[i] += 7;
            }
        }

        public double BasicFunctionsFi(int nodeIndex, double xi1, double xi2)
        {
            double result = 0;
            switch (nodeIndex)
            {
                case 0: result = -1.0 / 4 * (1 - xi1) * (1 - xi2) * (1 + xi1 + xi2); break;
                case 1: result = 1.0 / 2 * (1 - Math.Pow(xi1, 2)) * (1 - xi2); break;
                case 2: result = -1.0 / 4 * (1 + xi1) * (1 - xi2) * (1 - xi1 + xi2); break;
                case 3: result = 1.0 / 2 * (1 - Math.Pow(xi2, 2)) * (1 + xi1); break;
                case 4: result = -1.0 / 4 * (1 + xi1) * (1 + xi2) * (1 - xi1 - xi2); break;
                case 5: result = 1.0 / 2 * (1 - Math.Pow(xi1, 2)) * (1 + xi2); break;
                case 6: result = -1.0 / 4 * (1 - xi1) * (1 + xi2) * (1 + xi1 - xi2); break;
                case 7: result = 1.0 / 2 * (1 - Math.Pow(xi2, 2)) * (1 - xi1); break;
            }
            return result;
        }

        public double DerivativeOfFiOnXi(int nodeIndex, double xi1, double xi2, DerivativeOnXi derivativeOnXi)
        {
            double result = 0;
            switch (nodeIndex)
            {
                case 0: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? 1.0 / 4 * (1 - xi2) * (2 * xi1 + xi2) : 1.0 / 4 * (1 - xi1) * (xi1 + 2 * xi2); break;
                case 1: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? xi1 * (xi2 - 1)                       : 1.0 / 2 * (Math.Pow(xi1, 2) - 1); break;
                case 2: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? 1.0 / 4 * (1 - xi2) * (2 * xi1 - xi2) : 1.0 / 4 * (1 + xi1) * (2 * xi2 - xi1); break;
                case 3: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? 1.0 / 2 * (1 - Math.Pow(xi2, 2))      : -xi2 * (1 + xi1); break;
                case 4: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? 1.0 / 4 * (1 + xi2) * (2 * xi1 + xi2) : 1.0 / 4 * (1 + xi1) * (xi1 + 2 * xi2); break;
                case 5: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? -xi1 * (1 + xi2)                      : 1.0 / 2 * (1 - Math.Pow(xi1, 2)); break;
                case 6: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? 1.0 / 4 * (1 + xi2) * (2 * xi1 - xi2) : 1.0 / 4 * (1 - xi1) * (2 * xi2 - xi1); break;
                case 7: result = derivativeOnXi == DerivativeOnXi.DerivativeOnXi1 ? 1.0 / 2 * (Math.Pow(xi2, 2) - 1)      : xi2 * (xi1 - 1); break;
            }
            return result;
        }

        public double GetJacobian(int finiteElementIndex, double xi1, double xi2)
        {
            // Jacobian - Jacobi matrix determinant
            double sum = 0;
            for (int i = 0; i < SpecialData.NODES_NUMBER_IN_FINITE_ELEMENT; ++i) 
            {
                for (int j = 0; j < SpecialData.NODES_NUMBER_IN_FINITE_ELEMENT; ++j) 
                {
                    if (i != j) 
                    {
                        sum += DerivativeOfFiOnXi(i, xi1, xi2, DerivativeOnXi.DerivativeOnXi1) * DerivativeOfFiOnXi(j, xi1, xi2, DerivativeOnXi.DerivativeOnXi2)
                            * (NodesMatrix[(int)ConnectivityMatrix[finiteElementIndex, i]].X * NodesMatrix[(int)ConnectivityMatrix[finiteElementIndex, j]].Y
                            - NodesMatrix[(int)ConnectivityMatrix[finiteElementIndex, j]].X * NodesMatrix[(int)ConnectivityMatrix[finiteElementIndex, i]].Y);
                    }
                }
            }
            return sum;
        }

        public Matrix GetJacobiInverseMatrix(int finiteElementIndex, double xi1, double xi2)
        {
            int twoDimension = 2;
            Matrix result = new Matrix(twoDimension, twoDimension);
            for (int i = 0; i < SpecialData.NODES_NUMBER_IN_FINITE_ELEMENT; ++i)
            {
                double derivativeOfFiOnXi1 = DerivativeOfFiOnXi(i, xi1, xi2, DerivativeOnXi.DerivativeOnXi1);
                double derivativeOfFiOnXi2 = DerivativeOfFiOnXi(i, xi1, xi2, DerivativeOnXi.DerivativeOnXi2);

                Coordinate2D node = new Coordinate2D
                {
                    X = NodesMatrix[(int)ConnectivityMatrix[finiteElementIndex, i]].X,
                    Y = NodesMatrix[(int)ConnectivityMatrix[finiteElementIndex, i]].Y
                };

                Matrix nextMatrix = new Matrix(twoDimension, twoDimension);

                nextMatrix[0, 0] = derivativeOfFiOnXi2 * node.Y;
                nextMatrix[0, 1] = -derivativeOfFiOnXi1 * node.Y;
                nextMatrix[1, 0] = -derivativeOfFiOnXi2 * node.X;
                nextMatrix[1, 1] = derivativeOfFiOnXi1 * node.X;

                result += nextMatrix;
            }

            double reverseJacobian = Math.Pow(GetJacobian(finiteElementIndex, xi1, xi2), -1);
            for (int i = 0; i < result.RowsAmount; ++i)
            {
                for (int j = 0; j < result.ColumnsAmount; ++j) 
                {
                    result[i, j] *= reverseJacobian;
                }
            }

            return result;
        }

        public Matrix DerivativeOfFiOnAlpha(int finiteElementIndex, int nodeIndex, double xi1, double xi2)
        {
            Matrix derivativesOfFiOnXi1And2 = new Matrix(2, 1);
            derivativesOfFiOnXi1And2[0, 0] = DerivativeOfFiOnXi(nodeIndex, xi1, xi2, DerivativeOnXi.DerivativeOnXi1);
            derivativesOfFiOnXi1And2[1, 0] = DerivativeOfFiOnXi(nodeIndex, xi1, xi2, DerivativeOnXi.DerivativeOnXi2);

            return GetJacobiInverseMatrix(finiteElementIndex, xi1, xi2) * derivativesOfFiOnXi1And2;
        }
    }
}
