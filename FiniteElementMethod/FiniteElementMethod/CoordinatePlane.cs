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
        // FIELDS
        #region POINTS
        private Coordinate2D PointOnX;
        private Coordinate2D PointOnY;
        #endregion

        private CoordinateSystemConfig coordinateSystemConfig;
        private CylindricalShellConfig shellConfig;

        // PROPERTIES
        #region MATRICES
        public NodesMatrix NodesMatrix { get; private set; }
        public Matrix ConnectivityMatrix { get; private set; }
        public int[] BoundaryConditionMatrix { get; private set; }
        #endregion

        // CONSTRUCTORS
        public CoordinatePlane(CoordinateSystemConfig coordinateSystemConfig, CylindricalShellConfig shellConfig)
        {
            this.PointOnX = new Coordinate2D { X = coordinateSystemConfig.B, Y = coordinateSystemConfig.C };
            this.PointOnY = new Coordinate2D { X = coordinateSystemConfig.A, Y = coordinateSystemConfig.D };

            this.coordinateSystemConfig = coordinateSystemConfig;
            this.shellConfig = shellConfig;

            SplitCoordinatePlane(coordinateSystemConfig.N, coordinateSystemConfig.M);
        }

        // METHODS
        private void SplitCoordinatePlane(int nPartOnX, int mPartOnY)
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
            for (double i = nPartOnX < mPartOnY ? PointOnX.Y : PointOnY.X; i < (nPartOnX < mPartOnY ? PointOnY.Y + hOnY / 2 : PointOnX.X + hOnX / 2); i += nPartOnX < mPartOnY ? hOnY : hOnX) 
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

        public Matrix CreateMaterialElasticCharacteristicsMatrixBWithE()
        {
            double v = shellConfig.V, E = shellConfig.E, h = shellConfig.H;

            Matrix B = new Matrix(SpecialData.B_MATRIX_DIMENSION, SpecialData.B_MATRIX_DIMENSION);

            B[0, 0] = B[1, 1] = B[2, 2] = ((1 - v) * E * h) / ((1 + v) * (1 - 2 * v));
            B[0, 1] = B[0, 2] = B[1, 0] = B[1, 2] = B[2, 0] = B[2, 1] = (v * E * h) / ((1 + v) * (1 - 2 * v));

            B[3, 3] = B[4, 4] = B[5, 5] = 2 * E * h / (1 + v);

            B[6, 6] = B[7, 7] = ((1 - v) * E * h * h * h) / (12 * (1 + v) * (1 - 2 * v));
            B[6, 7] = B[7, 6] = (v * E * h * h * h) / (12 * (1 + v) * (1 - 2 * v));

            B[8, 8] = B[9, 9] = B[10, 10] = (2 * E * h * h * h) / (12 * (1 + v));
            return B;
        }

        public Matrix CreateDifferentialOperatorMatrixCl(int finiteElementIndex, int nodeIndex, int GausseNodeIndex)
        {
            Matrix Cl = new Matrix(11, 6);

            double firstGaussElem = SpecialData.GaussNodeMatrix()[GausseNodeIndex, 0],
                                    secondGaussElem = SpecialData.GaussNodeMatrix()[GausseNodeIndex, 1],
                                    R = shellConfig.R,
                                    k2 = 1 / R;

            Matrix derivativeOnAlpha = DerivativeOfFiOnAlpha(finiteElementIndex, nodeIndex, firstGaussElem, secondGaussElem);

            double d1 = derivativeOnAlpha[0, 0], d2 = derivativeOnAlpha[1, 0], phi = BasicFunctionsFi(nodeIndex, firstGaussElem, secondGaussElem);

            Cl[0, 0] = d1;
            Cl[1, 1] = d2 / R;
            Cl[1, 2] = phi / R;
            Cl[2, 5] = phi;
            Cl[3, 0] = d2 / (2 * R);
            Cl[3, 1] = d1 / 2;
            Cl[4, 2] = d1 / 2;
            Cl[4, 3] = phi / 2;
            Cl[5, 1] = -phi / (2 * R);
            Cl[5, 2] = d2 / (2 * R);
            Cl[5, 4] = phi / 2;
            Cl[6, 3] = d1;
            Cl[7, 4] = d2 / R;
            Cl[7, 5] = k2 * phi;
            Cl[8, 1] = d1 / (2 * R);
            Cl[8, 3] = d2 / (2 * R);
            Cl[8, 4] = d1 / 2;
            Cl[9, 5] = d1 / 2;
            Cl[10, 5] = d2 / (2 * R);

            return Cl;
        }

        public Matrix CreateGlobalMatrixForGaussNodeIndex(int finiteElementIndex, int gaussNodeIndex)
        {
            Matrix BE = CreateMaterialElasticCharacteristicsMatrixBWithE();

            Matrix[] Cl = new Matrix[8];

            for (int nodeIndex = 0; nodeIndex < SpecialData.NODES_NUMBER_IN_FINITE_ELEMENT; ++nodeIndex) 
            {
                Cl[nodeIndex] = CreateDifferentialOperatorMatrixCl(finiteElementIndex, nodeIndex, gaussNodeIndex);
            }

            // 8*8 of 6*6 matrices
            Matrix[,] matricesForGlobalMatrixForGaussNode = new Matrix[8, 8];
            for (int i = 0; i < matricesForGlobalMatrixForGaussNode.GetLength(0); ++i)
            {
                for (int j = 0; j < matricesForGlobalMatrixForGaussNode.GetLength(1); ++j)
                {
                    // 6*6
                    matricesForGlobalMatrixForGaussNode[i, j] = Cl[i].GetTransposeMatrix() * BE * Cl[j];
                }
            }

            // expand to 48*48 matrix
            Matrix globalMatrixForGaussNode = new Matrix(48, 48);
            for (int i = 0; i < 8; ++i)
            {
                for (int j = 0; j < 8; ++j)
                {
                    for (int k1 = 0; k1 < 6; ++k1)
                    {
                        for (int k2 = 0; k2 < 6; ++k2)
                        {
                            globalMatrixForGaussNode[i * 6 + k1, j * 6 + k2] = matricesForGlobalMatrixForGaussNode[i, j][k1, k2];
                        }
                    }
                }
            }

            return globalMatrixForGaussNode;
        }

        public Matrix CreateGlobalMatrix(int finiteElementIndex)
        {
            double R = shellConfig.R;

            double[] gaussWeights = SpecialData.GaussWeights();
            double[,] gaussNodes = SpecialData.GaussNodeMatrix();

            Matrix globalMatrix = new Matrix(48, 48);
            for (int gaussNodeIndex = 0; gaussNodeIndex < gaussNodes.GetLength(0); ++gaussNodeIndex) // 9
            {
                double detJ = GetJacobian(finiteElementIndex, gaussNodes[gaussNodeIndex, 0], gaussNodes[gaussNodeIndex, 1]);
                double M = gaussWeights[gaussNodeIndex / 3] * gaussWeights[gaussNodeIndex % 3] * detJ * R;

                globalMatrix += CreateGlobalMatrixForGaussNodeIndex(finiteElementIndex, gaussNodeIndex) * M;
            }
            return globalMatrix;
        }

        public double BasicFunctionsFi(int nodeIndex, double xi1, double xi2)
        {
            double result = 0;
            switch (nodeIndex)
            {
                case 0: result = -1.0 / 4 * (1 - xi1) * (1 - xi2) * (1 + xi1 + xi2); break;
                case 1: result = 1.0 / 2 * (1 - Math.Pow(xi1, 2)) * (1 - xi2);       break;
                case 2: result = -1.0 / 4 * (1 + xi1) * (1 - xi2) * (1 - xi1 + xi2); break;
                case 3: result = 1.0 / 2 * (1 - Math.Pow(xi2, 2)) * (1 + xi1);       break;
                case 4: result = -1.0 / 4 * (1 + xi1) * (1 + xi2) * (1 - xi1 - xi2); break;
                case 5: result = 1.0 / 2 * (1 - Math.Pow(xi1, 2)) * (1 + xi2);       break;
                case 6: result = -1.0 / 4 * (1 - xi1) * (1 + xi2) * (1 + xi1 - xi2); break;
                case 7: result = 1.0 / 2 * (1 - Math.Pow(xi2, 2)) * (1 - xi1);       break;
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

            return result * reverseJacobian;
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
