using System;
using System.IO;

using FiniteElementMethod;
using FiniteElementMethod.Matrices;
using FiniteElementMethod.Additions;

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


            Matrix[,] globalMatrix = coordinatePlane.CreateGlobalMatrix();
            Matrix left = ToMatrix(globalMatrix);
            Write(globalMatrix, "../../Global matrix.txt");

            Matrix[] globalVector = coordinatePlane.CreateGlobalVector();
            Matrix right = ToVectorInMatrixView(globalVector);
            Write(globalVector, "../../Global vector.txt");

            GaussSole gauss = new GaussSole(aMatrix: left, bVector: right, eps: 1e-50);
            Matrix result = gauss.Run();
            Matrix.Write(result, "../../RESULT.txt");

            Console.WriteLine("Done");


            Console.ReadLine();
        }

        static void Write(Matrix[,] m, string filePath)
        {
            using (StreamWriter writer = new StreamWriter(filePath))
            {
                for (int rows = 0; rows < m.GetLength(0); ++rows)
                {
                    for (int matrixRow = 0; matrixRow < 6; ++matrixRow)
                    {
                        for (int cols = 0; cols < m.GetLength(1); ++cols)
                        {
                            for (int matrixCols = 0; matrixCols < 6; ++matrixCols)
                            {
                                writer.Write(m[rows, cols][matrixRow, matrixCols].ToString().PadRight(25));
                            }
                        }
                        writer.WriteLine();
                    }
                }
            }
        }

        static void Write(Matrix[] v, string filePath)
        {
            using (StreamWriter writer = new StreamWriter(filePath))
            {
                for (int block = 0; block < v.GetLength(0); ++block)
                {
                    for (int i = 0; i < v[block].RowsAmount; ++i)
                    {
                        for (int j = 0; j < v[block].ColumnsAmount; ++j)
                        {
                            writer.Write(v[block][i, j].ToString().PadRight(25));
                        }
                        writer.WriteLine();
                    }
                    writer.WriteLine();
                }
            }
        }

        static Matrix ToMatrix(Matrix[,] m)
        {
            Matrix matrix = new Matrix(m.GetLength(0) * 6, m.GetLength(1) * 6);
            for (int rows = 0; rows < m.GetLength(0); ++rows)
            {
                for (int matrixRow = 0; matrixRow < 6; ++matrixRow)
                {
                    for (int cols = 0; cols < m.GetLength(1); ++cols)
                    {
                        for (int matrixCols = 0; matrixCols < 6; ++matrixCols)
                        {
                            matrix[6 * rows + matrixRow, 6 * cols + matrixCols] = m[rows, cols][matrixRow, matrixCols];
                        }
                    }
                }
            }
            return matrix;
        }

        static Matrix ToVectorInMatrixView(Matrix[] v)
        {
            Matrix vector = new Matrix(v.GetLength(0) * 6, 1);

            for (int rows = 0; rows < v.GetLength(0); ++rows)
            {
                for (int i = 0; i < v[rows].RowsAmount; ++i)
                {
                    for (int j = 0; j < v[rows].ColumnsAmount; ++j)
                    {
                        vector[rows * 6 + i, j] = v[rows][i, j];
                    }
                }
            }
            return vector;
        }
    }
}
