using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethod.Matrices
{
    public class Matrix
    {
        // FIELDS
        private double[,] matrix;

        // PROPERTIES
        public int RowsAmount { get; }
        public int ColumnsAmount { get; }

        // CONSTRUCTORS
        public Matrix(int n, int m)
        {
            RowsAmount = n;
            ColumnsAmount = m;
            matrix = new double[RowsAmount, ColumnsAmount];
        }

        // INDEXERS
        public double this[int row, int column]
        {
            get
            {
                return matrix[row, column];
            }
            set
            {
                matrix[row, column] = value;
            }
        }

        // METHODS
        public void Show()
        {
            for (int i = 0; i < RowsAmount; ++i) 
            {
                for (int j = 0; j < ColumnsAmount; ++j) 
                {
                    Console.Write($"{matrix[i, j],-20}");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        public Matrix Clone()
        {
            Matrix cloneMatrix = new Matrix(RowsAmount, ColumnsAmount);
            for (int i = 0; i < RowsAmount; ++i)
            {
                for (int j = 0; j < ColumnsAmount; ++j)
                {
                    cloneMatrix[i, j] = this[i, j];
                }
            }
            return cloneMatrix;
        }

        public Matrix GetTransposeMatrix()
        {
            Matrix transposeMatrix = new Matrix(ColumnsAmount, RowsAmount);

            for (int i = 0; i < RowsAmount; ++i)
            {
                for (int j = 0; j < ColumnsAmount; ++j)
                {
                    transposeMatrix[j, i] = this[i, j];
                }
            }
            return transposeMatrix;
        }

        // OPERATORS
        public static Matrix operator +(Matrix first, Matrix second)
        {
            if (first.RowsAmount == second.RowsAmount && first.ColumnsAmount == second.ColumnsAmount) 
            {
                Matrix sum = new Matrix(first.RowsAmount, first.ColumnsAmount);
                for (int i = 0; i < first.RowsAmount; ++i)
                {
                    for (int j = 0; j < first.ColumnsAmount; ++j)
                    {
                        sum[i, j] = first[i, j] + second[i, j];
                    }
                }
                return sum;
            }
            throw new ArgumentException("Different quantity of rows and columns is in matrices!");
        }

        public static Matrix operator *(Matrix first, Matrix second)
        {
            if (first.ColumnsAmount == second.RowsAmount)
            {
                Matrix product = new Matrix(first.RowsAmount, second.ColumnsAmount);
                for (int i = 0; i < first.RowsAmount; ++i)
                {
                    for (int j = 0; j < second.ColumnsAmount; ++j)
                    {
                        for (int k = 0; k < first.ColumnsAmount; ++k) 
                        {
                            product.matrix[i, j] += first.matrix[i, k] * second.matrix[k, j];
                        }
                    }
                }
                return product;
            }
            throw new ArgumentException("There is different quantity of columns in first matrix and rows in second matrix!");
        }

        public static Matrix operator *(Matrix matrix, double value)
        {
            Matrix result = matrix.Clone();
            for (int i = 0; i < result.RowsAmount; ++i)
            {
                for (int j = 0; j < result.ColumnsAmount; ++j)
                {
                    result[i, j] *= value;
                }
            }
            return result;
        }
    }
}
