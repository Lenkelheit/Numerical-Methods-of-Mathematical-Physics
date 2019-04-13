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
                    Console.Write($"{matrix[i, j],-4}");
                }
                Console.WriteLine();
            }
        }

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
    }
}
