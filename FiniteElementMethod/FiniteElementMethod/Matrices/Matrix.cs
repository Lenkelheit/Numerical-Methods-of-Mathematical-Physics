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
        private int[,] matrix;

        // PROPERTIES
        public int RowsAmount { get; }
        public int ColumnsAmount { get; }

        // CONSTRUCTORS
        public Matrix(int n, int m)
        {
            RowsAmount = n;
            ColumnsAmount = m;
            matrix = new int[RowsAmount, ColumnsAmount];
        }

        // INDEXERS
        public int this[int row, int column]
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
    }
}
