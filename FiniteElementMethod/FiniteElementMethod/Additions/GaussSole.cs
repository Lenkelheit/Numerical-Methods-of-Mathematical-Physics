using FiniteElementMethod.Matrices;

using static System.Math;

namespace FiniteElementMethod.Additions
{
    public class GaussSole
    {
        // FIELDS
        private Matrix aMatrix;
        private Matrix workAMatrix;

        private Matrix bVector;
        private Matrix workBVector;

        private Matrix xVector;

        private double eps; // accurancy order for real numbers
        private int size;

        // CONSTRUCTORS
        public GaussSole(Matrix aMatrix, Matrix bVector)
            : this(aMatrix, bVector, 0.0001) { }

        public GaussSole(Matrix aMatrix, Matrix bVector, double eps)
        {
            if (aMatrix == null) throw new System.ArgumentNullException("A matrix is null.");
            if (bVector == null) throw new System.ArgumentNullException("B vector is null.");
            if (!bVector.IsVector) throw new System.InvalidOperationException("B should be vector");
            
            int bSize = bVector.RowsAmount;
            int aSize = aMatrix.RowsAmount;

            if (aSize != bSize)
            {
                throw new System.ArgumentException("The amount of rows in A matrix should be equal to amount of rows in B Vector");
            }

            this.aMatrix = aMatrix;
            this.workAMatrix = aMatrix.Clone();

            this.bVector = bVector;
            this.workBVector = bVector.Clone();

            this.xVector = new Matrix(n: bSize, m: 1);

            this.size = bSize;
            this.eps = eps;
        }

        // METHODS        
        public Matrix Run()
        {
            int[] index = InitIndex();
            GaussForwardStroke(index);
            GaussBackwardStroke(index);
            return xVector;
        }

        // PRIVATE METHODS
        // initialize an array of row's indices
        private int[] InitIndex()
        {
            int[] index = new int[size];
            for (int i = 0; i < index.Length; ++i)
            {
                index[i] = i;
            }
            return index;
        }

        // finding a main element in matrix
        private double FindR(int row, int[] index)
        {
            int max_index = row;
            double max = workAMatrix[row, index[max_index]];
            double max_abs = Abs(max);

            // if (row < size - 1)
            for (int cur_index = row + 1; cur_index < size; ++cur_index)
            {
                double cur = workAMatrix[row, index[cur_index]];
                double cur_abs = Abs(cur);
                if (cur_abs > max_abs)
                {
                    max_index = cur_index;
                    max = cur;
                    max_abs = cur_abs;
                }
            }

            if (max_abs < eps)
            {
                if (Abs(workBVector[row, 0]) > eps) throw new System.ArithmeticException("The system of equations is incompatible.");
                else throw new System.ArithmeticException("The system of equations has many solutions.");
            }

            // swaping row's indices
            int temp = index[row];
            index[row] = index[max_index];
            index[max_index] = temp;

            return max;
        }

        private void GaussForwardStroke(int[] index)
        {
            // moving on each row from top to bottom
            for (int i = 0; i < size; ++i)
            {
                // 1) select main element
                double r = FindR(i, index);

                // 2) changing current A matrix row
                for (int j = 0; j < size; ++j)
                {
                    workAMatrix[i, j] /= r;
                }

                // 3) changing element at index I in B vector
                workBVector[i, 0] /= r;

                // 4) subtract the current row from all the rows below
                for (int k = i + 1; k < size; ++k)
                {
                    double p = workAMatrix[k, index[i]];
                    for (int j = i; j < size; ++j)
                    {
                        workAMatrix[k, index[j]] -= workAMatrix[i, index[j]] * p;
                    }
                    workBVector[k, 0] -= workBVector[i, 0] * p;
                    workAMatrix[k, index[i]] = 0D;
                }
            }
        }

        private void GaussBackwardStroke(int[] index)
        {
            // moving on each row from bottom to top
            for (int i = size - 1; i >= 0; --i)
            {
                // 1) initialize start value of x
                double x_i = workBVector[i, 0];

                // 3) updating x
                for (int j = i + 1; j < size; ++j)
                {
                    x_i -= xVector[index[j], 0] * workAMatrix[i, index[j]];
                }
                xVector[index[i], 0] = x_i;
            }
        }
    }
}
