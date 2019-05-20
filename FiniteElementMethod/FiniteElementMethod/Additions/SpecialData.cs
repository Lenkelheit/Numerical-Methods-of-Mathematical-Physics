using System;

namespace FiniteElementMethod.Additions
{
    public static class SpecialData
    {
        // CONSTANTS
        public const int NODES_NUMBER_IN_FINITE_ELEMENT = 8;

        public const int NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT = 3;

        public const int B_MATRIX_DIMENSION = 11;

        // METHODS
        public static double[,] GaussNodeMatrix()
        {
            return new double[9, 2]
            {
                { -0.77459,  -0.77459 },
                {  0,        -0.77459 },
                {  0.77459,  -0.77459 },
                { -0.77459,   0       },
                {  0,         0       },
                {  0.77459,   0       },
                { -0.77459,   0.77459 },
                {  0,         0.77459 },
                {  0.77459,   0.77459 }
            };
        }

        public static double[] GaussWeights()
        {
            return new double[3]
            {
                0.5555555555,
                0.8888888888,
                0.5555555555
            };
        }

        public static bool IsEven(int number)
        {
            return number % 2 == 0;
        }

        public static void ShowArrayWithIndices(int[] array)
        {
            for (int i = 0; i < array.Length; ++i) 
            {
                Console.WriteLine($"{i}: {array[i]}");
            }
            Console.WriteLine();
        }
    }
}
