using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethod.Additions
{
    public static class SpecialData
    {
        // CONSTANTS
        public const int NODES_NUMBER_IN_FINITE_ELEMENT = 8;

        public const int NODES_NUMBER_ON_EDGE_IN_FINITE_ELEMENT = 3;

        // METHODS
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
        }
    }
}
