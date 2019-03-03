using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethod.Additions
{
    public class NodesMatrix
    {
        // FIELDS
        private Coordinate2D[][] nodesCoordinate;

        // PROPERTIES
        public int NodesOnMAmount { get; }
        public int EvenNodesOnNAmount { get; }
        public int OddNodesOnNAmount { get; }

        // CONSTRUCTORS
        public NodesMatrix(int n, int m)
        {
            int NodesOnMAmount = 2 * m + 1;
            int EvenNodesOnNAmount = 2 * n + 1;
            int OddNodesOnNAmount = n + 1;

            nodesCoordinate = new Coordinate2D[NodesOnMAmount][];
            for (int i = 0; i < NodesOnMAmount; ++i) 
            {
                if (SpecialData.IsEven(i))
                {
                    nodesCoordinate[i] = new Coordinate2D[EvenNodesOnNAmount];
                }
                else
                {
                    nodesCoordinate[i] = new Coordinate2D[OddNodesOnNAmount];
                }
            }
        }

        // INDEXERS
        public Coordinate2D this[int index]
        {
            get
            {
                Pair indices = CreateIndicesOnNAndM(index);
                return nodesCoordinate[indices.First][indices.Second];
            }
            set
            {
                Pair indices = CreateIndicesOnNAndM(index);
                nodesCoordinate[indices.First][indices.Second] = value;
            }
        }

        public Coordinate2D this[int first, int second]
        {
            get
            {
                return nodesCoordinate[first][second];
            }
            set
            {
                nodesCoordinate[first][second] = value;
            }
        }

        // METHODS
        public void Show()
        {
            for (int i = 0; i < nodesCoordinate.Length; ++i) 
            {
                for (int j = 0; j < nodesCoordinate[i].Length; ++j)
                {
                    Console.Write($"({nodesCoordinate[i][j].X}, {nodesCoordinate[i][j].Y}) ");
                }
                Console.WriteLine();
            }
        }

        private Pair CreateIndicesOnNAndM(int index)
        {
            // converts index to two indices for nodesCoordinate array
            int i = 0;
            while (index >= nodesCoordinate[i].Length) 
            {
                index -= nodesCoordinate.GetLength(i);
                ++i;
            }
            return new Pair { First = i, Second = index };
        }
    }
}
