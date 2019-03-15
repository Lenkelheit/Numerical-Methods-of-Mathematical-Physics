using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FiniteElementMethod.Additions;

namespace FiniteElementMethod.Matrices
{
    public class NodesMatrix
    {
        // FIELDS
        private Coordinate2D[][] nodesCoordinate;

        // PROPERTIES
        public int NodesOnNAmount { get; }
        public int EvenNodesOnMAmount { get; }
        public int OddNodesOnMAmount { get; }

        // CONSTRUCTORS
        public NodesMatrix(int n, int m)
        {
            NodesOnNAmount = 2 * n + 1;
            EvenNodesOnMAmount = 2 * m + 1;
            OddNodesOnMAmount = m + 1;

            nodesCoordinate = new Coordinate2D[NodesOnNAmount][];
            for (int i = 0; i < NodesOnNAmount; ++i) 
            {
                if (SpecialData.IsEven(i))
                {
                    nodesCoordinate[i] = new Coordinate2D[EvenNodesOnMAmount];
                }
                else
                {
                    nodesCoordinate[i] = new Coordinate2D[OddNodesOnMAmount];
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

        public Coordinate2D this[int row, int column]
        {
            get
            {
                return nodesCoordinate[row][column];
            }
            set
            {
                nodesCoordinate[row][column] = value;
            }
        }

        // METHODS
        public void Show()
        {
            for (int i = 0; i < nodesCoordinate.Length; ++i) 
            {
                for (int j = 0; j < nodesCoordinate[i].Length; ++j)
                {
                    Console.Write($"{"(" + nodesCoordinate[i][j].X + ", " + nodesCoordinate[i][j].Y + ")",-15}");
                }
                Console.WriteLine();
            }
        }

        public int CreateGlobalIndex(int n, int m)
        {
            int globalIndex = 0;
            for (int i = 0; i < n; ++i) 
            {
                globalIndex += nodesCoordinate[i].Length;
            }
            return globalIndex + m;
        }

        private Pair CreateIndicesOnNAndM(int index)
        {
            // converts index to two indices for nodesCoordinate array
            int i = 0;
            while (index >= nodesCoordinate[i].Length) 
            {
                index -= nodesCoordinate[i].Length;
                ++i;
            }
            return new Pair { First = i, Second = index };
        }
    }
}
