using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethod.Additions
{
    public class FiniteElement
    {
        // CONSTANTS
        private const int NODES_NUMBER = 8;

        // FIELDS
        public double[] Nodes { get; set; }

        // CONSTRUCTORS
        public FiniteElement()
        {
            Nodes = new double[NODES_NUMBER];
        }
    }
}
