using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementMethod.Additions
{
    public enum NodeSide
    {
        Middle = 0,
        Left = 1,
        Top = 2,
        Right = 4,
        Bottom = 7,
        LeftTop = 3,
        LeftBottom = 8,
        RightTop = 6,
        RightBottom = 11
    }
}
