using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class BasisByMode
    {
        #region properties
        /// <summary>
        /// Number mode this is
        /// </summary>
        public int ModeN { get; private set; }

        /// <summary>
        /// Value of l, the vibrational angular momentum
        /// </summary>
        public int L { get; private set; }

        /// <summary>
        /// Value of v, the principal vibrational quantum number
        /// </summary>
        private int V;
        public int v
        {
            get { return V; }
            //set { nv = value; }
        }//end property nv

        /// <summary>
        /// Value of omega
        /// </summary>
        private double ModeOmega;
        public double modeOmega
        {
            get { return ModeOmega; }
            //set { nModeOmega = value; }
        }//end property nModeOmega

        /// <summary>
        /// Boolean indicating whether this is a nondegenerate mode (true) or degenerate mode (false)
        /// </summary>
        private bool SymmetryIsA;
        public bool symmetryIsA
        {
            get { return SymmetryIsA; }
            //set { nsymmetryIsA = value; }
        }//end property symmetry

        /// <summary>
        /// Value of the anharmonicity
        /// </summary>
        private double Anharmonicity;
        public double anharmonicity
        {
            get { return Anharmonicity; }
            //set { nAnharmonicity = value; }
        }//end property anharmonicity

        /// <summary>
        /// Value of D
        /// </summary>
        private double nDBasis;
        public double DBasis
        {
            get { return nDBasis; }
            //set { nDBasis = value; }
        }//end property DBasis

        /// <summary>
        /// Value of K
        /// </summary>
        private double nKBasis;
        public double KBasis
        {
            get { return nKBasis; }
            //set { nKBasis = value; }
        }//end property KBasis
        #endregion properties

        /// <summary>
        /// Constructor for Basis object
        /// </summary>
        /// <param name="symmetryIsA">
        /// Boolean value indicating if vector is A type (true) or not (false)
        /// </param>
        /// <param name="modeN">
        /// What number mode this is
        /// </param>
        /// <param name="l">
        /// The value of the vibrational angular momentum (always zero if symmetryIsA is true)
        /// </param>
        /// <param name="v">
        /// Value of the principal vibrational quantum number
        /// </param>
        /// <param name="modeOmega">
        /// Value of omega
        /// </param>
        /// <param name="anharmonicity">
        /// Value of the anharmonicity constant
        /// </param>
        /// <param name="D">
        /// Value of D
        /// </param>
        /// <param name="K">
        /// Value of K
        /// </param>
        public BasisByMode(bool symmetryIsA, int modeN, int l, int v, double modeOmega, double anharmonicity, double D, double K)
        {
            SymmetryIsA = symmetryIsA;
            ModeN = modeN;
            L = l;
            V = v;
            ModeOmega = modeOmega;
            Anharmonicity = anharmonicity;
            nDBasis = D;
            nKBasis = K;
        }

        /// <summary>
        /// Function to generate all v/l combinations for a given mode
        /// </summary>
        /// <param name="mode1">
        /// The ModeInfo object for the mode who's v/l combinations are being generated.
        /// </param>
        /// <param name="modeNumber">
        /// What number mode in the input file is it.
        /// </param>
        /// <returns>
        /// List of all possible BasisByMode objects for a given mode.
        /// </returns>
        public static List<BasisByMode> genVLCombinations(ModeInfo mode1, int modeNumber)
        {
            //List to return
            List<BasisByMode> basisVectors = new List<BasisByMode>();

            for (int b = 0; b <= mode1.modeVMax; b++)//loop to go through vMax
            {
                if (mode1.IsAType == false)//test to see if e type, if so then go to code here
                {
                    for (int c = b; c >= -b; c -= 2)//to iterate through all l values, changed from being Modes[i].vmax to b
                    {
                        basisVectors.Add(new BasisByMode(false, modeNumber, c, b, mode1.modeOmega, mode1.wExe, mode1.D, mode1.K));
                    }//end for all l
                }//end AType false if
                else//if not e type then it's A type and go here, only difference is l is always 0 here
                {
                    basisVectors.Add(new BasisByMode(true, modeNumber, 0, b, mode1.modeOmega, mode1.wExe, mode1.D, mode1.K));
                }//end else
            }//end for
            return basisVectors;
        }//end method genBasis
    }//end class BasisByMode
}
