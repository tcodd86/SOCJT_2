using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Basis
    {
        #region properties
        /// <summary>
        /// Number mode this is
        /// </summary>
        private int nModeN;
        public int modeN
        {
            get { return nModeN; }
            //set { nModeN = value; }
        }//end property nModeN

        /// <summary>
        /// Value of l, the vibrational angular momentum
        /// </summary>
        private int nl;
        public int l
        {
            get { return nl; }
            //set { nl = value; }
        }//end property l

        /// <summary>
        /// Value of v, the principal vibrational quantum number
        /// </summary>
        private int nv;
        public int v
        {
            get { return nv; }
            //set { nv = value; }
        }//end property nv

        /// <summary>
        /// Value of omega
        /// </summary>
        private double nModeOmega;
        public double modeOmega
        {
            get { return nModeOmega; }
            //set { nModeOmega = value; }
        }//end property nModeOmega

        /// <summary>
        /// Boolean indicating whether this is a nondegenerate mode (true) or degenerate mode (false)
        /// </summary>
        private bool nsymmetryIsA;
        public bool symmetryIsA
        {
            get { return nsymmetryIsA; }
            //set { nsymmetryIsA = value; }
        }//end property symmetry

        /// <summary>
        /// Value of the anharmonicity
        /// </summary>
        private double nAnharmonicity;
        public double anharmonicity
        {
            get { return nAnharmonicity; }
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
        /// <param name="DBasis">
        /// Value of D
        /// </param>
        /// <param name="KBasis">
        /// Value of K
        /// </param>
        public Basis(bool symmetryIsA, int modeN, int l, int v, double modeOmega, double anharmonicity, double DBasis, double KBasis)
        {
            nsymmetryIsA = symmetryIsA;
            nModeN = modeN;
            nl = l;
            nv = v;
            nModeOmega = modeOmega;
            nAnharmonicity = anharmonicity;
            nDBasis = DBasis;
            nKBasis = KBasis;
        }

        public static List<Basis> genBasisVectors(Mode mode1, int modeNumber)
        {
            List<Basis> basisVectors = new List<Basis>();
            //int basisNumber = 0;
            for (int b = 0; b <= mode1.modeVMax; b++)//to go through vMax
            {
                if (mode1.IsAType == false)//test to see if e type, if so then go to code here
                {
                    for (int c = b; c >= -b; c -= 2)//to iterate through all l values, changed from being Modes[i].vmax to b, this should be right
                    {
                        basisVectors.Add(new Basis(false, modeNumber, c, b, mode1.modeOmega, mode1.wExe, mode1.D, mode1.K));
                        /*
                        basisVectors.Add(new Basis());
                        basisVectors[basisNumber].symmetryIsA = false;
                        basisVectors[basisNumber].modeN = modeNumber;                                    
                        basisVectors[basisNumber].l = c;
                        basisVectors[basisNumber].v = b;
                        basisVectors[basisNumber].modeOmega = mode1.modeOmega;
                        basisVectors[basisNumber].anharmonicity = mode1.wExe;
                        basisVectors[basisNumber].DBasis = mode1.D;
                        basisVectors[basisNumber].KBasis = mode1.K;
                        basisNumber++;
                        */
                    }//end for all l
                }//end AType false if
                else//if not e type then it's A type and go here, only difference is l is always 0 here
                {
                    basisVectors.Add(new Basis(true, modeNumber, 0, b, mode1.modeOmega, mode1.wExe, mode1.D, mode1.K));
                    /*
                    basisVectors.Add(new Basis());
                    basisVectors[basisNumber].symmetryIsA = true;
                    basisVectors[basisNumber].modeN = modeNumber;
                    basisVectors[basisNumber].l = 0;
                    basisVectors[basisNumber].v = b;
                    basisVectors[basisNumber].modeOmega = mode1.modeOmega;
                    basisVectors[basisNumber].anharmonicity = mode1.wExe;
                    basisVectors[basisNumber].DBasis = mode1.D;
                    basisVectors[basisNumber].KBasis = mode1.K;
                    basisNumber++;
                    */
                }//end else
            }//end for
            return basisVectors;
        }//end method genBasis
    }//end class Basis
}
