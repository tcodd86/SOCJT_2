using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Basis
    {
        #region properties
        private int nModeN;
        public int modeN
        {
            get { return nModeN; }
            set { nModeN = value; }
        }//end property nModeN

        private int nl;
        public int l
        {
            get { return nl; }
            set { nl = value; }
        }//end property l

        private int nv;
        public int v
        {
            get { return nv; }
            set { nv = value; }
        }//end property nv

        private double nModeOmega;
        public double modeOmega
        {
            get { return nModeOmega; }
            set { nModeOmega = value; }
        }//end property nModeOmega

        private bool nsymmetryIsA;
        public bool symmetryIsA
        {
            get { return nsymmetryIsA; }
            set { nsymmetryIsA = value; }
        }//end property symmetry

        private double nAnharmonicity;
        public double anharmonicity
        {
            get { return nAnharmonicity; }
            set { nAnharmonicity = value; }
        }//end property anharmonicity

        private double nDBasis;
        public double DBasis
        {
            get { return nDBasis; }
            set { nDBasis = value; }
        }//end property DBasis

        private double nKBasis;
        public double KBasis
        {
            get { return nKBasis; }
            set { nKBasis = value; }
        }//end property KBasis
        #endregion properties

        public static List<Basis> genBasisVectors(Mode mode1, int modeNumber)
        {
            List<Basis> basisVectors = new List<Basis>();
            int basisNumber = 0;
            for (int b = 0; b <= mode1.modeVMax; b++)//to go through vMax
            {
                if (mode1.IsAType == false)//test to see if e type, if so then go to code here
                {
                    for (int c = b; c >= -b; c -= 2)//to iterate through all l values, changed from being Modes[i].vmax to b, this should be right
                    {
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
                    }//end for all l
                }//end AType false if
                else//if not e type then it's A type and go here, only difference is l is always 0 here
                {
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
                }//end else
            }//end for
            return basisVectors;
        }//end method genBasis
    }//end class Basis
}
