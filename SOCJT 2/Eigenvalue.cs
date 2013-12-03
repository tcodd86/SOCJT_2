using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Eigenvalue
    {
        #region Properties
        private decimal npJ;
        public decimal pJ
        {
            get { return npJ; }
            set { npJ = value; }
        }

        private int nnJ;
        public int nJ
        {
            get { return nnJ; }
            set { nnJ = value; }
        }

        private decimal nSig;
        public decimal Sig
        {
            get { return nSig; }
            set { nSig = value; }
        }

        private double nEv;
        public double Ev
        {
            get { return nEv; }
            set { nEv = value; }
        }

        public bool isa1 { get; private set; }

        #endregion

        public Eigenvalue(decimal pJa, int nJa, decimal Siga, double Eva, bool isa1)
        {
            pJ = pJa;
            nJ = nJa;
            Sig = Siga;
            Ev = Eva;
            this.isa1 = isa1;
        }//end constructor
    }//end class Eigenvalue
}
