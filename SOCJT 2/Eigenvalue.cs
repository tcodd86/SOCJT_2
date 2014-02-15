using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Eigenvalue
    {
        #region Properties
        /// <summary>
        /// J block that this eigenvalue is in
        /// </summary>
        public decimal pJ { get; private set; }

        /// <summary>
        /// Which eigenvalue this is in the j-block
        /// </summary>
        public int nJ { get; set; }

        /// <summary>
        /// Value of Sigma (total spin)
        /// </summary>
        public decimal Sig { get; private set; }

        /// <summary>
        /// Eigenvalue
        /// </summary>
        public double Ev { get; set; }

        /// <summary>
        /// True if symmetric (a1)
        /// </summary>
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
