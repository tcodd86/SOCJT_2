using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class JBasisVector
    {
        private decimal nJ;
        public decimal J
        {
            get { return nJ; }
            set { nJ = value; }
        }//end property J

        private int nlambda;
        public int Lambda
        {
            get { return nlambda; }
            set { nlambda = value; }
        }//end property lambda

        public List<Basis> modesInVec = new List<Basis>();

        public void setVec(List<Basis> basis, int lam)
        {
            int l = 0;
            modesInVec.Clear();
            modesInVec.AddRange(basis);
            for (int i = 0; i < basis.Count; i++)
            {
                l += basis[i].l;
            }
            nJ = (decimal) l + ((decimal) lam) / 2M;
            Lambda = lam;
        }//end method setVec


        /// <summary>
        /// Function to generate a list of all possible JBasisVectors.
        /// </summary>
        /// <param name="basisByMode">
        /// List of Basis objects sorted by mode.
        /// </param>
        /// <param name="numModes">
        /// Number of modes.
        /// </param>
        /// <param name="minJ">
        /// Minimum j value.
        /// </param>
        /// <param name="maxJ">
        /// Maximum j value.
        /// </param>
        /// <returns>
        /// List of all possible JBasisVector objects.
        /// </returns>
        static public List<JBasisVector> genJVecs(List<List<Basis>> basisByMode, int numModes, decimal minJ, decimal maxJ)
        {
            List<JBasisVector> hamBasisSet = new List<JBasisVector>();
            List<Basis> modesInVec = new List<Basis>();
            List<int> ints = new List<int>();
            int[] count;
            int nModes;
            bool keepGoing = true;
            for (int i = 0; i < numModes; i++)
            {
                ints.Add(new int());
                ints[i] = basisByMode[i].Count;
            }//end for
            nModes = numModes;
            count = ints.ToArray();
            for (int i = 0; i < count.Length; i++)
            {
                count[i] = 0;
            }//end for
            while (keepGoing == true)
            {
                //recursiveMethod(basisByMode, numModes);
                modesInVec.Clear();
                for (int n = 0; n < count.Length; n++)
                {
                    int temp = count[n];
                    modesInVec.Add(basisByMode[n][temp]);
                }//end for
                countKeeper(ref count, 0, ints, ref keepGoing);

                for (int i = -1; i < 2; i += 2)
                {
                    JBasisVector vector = new JBasisVector();//moved this into the for loop so that I'm generating a new vector object for each value of lambda
                    vector.setVec(modesInVec, i);
                    if (vector.J <= maxJ & vector.J >= minJ)//flesh out this method to make it work THIS IS WHERE I PUT A BREAKPOINT
                    {
                        hamBasisSet.Add(vector);
                    }//end if
                }
            }//end while
            return hamBasisSet;
        }//end method genJVecs

        /// <summary>
        /// Recursive method to allow user to count through all possible combinations of basis functions.
        /// </summary>
        /// <param name="count">
        /// Array of ints that correspond to current 'count' or combination of Basis objects.
        /// </param>
        /// <param name="n">
        /// Gives index of array to be incremented (corrsponds to which mode).
        /// </param>
        static private void countKeeper(ref int[] count, int n, List<int> ints, ref bool keepGoing)
        {
            count[n] += 1;
            if (count[n] == ints[n])
            {
                if (n + 1 != count.Length)
                {
                    count[n] = 0;
                    n += 1;
                    countKeeper(ref count, n, ints, ref keepGoing);
                }//end if
                else
                {
                    keepGoing = false;
                }//end else
            }//end if
        }//end method countKeeper
    }//end class JBasisVector
}
