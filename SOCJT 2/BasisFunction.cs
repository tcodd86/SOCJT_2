using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace ConsoleApplication1
{
    class BasisFunction
    {
        /// <summary>
        /// Value of j
        /// </summary>
        private decimal nJ;
        public decimal J
        {
            get { return nJ; }
            set { nJ = value; }
        }//end property J

        /// <summary>
        /// Value of Lambda, label used to distinguish betwen two components of degenerate electronic state
        /// </summary>
        private int nlambda;
        public int Lambda
        {
            get { return nlambda; }
            set { nlambda = value; }
        }//end property lambda

        /// <summary>
        /// List of Basis objects to store information for each mode
        /// </summary>
        public List<BasisByMode> modesInVec = new List<BasisByMode>();

        /// <summary>
        /// Constructor for BasisFunction
        /// </summary>
        /// <param name="basis">
        /// List of BasisByMode objects
        /// </param>
        /// <param name="lam">
        /// Lambda, the label distinguishing between the two components of the degenerate electronic state
        /// </param>
        public BasisFunction(List<BasisByMode> basis, int lam)
        {
            int l = 0;
            modesInVec.Clear();
            modesInVec.AddRange(basis);
            for (int i = 0; i < basis.Count; i++)
            {
                l += basis[i].l;
            }
            nJ = (decimal)l + ((decimal)lam) / 2M;
            Lambda = lam;
        }//end constructor

        /// <summary>
        /// Class to implement IComparer interface for BasisFunctions for sorting
        /// </summary>
        private class sortBasisFunctionsHelper : IComparer<BasisFunction>
        {
            int IComparer<BasisFunction>.Compare(BasisFunction a, BasisFunction b)
            {
                int val = 0;
                //start at last mode and check it, if it matches move to next mode
                //try leaving lambda out of it completely
                /*
                //start with lambda values
                if (a.Lambda > b.Lambda)
                {
                    val = 1;
                    return val;
                }
                if (a.Lambda < b.Lambda)
                {
                    val = -1;
                    return val;
                }
                //if lambdas are the same then check v and l for each mode until a difference is found
                */
                for(int place = a.modesInVec.Count - 1; place >= 0; place--)
                {
                    //check v
                    if (a.modesInVec[place].v > b.modesInVec[place].v)
                    {
                        val = 1;
                        break;
                    }
                    if (a.modesInVec[place].v < b.modesInVec[place].v)
                    {
                        val = -1;
                        break;
                    }
                    //if neither of those is true then the two v values are equal, then check l values
                    if (a.modesInVec[place].l > b.modesInVec[place].l)
                    {
                        val = 1;
                        break;
                    }
                    if (a.modesInVec[place].l < b.modesInVec[place].l)
                    {
                        val = -1;
                        break;
                    }
                    //if neither of these is true then the two l values are identical, check the next mode in the list
                }//end for loop
                return val;
            }//end comparer            
        }//end class SortBasisFunctions

        /// <summary>
        /// Functiont IComparer to pass in sort functions
        /// </summary>
        /// <returns>
        /// IComparer for custom sorting of BasisFunction objects
        /// </returns>
        public static IComparer<BasisFunction> sortBasisFunctions()
        {
            return (IComparer<BasisFunction>) new BasisFunction.sortBasisFunctionsHelper();
        }//and sortBasisFunctions

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
        static public List<BasisFunction> genJVecs(List<List<BasisByMode>> basisByMode, int numModes, decimal minJ, decimal maxJ)
        {
            //List of basis functions, to be returned
            List<BasisFunction> hamBasisSet = new List<BasisFunction>();

            //List of BasisByMode objects, this list is n long where n is the number of modes.  Contains
            //the info for each mode's v and l values
            List<BasisByMode> modesInVec = new List<BasisByMode>();

            //this array stores the upper bounds for the countkeeper function
            int[] maxValues = new int[numModes];
            for (int i = 0; i < numModes; i++)
            {
                maxValues[i] = basisByMode[i].Count;
            }//end for
            
            //this array stores the count values for comparison agains the maxValues
            int[] count = new int[numModes];
            for (int i = 0; i < count.Length; i++)
            {
                count[i] = 0;
            }//end for
            
            //boolean value used to terminate generating BasisFunctions when all possible combinations have been generated
            bool keepGoing = true;

            //loop that runs until all possible combinations of BasisByMode values and Lambda have been generated
            while (keepGoing == true)
            {
                //resets the list each iteration
                modesInVec.Clear();

                //adds the BasisByMode objects for each mode for a given set of count[] values
                for (int n = 0; n < count.Length; n++)
                {
                    modesInVec.Add(basisByMode[n][count[n]]);
                }//end for

                //call to function to increase count by one
                countKeeper(ref count, 0, maxValues, ref keepGoing);

                //This generates the same basis vector with Lambda = +1 and -1
                for (int i = -1; i < 2; i += 2)
                {
                    BasisFunction vector = new BasisFunction(modesInVec, i);

                    //checks that the value of J for 'vector' is within the bounds set in the input file.  If so the vector is added to the basis set
                    if (vector.J <= maxJ & vector.J >= minJ)
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
        static private void countKeeper(ref int[] count, int n, int[] ints, ref bool keepGoing)
        {
            count[n] += 1;
            if (count[n] == ints[n])
            {
                //checks if the count is now done.
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
