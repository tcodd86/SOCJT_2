﻿using System;
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
        public decimal J { get; private set; }

        /// <summary>
        /// Value of Lambda, label used to distinguish betwen two components of degenerate electronic state
        /// </summary>
        public int Lambda { get; private set; }

        /// <summary>
        /// HashCode for storing basis function positions.
        /// </summary>
        public string HashCode { get; private set; }

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
                l += basis[i].L;
            }
            J = (decimal)l + ((decimal)lam) / 2M;
            Lambda = lam;
            HashCode = GenerateHashCode(basis, lam);
        }//end constructor

        /// <summary>
        /// Class to implement IComparer interface for BasisFunctions for sorting
        /// </summary>
        private class SortBasisFunctionsHelper : IComparer<BasisFunction>
        {
            int IComparer<BasisFunction>.Compare(BasisFunction a, BasisFunction b)
            {
                int val = 0;
                for(int place = a.modesInVec.Count - 1; place >= 0; place--)
                {
                    //check v
                    if (a.modesInVec[place].V > b.modesInVec[place].V)
                    {
                        val = 1;
                        break;
                    }
                    if (a.modesInVec[place].V < b.modesInVec[place].V)
                    {
                        val = -1;
                        break;
                    }
                    //if neither of those is true then the two v values are equal, then check l values
                    if (a.modesInVec[place].L > b.modesInVec[place].L)
                    {
                        val = 1;
                        break;
                    }
                    if (a.modesInVec[place].L < b.modesInVec[place].L)
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
        public static IComparer<BasisFunction> SortBasisFunctions()
        {
            return (IComparer<BasisFunction>) new BasisFunction.SortBasisFunctionsHelper();
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
        static public List<BasisFunction> GenJVecs(List<List<BasisByMode>> basisByMode, int numModes, decimal minJ, decimal maxJ)
        {
            //List of basis functions, to be returned
            List<BasisFunction> hamBasisSet = new List<BasisFunction>();
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
            
            bool keepGoing = true;

            //runs until all possible combinations of BasisByMode values and Lambda have been generated
            while (keepGoing == true)
            {
                modesInVec.Clear();

                //adds the BasisByMode objects for each mode for a given set of count[] values
                for (int n = 0; n < count.Length; n++)
                {
                    modesInVec.Add(basisByMode[n][count[n]]);
                }//end for

                //call to function to increase count by one
                CountKeeper(ref count, 0, maxValues, ref keepGoing);

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
        /// <param name="ints">
        /// Array of max values for the CountKeeper
        /// </param>
        /// <param name="keepGoing">
        /// Boolean that's true while the count needs to keep going and false when it's done.
        /// </param>
        static private void CountKeeper(ref int[] count, int n, int[] ints, ref bool keepGoing)
        {
            count[n] += 1;
            if (count[n] == ints[n])
            {
                //checks if the count is now done.
                if (n + 1 != count.Length)
                {
                    count[n] = 0;
                    n += 1;
                    CountKeeper(ref count, n, ints, ref keepGoing);
                }//end if
                else
                {
                    keepGoing = false;
                }//end else
            }//end if
        }//end method countKeeper

        /// <summary>
        /// Function to generate hashcode string for BasisFunction object
        /// </summary>
        /// <param name="Modes">
        /// List of BasisByMode objects representing modes in the BasisFunction
        /// </param>
        /// <param name="Lambda">
        /// Value of Lambda
        /// </param>
        /// <returns>
        /// String to be used as key for position lookup in dictionary
        /// </returns>
        private static string GenerateHashCode(List<BasisByMode> Modes, int Lambda)
        {
            string s = "";
            for (int i = 0; i < Modes.Count; i++)
            {
                s += Modes[i].V + " ";
                s += Modes[i].L + " ";
            }
            s += Lambda;
            return s;
        }

        /// <summary>
        /// Public interface to find hashcode of a given basis function
        /// </summary>
        /// <param name="a">
        /// Basis function for which the hashcode is desired
        /// </param>
        /// <returns>
        /// Key of basisfunction as a string
        /// </returns>
        public static string GenerateHashCode(BasisFunction a)
        {
            return GenerateHashCode(a.modesInVec, a.Lambda);
        }

        /// <summary>
        /// Function to generate necessary HashCode based on a vlLambda array and the number of modes
        /// </summary>
        /// <param name="vlLambda">
        /// Array with values of v and l for each mode and lambda and J
        /// </param>
        /// <param name="nModes">
        /// Number of modes in the calculation
        /// </param>
        /// <returns>
        /// String representing the needed key.
        /// </returns>
        /// <order>
        /// Boolean. True (or blank) if vlLambda is in v1, v2, ....., l1, l2, ....., Lambda order. False if in v1, l1, v2, l2, ....., Lambda order.
        /// </order>
        public static string GenerateHashCode(int[] vlLambda, int nModes, bool order = true)
        {
            string s = "";
            if (order)
            {
                for (int i = 0; i < nModes; i++)
                {
                    s += vlLambda[i] + " ";
                    s += vlLambda[i + nModes] + " ";
                }                
            }
            else
            {
                for (int i = 0; i < 2 * nModes; i += 2)
                {
                    s += vlLambda[i] + " "; //v
                    s += vlLambda[i + 1] + " "; //l
                }
            }
            s += vlLambda[nModes * 2]; //lambda
            return s;
        }
    }//end class JBasisVector
}
