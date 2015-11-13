using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Seed
    {
        /// <summary>
        /// List of Arrays. Each array stores v, l, Lambda for a basis function in the seed in the order v1, l1, ... , Lambda. Each list is a new basis function in the expansion.
        /// </summary>
        public List<int[]> vlLambdaSeed = new List<int[]>();

        /// <summary>
        /// Number of positions to set to nonzero values in seed.
        /// </summary>
        public int SeedIndex { get; private set; }

        /// <summary>
        /// List that stores the coefficient on each basis function in the seed vector. 
        /// </summary>
        public List<double> SeedValue = new List<double>();

        private Object thisLock = new Object();

        public Seed(string SeedFile, int nModes) // This function reads the seed file and stores everything needed to generate the basis positions.
        {
            lock (thisLock) // I don't know if I need this
            {
                SeedIndex = 0;
                int ii = 0;
                string[] FullSeedFile = FileInfo.FileRead(SeedFile);
                int[] tmpArray = new int[2 * nModes + 1];
                bool ValueSkip = false;

                for(int i = 0; i < FullSeedFile.Length; i++)
                {
                    if (FullSeedFile[i].ToUpper() == "VALUE")
                    {
                        SeedValue.Add(Convert.ToDouble(FullSeedFile[i + 1])); // Adds the coefficient to the list
                        ValueSkip = true; // Tells program to skip the "/" line
                        continue;
                    }
                    if (ValueSkip) // Skips the "/" Line
                    {
                        ValueSkip = false;
                        continue;
                    }
                    if(FullSeedFile[i] == "/")
                    {
                        ii = 0; // Reset index
                        vlLambdaSeed.Add(tmpArray); // Add array to list
                        SeedIndex++; // Count the basis function
                        tmpArray = new int[2 * nModes + 1]; // Free the array to reuse
                        continue;
                    }
                    tmpArray[ii] = int.Parse(FullSeedFile[i]); // Adds each v or l or Lambda to array
                    ii++; // Raises index
                }
            } // End lock
        } // end Seed
    } // end Class
}