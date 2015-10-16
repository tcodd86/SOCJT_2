using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Seed
    {
        /// <summary>
        /// Array that stores v, l, Lambda for the seed in the order v1, l1, ... , Lambda
        /// </summary>
        public List<int[]> vlLambdaSeed = new List<int[]>();

        /// <summary>
        /// Number of positions to set to nonzero values in seed.
        /// </summary>
        public int SeedIndex { get; private set; }

        private Object thisLock = new Object();

        public Seed(string SeedFile, int nModes)
        {
            lock (thisLock)
            {
                SeedIndex = 0;
                int ii = 0;
                string line;
                int[] tmpArray = new int[2 * nModes + 1];

                System.IO.StreamReader SeedVectorFile = new System.IO.StreamReader(SeedFile);

                while ((line = SeedVectorFile.ReadLine()) != null)
                {
                    if (line == "/")
                    {
                        ii = 0;
                        vlLambdaSeed.Add(tmpArray);
                        SeedIndex++;
                        tmpArray = new int[2 * nModes + 1]; // Recreate array so that a new array is referenced.
                        continue;
                    }
                    tmpArray[ii] = int.Parse(line); // Array of nonzero elements in the seed vector. vlLambdaSeed[0] is one position, vlLambdaSeed[1] is another, and so on...
                    ii++;
                }
            } // End lock
        } // end Seed
    } // end Class
}