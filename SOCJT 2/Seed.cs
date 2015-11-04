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

        /// <summary>
        /// Array that stores the relative value of each seed.
        /// </summary>
        public List<double> SeedValue = new List<double>();

        private Object thisLock = new Object();

        public Seed(string SeedFile, int nModes)
        {
            lock (thisLock)
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
                        SeedValue.Add(Convert.ToDouble(FullSeedFile[i + 1]));
                        ValueSkip = true;
                        continue;
                    }
                    if (ValueSkip)
                    {
                        ValueSkip = false;
                        continue;
                    }
                    if(FullSeedFile[i] == "/")
                    {
                        ii = 0;
                        vlLambdaSeed.Add(tmpArray);
                        SeedIndex++;
                        tmpArray = new int[2 * nModes + 1];
                        continue;
                    }
                    tmpArray[ii] = int.Parse(FullSeedFile[i]);
                    ii++;
                }

                //SeedIndex = 0;
                //int ii = 0;
                //string line;
                //int[] tmpArray = new int[2 * nModes + 1];

                //System.IO.StreamReader SeedVectorFile = new System.IO.StreamReader(SeedFile);

                //while ((line = SeedVectorFile.ReadLine()) != null)
                //{
                //    if (line == "/")
                //    {
                //        ii = 0;
                //        vlLambdaSeed.Add(tmpArray); // Array of nonzero elements in the seed vector. vlLambdaSeed[0] is one position, vlLambdaSeed[1] is another, and so on...
                //        SeedIndex++;
                //        tmpArray = new int[2 * nModes + 1]; // Recreate array so that a new array is referenced.
                //        continue;
                //    }
                //    try
                //    {
                //        tmpArray[ii] = int.Parse(line);
                //    }
                //    catch (IndexOutOfRangeException)
                //    {
                //        throw new Exception("Improper number of quantum numbers in seed vector.");
                //    }
                //    ii++;
                //}
            } // End lock
        } // end Seed
    } // end Class
}