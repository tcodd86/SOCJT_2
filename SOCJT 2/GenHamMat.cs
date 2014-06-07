using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections.Concurrent;

namespace ConsoleApplication1
{
    static class GenHamMat
    {
        public static List<Dictionary<string, int>> basisPositions;

        /// <summary>
        /// Generates a list of alglib sparsematrix objects where the first is the diagonal elements and the remaining ones are for D and K for each mode and then the matrices for the cross terms.
        /// </summary>
        /// <param name="basisVectorsByJ">
        /// List of JBasisVectors sorted by J.
        /// </param>
        /// <param name="isQuad">
        /// True means a quadratic basis set is used and quadratic terms will be calculated.
        /// </param>
        /// <param name="input">
        /// FileInfo object containing input information.
        /// </param>
        /// <param name="nColumns">
        /// Integer with order of A.
        /// </param>
        /// <param name="par">
        /// Value indicating max degree of parallelism for matrix generation loop.
        /// </param>
        /// <param name="diagOnly">
        /// Boolean to see if only the diagonal elements should be generated
        /// </param>
        /// <returns>
        /// List of sparsematrix objects corresponding to different parameters.
        /// </returns>
        public static List<alglib.sparsematrix> GenMatrixHash(List<BasisFunction> basisVectorsByJ, bool isQuad, FileInfo input, out int nColumns, int par, bool diagOnly, int position)
        {
            int matSize = basisVectorsByJ.Count;
            nColumns = matSize;
            bool bilinear = false;
            int nModes = input.nModes;

            //Lists to store the positions of the A and E vecs with cross-terms coupling
            List<int> biAVecPos = new List<int>();
            List<int> biEVecPos = new List<int>();
            List<int> eVecPos = new List<int>();
            List<int> aVecPos = new List<int>();
            //List to store the alglib sparse matrices for each variable
            var matList = new List<alglib.sparsematrix>();
            //The Tuple stores the row and column of the matrix element in the two int values and the value at that matrix element as the double
            //order of row and column is unimportant since the matrix is symmetric
            List<ConcurrentBag<Tuple<int, int, double>>> matrixPos = new List<ConcurrentBag<Tuple<int, int, double>>>();

            //this array stores the v and l values for each mode for each basis function as well as Lambda and J
            //all v values are stored in elements 0 through nmodes - 1 and l is in nmodes through 2*nmodes - 1
            //Lambda is stored in element 2*nmodes and J is stored as an int as (J - 0.5) in 2 * nmodes + 1
            //the nested loops initialize the vlLambda arrays for each basis function
            //This is ugly but makes the matrix generation much faster than accessing the objects in the matrix generation and saves me a ton of re-coding
            int[,] vlLambda = new int[matSize, nModes * 2 + 2];
            int[] hashStorage = new int[nModes * 2 + 1];
            //Dictionary<string, int> basisPositions = new Dictionary<string, int>();

            //basisPositions[position] = new Dictionary<string, int>();
            for (int i = 0; i < matSize; i++)
            {
                for (int j = 0; j < nModes; j++)
                {
                    vlLambda[i, j] = hashStorage[j] = basisVectorsByJ[i].modesInVec[j].V;
                    vlLambda[i, j + nModes] = hashStorage[j + nModes] = basisVectorsByJ[i].modesInVec[j].L;
                }
                vlLambda[i, nModes * 2] = hashStorage[nModes * 2] = basisVectorsByJ[i].Lambda;
                vlLambda[i, nModes * 2 + 1] = (int)(basisVectorsByJ[i].J - 0.5M);
                basisPositions[position].Add(BasisFunction.GenerateHashCode(hashStorage, nModes), i);
            }//end loop to make vlLambda

            //generate an array to store omega, omegaExe, D, K and degeneracy of each mode
            //only one copy of this is needed since those values are the same across all basis functions
            var modeVals = new double[nModes, 5];
            for (int i = 0; i < nModes; i++)
            {
                modeVals[i, 0] = basisVectorsByJ[0].modesInVec[i].ModeOmega;
                modeVals[i, 1] = basisVectorsByJ[0].modesInVec[i].Anharmonicity;
                modeVals[i, 2] = basisVectorsByJ[0].modesInVec[i].DBasis;
                modeVals[i, 3] = basisVectorsByJ[0].modesInVec[i].KBasis;
                modeVals[i, 4] = 2.0;//degeneracy = 2 by default
                if (basisVectorsByJ[0].modesInVec[i].SymmetryIsA)
                {
                    modeVals[i, 4] = 1.0;
                    aVecPos.Add(i);
                }
                else
                {
                    modeVals[i, 4] = 2.0;
                    eVecPos.Add(i);
                }
            }//end loop to make modeVals[,] array

            #region HO Terms
            //just the harmonic oscillator terms.  Here only the diagonal elements are generated
            var diag = new alglib.sparsematrix();
            alglib.sparsecreate(matSize, matSize, matSize, out diag);

            //loop through all diagonal elements
            for (int n = 0; n < matSize; n++)
            {
                //one mode harmonic and anharmonic terms
                for (int i = 0; i < input.nModes; i++)
                {
                    double temp = modeVals[i, 0] * ((double)vlLambda[n, i] + (double)modeVals[i, 4] / 2.0) - modeVals[i, 1] * Math.Pow((vlLambda[n, i] + (double)modeVals[i, 4] / 2.0), 2.0);
                    alglib.sparseadd(diag, n, n, temp);
                }//end loop over modes
                continue;
            }//end HO loop
            //add diagonal matrix elements to matList as element [0]
            matList.Add(diag);
            #endregion

            //if this is not the first call to the Hamiltonian generation function then only the diagonal elements are needed
            //therefore, if diagOnly is true, just return matList with only the diagonal portion of the matrix.
            if (diagOnly)
            {
                return matList;
            }

            //loop through each mode to see if a sparse matrix is needed for D and K
            for (int i = 0; i < nModes; i++)
            {
                if (!basisVectorsByJ[0].modesInVec[i].SymmetryIsA || basisVectorsByJ[0].modesInVec[i].DBasis != 0.0)
                {
                    //means this mode is degenerate, add matrices for D and K
                    alglib.sparsematrix tempMat = new alglib.sparsematrix();
                    alglib.sparsecreate(matSize, matSize, matSize * (nModes + 1), out tempMat);
                    matList.Add(tempMat);
                    alglib.sparsematrix tempMat2 = new alglib.sparsematrix();
                    alglib.sparsecreate(matSize, matSize, matSize * (nModes + 1), out tempMat2);
                    matList.Add(tempMat2);
                }
                else
                {
                    //add empty matrices for D and K for non degenerate modes. Probably not a huge deal but saves some memory
                    alglib.sparsematrix tempMat = new alglib.sparsematrix();
                    alglib.sparsecreate(matSize, matSize, 0, out tempMat);
                    matList.Add(tempMat);
                    alglib.sparsematrix tempMat2 = new alglib.sparsematrix();
                    alglib.sparsecreate(matSize, matSize, 0, out tempMat2);
                    matList.Add(tempMat2);
                }
            }

            //initialize cross-terms and generate biAVecPos and biEVecPos lists
            BilinearInitialization(basisVectorsByJ[0].modesInVec, nModes, out bilinear, out biAVecPos, out biEVecPos, input.CrossTermMatrix);
            bool crossQuad = false;
            var crossQuadPos = CrossQuadraticInitialization(eVecPos, out crossQuad, input.CrossTermMatrix);

            //add any matrices needed for cross-terms
            if (input.CrossTermMatrix != null)
            {
                for (int i = 0; i < nModes; i++)
                {
                    for (int j = 0; j < nModes; j++)
                    {
                        //add a new sparsematrix for each nonzero cross-term element
                        if (input.CrossTermMatrix[i, j] != 0.0)
                        {
                            alglib.sparsematrix tempMat = new alglib.sparsematrix();
                            alglib.sparsecreate(matSize, matSize, matSize * (nModes + 1), out tempMat);
                            matList.Add(tempMat);
                        }
                    }
                }
            }//enf if crossTermMatrix == null

            //initialize the concurrentBags for each matrix being generated. -1 in loop bounds because diagonal already done
            for (int n = 0; n < matList.Count - 1; n++)
            {
                matrixPos.Add(new ConcurrentBag<Tuple<int, int, double>>());
            }

            //set up the settings for the parallel foreach loop
            var rangePartitioner = Partitioner.Create(0, matSize);
            ParallelOptions parOp = new ParallelOptions();
            parOp.MaxDegreeOfParallelism = par;
            Parallel.ForEach(rangePartitioner, parOp, (range, loopState) =>
            {
                //indexes n and m are for the rows and columns of the matrix respectively
                for (int n = range.Item1; n < range.Item2; n++)
                {
                    int[] hashInt = new int[nModes * 2 + 1];
                    for (int a = 0; a < hashInt.Length; a++)
                    {
                        hashInt[a] = vlLambda[n, a];
                    }
                    int[] tempInt;
                    int m;
                    double temp;
                    string hashCode;
                    for (int aModes = 0; aModes < aVecPos.Count; aModes++)
                    {
                        #region SymmetricGradient
                        for (int deltaV = -1; deltaV < 2; deltaV += 2)
                        {
                            tempInt = (int[])hashInt.Clone();
                            DeltaVL(ref tempInt, aVecPos[aModes], deltaV, 0, nModes, false);
                            hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                            if (basisPositions[position].TryGetValue(hashCode, out m))
                            { 
                                //assign temp here.
                                temp = (double)vlLambda[n, aVecPos[aModes]];
                                if (deltaV == 1)
                                {
                                    temp += 1.0;
                                }
                                temp = Math.Sqrt(temp);
                                Tuple<int, int, double> ttTemp = new Tuple<int, int, double>(n, m, temp);// basisVectorsByJ[n].modesInVec[mode].v     basisVectorsByJ[n].modesInVec[mode].l
                                matrixPos[2 * aVecPos[aModes]].Add(ttTemp);
                            }
                        }
                        #endregion
                    }
                    for (int eModes = 0; eModes < eVecPos.Count; eModes++)
                    {
                        #region Linear
                        //linear portion
                        //now check each EVecPos mode for a linear element by taking vlLambda[n, ] and changing the values appropriately to find a linear element. Store in tempint
                        //tempInt = (int[])hashInt.Clone();                         
                        for(int deltal = -1; deltal < 2; deltal += 2)
                        {
                            for (int deltaV = -1; deltaV < 2; deltaV += 2)
                            {
                                tempInt = (int[])hashInt.Clone();
                                DeltaVL(ref tempInt, eVecPos[eModes], deltaV, deltal, nModes);
                                hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                                //if it exists, then assign it to the linear value                            
                                if (basisPositions[position].TryGetValue(hashCode, out m))
                                {
                                    if (m > n)
                                    {
                                        temp = LinearMatrixElement(vlLambda, n, eVecPos[eModes], nModes, deltal, deltaV);
                                        Tuple<int, int, double> ttTemp = new Tuple<int, int, double>(n, m, temp);// basisVectorsByJ[n].modesInVec[mode].v     basisVectorsByJ[n].modesInVec[mode].l
                                        matrixPos[2 * eVecPos[eModes]].Add(ttTemp);
                                    }
                                }
                            }//end loop over deltaV
                        }//end loop over linear l values
                        #endregion

                        #region Quadratic
                        //check for all quadratic elements
                        //reset tempInt values for Quadratic elements
                        for (int ll = -2; ll < 3; ll += 4)
                        {
                            tempInt = (int[])hashInt.Clone();
                            //for bottom matrix element on page
                            DeltaVL(ref tempInt, eVecPos[eModes], 2, ll, nModes);
                            hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                            if (basisPositions[position].TryGetValue(hashCode, out m))
                            {
                                if (m > n)
                                {
                                    temp = (1 / 4D * Math.Sqrt((vlLambda[n, eVecPos[eModes]] + (ll / 2) * vlLambda[n, nModes + eVecPos[eModes]] + 4D) * (vlLambda[n, eVecPos[eModes]] + (ll / 2) * vlLambda[n, nModes + eVecPos[eModes]] + 2)));
                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                    matrixPos[eVecPos[eModes] * 2 + 1].Add(tTemp);
                                }
                            }

                            //reset tempInt values for Quadratic elements
                            tempInt = (int[])hashInt.Clone();
                            //for middle matrix element on page
                            DeltaVL(ref tempInt, eVecPos[eModes], 0, ll, nModes);
                            hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                            if (basisPositions[position].TryGetValue(hashCode, out m))
                            {
                                if (m > n)
                                {
                                    temp = (1D / 2D * Math.Sqrt((vlLambda[n, eVecPos[eModes]] + (ll / 2) * vlLambda[n, nModes + eVecPos[eModes]] + 2) * (vlLambda[n, eVecPos[eModes]] - (ll / 2) * vlLambda[n, eVecPos[eModes] + nModes])));
                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                    matrixPos[eVecPos[eModes] * 2 + 1].Add(tTemp);
                                }
                            }

                            tempInt = (int[])hashInt.Clone();
                            //for top matrix element on page
                            DeltaVL(ref tempInt, eVecPos[eModes], -2, ll, nModes);
                            hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                            if (basisPositions[position].TryGetValue(hashCode, out m))
                            {
                                if (m > n)
                                {
                                    temp = (1 / 4D * Math.Sqrt((vlLambda[n, eVecPos[eModes]] - (ll / 2) * vlLambda[n, nModes + eVecPos[eModes]]) * (vlLambda[n, eVecPos[eModes]] - (ll / 2) * vlLambda[n, nModes + eVecPos[eModes]] - 2)));
                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                    matrixPos[eVecPos[eModes] * 2 + 1].Add(tTemp);
                                }
                            }
                        }//end quadratic loop over l values
                        #endregion
                    }//end loop over E-modes

                    int crossCount = 0;
                    #region Bilinear
                    if (bilinear)
                    {
                        //loop through A and E vecs with possible coupling
                        for (int a = 0; a < biAVecPos.Count; a++)
                        {
                            for (int e = 0; e < biEVecPos.Count; e++)
                            {
                                //value to keep track of which cross-term matrix we're on.
                                crossCount = a + e;

                                //this is because the cross-term couplings for JT terms are stored in the upper-diagonal of the cross-term matrix
                                int row;
                                int column;
                                if (biAVecPos[a] > biEVecPos[e])
                                {
                                    column = biAVecPos[a];
                                    row = biEVecPos[e];
                                }
                                else
                                {
                                    column = biEVecPos[e];
                                    row = biAVecPos[a];
                                }
                                for (int ll = -1; ll < 2; ll += 2)
                                {
                                    for (int dve = -1; dve < 2; dve += 2)
                                    {
                                        for (int dva = -1; dva < 2; dva += 2)
                                        {
                                            double oneORnone = 0.0;
                                            if (dva == 1)
                                            {
                                                oneORnone = 1.0;
                                            }
                                            double twoORnone = 0.0;
                                            int pre = -1;
                                            if (dve == 1)
                                            {
                                                twoORnone = 2.0;
                                                pre = 1;
                                            }
                                            tempInt = (int[])hashInt.Clone();
                                            //delta vl for e mode
                                            DeltaVL(ref tempInt, biEVecPos[e], dve, ll, nModes);
                                            //delta vl for a mode
                                            DeltaVL(ref tempInt, biAVecPos[a], dva, 0, nModes, false);
                                            hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                                            //formula for bilinear matrix element
                                            if (basisPositions[position].TryGetValue(hashCode, out m))
                                            {
                                                if (m > n)
                                                {
                                                    temp = 0.5 * Math.Sqrt(((double)vlLambda[n, biAVecPos[a]] + oneORnone) * ((double)vlLambda[n, biEVecPos[e]] + ll * pre * vlLambda[n, nModes + biEVecPos[e]] + twoORnone));//ll * (double)vlLambda[n, EVecPos[eModes] + nModes] + 2D))
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matrixPos[2 * nModes + crossCount].Add(tTemp);
                                                }
                                            }
                                        }//end loop over delta v values for a mode
                                    }//end loop over delta v values for e mode
                                }//end loop over l values for e mode
                            }//end for loop over evec positions
                        }//end for loop over a vec positions
                    }//end bilinear if
                    #endregion

                    #region Cross-Quadratic
                    if (crossQuad)
                    {
                        for (int crossTerm = 0; crossTerm < crossQuadPos.Count; crossTerm += 2)
                        {
                            for (int deltal = -1; deltal < 2; deltal += 2)
                            {
                                for (int deltaV = -1; deltaV < 2; deltaV += 2)
                                {                                    
                                    for (int deltal2 = -1; deltal2 < 2; deltal2 += 2)
                                    {
                                        for (int deltaV2 = -1; deltaV2 < 2; deltaV2 += 2)
                                        {
                                            for (int changeLambda = 0; changeLambda < 2; changeLambda++)
                                            {
                                                //tbool makes it so that every other iteration of changeLambda is electronically off-diagonal
                                                bool tbool = false;
                                                if (changeLambda == 0)
                                                {
                                                    tbool = true;
                                                }
                                                tempInt = (int[])hashInt.Clone();
                                                //get DeltaVL for the first mode
                                                DeltaVL(ref tempInt, eVecPos[crossQuadPos[crossTerm]], deltaV, deltal, nModes, tbool);
                                                //get DeltaVL for the second mode
                                                DeltaVL(ref tempInt, eVecPos[crossQuadPos[crossTerm + 1]], deltaV2, deltal2, nModes);
                                                //generate the hashcode
                                                hashCode = BasisFunction.GenerateHashCode(tempInt, nModes);
                                                //if it exists, then assign it to the linear value                            
                                                if (basisPositions[position].TryGetValue(hashCode, out m))
                                                {
                                                    if (m > n)
                                                    {
                                                        temp = LinearMatrixElement(vlLambda, n, eVecPos[crossQuadPos[crossTerm]], nModes, deltal, deltaV);
                                                        temp *= LinearMatrixElement(vlLambda, n, eVecPos[crossQuadPos[crossTerm + 1]], nModes, deltal, deltaV);
                                                        Tuple<int, int, double> ttTemp = new Tuple<int, int, double>(n, m, temp);// basisVectorsByJ[n].modesInVec[mode].v     basisVectorsByJ[n].modesInVec[mode].l
                                                        matrixPos[2 * nModes + crossCount].Add(ttTemp);
                                                    }
                                                }
                                            }//end loop over change lambda
                                        }//end loop over deltaV
                                    }//end loop over linear l values
                                }//end loop over deltaV
                            }//end loop over linear l values
                            crossCount++;
                        }//end loop over cross quadratic terms
                    }
                    #endregion
                }//row for loop
            }//end anonymous function in parallel for loop
            );//end parallel for
            
            //actually add all of the matrix elements to the matrices
            //Start at 0 because matrixPos only has off diagonal elements.  Add elements to matList[i + 1] because matList already has diagonal elements in position 0.
            for (int i = 0; i < matrixPos.Count; i++)
            {
                //add all calculated matrix elements to the appropriate sparsematrices
                foreach (Tuple<int, int, double> spot in matrixPos[i])
                {
                    alglib.sparseadd(matList[i + 1], spot.Item1, spot.Item2, spot.Item3);
                }
            }
            return matList;
        }//end method genMatrix

        /// <summary>
        /// Generates the value for a linear matrix element.
        /// </summary>
        /// <param name="vlLambda">
        /// v, l, and Lambda for this mode
        /// </param>
        /// <param name="n">
        /// Which basis function is being used
        /// </param>
        /// <param name="index">
        /// Which mode is being used.
        /// </param>
        /// <param name="nModes">
        /// Number of modes
        /// </param>
        /// <param name="deltal">
        /// What the delta l is
        /// </param>
        /// <param name="deltaV">
        /// What the delta v is
        /// </param>
        /// <returns>
        /// Double value for everything inside the radical
        /// </returns>
        private static double LinearMatrixElement(int[,] vlLambda, int n, int index, int nModes, int deltal, int deltaV)
        {
            int plus = 0;
            if (deltaV == 1)
            {
                plus = 2;
            }
            return Math.Sqrt((double)vlLambda[n, index] + deltaV * deltal * (double)vlLambda[n, index + nModes] + plus);
        }

        /// <summary>
        /// Changes appropriate values for finding a given matrix element
        /// </summary>
        /// <param name="vlArray">
        /// Integer array of basis function values
        /// </param>
        /// <param name="Mode">
        /// Which mode is being changed
        /// </param>
        /// <param name="deltaV">
        /// Delta v
        /// </param>
        /// <param name="deltaL">
        /// Delta l
        /// </param>
        /// <param name="nModes">
        /// Number of Modes.
        /// </param>
        /// <param name="switchLambda">
        /// Boolean to say whether the sign of lambda should be changed or not.
        /// </param>
        private static void DeltaVL(ref int[] vlArray, int Mode, int deltaV, int deltaL, int nModes, bool switchLambda = true)
        {
            vlArray[Mode] += deltaV;
            vlArray[Mode + nModes] += deltaL;
            if (switchLambda)
            {
                vlArray[2 * nModes] *= -1;
            }
        }

        /// <summary>
        /// Generates a list containing the positions of modes with cross-quadratic coupling.
        /// </summary>
        /// <param name="evecPosition">
        /// List containing the positions of degenerate (e type) modes
        /// </param>
        /// <param name="crossQuadratic">
        /// Boolean to be set if there are cross-quadratic terms
        /// </param>
        /// <param name="crossTermMatrix">
        /// Cross-term matrix
        /// </param>
        /// <returns>
        /// List whose elements are the modes with cross-quadratic coupling.
        /// </returns>
        public static List<int> CrossQuadraticInitialization(List<int> evecPosition, out bool crossQuadratic, double[,] crossTermMatrix)
        {
            crossQuadratic = false;
            var crossQuadModes = new List<int>();
            //check the upper diagonal of the crossterm matrix to see if any emodes are coupled to one another
            for (int i = 0; i < evecPosition.Count; i++)
            {
                for (int j = i + 1; j < evecPosition.Count; j++)
                { 
                    //check here to see if the cross term matrix element between evecPosition[i] and evecPosition[j] is 0
                    if (crossTermMatrix[evecPosition[i], evecPosition[j]] != 0.0)
                    {
                        crossQuadModes.Add(evecPosition[i]);
                        crossQuadModes.Add(evecPosition[j]);
                        crossQuadratic = true;
                    }
                }
            }
            return crossQuadModes;
        }//end method CrossQuadraticInitialization

        /// <summary>
        /// This function initializes the biAVecPos and biEVecPos lists which tell which A and E vecs have cross-term coupling.  Also finds if there is any bilinear coupling.
        /// </summary>
        /// <param name="modesInVec">
        /// List BasisByMode which is used to find the symmetry of all of the modes.
        /// </param>
        /// <param name="nModes">
        /// Int indicating how many modes are in the input file.
        /// </param>
        /// <param name="bilinear">
        /// Boolean used to indicate if there is any bilinear coupling
        /// </param>
        /// <param name="biAVecPos">
        /// List which will contain the positions (in the mode list) of any A modes with bilinear coupling
        /// </param>
        /// <param name="biEVecPos">
        /// List which will contain the positions (in the mode list) of any E modes with bilinear coupling
        /// </param>
        /// <param name="crossTermMatrix">
        /// Cross term matrix containing coupling the actual coupling terms.
        /// </param>
        public static void BilinearInitialization(List<BasisByMode> modesInVec, int nModes, out bool bilinear, out List<int> biAVecPos, out List<int> biEVecPos, double[,] crossTermMatrix)
        {
            bool containsAVecs = false;
            bilinear = false;
            List<int> AVecPos = new List<int>();
            List<int> EVecPos = new List<int>();
            biAVecPos = new List<int>();
            biEVecPos = new List<int>();
            for (int i = 0; i < nModes; i++)
            {
                if (modesInVec[i].SymmetryIsA == true)
                {
                    AVecPos.Add(i);
                    containsAVecs = true;
                    bilinear = true;
                    continue;
                }
                else
                {
                    EVecPos.Add(i);
                }
            }//end for

            //for bilinear coupling
            if (containsAVecs == true)
            {
                //these lists contain all postions of A and E modes
                for (int i = 0; i < AVecPos.Count; i++)
                {
                    biAVecPos.Add(AVecPos[i]);
                }
                for (int i = 0; i < EVecPos.Count; i++)
                {
                    biEVecPos.Add(EVecPos[i]);
                }

                //loop to eliminate any A modes that have no cross-coupling from the lists
                //after this and the following loops there should be no evecs or avecs in these lists which do not have bilinear coupling
                for (int i = 0; i < biAVecPos.Count; i++)
                {
                    for (int j = 0; j < biEVecPos.Count; j++)
                    {
                        //to make sure to only use the correct positions in the cross-term matrix
                        if (biAVecPos[i] > biEVecPos[j])
                        {
                            //if the coupling element is nonzero then no need to remove position so break and skip
                            if (crossTermMatrix[biEVecPos[j], biAVecPos[i]] != 0)
                            {
                                break;
                            }
                        }
                        else
                        {
                            if (crossTermMatrix[biAVecPos[i], biEVecPos[j]] != 0)
                            {
                                break;
                            }
                        }
                        //means it's reached the end of the evec list and has not found a coupling term between it and an A vector
                        if (j == biEVecPos.Count - 1)
                        {
                            //remove the A vector from the list if there's no coupling term
                            biAVecPos.RemoveAt(i);
                            i--;
                            break;
                        }
                    }
                }

                //loop to eliminate any E vectors that have no cross-coupling
                for (int i = 0; i < biEVecPos.Count; i++)
                {
                    for (int j = 0; j < biAVecPos.Count; j++)
                    {
                        if (biAVecPos[j] > biEVecPos[i])
                        {
                            if (crossTermMatrix[biEVecPos[i], biAVecPos[j]] != 0)
                            {
                                break;
                            }
                        }
                        else
                        {
                            if (crossTermMatrix[biAVecPos[j], biEVecPos[i]] != 0)
                            {
                                break;
                            }
                        }
                        if (j == biAVecPos.Count - 1)
                        {
                            biEVecPos.RemoveAt(i);
                            i--;
                            break;
                        }
                    }
                }
                //if there's no A or E modes left then bilinear = false;
                if (biAVecPos.Count == 0 || biEVecPos.Count == 0)
                {
                    bilinear = false;
                    return;
                }
            }
        }//end crossTermInitialization

        /// <summary>
        /// Reads through a List of JBasisVectors and returns those with a specified j value in a new List.
        /// </summary>
        /// <param name="jVecs">
        /// List of all JBasisVectors.
        /// </param>
        /// <param name="J">
        /// Desired j value.
        /// </param>
        /// <returns>
        /// List of JBasisVectors with the specified value of j.
        /// </returns>
        public static List<BasisFunction> SortByJ(List<BasisFunction> jVecs, decimal J)
        {
            List<BasisFunction> outList = new List<BasisFunction>();
            //count = 0;
            foreach (BasisFunction vector in jVecs)
            {
                if (vector.J == J)
                {
                    outList.Add(vector);
                }//end if
            }//end foreach
            return outList;
        }//end method sortByJ
    }//end class GenHamMat
}
