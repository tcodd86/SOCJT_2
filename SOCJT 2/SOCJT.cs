using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Collections;

namespace ConsoleApplication1
{
    class SOCJT
    {
        private static bool matricesMade = false;

        //keep this so that on concurrent calls the matrix does not need to be regenerated
        public static List<List<alglib.sparsematrix>> fitHamList { get; private set; }

        public List<List<BasisFunction>> basisSet { get; private set; }

        private Eigenvalue[] nfinalList;
        public Eigenvalue[] finalList
        {
            get { return nfinalList; }
            set { nfinalList = value; }
        }

        private List<string> nOutput;
        public List<string> outp
        {
            get { return nOutput; }
            set { nOutput = value; }
        }

        //this list stores the untransformed eigenvectors
        public List<double[,]> lanczosEVectors { get; private set; }

        public List<string> SOCJTroutine(List<ModeInfo> Modes, bool isQuad, string[] inputFile, FileInfo input)
        {
            //Sets minimum and maximum j values.
            Stopwatch measurer = new Stopwatch();
            long howMuchTime;
            decimal jMin;
            decimal jMax;
            //If is quadratic makes sure that the maxJ is at least 7.5
            if (isQuad == true)
            {
                if (input.maxJ < 7.5M)
                {
                    input.maxJ = 7.5M;
                }
            }
            if (isQuad == true)
            {
                jMax = input.maxJ;
                jMin = -input.maxJ;
            }//end if
            else
            {
                jMax = input.maxJ;
                if (input.MinJBool == true)
                {
                    jMin = input.minJ;
                }
                else
                {
                    jMin = 0.5M;
                }
            }//end if

            //Makes a List of Lists of Basis objects with each List of Basis objects being for one mode.
            List<List<BasisByMode>> basisByMode = new List<List<BasisByMode>>();
            for (int i = 0; i < input.nModes; i++)
            {
                basisByMode.Add(new List<BasisByMode>());
                basisByMode[i] = BasisByMode.GenVLCombinations(Modes[i], i);
            }//end for            

            //Generates all of the JBasisVectors to be used in calculation.
            List<BasisFunction> hamiltonianVecs = BasisFunction.GenJVecs(basisByMode, input.nModes, jMin, jMax);
            //hamiltonianVecs are the total list of all J vectors in the Hamiltonian

            //Sorts the hamiltonianVecs by J and puts them into a List of Lists of JBasisVectors.
            List<List<BasisFunction>> jBasisVecsByJ = new List<List<BasisFunction>>();
            for (decimal i = jMin; i <= jMax; i++)
            {
                jBasisVecsByJ.Add(GenHamMat.SortByJ(hamiltonianVecs, i));
            }//end for loop

            //Initializes Lists to hold the Hamiltonian matrices, eigenvectors, eigenvalues, basis vectors for the output file and number of columns for each j matrix respectively.
            //Also initializes the Dictionary list for storing the positions of the basis functions
            List<double[,]> zMatrices = new List<double[,]>();
            List<double[]> eigenvalues = new List<double[]>();
            List<List<BasisFunction>> JvecsForOutuput = new List<List<BasisFunction>>();
            List<BasisFunction>[] jbasisoutA = new List<BasisFunction>[0];
            List<int> numColumns = new List<int>();
            GenHamMat.basisPositions = new List<Dictionary<string, int>>();

            #region Hamiltonian
            List<alglib.sparsematrix> sHamMatrix = new List<alglib.sparsematrix>();
            alglib.sparsematrix[] array1;
            int[] numcolumnsA;
            //smallMat = false;            
            //Creates the Hamiltonian matrices for linear cases            
            int numQuadMatrix = 0;
            List<int> a = new List<int>();

            if (input.M > input.NumberOfIts)
            {
                throw new BasisSetTooSmallException(false);
            }

            if (isQuad == false)
            {
                #region LinearHamiltonian
                measurer.Reset();
                measurer.Start();
                int h = 0;
                array1 = new alglib.sparsematrix[jBasisVecsByJ.Count];
                List<int> basisSize = new List<int>();
                basisSet = jBasisVecsByJ;
                if (!matricesMade)
                {
                    fitHamList = new List<List<alglib.sparsematrix>>();
                    for (int m = 0; m < jBasisVecsByJ.Count; m++)
                    {
                        fitHamList.Add(new List<alglib.sparsematrix>());
                    }
                    //check here to see if matrix file should be used, if so then read from file and make them, set matricesMade to true
                    //matricesMade = matReadFunction(input);
                    basisSize = matReadFunction(input, ref matricesMade);
                }
                numcolumnsA = new int[jBasisVecsByJ.Count];


                //put items up to length in here
                for (int i = (int)(jMin - 0.5M); i < (int)(jMax + 0.5M); i++)
                {
                    GenHamMat.basisPositions.Add(new Dictionary<string, int>());
                }

                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.ParJ;
                try
                {
                    Parallel.For((int)(jMin - 0.5M), (int)(jMax + 0.5M), options, i =>
                    {
                        int nColumns;
                        if (jBasisVecsByJ[i].Count != 0)//changed from h to i                    
                        {
                            //this checks if the matrix was read from file, and if so if the basis set is the correct size
                            if (basisSize.Count != 0)
                            {
                                if (basisSize[i] != jBasisVecsByJ[i].Count)
                                {
                                    throw new MatrixFileError();
                                }
                            }
                            if (!matricesMade)//if matrices not made then generate all matrices
                            {
                                fitHamList[i] = GenHamMat.GenMatrixHash(jBasisVecsByJ[i], isQuad, input, out nColumns, input.ParMatrix, false, i);
                            }
                            else//this makes sure that the diagonal portion is regenerated on each call.
                            {
                                fitHamList[i][0] = GenHamMat.GenMatrixHash(jBasisVecsByJ[i], isQuad, input, out nColumns, input.ParMatrix, true, i)[0];
                            }
                            numcolumnsA[i] = nColumns;
                            if (numcolumnsA[i] < input.M)
                            {
                                a.Add(0);
                            }
                            h++;
                        }
                        else
                        {
                            a.Add(0);
                        }
                    }//end for loop                                    
                    );
                }
                catch (AggregateException ae)
                {
                    foreach (var e in ae.InnerExceptions)
                    {
                        if (e is MatrixFileError)
                        {
                            throw new MatrixFileError();
                        }
                        else
                        {
                            throw;
                        }
                    }

                }
                measurer.Stop();
                howMuchTime = measurer.ElapsedMilliseconds;
                input.MatrixGenerationTime = (double)howMuchTime / 1000D;
                for (int i = (int)(jMin - 0.5M); i < (int)(jMax + 0.5M); i++)
                {
                    zMatrices.Add(new double[0, 0]);
                    eigenvalues.Add(new double[0]);
                }
                matricesMade = true;
                //handles errors where the basis set is too small
                if (a.Count > 0)
                {
                    throw new BasisSetTooSmallException(true);
                }
                #endregion
            }//end if  

            //Creates the Hamiltonian matrices for quadratic cases.            
            else
            {
                #region QuadraticHamiltonian
                measurer.Reset();
                measurer.Start();
                int dynVar1 = (int)(jMax - 1.5M);
                int dynVar2 = dynVar1 / 3;
                List<int> basisSize = new List<int>();
                array1 = new alglib.sparsematrix[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                if (!matricesMade)
                {
                    fitHamList = new List<List<alglib.sparsematrix>>();
                    for (int m = 0; m < jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2; m++)
                    {
                        fitHamList.Add(new List<alglib.sparsematrix>());
                    }
                    //check here to see if matrix file should be used, if so then read from file and make them, set matricesMade to true
                    basisSize = matReadFunction(input, ref matricesMade);
                }
                numcolumnsA = new int[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                jbasisoutA = new List<BasisFunction>[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6

                for (int i = jBasisVecsByJ.Count / 2; i < jBasisVecsByJ.Count - dynVar1; i++)
                {
                    GenHamMat.basisPositions.Add(new Dictionary<string, int>());
                }

                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.ParJ;
                try
                {
                    Parallel.For(jBasisVecsByJ.Count / 2, jBasisVecsByJ.Count - dynVar1, options, i =>//changed to dynVar1 from 6
                    {
                        List<BasisFunction> quadVecs = new List<BasisFunction>();
                        int nColumns;
                        for (int v = -dynVar2 - 1; v <= dynVar2; v++)
                        {
                            if (i + v * 3 >= 0)
                            {
                                quadVecs.AddRange(jBasisVecsByJ[i + v * 3]);
                            }
                        }
                        //this checks if the matrix was read from file, and if so if the basis set is the correct size
                        if (basisSize.Count != 0)
                        {
                            if (basisSize[i - jBasisVecsByJ.Count / 2] != quadVecs.Count)
                            {
                                throw new MatrixFileError();
                            }
                        }
                        //if matrices aren't made then generate all of them
                        if (!matricesMade)
                        {
                            fitHamList[i - jBasisVecsByJ.Count / 2] = GenHamMat.GenMatrixHash(quadVecs, isQuad, input, out nColumns, input.ParMatrix, false, i - jBasisVecsByJ.Count / 2);
                        }
                        else//If they are made then just generate the diagonal elements.
                        {
                            fitHamList[i - jBasisVecsByJ.Count / 2][0] = GenHamMat.GenMatrixHash(quadVecs, isQuad, input, out nColumns, input.ParMatrix, true, i - jBasisVecsByJ.Count / 2)[0];
                        }

                        jbasisoutA[i - jBasisVecsByJ.Count / 2] = quadVecs;
                        numcolumnsA[i - jBasisVecsByJ.Count / 2] = nColumns;
                        if (numcolumnsA[i - jBasisVecsByJ.Count / 2] < input.M)
                        {
                            a.Add(0);
                        }
                        numQuadMatrix++;
                    }
                    );
                }
                catch (AggregateException ae)
                {
                    foreach (var e in ae.InnerExceptions)
                    {
                        if (e is MatrixFileError)
                        {
                            throw new MatrixFileError();
                        }
                        else
                        {
                            throw;
                        }
                    }

                }
                basisSet = new List<List<BasisFunction>>();
                for (int jj = 0; jj < jbasisoutA.Length; jj++)
                {
                    basisSet.Add(jbasisoutA[jj]);
                }
                matricesMade = true;
                measurer.Stop();
                howMuchTime = measurer.ElapsedMilliseconds;
                input.MatrixGenerationTime = (double)howMuchTime / 1000D;
                if (a.Count > 0)
                {
                    throw new BasisSetTooSmallException(true);
                }
                for (int i = 0; i < 2; i++)
                {
                    zMatrices.Add(new double[0, 0]);
                    eigenvalues.Add(new double[0]);
                }
                #endregion
            }//end else
            #endregion

            //list where each element of mat is a list of alglib.sparsematrix objects.  One for each off-diagonal parameter (D, K, B)
            var mat = new List<List<alglib.sparsematrix>>();
            bool bilinear = false;
            var biAVecPos = new List<int>();
            var biEVecPos = new List<int>();
            GenHamMat.BilinearInitialization(jBasisVecsByJ[0][0].modesInVec, input.nModes, out bilinear, out biAVecPos, out biEVecPos, input.CrossTermMatrix);
            //code here to convert the alglib matrices to matrices for each j block


            for (int i = 0; i < fitHamList.Count; i++)
            {
                int whichCrossMatrix = 0;
                mat.Add(new List<alglib.sparsematrix>());
                int count;
                int DorK;
                double val = 0.0;
                //this adds the diagonal elements to the list
                mat[i].Add(fitHamList[i][0]);
                //go through each member of the list and multiply it by the appropriate value, combine them all
                //skip the first one which is the diagonal elements and isn't multiplied by anything.
                for (int j = 1; j < fitHamList[i].Count; j++)
                {
                    count = (j - 1) / 2;
                    DorK = (j - 1) % 2;
                    if (count < input.nModes)
                    {
                        if (DorK == 0)
                        {
                            if (!Modes[count].IsAType)
                            {
                                val = Math.Sqrt(Modes[count].D) * Modes[count].modeOmega;
                            }
                            else
                            {
                                //val = Math.Sqrt(Modes[count].D);
                                val = Modes[count].D;
                            }
                        }
                        else
                        {
                            val = Modes[count].K * Modes[count].modeOmega;
                        }
                    }
                    else//means it's a cross term. loop over relevant E and A terms in same order as in genFitMatrix function
                    {
                        int crossMatrixCounter = 0;
                        for (int blOrNot = 0; blOrNot < 2; blOrNot++)
                        {
                            for (int row = 0; row < input.CrossTermMatrix.GetLength(0); row++)
                            {
                                for (int column = row + 1; column < input.CrossTermMatrix.GetLength(0); column++)
                                {
                                    if (input.CrossTermMatrix[row, column] != 0.0)
                                    {
                                        //now test to see if this is a bilinear or cross quadratic term. do bilinear first, then cross quadratic.
                                        //then can get the value
                                        //use blOrNot to see if we should be doing bilinear or not
                                        //if this is a bilinear term, only add it if blOrNot == 0
                                        //crossMatrixCounter will keep going 
                                        if ((biAVecPos.Exists(x => x == row) || biAVecPos.Exists(x => x == column)) && blOrNot == 0)
                                        {
                                            //means this is a bilinear term                                            
                                            if (crossMatrixCounter == whichCrossMatrix)
                                            {
                                                val = input.CrossTermMatrix[row, column];
                                            }//end conditional to see if this is the cross-term element we want    
                                            crossMatrixCounter++;
                                        }
                                        if (!(biAVecPos.Exists(x => x == row) || biAVecPos.Exists(x => x == column)) && blOrNot == 1)
                                        {
                                            //means this is a cross-quadratic term
                                            if (crossMatrixCounter == whichCrossMatrix)
                                            {
                                                val = input.CrossTermMatrix[row, column];
                                            }//end conditional to see if this is the cross-term element we want  
                                            crossMatrixCounter++;
                                        }
                                    }//end conditional to see if this matrix element is 0
                                }//end loop over columns of cross-term matrix
                            }//end loop over rows of cross-term matrix
                        }//end loop to go through the cross-term matrix twice
                        whichCrossMatrix++;
                    }//end else for counting if it's D / K or cross-Term                    
                    mat[i].Add(cTimesSparse(fitHamList[i][j], val));
                }//end loop over fitHamList
            }
            //now need to convert the lists of matrices in each mat element into a single object
            //here convert the alglib matrices to the appropriate things.
            for (int i = 0; i < array1.Length; i++)
            {
                array1[i] = aggregator(mat[i]);//some functio to put in aggregate of mat[i] sparsematrices
            }

            //move this to after genFitMatricss are treated so that SO code does not need to be changed
            #region Spin Orbit
            //add SO stuff here
            if (input.IncludeSO == true)
            {
                //Jvecs for output
                //
                //var tempBasisList = new List<List<BasisFunction>>();
                decimal minS = input.S * -1M;
                List<alglib.sparsematrix> tempMatList = new List<alglib.sparsematrix>();
                List<int> tempNumbColumns = new List<int>();

                for (int i = 0; i < array1.Length; i++)
                {
                    for (decimal j = minS; j <= input.S; j++)
                    {
                        alglib.sparsematrix tempMat = new alglib.sparsematrix();
                        alglib.sparsecopy(array1[i], out tempMat);
                        tempMatList.Add(tempMat);
                        if (isQuad)
                        {
                            JvecsForOutuput.Add(jbasisoutA[i]);
                        }
                        else
                        {
                            JvecsForOutuput.Add(jBasisVecsByJ[i]);
                        }
                        //jbasisoutA for quadratic
                        //jBasisVecsByJ for linear
                        //tempBasisList.Add(JvecsForOutuput[i]);
                        for (int k = 0; k < numcolumnsA[i]; k++)
                        {
                            double temp = input.Azeta * (double)JvecsForOutuput[JvecsForOutuput.Count - 1][k].Lambda * (double)j;
                            alglib.sparseadd(tempMatList[tempMatList.Count - 1], k, k, temp);
                        }//end loop over diagonal matrix elements
                        tempNumbColumns.Add(numcolumnsA[i]);
                        //means SO only in j = 0.5 block for quadratic cases
                        //if (i > 0 && isQuad)
                        //means SO only in degenerate blocks
                        if ((i - 1) % 3 == 0)
                        {
                            break;
                        }
                    }//end loop over values of S
                }//end loop over all previouly made sparseMatrices

                numcolumnsA = null;
                numcolumnsA = tempNumbColumns.ToArray();
                array1 = null;
                array1 = tempMatList.ToArray();
                //JvecsForOutuput.Clear();
                //JvecsForOutuput = tempBasisList;
                zMatrices.Clear();
                eigenvalues.Clear();
                for (int i = 0; i < array1.Length; i++)
                {
                    zMatrices.Add(new double[0, 0]);
                    eigenvalues.Add(new double[0]);
                }
            }//end if inclSO == true
            #endregion

            for (int i = 0; i < array1.Length; i++)
            {
                alglib.sparseconverttocrs(array1[i]);
            }

            #region Seed
            // Makes an Array of List where the first index is Floor(j) and the second index are all elements in the seed vector to be non zero.
            // The Lanczos routine is parallelized by j, and the routine takes the list at [j] to implement the seed vector. Empty list will indicate that no seed is being used.
            var SeedPositionsByJ = new List<int>[array1.Length];
            var SeedCoefficientsByJ = new List<double>[array1.Length];
            for (int i = 0; i < array1.Length; i++)
            {
                SeedPositionsByJ[i] = new List<int>(); // Initialize seed position list for each j
                SeedCoefficientsByJ[i] = new List<double>();
            }
            if (input.useSeed)
            {
                SeedPositionsByJ = GenerateSeedPositions(input.SeedFile, input.nModes, isQuad, ref SeedCoefficientsByJ);
            }
            #endregion

            #region Lanczos
            int[] IECODE = new int[array1.Length];
            int[] ITER = new int[array1.Length];
            //actually diagonalizes the Hamiltonian matrix
            measurer.Reset();
            measurer.Start();
            //if the evecs of the lanczos matrices will need to be stored then save the lanczos matrices.
            if (input.PrintVector && !input.BlockLanczos && array1[0].innerobj.m >= Lanczos.basisSetLimit)
            {
                lanczosEVectors = new List<double[,]>();
                for (int i = 0; i < array1.Length; i++)
                {
                    lanczosEVectors.Add(new double[0, 0]);
                }
            }
            ParallelOptions options2 = new ParallelOptions();
            options2.MaxDegreeOfParallelism = input.ParJ;
            try
            {
                Parallel.For(0, array1.Length, options2, i =>//changed to array1.count from sHamMatrix.count
                {
                    double[] evs;
                    double[,] temp;
                    IECODE[i] = -1;

                    //add a parameter to count Lanczos iterations to set possible stopping criteria that way
                    //call MINVAL from here
                    if (!input.BlockLanczos)//means use naiveLanczos routine
                    {
                        ITER[i] = input.NumberOfIts;
                        evs = new double[input.M + 1];
                        temp = new double[numcolumnsA[i], input.M + 1];
                        Lanczos.NaiveLanczos(ref evs, ref temp, array1[i], input.NumberOfIts, input.Tolerance, input.PrintVector, SeedPositionsByJ[i], SeedCoefficientsByJ[i], i, input.FilePath);
                    }
                    else//means use block Lanczos from SOCJT
                    {
                        evs = new double[input.M + 1];
                        temp = new double[numcolumnsA[i], input.M + 1];//changed here to numcolumnsA
                        IECODE[i] = -1;
                        ITER[i] = Lanczos.MINVAL(numcolumnsA[i], input.M + 1, input.kFactor, input.M, input.NumberOfIts, input.Tolerance, 0, ref evs, ref temp, ref IECODE[i], array1[i], input.ParVectorMultiplication);
                    }

                    //initialize eigenvalues to have a length.                    
                    eigenvalues[i] = new double[evs.Length - 1];
                    for (int j = 0; j < evs.Length - 1; j++)
                    {
                        eigenvalues[i][j] = evs[j];
                    }
                    //I think this should be only for if block lanczos or naive lanczos with already calculated eigenvectors
                    if (input.BlockLanczos) // || (!input.BlockLanczos && array1[i].innerobj.m * input.NumberOfIts < Lanczos.basisSetLimit)) // I don't know what this second condition is for.
                    {
                        zMatrices[i] = new double[numcolumnsA[i], evs.Length - 1];//changed input.M to evs.Length - 1
                        for (int j = 0; j < numcolumnsA[i]; j++)
                        {
                            for (int k = 0; k < evs.Length - 1; k++)//changed input.M to evs.Length - 1                            
                            {
                                zMatrices[i][j, k] = temp[j, k];
                            }
                        }
                    }
                    //here if evectors are needed and hamiltonian is too large assign the lanczosEVectors to the 
                    if (!input.BlockLanczos && array1[0].innerobj.m >= Lanczos.basisSetLimit && input.PrintVector)
                    {
                        //assign the evecs of the lanczos matrices to the lanczosEVectors list.
                        lanczosEVectors[i] = new double[temp.GetLength(0), temp.GetLength(1)];
                        for (int j = 0; j < temp.GetLength(0); j++)
                        {
                            for (int k = 0; k < temp.GetLength(1); k++)
                            {
                                lanczosEVectors[i][j, k] = temp[j, k];
                            }
                        }
                    }
                    temp = null;
                    evs = null;
                }//end for
                );
            }//end try
            catch (AggregateException ae)
            {
                foreach (var e in ae.InnerExceptions)
                {
                    if (e is RepeaterError)
                    {
                        throw new RepeaterError();
                    }
                    else
                    {
                        throw;
                    }
                }

            }//end catch

            if (input.CheckEigenvector)
            {
                //LanczosChecker(input, array1, zMatrices, eigenvalues);
                int jBlock = (int)(input.JBlockEigenvector.Item1 - 0.5M);
                CheckJohnEigenvector(input, array1[jBlock], eigenvalues[jBlock][input.JBlockEigenvector.Item2 - 1]);
            }

            if (!input.BlockLanczos && array1[0].innerobj.m * input.NumberOfIts >= Lanczos.basisSetLimit && input.PrintVector)
            {
                //make it so that the output file generator does not try to print the values in the zmatrices which will be the eigenvectors of the lanczos matrix, not the hamiltonian
                input.PrintVector = false;
            }

            measurer.Stop();
            howMuchTime = measurer.ElapsedMilliseconds;
            input.DiagonalizationTime = (double)howMuchTime / 1000D;
            #endregion

            if (input.IncludeSO == false)
            {
                if (isQuad == false)
                {
                    JvecsForOutuput = jBasisVecsByJ;
                }//end if
                else
                {
                    for (int i = 0; i < jbasisoutA.Length; i++)
                    {
                        JvecsForOutuput.Add(jbasisoutA[i]);
                    }
                }//end else
            }

            //writes the eigenvectors to disk if that has been requested
            //if NTooBig == true then the separate eigenvector file will be made after the eigenvectors are calculated
            if (input.EVectorFile && (input.BlockLanczos || input.PrintVector))
            {
                writeVecFile(input, zMatrices, JvecsForOutuput, 0.0);
            }//end if

            //put code here to write eigenvectors to file with entire basis set.
            //do this by writing function to use jbasisvecsbyj which has all basis functions and
            //go through that list, where a basis function is in it and JvecsForOutput, pull the appropriate coefficient
            //from the zmatrices, otherwise just put in a 0.
            if (input.VectorFileComplete)
            {
                writeVecComplete(jBasisVecsByJ, JvecsForOutuput, zMatrices, input, eigenvalues, isQuad);
            }

            //this is where the intensity or overlap will be checked if necessary
            var overlaps = new List<double[]>(eigenvalues.Count);
            for (int count = 0; count < eigenvalues.Count; count++)
            {
                overlaps.Add(new double[eigenvalues[count].Length]);
            }

            if (input.Intensity && input.PrintVector)
            {
                //code here to read vector and take dot product
                double[] vector;
                if (input.JSInten)
                {
                    for (int jIndex = 0; jIndex < eigenvalues.Count(); jIndex++)
                    {
                        vector = EigenvectorReader(input.FilePath + input.VectorName, jIndex);
                        Lanczos.normalize(vector);
                        overlaps[jIndex] = Overlap(zMatrices[jIndex], vector);
                    }
                }
                else
                {
                    vector = ReadSOCJT2Vector(input.VectorName, input.VectorIndex, input.VectorJBlock, input.Special);
                    Lanczos.normalize(vector);
                    for (int jIndex = 0; jIndex < eigenvalues.Count(); jIndex++)
                    {
                        if (jIndex == (int)(input.VectorJBlock - 0.5M))
                        {
                            overlaps[jIndex] = Overlap(zMatrices[jIndex], vector);
                        }
                        else
                        {
                            overlaps[jIndex] = new double[zMatrices[jIndex].GetLength(0)];
                        }
                    }
                }
                foreach (double[] overlap in overlaps)
                {
                    Lanczos.normalize(overlap);
                }
            }

            List<string> linesToWrite = new List<string>();
            finalList = setAndSortEVs(eigenvalues, input.S, input.IncludeSO, zMatrices, JvecsForOutuput, input, overlaps);//add the eigenvectors so that the symmetry can be included as well
            linesToWrite = OutputFile.makeOutput(input, zMatrices, array1, JvecsForOutuput, eigenvalues, isQuad, finalList, IECODE, ITER);
            outp = linesToWrite;
            return linesToWrite;
        }//end SOCJT Routine

        /// <summary>
        /// Parses a SOCJT 2 vector file and returns the specified vector
        /// </summary>
        /// <param name="fileName">The name of vector file</param>
        /// <param name="index">Which eigenvector</param>
        /// <param name="jBlock">Which j-block</param>
        /// <param name="special">A parameter for a one off calculation. Leave at default.</param>
        /// <returns>
        /// The eigenvector read from file in the right order.
        /// </returns>
        private static double[] ReadSOCJT2Vector(string fileName, int index, decimal jBlock, bool special = false)
        {
            var evec = new List<string>();
            StreamReader vec = new StreamReader(fileName);
            string line;

            while ((line = vec.ReadLine()) != ("j-block " + Convert.ToString(jBlock)))
            {
                continue;
            }
            while ((line = vec.ReadLine()) != ("Eigenvector: " + Convert.ToString(index)))
            {
                continue;
            }
            while ((line = vec.ReadLine()) != ("Eigenvector: " + Convert.ToString(index + 1)) && vec.Peek() != -1)
            {
                evec.Add(line);
            }
            string[] parsedLine;
            int jPos = (int)(jBlock - 0.5M);
            var vecToReturn = new double[GenHamMat.basisPositions[jPos].Count()];
            char[] delimiters = new char[] { '\t', '\r', '=', ' ' };
            string hash;
            int[] vlLambda;
            int pos;
            double dummy;
            for (int i = 0; i < evec.Count(); i++)
            {
                parsedLine = evec[i].Split(delimiters, StringSplitOptions.RemoveEmptyEntries);
                //have this to ignore entries where there is only a space or nothing
                if (parsedLine.Length == 0 || !double.TryParse(parsedLine[0], out dummy))
                {
                    continue;
                }
                int nmodes = parsedLine.Length / 2 - 1;
                vlLambda = new int[nmodes * 2 + 1];
                for (int j = 1; j < parsedLine.Length; j++)
                {
                    vlLambda[j - 1] = Convert.ToInt32(parsedLine[j]);
                }
                hash = BasisFunction.GenerateHashCode(vlLambda, nmodes, false);
                if (GenHamMat.basisPositions[jPos].TryGetValue(hash, out pos))
                {
                    vecToReturn[pos] = FileInfo.ParseDouble(parsedLine[0]);
                }
            }
            //This is a special function used only in calculating intensities using special wavefunctions from the B state.
            if (special)
            {
                var zeroer = new double[vecToReturn.Length];
                for (int Lambda = -1; Lambda < 2; Lambda += 2)
                {
                    for (int v3v = 0; v3v < 2; v3v++)
                    {
                        for (int v4v = 0; v4v < 2; v4v++)
                        {
                            for (int v3l = -1; v3l < 2; v3l += 2)
                            {
                                for (int v4l = -1; v4l < 2; v4l += 2)
                                {
                                    if (v4v == 0 && v3v == 0)
                                    {
                                        continue;
                                    }
                                    if (v3v == 0)
                                    {
                                        v3l = 0;
                                    }
                                    if (v4v == 0)
                                    {
                                        v4l = 0;
                                    }
                                    vlLambda = new int[] { 0, 0, v3v, v3l, 1, v4l, Lambda };
                                    hash = BasisFunction.GenerateHashCode(vlLambda, 3, false);
                                    if (GenHamMat.basisPositions[jPos].TryGetValue(hash, out pos))
                                    {
                                        zeroer[pos] = 1.0;
                                    }
                                }//end v4l
                            }//end v3l
                        }//end v4v
                    }//end v3v
                }//end lambda
                for (int i = 0; i < zeroer.Length; i++)
                {
                    vecToReturn[i] = vecToReturn[i] * zeroer[i];
                }
            }
            return vecToReturn;
        }

        /// <summary>
        /// Calculates the overlap integral (dot product) of a provided vector and the eigenvectors
        /// </summary>
        /// <param name="eigenvectors">
        /// Eigenvectors of the Hamiltonian
        /// </param>
        /// <param name="overlapVector">
        /// Vector from John or from file reading.
        /// </param>
        /// <returns>
        /// Array of doubles which correspond to the dot products of the overlapVector with each of the eigenvectors.
        /// </returns>
        private static double[] Overlap(double[,] eigenvectors, double[] overlapVector)
        {
            //vector to store the overlaps
            var overlaps = new double[eigenvectors.GetLength(1)];
            for (int i = 0; i < overlaps.Length; i++)
            {
                for (int j = 0; j < overlapVector.Length; j++)
                {
                    overlaps[i] += overlapVector[j] * eigenvectors[j, i];
                }
                overlaps[i] = Math.Pow(overlaps[i], 2.0);
            }
            return overlaps;
        }

        /// <summary>
        /// Functiont to see if a vector from the Lanczos diagonalization is an eigenvector of the Hamiltonian.
        /// </summary>
        /// <param name="input">
        /// Fileinfo object
        /// </param>
        /// <param name="hamiltonianArray">
        /// Array witht the Hamiltonian matrices in them
        /// </param>
        /// <param name="zmatrices">
        /// Eigenvectors of the hamiltonians
        /// </param>
        /// <param name="eigenvalues">
        /// Eigenvalues of the Hamiltonians
        /// </param>
        private static void LanczosChecker(FileInfo input, alglib.sparsematrix[] hamiltonianArray, List<double[,]> zmatrices, List<double[]> eigenvalues)
        {
            int jBlock = (int)(input.JBlockEigenvector.Item1 - 0.5M);
            int whichEigenvalue = input.JBlockEigenvector.Item2 - 1;

            var temp = new double[zmatrices[jBlock].GetLength(0)];
            for (int row = 0; row < temp.Length; row++)
            {
                temp[row] = zmatrices[jBlock][row, whichEigenvalue];
            }
            EigenvectorCheck(hamiltonianArray[jBlock], eigenvalues[jBlock][whichEigenvalue], temp);
        }//end function EvalChecker

        private static void CheckJohnEigenvector(FileInfo input, alglib.sparsematrix hamiltonian, double eigenvalue)
        {
            var johnsVector = EigenvectorReader(input.FilePath + input.EigenvectorFileName, (int)(input.JBlockEigenvector.Item1 - 0.5M));
            EigenvectorCheck(hamiltonian, eigenvalue, johnsVector);
        }

        /// <summary>
        /// Checks if a vector is an eigenvector
        /// </summary>
        /// <param name="hamiltonianArray">
        /// Hamiltonian to be checked
        /// </param>
        /// <param name="eigenvalue">
        /// Eigenvalue to be compared against
        /// </param>
        /// <param name="temp">
        /// Vector being tested as an eigenvector of hamiltonianArray
        /// </param>
        private static void EigenvectorCheck(alglib.sparsematrix hamiltonianArray, double eigenvalue, double[] temp)
        {
            var V = new double[temp.Length];
            Lanczos.normalize(temp);
            alglib.sparsemv(hamiltonianArray, temp, ref V);
            double sum = 0.0;
            double magnitude = Lanczos.Magnitude(V);
            for (int row = 0; row < temp.Length; row++)
            {
                sum += Math.Pow(V[row] / magnitude - temp[row], 2.0);
            }
            sum /= V.Length;
            if (eigenvalue < 0.0)
            {
                magnitude *= -1.0;
            }
            Console.WriteLine("Calculated Eigenvalue = " + Convert.ToString(magnitude));
            Console.WriteLine("Eigenvalue from Lanczos = " + Convert.ToString(eigenvalue));
            double ratio = eigenvalue / magnitude;
            Console.WriteLine("Ratio of Lanczos/Calculated = " + Convert.ToString(ratio));

            //then find the residuals of A*x - c*x where c is the magnitude from above and x is temp
            //means keep a running sum of V[i] - magnitude*temp[i]
            sum = 0.0;
            for (int row = 0; row < temp.Length; row++)
            {
                sum += Math.Pow(V[row] - magnitude * temp[row], 2.0);
            }
            sum /= V.Length;
            sum = Math.Sqrt(sum);
            Console.WriteLine("RMS error between vectors = " + sum);
            Console.ReadLine();
        }//end EigenvectorCheck

        /// <summary>
        /// Reads in and parses an eigenvector from John
        /// </summary>
        /// <param name="filepath">
        /// What the filename is of the eigenvector including filepath.
        /// </param>
        /// <param name="jIndex">
        /// Which j block this vector corresponds to.
        /// </param>
        /// <returns>
        /// The eigevnector with the coefficients from John's file in the appropriate spot.
        /// </returns>
        private static double[] EigenvectorReader(string filepath, int jIndex)
        {
            double[] eigenvector = new double[GenHamMat.basisPositions[jIndex].Count()];
            string[] vecFile = FileInfo.FileRead(filepath);
            double temp = 0.0;
            for (int position = 0; position < vecFile.Length; position += 8)
            {
                int[] vlLambda = new int[7];
                vlLambda[0] = Convert.ToInt32(vecFile[position]);//v for v1
                vlLambda[3] = 0;//l for v1
                vlLambda[1] = Convert.ToInt32(vecFile[position + 4]);//v for v3
                vlLambda[4] = Convert.ToInt32(vecFile[position + 2]);//l for v3
                vlLambda[2] = Convert.ToInt32(vecFile[position + 5]);//v for v4
                vlLambda[5] = Convert.ToInt32(vecFile[position + 3]);//l for v4
                if (vecFile[position + 7][vecFile[position + 7].Length - 1] == '+')//for Lambda
                {
                    vlLambda[6] = 1;
                }
                else
                {
                    vlLambda[6] = -1;
                }
                string hashCode = BasisFunction.GenerateHashCode(vlLambda, 3);
                int index = 0;
                if (GenHamMat.basisPositions[jIndex].TryGetValue(hashCode, out index))
                {
                    temp = Convert.ToDouble(vecFile[position + 6]);
                    if (temp != 0.0)
                    {
                        eigenvector[index] = Convert.ToDouble(vecFile[position + 6]);
                    }
                    else
                    {
                        string s = vecFile[position + 7].Substring(0, vecFile[position + 7].Length - 1);
                        eigenvector[index] = Convert.ToDouble(s);
                    }
                }
            }
            return eigenvector;
        }

        /// <summary>
        /// Writes the eigenvectors to a separate file including basis functions with a 0.0 coefficient.
        /// </summary>
        /// <param name="input">
        /// FileInfo object.
        /// </param>
        /// <param name="zMatrices">
        /// Eigenvectors.
        /// </param>
        /// <param name="JvecsForOutuput">
        /// Basis Set.
        /// </param>
        public static void writeVecFile(FileInfo input, List<double[,]> zMatrices, List<List<BasisFunction>> JvecsForOutuput, double evMin, string file = "")
        {
            StringBuilder vecFile = new StringBuilder();
            vecFile.AppendLine("VecFile " + input.Title);
            vecFile.AppendLine(" ");
            decimal S = input.S;
            for (int m = 0; m < zMatrices.Count; m++)
            {
                vecFile.AppendLine("*************************************");
                vecFile.AppendLine(" ");
                vecFile.AppendLine("j-block " + ((decimal)m + 0.5M));
                vecFile.AppendLine(" ");
                vecFile.AppendLine("*************************************");
                vecFile.AppendLine(" ");
                for (int o = 0; o < zMatrices[m].GetLength(1); o++)
                {
                    vecFile.AppendLine("Eigenvector: " + (o + 1));
                    vecFile.AppendLine(" ");
                    OutputFile.vecBuilder(input, JvecsForOutuput[m], vecFile, zMatrices[m], o, evMin, S);
                    vecFile.AppendLine(" ");
                }
                vecFile.AppendLine(" ");
                if (input.IncludeSO)
                {
                    S++;
                    if (S > -1.0M * input.S)
                    {
                        S = input.S;
                    }
                }
            }
            List<string> vecFileOut = new List<string>();
            vecFileOut.Add(vecFile.ToString());
            string title;
            if (file == "")
            {
                title = input.FilePath + input.Title + "_vec.out";
            }
            else
            {
                title = file;
            }
            File.WriteAllLines(title, vecFileOut);
        }

        /// <summary>
        /// Writes the eigenvectors to file using the complete basis set for all eigenvectors.
        /// This means that (in a quadratic calculation) the eigenvectors for an e level
        /// would include the basis functions with a symmetry and those with neither symmetry.
        /// This is useful for comparison to calculations where the complete basis set is used.
        /// </summary>
        /// <param name="jBasisVecsByJ">
        /// All basis functions sorted by J value.
        /// </param>
        /// <param name="JvecsForOutput">
        /// Only symmetry allowed basis functions sorted by J (linear) or symmetry (quadratic)
        /// </param>
        /// <param name="tempMat">
        /// The coefficients of the eigenvectors using the JvecsForOutput basis.
        /// </param>
        /// <param name="input">
        /// FileInfo object.
        /// </param>
        /// <param name="eigenvalues">
        /// Calculated eigenvalues.
        /// </param>
        /// <param name="isQuad">
        /// True if this is a quadrtic calculation, false if not.
        /// </param>
        private static void writeVecComplete(List<List<BasisFunction>> jBasisVecsByJ, List<List<BasisFunction>> JvecsForOutput, List<double[,]> tempMat, FileInfo input, List<double[]> eigenvalues, bool isQuad)
        {
            double ZPE = eigenvalues[0][0];
            if (ZPE < 0.0)
            {
                ZPE *= -1.0;
            }
            var tempEvalue = new List<double[]>(eigenvalues.Count);

            for (int i = 0; i < eigenvalues.Count; i++)
            {
                tempEvalue.Add(new double[eigenvalues[i].Length]);
                //tempEvalue[i] = new double[eigenvalues[i].Length];
                for (int j = 0; j < eigenvalues[i].Length; j++)
                {
                    tempEvalue[i][j] = eigenvalues[i][j] + ZPE;
                }
            }
            StringBuilder file = new StringBuilder();
            file.AppendLine(" ");
            file.AppendLine("Eigenvectors in complete basis set for " + input.Title);
            file.AppendLine(" ");
            //loop over jBlocks
            for (int jBlockIndex = 0; jBlockIndex < eigenvalues.Count; jBlockIndex++)
            {
                file.AppendLine(" ");
                file.AppendLine("J-Block " + ((decimal)jBlockIndex + 0.5M));
                file.AppendLine(" ");
                //loop over all of the eigenvalues found
                for (int eigenvectorIndex = 0; eigenvectorIndex < eigenvalues[jBlockIndex].Length; eigenvectorIndex++)
                {
                    file.AppendLine(" " + "\r");
                    file.AppendLine("Eigenvalue" + "\t" + Convert.ToString(eigenvectorIndex + 1) + " = " + String.Format("{0,10:0.0000}", tempEvalue[jBlockIndex][eigenvectorIndex]));
                    bool a1 = SOCJT.isA(JvecsForOutput[jBlockIndex], tempMat[jBlockIndex], eigenvectorIndex, input, false);
                    if (a1)
                    {
                        file.AppendLine("Vector is Type 1");
                    }
                    else
                    {
                        file.AppendLine("Vector is Type 2");
                    }
                    file.AppendLine(" " + "\r");
                    file.Append(" Coefficient  " + "\t");
                    for (int h = 0; h < input.nModes; h++)
                    {
                        file.Append("v(" + Convert.ToString(h + 1) + ")" + "\t" + "l(" + Convert.ToString(h + 1) + ")" + "\t");
                    }
                    file.Append("lambda");
                    decimal S = input.S;
                    for (int j = 0; j < jBasisVecsByJ.Count; j++)//goes through basis vectors
                    {
                        //first split into quadratic and linear
                        if (isQuad)
                        {
                            for (int k = 0; k < jBasisVecsByJ[j].Count; k++)
                            {
                                string hashCode = BasisFunction.GenerateHashCode(jBasisVecsByJ[j][k]);
                                int m = 0;
                                if (GenHamMat.basisPositions[jBlockIndex].TryGetValue(hashCode, out m))
                                {
                                    writeVec(tempMat[jBlockIndex][m, eigenvectorIndex], jBasisVecsByJ[j][k], file, S);
                                }
                                else
                                {
                                    writeVec(0.0, jBasisVecsByJ[j][k], file, S);
                                }
                            }//end for loop over J
                        }//end if isQuad
                        else
                        {
                            if (jBasisVecsByJ[j][0].J == (decimal)jBlockIndex + 0.5M)
                            {
                                //write actual values to the eigenvector
                                for (int k = 0; k < jBasisVecsByJ[j].Count; k++)
                                {
                                    writeVec(tempMat[jBlockIndex][k, eigenvectorIndex], jBasisVecsByJ[j][k], file, S);
                                }//end loop over BasisFunctions in j value for jBasisVecsByJ
                            }
                            else
                            {
                                //just write a zero for all of the coefficients for this j value
                                for (int k = 0; k < jBasisVecsByJ[j].Count; k++)
                                {
                                    writeVec(0.0, jBasisVecsByJ[j][k], file, S);
                                }//end loop over BasisFunctions in j value for jBasisVecsByJ
                            }
                        }//end linear option
                        if (input.IncludeSO)
                        {
                            S++;
                            if (S > -1.0M * input.S)
                            {
                                S = input.S;
                            }
                        }
                    }//end loop over jValues in jbasisVecsByJ
                    file.AppendLine("\r");
                }//end loop over the number of eigenvalues
            }//end loop over j blocks

            List<string> vecFileOut = new List<string>();
            vecFileOut.Add(file.ToString());
            File.WriteAllLines((input.FilePath + input.Title + "_CompleteVector.out"), vecFileOut);
        }//end function writeCompleteVec

        /// <summary>
        /// Addes a basis function to an eigenvector including its coefficient
        /// </summary>
        /// <param name="coefficient">
        /// Coefficient of the basis function
        /// </param>
        /// <param name="func">
        /// Basis function being added
        /// </param>
        /// <param name="file">
        /// StringBuilder to add the function to.
        /// </param>
        public static void writeVec(double coefficient, BasisFunction func, StringBuilder file, decimal S)
        {
            file.AppendLine("\t");
            file.Append(String.Format("{0,14:0.0000000000}", coefficient));
            for (int m = 0; m < func.modesInVec.Count; m++)//goes through each mode
            {
                file.Append("\t" + "  " + Convert.ToString(func.modesInVec[m].V) + "\t" + String.Format("{0,3}", func.modesInVec[m].L));//  "  " + Convert.ToString(jBasisVecsByJ[i][h].modesInVec[m].l));
            }
            file.Append("\t" + String.Format("{0,4}", func.Lambda) + "\t" + String.Format("{0,5}", S));
        }//end of function writeVec

        /// <summary>
        /// Function to determine if a given basis function is found in an eigenvector.
        /// </summary>
        /// <param name="jBasisVec">
        /// The basis function to be looked for
        /// </param>
        /// <param name="eigenvector">
        /// The eigenvector to check for jBasisVec
        /// </param>
        /// <param name="place">
        /// This will be the index of the jBasisVec if it is found.
        /// </param>
        /// <param name="start">
        /// The place to start in eigenvector. Used to avoid double checking basis functions which have already been found.
        /// </param>
        /// <returns>
        /// Boolean to indicate whether or not the basisfunction is in the eigenvector.
        /// </returns>
        private static bool isInBasis(BasisFunction jBasisVec, List<BasisFunction> eigenvector, ref int place, int start)
        {
            bool In = false;
            for (int i = start; i < eigenvector.Count; i++)
            {
                if (jBasisVec.J == eigenvector[i].J && jBasisVec.Lambda == eigenvector[i].Lambda)
                {
                    for (int j = 0; j < jBasisVec.modesInVec.Count; j++)
                    {

                        if (jBasisVec.modesInVec[j].V == eigenvector[i].modesInVec[j].V && jBasisVec.modesInVec[j].L == eigenvector[i].modesInVec[j].L)
                        {
                            if (j == jBasisVec.modesInVec.Count - 1)
                            {
                                In = true;
                                place = i;
                            }
                            continue;
                        }//end if for v and l
                        else
                        {
                            break;
                        }
                    }//end for loop
                }//end if for J and Lambda
            }//end loop over all eigenvectors
            return In;
        }//end isInBasis

        /// <summary>
        /// Function to generate array of Eigenvalues for final output file.
        /// </summary>
        /// <param name="evs">
        /// Eigenvalues for all j blocks
        /// </param>
        /// <param name="S">
        /// Spin
        /// </param>
        /// <param name="inclSO">
        /// True if there is SO coupling, false if not.
        /// </param>
        /// <param name="zMatrices">
        /// Eigenvectors for each j-block
        /// </param>
        /// <param name="jvecs">
        /// Basis set for each j-block
        /// </param>
        /// <param name="input">
        /// FileInfo object
        /// </param>
        /// <returns>
        /// Eigenvalue array with eigenvalue objects all initialized and sorted by value.
        /// </returns>
        public static Eigenvalue[] setAndSortEVs(List<double[]> evs, decimal S, bool inclSO, List<double[,]> zMatrices, List<List<BasisFunction>> jvecs, FileInfo input, List<double[]> overlap)
        {
            List<Eigenvalue> eigen = new List<Eigenvalue>();
            int counter = 0;
            decimal J = 0.5M;
            S = S * -1M;
            decimal tempS = S;
            decimal maxS = S;
            if (inclSO == true)
            {
                maxS = maxS * -1M;
            }
            for (int i = 0; i < evs.Count; i++)
            {
                for (int j = 0; j < evs[i].Length; j++)
                {
                    //add call to symmetry checker function here.
                    bool tbool = isA(jvecs[i], zMatrices[i], j, input, false);
                    if (input.Intensity && input.PrintVector)
                    {
                        eigen.Add(new Eigenvalue(J, j + 1, tempS, evs[i][j], tbool, overlap[i][j]));
                    }
                    else
                    {
                        eigen.Add(new Eigenvalue(J, j + 1, tempS, evs[i][j], tbool));
                    }
                }
                if (tempS < maxS && (J - 1.5M) % 3M != 0M)
                {
                    tempS++;
                }
                else
                {
                    tempS = S;
                    J++;
                }
                counter += evs[i].Length;
            }
            Eigenvalue[] eigenarray = eigen.ToArray();
            bubbleSort(ref eigenarray);
            if (input.useAbsoluteEV == false)
            {
                double ZPE = eigenarray[0].Evalue;
                int[] temp = new int[evs.Count];
                //for (int i = 0; i < evs.Count; i++)
                //{
                //    temp[i] = 1;
                //}
                for (int i = 0; i < eigenarray.Length; i++)
                {
                    eigenarray[i].Evalue = eigenarray[i].Evalue - ZPE;
                }
            }//end ZPE loop
            int SOnumb = (int)(-2M * S) + 1;
            if (inclSO == false)
            {
                SOnumb = 1;
            }
            //for (int i = 0; i < eigenarray.Length; i++)
            //{
            //    int Snumb = (int)(eigenarray[i].Sigma - S);
            //    int j = (int)(eigenarray[i].JBlock - 0.5M);
            //    int place = j * SOnumb + Snumb;
            //    eigenarray[i].Number = temp[place];
            //    temp[place]++;
            //}
            return eigenarray;
        }

        /// <summary>
        /// Method to determine if a vector is type 1 (symmetric) or type 2 (antisymmetric).
        /// </summary>
        /// <param name="jBasisVecsByJ">
        /// Basis set of vector being looked at
        /// </param>
        /// <param name="tempMat">
        /// Matrix containing the eigenvectors.
        /// </param>
        /// <param name="j">
        /// Which eigenvector to look at.
        /// </param>
        /// <param name="input">
        /// FileInfo object
        /// </param>
        /// <param name="overRide">
        /// Boolean to simply return type 2 for all eigenvectors because this function is being called with the lanczos eigenvectors, not the true eigenvectors.
        /// </param>
        /// <returns>
        /// True if symmetric, false if antisymmetric.
        /// </returns>
        public static bool isA(List<BasisFunction> jBasisVecsByJ, double[,] tempMat, int j, FileInfo input, bool overRide)
        {
            bool a1 = false;
            double temp = 0.0;
            int[] tempVL = new int[input.nModes * 2 + 1];
            //this conditional is because if naive lanczos is used and eigenvectors are not calculated the tempMat actually contains the eigenvectors of the the lanczos matrix, not the Hamiltonian
            //boolean overRide lets me skip this when it's called to do separate eigenvector calculations
            if (!input.BlockLanczos && !input.PrintVector && !overRide)
            {
                return false;
            }
            for (int m = 0; m < jBasisVecsByJ.Count; m++)
            {
                if (Math.Abs(tempMat[m, j]) > temp)
                {
                    for (int n = 0; n < input.nModes; n++)
                    {
                        tempVL[n * 2] = jBasisVecsByJ[m].modesInVec[n].V;
                        tempVL[n * 2 + 1] = jBasisVecsByJ[m].modesInVec[n].L;
                    }
                    tempVL[input.nModes * 2] = jBasisVecsByJ[m].Lambda;
                    temp = tempMat[m, j];
                }
            }

            for (int m = 0; m < jBasisVecsByJ.Count; m++)
            {
                if (jBasisVecsByJ[m].Lambda == -1 * tempVL[input.nModes * 2])
                {
                    int tempInt = 0;
                    for (int v = 0; v < input.nModes; v++)
                    {
                        if (jBasisVecsByJ[m].modesInVec[v].V == tempVL[v * 2])
                        {
                            if (jBasisVecsByJ[m].modesInVec[v].L == -1 * tempVL[v * 2 + 1])
                            {
                                tempInt++;
                            }
                        }
                        if (tempInt == input.nModes)
                        {
                            if (temp / tempMat[m, j] > 0)
                            {
                                a1 = true;
                            }
                            m = jBasisVecsByJ.Count;
                        }
                    }
                }
            }
            return a1;
        }

        /// <summary>
        /// Method to sort eigenvalues in increasing order
        /// </summary>
        /// <param name="arr">
        /// Eigenvalues to be sorted
        /// </param>
        private static void bubbleSort(ref Eigenvalue[] arr)
        {
            bool swapped = true;
            int j = 0;
            Eigenvalue tmp;
            while (swapped == true)
            {
                swapped = false;
                j++;
                for (int i = 0; i < arr.Length - j; i++)
                {
                    if (arr[i].Evalue > arr[i + 1].Evalue)
                    {
                        tmp = arr[i];
                        arr[i] = arr[i + 1];
                        arr[i + 1] = tmp;
                        swapped = true;
                    }//end if
                }//end for
            }//end while           
        }//end method bubbleSort

        /// <summary>
        /// Used to multiply all elements in a sparse matrix A by the value val
        /// </summary>
        /// <param name="A">
        /// Sparse matrix to be multiplied.
        /// </param>
        /// <param name="val">
        /// Value to be multiplied
        /// </param>
        /// <returns>
        /// Sparesmatrix containing the original matrix A times the double val.
        /// </returns>
        private static alglib.sparsematrix cTimesSparse(alglib.sparsematrix A, double val)
        {
            int i;
            int j;
            double oldVal;
            int t0 = 0;
            int t1 = 0;
            alglib.sparsematrix B = new alglib.sparsematrix();
            alglib.sparsecreate(A.innerobj.m, A.innerobj.n, out B);
            while (alglib.sparseenumerate(A, ref t0, ref t1, out i, out j, out oldVal))
            {
                alglib.sparseadd(B, i, j, oldVal * val);
            }
            return B;
        }//end method cTimesSparse

        /// <summary>
        /// Takes a list of alglib.sparsematrix objects and combines them.
        /// </summary>
        /// <param name="mat">
        /// List of alglib.sparsematrix objects to be combined
        /// </param>
        /// <returns>
        /// alglib.sparsematrix object which contains all elements of all alglib.sparsematrix objects in the List mat.
        /// </returns>
        private static alglib.sparsematrix aggregator(List<alglib.sparsematrix> mat)
        {
            int i;
            int j;
            double oldVal;
            alglib.sparsematrix B = new alglib.sparsematrix();
            alglib.sparsecreate(mat[0].innerobj.m, mat[0].innerobj.n, out B);
            for (int m = 0; m < mat.Count; m++)
            {
                int t0 = 0;
                int t1 = 0;
                while (alglib.sparseenumerate(mat[m], ref t0, ref t1, out i, out j, out oldVal))
                {
                    alglib.sparseadd(B, i, j, oldVal);
                    //adds lower diagonal matrix elements for all off-diagonal matrices (first element in mat is the diagonal elements so don't do this for it). 
                    if (m != 0)
                    {
                        alglib.sparseadd(B, j, i, oldVal);
                    }
                }
            }
            return B;
        }

        /// <summary>
        /// Function that reads a matrix file from file and initializes hamiltonian matrices.
        /// </summary>
        /// <param name="input">
        /// initialized FileInfo object.
        /// </param>
        /// <param name="matricesMade">
        /// Set to true if the matrices are generated after this function call.
        /// </param>
        /// <returns>
        /// Boolean indicating whether or not the fitHamList matrices have been initialized from matFile.
        /// </returns>
        private static List<int> matReadFunction(FileInfo input, ref bool matricesMade)
        {
            matricesMade = false;
            List<int> basisSizeList = new List<int>();
            if (input.UseMatrixFile && input.MatrixMade)
            {
                string[] matFile = { };
                try
                {
                    matFile = FileInfo.FileRead(input.MatrixFilePath);
                }
                catch (FileNotFoundException)
                {
                    throw new FileNotFoundException("The matrix file does not exist.");
                }
                catch
                {
                    throw new Exception("Matrix File Error." + "\r" + "Please check the matrix file and try again.");
                }
                matricesMade = true;
                //read the matrix into memory
                int index = 0, basisSize = 0;
                for (int i = 0; i < matFile.Length; i++)
                {
                    if (matFile[i] == "List")
                    {
                        index = Convert.ToInt32(matFile[i + 1]);
                        basisSize = Convert.ToInt32(matFile[i + 2]);
                        basisSizeList.Add(basisSize);
                        //this is the matrix for the diagonal elements.  These are added seperately.
                        var B = new alglib.sparsematrix();
                        alglib.sparsecreate(basisSize, basisSize, out B);
                        fitHamList[index].Add(B);
                        continue;
                    }//end if matFile[i] == List
                    if (matFile[i] == "Matrix")
                    {
                        var B = new alglib.sparsematrix();
                        alglib.sparsecreate(basisSize, basisSize, out B);
                        i += 2;
                        while (matFile[i] != "List" && matFile[i] != "Matrix")
                        {
                            alglib.sparseadd(B, Convert.ToInt32(matFile[i]), Convert.ToInt32(matFile[i + 1]), FileInfo.ParseDouble(matFile[i + 2]));
                            i += 3;
                            if (i >= matFile.Length)
                            {
                                break;
                            }
                        }
                        i--;
                        fitHamList[index].Add(B);
                    }//end if matFile[i] == Matrix
                }//end for loop over entire matrix file
            }//end if
            //return matricesMade;
            return basisSizeList;
        }//end matReadFunction

        /// <summary>
        /// Function that takes the inputted seed vector information and creates an array of lists of integers and an array of lists of doubles. Each array element corresponds to a different j value. 
        /// Each list holds the positions of basis function to be nonzero in the seed vector, or the coefficient of that seed vector. They have to be partitioned by j because SOCJT handles each j block
        /// separately. The function specifically goes through each basis function and calculates the j value. If we are in the nonquadratic case, then the index for the jth block is just Floor(j). 
        /// If we are in the quadratic case, then the index is 0 for j = 1/2 +- 3n (these are degenerate with j = 5/2 +- 3n. so j = 5/2 +- 3n should not be included and the function skips it) and 1 for
        /// j = 3/2 +- 3n, e and a1/a2 respectively. The position and coefficient is added onto the list stored in the array at this point.
        /// </summary>
        /// <param name="SeedFile"></param>
        /// String for the name of the seed file.
        /// <param name="nModes"></param>
        /// Integer for the number of modes.
        /// <param name="isQuad"></param>
        /// Boolean for whether or not the Hamiltonian is in the quadratic or nonquadratic case.
        /// <param name="SeedCoefficientsByJ"></param>
        /// Array of lists that stores the coefficients for each basis function in the seed vector in order, sorted by j.
        /// <returns></returns>
        private static List<int>[] GenerateSeedPositions(string SeedFile, int nModes, bool isQuad, ref List<double>[] SeedCoefficientsByJ) // This generates a List of integers in the second dimension which holds all basis functions to be nonzero. The first dimension holds the Floor(j) value. - HT
        {
            Seed SeedVector = new Seed(SeedFile, nModes);
            var SeedPositionsByJ = new List<int>[GenHamMat.basisPositions.Count]; // An array of lists. The first index points to the j list and in the j list are all positions for that j.
            SeedCoefficientsByJ = new List<double>[GenHamMat.basisPositions.Count]; // This holds the coefficients of the basis function in the seed vector for each j.
            for (int i = 0; i < GenHamMat.basisPositions.Count; i++)
            {
                SeedPositionsByJ[i] = new List<int>(); // Initialize a new list for each j.
                SeedCoefficientsByJ[i] = new List<double>();
            }

            for (int i = 0; i < SeedVector.SeedIndex; i++) // Goes through each position on the seed vector and sorts them by j.
            {
                int position; // Stores position in the basis.
                string tmpHash = BasisFunction.GenerateHashCode(SeedVector.vlLambdaSeed[i], nModes, false); // hashcode to get the position
                decimal jBlock = 0; // Will store j
                int jIndex = 0; // Which index corresponds to that value of j
                for (int j = 0; j < nModes; j++)
                {
                    jBlock += SeedVector.vlLambdaSeed[i][1 + 2 * j]; // Sum of l
                }
                jBlock += (decimal)SeedVector.vlLambdaSeed[i][2 * nModes] / 2; // Total, jBlock = \sum_i l_i + Lambda / 2 (the j value)
                if (isQuad == false)
                {
                    jIndex = (int)(Math.Abs(jBlock) - 0.5M); // Floor(j) is the index for each j value in the nonquadratic case
                }
                if (isQuad == true)
                {
                    if ((int)(jBlock - 1.5M) % 3 == 0) // j = 3/2 + 3n means a1/a2 block
                    {
                        jIndex = 1; // Index for a1/a2 block
                    }
                    else
                    {
                        //if ((int)(jBlock - 2.5M) % 3 == 0) // Means j = 5/2 + 3n, which is degerate with 1/2 + 3n and not included.
                        //{
                        //    continue;
                        //}
                        jIndex = 0; // Index for e block
                    }
                }
                GenHamMat.basisPositions[jIndex].TryGetValue(tmpHash, out position); // Generates position of the basis function. Requires the j index and hashcode. Stores into position.
                if (position == 0)
                {
                    throw new System.ArgumentException("Error in seed vector. Recheck seed file. \nReminder that j = 5/2 + 3n should not be included."); // It is highly unlikely that position 0 will be used, and this is the default for improper vectors so I throw an error.
                }
                SeedPositionsByJ[jIndex].Add(position); // Adds position
                SeedCoefficientsByJ[jIndex].Add(SeedVector.SeedCoefficient[i]); // Add coefficient
            }
            return SeedPositionsByJ;
        }//end GenerateSeedPositions function
    }//end class SOCJT
}