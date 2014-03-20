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
                if (input.minJBool == true)
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
                basisByMode[i] = BasisByMode.genVLCombinations(Modes[i], i);
            }//end for            

            //Generates all of the JBasisVectors to be used in calculation.
            List<BasisFunction> hamiltonianVecs = BasisFunction.genJVecs(basisByMode, input.nModes, jMin, jMax);
            //hamiltonianVecs are the total list of all J vectors in the Hamiltonian

            //Sorts the hamiltonianVecs by J and puts them into a List of Lists of JBasisVectors.
            List<List<BasisFunction>> jBasisVecsByJ = new List<List<BasisFunction>>();            
            for (decimal i = jMin; i <= jMax; i++)
            {
                jBasisVecsByJ.Add(GenHamMat.SortByJ(hamiltonianVecs, i));
            }//end for loop

            //Initializes Lists to hold the Hamiltonian matrices, eigenvectors, eigenvalues, basis vectors for the output file and number of columns for each j matrix respectively.
            //List<double[,]> hamMatrices = new List<double[,]>(); ---moved into if statement so it's only made if it's needed.
            List<double[,]> zMatrices = new List<double[,]>();
            List<double[]> eigenvalues = new List<double[]>();
            List<List<BasisFunction>> JvecsForOutuput = new List<List<BasisFunction>>();
            List<BasisFunction>[] jbasisoutA = new List<BasisFunction>[0];
            List<int> numColumns = new List<int>();


            #region Hamiltonian
            List<alglib.sparsematrix> sHamMatrix = new List<alglib.sparsematrix>();            
            alglib.sparsematrix[] array1;            
            int[] numcolumnsA;
            //smallMat = false;            
            //Creates the Hamiltonian matrices for linear cases            
            int numQuadMatrix = 0;            
            List<int> a = new List<int>();

            if (input.M > input.noIts)
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

                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.parJ;
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
                                //fitHamList[i] = GenHamMat.genFitMatrix(jBasisVecsByJ[i], isQuad, input, out nColumns, input.parMat, false);
                                fitHamList[i] = GenHamMat.GenMatrixHash(jBasisVecsByJ[i], isQuad, input, out nColumns, input.parMat, false);                                
                            }
                            else//this makes sure that the diagonal portion is regenerated on each call.
                            {
                                //fitHamList[i][0] = GenHamMat.genFitMatrix(jBasisVecsByJ[i], isQuad, input, out nColumns, input.parMat, true)[0];
                                fitHamList[i][0] = GenHamMat.GenMatrixHash(jBasisVecsByJ[i], isQuad, input, out nColumns, input.parMat, true)[0];
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
                input.matGenTime = (double)howMuchTime / 1000D;
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

                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.parJ;
                try
                {
                    Parallel.For(jBasisVecsByJ.Count / 2, jBasisVecsByJ.Count - dynVar1, options, i =>//changed to dynVar1 from 6
                    {
                        List<BasisFunction> quadVecs = new List<BasisFunction>();
                        int nColumns;

                        //for (int v = tempDynVar2; v <= dynVar2; v++)
                        for (int v = -dynVar2; v <= dynVar2; v++)                        
                        {
                            quadVecs.AddRange(jBasisVecsByJ[i + v * 3]);
                        }
                        //this checks if the matrix was read from file, and if so if the basis set is the correct size
                        if (basisSize.Count != 0)                        
                        { 
                            if(basisSize[i - jBasisVecsByJ.Count / 2] != quadVecs.Count)
                            {
                                throw new MatrixFileError();
                            }
                        }
                        //if matrices aren't made then generate all of them
                        if (!matricesMade)
                        {
                            //fitHamList[i - jBasisVecsByJ.Count / 2] = GenHamMat.genFitMatrix(quadVecs, isQuad, input, out nColumns, input.parMat, false);
                            fitHamList[i - jBasisVecsByJ.Count / 2] = GenHamMat.GenMatrixHash(quadVecs, isQuad, input, out nColumns, input.parMat, false);       
                        }                        
                        else//If they are made then just generate the diagonal elements.
	                    {
                            //fitHamList[i - jBasisVecsByJ.Count / 2][0] = GenHamMat.genFitMatrix(quadVecs, isQuad, input, out nColumns, input.parMat, true)[0];
                            fitHamList[i - jBasisVecsByJ.Count / 2][0] = GenHamMat.GenMatrixHash(quadVecs, isQuad, input, out nColumns, input.parMat, true)[0];		 
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
                input.matGenTime = (double)howMuchTime / 1000D;
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
            GenHamMat.CrossTermInitialization(jBasisVecsByJ[0][0].modesInVec, input.nModes, out bilinear, out biAVecPos, out biEVecPos, input.crossTermMatrix);
            //code here to convert the alglib matrices to matrices for each j block
            for (int i = 0; i < fitHamList.Count; i++)
            {
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
                            val = Math.Sqrt(Modes[count].D) * Modes[count].modeOmega;
                        }
                        else
                        {
                            val = Modes[count].K * Modes[count].modeOmega;
                        }
                    }
                    else//means it's a cross term. loop over relevant E and A terms in same order as in genFitMatrix function
                    {
                        for (int aa = 0; aa < biAVecPos.Count; aa++)
                        {
                            for (int e = 0; e < biEVecPos.Count; e++)
                            {
                                int crossCount = aa + e;
                                int row;
                                int column;
                                if (biAVecPos[aa] > biEVecPos[e])
                                {
                                    column = biAVecPos[aa];
                                    row = biEVecPos[e];
                                }
                                else
                                {
                                    column = biEVecPos[e];
                                    row = biAVecPos[aa];
                                }
                                val = input.crossTermMatrix[row, column];
                            }//end loop over e elements
                        }//end loop over a elements
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
            if (input.inclSO == true)
            {
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
                        for (int k = 0; k < numcolumnsA[i]; k++)
                        {
                            double temp = input.Azeta * (double)jbasisoutA[i][k].Lambda * (double)j;
                            alglib.sparseadd(tempMatList[tempMatList.Count - 1], k, k, temp);
                        }//end loop over diagonal matrix elements
                        tempNumbColumns.Add(numcolumnsA[i]);
                        //means SO only in j = 0.5 block for quadratic cases
                        if (i > 0 && isQuad)
                        {
                            break;
                        }
                    }//end loop over values of S
                }//end loop over all previouly made sparseMatrices

                numcolumnsA = null;
                numcolumnsA = tempNumbColumns.ToArray();
                array1 = null;
                array1 = tempMatList.ToArray();
                zMatrices.Clear();
                eigenvalues.Clear();
                for (int i = 0; i < array1.Length; i++)
                {
                    zMatrices.Add(new double[0, 0]);
                    eigenvalues.Add(new double[0]);
                }
            }//end if inclSO == true
            #endregion

            else
            {
                for (int i = 0; i < array1.Length; i++)
                {
                    alglib.sparseconverttocrs(array1[i]);
                }
            }

            #region Lanczos
            int[] IECODE = new int[array1.Length];
            int[] ITER = new int[array1.Length];
            //actually diagonalizes the Hamiltonian matrix
            measurer.Reset();
            measurer.Start();
            //if the evecs of the lanczos matrices will need to be stored then save the lanczos matrices.
            if (input.pVector && !input.blockLanczos && array1[0].innerobj.m >= Lanczos.basisSetLimit)
            {
                lanczosEVectors = new List<double[,]>();
                for (int i = 0; i < array1.Length; i++)
                { 
                    lanczosEVectors.Add(new double[0,0]);
                }
            }
            ParallelOptions options2 = new ParallelOptions();
            options2.MaxDegreeOfParallelism = input.parJ;
            try
            {
                Parallel.For(0, array1.Length, options2, i =>//changed to array1.count from sHamMatrix.count
                {
                    double[] evs;
                    double[,] temp;
                    IECODE[i] = -1;

                    //add a parameter to count Lanczos iterations to set possible stopping criteria that way
                    //call MINVAL from here
                    if (!input.blockLanczos)//means use naiveLanczos routine
                    {
                        ITER[i] = input.noIts;
                        evs = new double[input.M + 1];
                        temp = new double[numcolumnsA[i], input.M + 1];
                        Lanczos.NaiveLanczos(ref evs, ref temp, array1[i], input.noIts, input.tol, input.pVector, i, input.filePath);
                    }
                    else//means use block Lanczos from SOCJT
                    {
                        evs = new double[input.M + 1];
                        temp = new double[numcolumnsA[i], input.M + 1];//changed here to numcolumnsA
                        IECODE[i] = -1;
                        ITER[i] = Lanczos.MINVAL(numcolumnsA[i], input.M + 1, input.kFactor, input.M, input.noIts, input.tol, 0, ref evs, ref temp, ref IECODE[i], array1[i], input.parVec);
                    }
                 
                    //initialize eigenvalues to have a length.                    
                    eigenvalues[i] = new double[evs.Length - 1];                    
                    for (int j = 0; j < evs.Length - 1; j++)                    
                    {                    
                        eigenvalues[i][j] = evs[j];                        
                    }
                    //I think this should be only for if block lanczos or naive lanczos with already calculated eigenvectors
                    if (input.blockLanczos || (!input.blockLanczos && array1[i].innerobj.m < Lanczos.basisSetLimit))
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
                    if (!input.blockLanczos && array1[0].innerobj.m >= Lanczos.basisSetLimit && input.pVector)
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
            if (!input.blockLanczos && array1[0].innerobj.m >= Lanczos.basisSetLimit && input.pVector)
            {
                //make it so that the output file generator does not try to print the values in the zmatrices which will be the eigenvectors of the lanczos matrix, not the hamiltonian
                input.pVector = false;
            }

            measurer.Stop();
            howMuchTime = measurer.ElapsedMilliseconds;
            input.diagTime = (double)howMuchTime / 1000D;
            #endregion

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

            //writes the eigenvectors to disk if that has been requested
            //if NTooBig == true then the separate eigenvector file will be made after the eigenvectors are calculated
            if (input.vecFile && (input.blockLanczos || input.pVector))
            {
                writeVecFile(input, zMatrices, JvecsForOutuput);
            }//end if

            //put code here to write eigenvectors to file with entire basis set.
            //do this by writing function to use jbasisvecsbyj which has all basis functions and
            //go through that list, where a basis function is in it and JvecsForOutput, pull the appropriate coefficient
            //from the zmatrices, otherwise just put in a 0.
            if (input.vecFileComplete)
            {
                writeVecComplete(jBasisVecsByJ, JvecsForOutuput, zMatrices, input, eigenvalues, isQuad);
            }

            List<string> linesToWrite = new List<string>();
            finalList = setAndSortEVs(eigenvalues, input.S, input.inclSO, zMatrices, JvecsForOutuput, input);//add the eigenvectors so that the symmetry can be included as well
            linesToWrite = OutputFile.makeOutput(input, zMatrices, array1, JvecsForOutuput, eigenvalues, isQuad, finalList, IECODE, ITER);                
            outp = linesToWrite;                
            return linesToWrite;   
        }//end SOCJT Routine

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
        public static void writeVecFile(FileInfo input, List<double[,]> zMatrices, List<List<BasisFunction>> JvecsForOutuput)
        {
            StringBuilder vecFile = new StringBuilder();
            vecFile.AppendLine("VecFile " + input.title);
            vecFile.AppendLine(" ");
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
                    OutputFile.vecBuilder(input, JvecsForOutuput[m], vecFile, zMatrices[m], o, 0.0);
                    vecFile.AppendLine(" ");
                }
                vecFile.AppendLine(" ");
            }
            List<string> vecFileOut = new List<string>();
            vecFileOut.Add(vecFile.ToString());
            File.WriteAllLines((input.filePath + input.title + "_vec.out"), vecFileOut);
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
            for (int i = 0; i < eigenvalues.Count; i++)
            {
                for (int j = 0; j < eigenvalues[i].Length; j++)
                {
                    eigenvalues[i][j] += ZPE;
                }
            }
            StringBuilder file = new StringBuilder();
            file.AppendLine(" ");
            file.AppendLine("Eigenvectors in complete basis set for " + input.title);
            file.AppendLine(" ");            
            //loop over jBlocks
            for (int i = 0; i < eigenvalues.Count; i++)
            {
                file.AppendLine(" ");
                file.AppendLine("J-Block " + ((decimal)i + 0.5M));
                file.AppendLine(" ");                

                //here make a list of possible j values in this eigenvector to save time
                List<decimal> possibleJVals = new List<decimal>();
                for (int n = i; n < input.maxJ; n += 3)
                {
                    possibleJVals.Add((decimal)n + 0.5M);
                    if (n == i)
                    {
                        continue;
                    }
                    possibleJVals.Add((decimal)n * -1M + (decimal)i * 2M + 0.5M);
                }//end loop over possible j values
                possibleJVals.Sort();

                //loop over all of the eigenvalues found
                for (int l = 0; l < eigenvalues[i].Length; l++)
                {                    
                    file.AppendLine(" " + "\r");
                    file.AppendLine("Eigenvalue" + "\t" + Convert.ToString(l + 1) + " = " + String.Format("{0,10:0.0000}", eigenvalues[i][l]));
                    bool a1 = SOCJT.isA(JvecsForOutput[i], tempMat[i], l, input, false);
                    if (a1)
                    {
                        file.AppendLine("Vector is Type 1");
                    }
                    else
                    {
                        file.AppendLine("Vector is Type 2");
                    }
                    file.AppendLine(" " + "\r");
                    file.Append("Coefficient" + "\t");
                    for (int h = 0; h < input.nModes; h++)
                    {
                        file.Append("v(" + Convert.ToString(h + 1) + ")" + "\t" + "l(" + Convert.ToString(h + 1) + ")" + "\t");
                    }
                    file.Append("lambda");
                    
                    for (int j = 0; j < jBasisVecsByJ.Count; j++)//goes through basis vectors
                    {
                        int place = 0;
                        
                        //first split into quadratic and linear
                        if (isQuad)
                        {
                            //first check to see if this J value will have nonzero coefficients
                            bool rightJ = false;
                            for (int b = 0; b < possibleJVals.Count; b++)
                            {
                                if (jBasisVecsByJ[j][0].J == possibleJVals[b])
                                {
                                    rightJ = true;
                                    break;
                                }
                            }//end loop over possible j values

                            //if this j value is not a j value with a nonzero coefficient just put 0 in
                            if (!rightJ)
                            { 
                                for(int k = 0; k < jBasisVecsByJ[j].Count; k++)
                                {
                                    writeVec(0.0, jBasisVecsByJ[j][k], file);
                                }
                            }//end if
                            else//means this j has possible nonzero matrix elements
                            {
                                for (int k = 0; k < jBasisVecsByJ[j].Count; k++)
                                {
                                    bool temp = isInBasis(jBasisVecsByJ[j][k], JvecsForOutput[i], ref place, 0);
                                    if (!temp)
                                    {
                                        writeVec(0.0, jBasisVecsByJ[j][k], file);
                                    }
                                    else
                                    {
                                        writeVec(tempMat[i][place, l], jBasisVecsByJ[j][k], file);
                                    }
                                }//end for loop over J
                            }//end else
                        }//end if isQuad
                        else
                        {
                            if (jBasisVecsByJ[j][0].J == (decimal)i + 0.5M)
                            {
                                //write actual values to the eigenvector
                                for (int k = 0; k < jBasisVecsByJ[j].Count; k++)
                                {
                                    writeVec(tempMat[i][k, l], jBasisVecsByJ[j][k], file);
                                }//end loop over BasisFunctions in j value for jBasisVecsByJ
                            }
                            else
                            {
                                //just write a zero for all of the coefficients for this j value
                                for (int k = 0; k < jBasisVecsByJ[j].Count; k++)
                                {
                                    writeVec(0.0, jBasisVecsByJ[j][k], file);
                                }//end loop over BasisFunctions in j value for jBasisVecsByJ
                            }
                        }//end linear option
                    }//end loop over jValues in jbasisVecsByJ
                    file.AppendLine("\r");
                }//end loop over the number of eigenvalues
            }//end loop over j blocks

            List<string> vecFileOut = new List<string>();
            vecFileOut.Add(file.ToString());
            File.WriteAllLines((input.filePath + input.title + "_CompleteVector.out"), vecFileOut);
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
        public static void writeVec(double coefficient, BasisFunction func, StringBuilder file)
        {
            file.AppendLine("\t");
            file.Append(String.Format("{0,10:0.000000}", coefficient));
            for (int m = 0; m < func.modesInVec.Count; m++)//goes through each mode
            {
                file.Append("\t" + "  " + Convert.ToString(func.modesInVec[m].v) + "\t" + String.Format("{0,3}", func.modesInVec[m].l));//  "  " + Convert.ToString(jBasisVecsByJ[i][h].modesInVec[m].l));
            }
            file.Append("\t" + String.Format("{0,4}", func.Lambda));
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
                    for(int j = 0; j < jBasisVec.modesInVec.Count; j++)
                    {
                    
                        if (jBasisVec.modesInVec[j].v == eigenvector[i].modesInVec[j].v && jBasisVec.modesInVec[j].l == eigenvector[i].modesInVec[j].l)
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
        public static Eigenvalue[] setAndSortEVs(List<double[]> evs, decimal S, bool inclSO, List<double[,]> zMatrices, List<List<BasisFunction>>jvecs, FileInfo input)
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
                    eigen.Add(new Eigenvalue(J, j + 1, tempS, evs[i][j], tbool));
                }
                if (tempS < maxS)
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
            double ZPE = eigenarray[0].Ev;
            int[] temp = new int[evs.Count];
            for (int i = 0; i < evs.Count; i++)
            {
                temp[i] = 1;
            }
            for (int i = 0; i < eigenarray.Length; i++)
            {
                eigenarray[i].Ev = eigenarray[i].Ev - ZPE;
            }
            int SOnumb = (int)(-2M * S) + 1;
            if (inclSO == false)
            {
                SOnumb = 1;
            }
            for (int i = 0; i < eigenarray.Length; i++)
            {
                int Snumb = (int)(eigenarray[i].Sig - S);
                int j = (int)(eigenarray[i].pJ - 0.5M);
                int place = j * SOnumb + Snumb;
                eigenarray[i].nJ = temp[place];
                temp[place]++;
            }
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
            if (!input.blockLanczos && !input.pVector && !overRide)
            {
                return false;
            }
            for (int m = 0; m < jBasisVecsByJ.Count; m++)
            {
                if (Math.Abs(tempMat[m, j]) > temp)
                {
                    for (int n = 0; n < input.nModes; n++)
                    {
                        tempVL[n * 2] = jBasisVecsByJ[m].modesInVec[n].v;
                        tempVL[n * 2 + 1] = jBasisVecsByJ[m].modesInVec[n].l;
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
                        if (jBasisVecsByJ[m].modesInVec[v].v == tempVL[v * 2])
                        {
                            if (jBasisVecsByJ[m].modesInVec[v].l == -1 * tempVL[v * 2 + 1])
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
                    if (arr[i].Ev > arr[i + 1].Ev)
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
            while(alglib.sparseenumerate(A, ref t0, ref t1, out i, out j, out oldVal))
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
        /// <param name="fitHamList">
        /// List of matrices. If matFile is to be used and exists this list will contain initialized off-diagonal matrices after function call.
        /// </param>
        /// <returns>
        /// Boolean indicating whether or not the fitHamList matrices have been initialized from matFile.
        /// </returns>
        private static List<int> matReadFunction(FileInfo input, ref bool matricesMade)
        {
            matricesMade = false;
            List<int> basisSizeList = new List<int>();
            if (input.useMatFile && input.matMade)
            {
                string[] matFile = { };
                try
                {
                    matFile = FileInfo.fileRead(input.matFilePath);
                }
                catch(FileNotFoundException)
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
                            alglib.sparseadd(B, Convert.ToInt32(matFile[i]), Convert.ToInt32(matFile[i + 1]), FileInfo.parseDouble(matFile[i + 2]));
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
    }//end class SOCJT
}
