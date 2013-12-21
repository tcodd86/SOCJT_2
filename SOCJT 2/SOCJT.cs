using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;

namespace ConsoleApplication1
{
    class SOCJT
    {
        private static bool matricesMade = false;

        //keep this so that on concurrent calls the matrix does not need to be regenerated
        public static List<List<alglib.sparsematrix>> fitHamList { get; private set; }

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
                jBasisVecsByJ.Add(GenHamMat.sortByJ(hamiltonianVecs, i));
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
                measurer.Reset();
                measurer.Start();
                int h = 0;                
                array1 = new alglib.sparsematrix[jBasisVecsByJ.Count];
                if (!matricesMade)
                {
                    fitHamList = new List<List<alglib.sparsematrix>>();
                    for (int m = 0; m < jBasisVecsByJ.Count; m++)
                    {
                        fitHamList.Add(new List<alglib.sparsematrix>());
                    }
                    //check here to see if matrix file should be used, if so then read from file and make them, set matricesMade to true
                    matricesMade = matReadFunction(input);
                }
                numcolumnsA = new int[jBasisVecsByJ.Count];

                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.parJ;          
                Parallel.For((int)(jMin - 0.5M), (int)(jMax + 0.5M), options, i =>                
                {                
                    int nColumns;                    
                    if (jBasisVecsByJ[i].Count != 0)//changed from h to i                    
                    {
                        if (!matricesMade)//if matrices not made then generate all matrices
                        {
                            fitHamList[i] = GenHamMat.genFitMatrix(jBasisVecsByJ[i], isQuad, input, out nColumns, input.parMat, false);
                            matricesMade = true;
                        }
                        else//this makes sure that the diagonal portion is regenerated on each call.
	                    {
                            fitHamList[i][0] = GenHamMat.genFitMatrix(jBasisVecsByJ[i], isQuad, input, out nColumns, input.parMat, true)[0];		 
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
                measurer.Stop();
                howMuchTime = measurer.ElapsedMilliseconds;
                input.matGenTime = (double)howMuchTime / 1000D;
                for (int i = (int)(jMin - 0.5M); i < (int)(jMax + 0.5M); i++)
                {
                    zMatrices.Add(new double[0, 0]);
                    eigenvalues.Add(new double[0]);
                }

                //handles errors where the basis set is too small
                if (a.Count > 0)
                {
                    throw new BasisSetTooSmallException(true);
                }
            }//end if  

            //Creates the Hamiltonian matrices for quadratic cases.            
            else            
            {
                measurer.Reset();
                measurer.Start();
                int dynVar1 = (int)(jMax - 1.5M);
                int dynVar2 = dynVar1 / 3;
                array1 = new alglib.sparsematrix[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                if (!matricesMade)
                {
                    fitHamList = new List<List<alglib.sparsematrix>>();
                    for (int m = 0; m < jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2; m++)
                    {
                        fitHamList.Add(new List<alglib.sparsematrix>());
                    }
                    //check here to see if matrix file should be used, if so then read from file and make them, set matricesMade to true
                    matricesMade = matReadFunction(input);
                }
                numcolumnsA = new int[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                jbasisoutA = new List<BasisFunction>[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                //this tells how much time has passed this could be used to time out different parts of code      

                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.parJ;
                
                Parallel.For(jBasisVecsByJ.Count / 2, jBasisVecsByJ.Count - dynVar1, options, i =>//changed to dynVar1 from 6
                {
                    List<BasisFunction> quadVecs = new List<BasisFunction>();
                    int nColumns;
                    for (int v = -dynVar2; v <= dynVar2; v++)
                    {
                        quadVecs.AddRange(jBasisVecsByJ[i + v * 3]);
                    }

                    //made specialHam matrix the default and not optional
                    if (!matricesMade)
                    {
                        fitHamList[i - jBasisVecsByJ.Count / 2] = GenHamMat.genFitMatrix(quadVecs, isQuad, input, out nColumns, input.parMat, false);                        
                    }                        
                    else//this makes sure that the diagonal portion is regenerated on each call.
	                {
                        fitHamList[i - jBasisVecsByJ.Count / 2][0] = GenHamMat.genFitMatrix(quadVecs, isQuad, input, out nColumns, input.parMat, true)[0];		 
	                }

                    jbasisoutA[i - jBasisVecsByJ.Count / 2] = quadVecs;
                    numcolumnsA[i - jBasisVecsByJ.Count / 2] = nColumns;
                    //checks to make sure that
                    if (numcolumnsA[i - jBasisVecsByJ.Count / 2] < input.M)
                    {
                        a.Add(0);
                    }
                    numQuadMatrix++;
                }
                );
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
            }//end else
            #endregion            

            //list where each element of mat is a list of alglib.sparsematrix objects.  One for each off-diagonal parameter (D, K, B)
            var mat = new List<List<alglib.sparsematrix>>();
            bool bilinear = false;
            var biAVecPos = new List<int>();
            var biEVecPos = new List<int>();
            GenHamMat.crossTermInitialization(jBasisVecsByJ[0][0].modesInVec, input.nModes, out bilinear, out biAVecPos, out biEVecPos, input.crossTermMatrix);
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
            ParallelOptions options2 = new ParallelOptions();
            options2.MaxDegreeOfParallelism = input.parJ;
            try
            {
                Parallel.For(0, array1.Length, options2, i =>//changed to array1.count from sHamMatrix.count
                {
                    //this is where multithreading is needed
                    double[] evs;
                    double[,] temp;//changed here to numcolumnsA
                    IECODE[i] = -1;

                    //add a parameter to count Lanczos iterations to set possible stopping criteria that way
                    //call MINVAL from here
                    if (!input.blockLanczos)//means use naiveLanczos routine
                    {
                        ITER[i] = input.noIts;
                        evs = new double[input.M];
                        temp = new double[numcolumnsA[i], input.M];
                        Lanczos.NaiveLanczos(ref evs, ref temp, array1[i], input.noIts, input.tol, input.oldRandom, input.pVector);
                    }
                    else//means use block Lanczos from SOCJT
                    {
                        evs = new double[input.M + 1];
                        temp = new double[numcolumnsA[i], input.M + 1];//changed here to numcolumnsA
                        IECODE[i] = -1;
                        ITER[i] = Lanczos.MINVAL(numcolumnsA[i], input.M + 1, input.kFactor, input.M, input.noIts, input.tol, 0, ref evs, ref temp, ref IECODE[i], array1[i], input.parVec, input.oldRandom);
                    }
                 
                    //initialize eigenvalues to have a length.                    
                    eigenvalues[i] = new double[evs.Length - 1];                    
                    for (int j = 0; j < evs.Length - 1; j++)                    
                    {                    
                        eigenvalues[i][j] = evs[j];                        
                    }
                    zMatrices[i] = new double[numcolumnsA[i], evs.Length - 1];//changed input.M to evs.Length - 1
                    for (int j = 0; j < numcolumnsA[i]; j++)                        
                    {
                        for (int k = 0; k < evs.Length - 1; k++)//changed input.M to evs.Length - 1                            
                        {                            
                            zMatrices[i][j, k] = temp[j, k];                                
                        }                            
                    }
                    
                    temp = null;                    
                    evs = null;                    
                }//end for
                );
            }
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
                    
            List<string> linesToWrite = new List<string>();
            finalList = setAndSortEVs(eigenvalues, input.S, input.inclSO, zMatrices, JvecsForOutuput, input);//add the eigenvectors so that the symmetry can be included as well
            linesToWrite = OutputFile.makeOutput(input, zMatrices, array1, JvecsForOutuput, eigenvalues, isQuad, finalList, IECODE, ITER);                
            outp = linesToWrite;                
            return linesToWrite;   
        }//end SOCJT Routine

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
                    bool tbool = isA(jvecs[i], zMatrices[i], j, input);
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

        public static bool isA(List<BasisFunction> jBasisVecsByJ, double[,] tempMat, int j, FileInfo input)
        {
            bool a1 = false;
            double temp = 0.0;
            int[] tempVL = new int[input.nModes * 2 + 1];
            //this conditional is because if naive lanczos is used and eigenvectors are not calculated the tempMat actually contains the eigenvectors of the the lanczos matrix, not the Hamiltonian
            if (!input.blockLanczos && !input.pVector)
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
        }//end method bublleSort

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
        private static bool matReadFunction(FileInfo input)
        {
            bool matricesMade = false;
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
                //string[] matFile = FileInfo.fileRead(input.matFilePath);
                int index = 0, basisSize = 0;
                for (int i = 0; i < matFile.Length; i++)
                {
                    if (matFile[i] == "List")
                    {
                        index = Convert.ToInt32(matFile[i + 1]);
                        basisSize = Convert.ToInt32(matFile[i + 2]);
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
            return matricesMade;
        }//end matReadFunction
    }//end class SOCJT
}
