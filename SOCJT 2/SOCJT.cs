using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace ConsoleApplication1
{
    class SOCJT
    {
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
            bool smallMat = true;
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

 
            #region LargeMatrix                
            List<alglib.sparsematrix> sHamMatrix = new List<alglib.sparsematrix>();            
            alglib.sparsematrix[] array1;            
            int[] numcolumnsA;
            smallMat = false;            
            //Creates the Hamiltonian matrices for linear cases            
            int numQuadMatrix = 0;            
            List<int> a = new List<int>();
            
          
            if (isQuad == false)            
            {            
                int h = 0;                
                array1 = new alglib.sparsematrix[jBasisVecsByJ.Count];                
                numcolumnsA = new int[jBasisVecsByJ.Count];
                measurer.Reset();
                measurer.Start();
                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = input.parJ;
                //for (decimal i = jMin; i <= jMax; i++)                
                Parallel.For((int)(jMin - 0.5M), (int)(jMax + 0.5M), options, i =>                
                {                
                    int nColumns;                    
                    if (jBasisVecsByJ[i].Count != 0)//changed from h to i                    
                    {   
                        //array1[i] = GenHamMat.genMatrix(jBasisVecsByJ[i], isQuad, input, out nColumns, true, input.parMat);
                        //replaced line above with conditionals below
                        if (input.specialHam)
                        {
                            array1[i] = GenHamMat.genMatrix2(jBasisVecsByJ[i], isQuad, input, out nColumns, true, input.parMat);
                        }
                        else
                        {
                            array1[i] = GenHamMat.genMatrix(jBasisVecsByJ[i], isQuad, input, out nColumns, true, input.parMat);
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
                    throw new BasisSetTooSmallException();
                }
            }//end if  

            //Creates the Hamiltonian matrices for quadratic cases.            
            else            
            { 
                int dynVar1 = (int)(jMax - 1.5M);
                int dynVar2 = dynVar1 / 3;
                    array1 = new alglib.sparsematrix[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                    numcolumnsA = new int[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                    jbasisoutA = new List<BasisFunction>[jBasisVecsByJ.Count - dynVar1 - jBasisVecsByJ.Count / 2];//changed to dynVar1 from 6
                    //this tells how much time has passed this could be used to time out different parts of code                
                measurer.Reset();
                measurer.Start();
                
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


                        //array1[i - jBasisVecsByJ.Count / 2] = GenHamMat.genMatrix(quadVecs, isQuad, input, out nColumns, true, input.parMat);
                        //replaced the line above with the conditionals below
                        if (input.specialHam)
                        {
                            array1[i - jBasisVecsByJ.Count / 2] = GenHamMat.genMatrix2(quadVecs, isQuad, input, out nColumns, true, input.parMat);
                        }
                        else
                        {
                            array1[i - jBasisVecsByJ.Count / 2] = GenHamMat.genMatrix(quadVecs, isQuad, input, out nColumns, true, input.parMat);
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
                    measurer.Stop();                    
                    howMuchTime = measurer.ElapsedMilliseconds;
                    input.matGenTime = (double)howMuchTime / 1000D;
                    if (a.Count > 0)
                    {
                        throw new BasisSetTooSmallException();
                    }
                    for (int i = 0; i < 2; i++)
                    {
                        zMatrices.Add(new double[0, 0]);
                        eigenvalues.Add(new double[0]);
                    }
            }//end else

            //add SO stuff here
            if (input.inclSO == true)
            {
                decimal minS = input.S * -1M;
                List<alglib.sparsematrix> tempMatList = new List<alglib.sparsematrix>();
                List<int> tempNumbColumns = new List<int>();
                if (isQuad == true)
                {
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
                                //double temp2 = alglib.sparseget(array1[i], k, k);
                                //alglib.sparseset(tempMatList[tempMatList.Count - 1], k, k, temp + temp2);
                                alglib.sparseadd(tempMatList[tempMatList.Count - 1], k, k, temp);
                            }//end loop over diagonal matrix elements
                            tempNumbColumns.Add(numcolumnsA[i]);
                            if (i > 0)
                            {
                                break;
                            }
                        }//end loop over values of S
                    }//end loop over all previouly made sparseMatrices               
                }//end if
                else
                {
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
                        }//end loop over values of S
                    }//end loop over all previouly made sparseMatrices
                }//end else

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

            for (int i = 0; i < array1.Length; i++)
            {
                alglib.sparseconverttocrs(array1[i]);
            }

                int[] IECODE = new int[array1.Length];
                int[] ITER = new int[array1.Length];
                //actually diagonalizes the Hamiltonian matrix
                measurer.Reset();
                measurer.Start();
                ParallelOptions options2 = new ParallelOptions();
                options2.MaxDegreeOfParallelism = input.parJ;
                Parallel.For(0, array1.Length, options2, i =>//changed to array1.count from sHamMatrix.count
                {
                    //this is where multithreading is needed
                    double[] evs = new double[input.M + 1];
                    double[,] temp = new double[numcolumnsA[i], input.M + 1];//changed here to numcolumnsA
                    IECODE[i] = -1;

                    //call MINVAL from here                    
                    ITER[i] = Lanczos.MINVAL(numcolumnsA[i], input.M + 1, input.kFactor, input.M, input.noIts, input.tol, 0, ref evs, ref temp, ref IECODE[i], array1[i], input.parVec);
                                        
                    //initialize eigenvalues to have a length.                    
                    eigenvalues[i] = new double[evs.Length - 1];                    
                    for (int j = 0; j < evs.Length - 1; j++)                    
                    {                    
                        eigenvalues[i][j] = evs[j];                        
                    }                    
                    zMatrices[i] = new double[numcolumnsA[i], input.M];                    
                    if (input.pVector == true)                    
                    {                    
                        for (int j = 0; j < numcolumnsA[i]; j++)                        
                        {                        
                            for (int k = 0; k < input.M; k++)                            
                            {                            
                                zMatrices[i][j, k] = temp[j, k];                                
                            }                            
                        }                        
                    }                    
                    else                    
                    {                    
                        zMatrices[i] = null;                        
                    }
                    
                    temp = null;                    
                    evs = null;                    
                }//end for
                );
                measurer.Stop();
                howMuchTime = measurer.ElapsedMilliseconds;
                input.diagTime = (double)howMuchTime / 1000D;
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
                finalList = Eigenvalue.setAndSortEVs(eigenvalues, input.S, input.inclSO);                
                //dummy hamiltonian matrix list for outuput file generator                
                List<double[,]> hamMatrices = new List<double[,]>();                
                sHamMatrix = array1.ToList();                
                linesToWrite = OutputFile.makeOutput(input, zMatrices, hamMatrices, sHamMatrix, JvecsForOutuput, eigenvalues, isQuad, numColumns, finalList, true, IECODE, ITER);                
                outp = linesToWrite;                
                return linesToWrite;                
                #endregion
        }
    }
}
