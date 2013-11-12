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
        /// <summary>
        /// Generates an alglib sparsematrix object to be used in Lanczos diagonalization.
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
        /// <param name="Large">
        /// boolean value so that genMatrix method could be overloaded.
        /// </param>
        /// <returns></returns>
        public static alglib.sparsematrix genMatrix(List<BasisFunction> basisVectorsByJ, bool isQuad, FileInfo input, out int nColumns, bool Large, int par)
        {
            int matSize = basisVectorsByJ.Count;
            nColumns = matSize;
            int r = 0;
            //double temp = 0;

            alglib.sparsematrix A = new alglib.sparsematrix();
            alglib.sparsecreate(matSize, matSize, 10, out A);
            bool containsAVecs = false;
            bool bilinear = true;
            List<int> AVecPos = new List<int>();
            List<int> EVecPos = new List<int>();
            List<int> biAVecPos = new List<int>();
            List<int> biEVecPos = new List<int>();
            //List<Tuple<int, int, double>> matPos = new List<Tuple<int, int, double>>();
            ConcurrentBag<Tuple<int, int, double>> matPos = new ConcurrentBag<Tuple<int, int, double>>();
            int[] change = new int[3];

            #region CrossTermInitializationStuff
            if (input.includeCrossTerms == true)
            {

                for (int i = 0; i < input.nModes; i++)
                {
                    if (basisVectorsByJ[0].modesInVec[i].symmetryIsA == true)
                    {
                        AVecPos.Add(i);
                        containsAVecs = true;
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
                    //new lists just for bilinear coupling
                    //biAVecPos = AVecPos;
                    //biEVecPos = EVecPos;
                    for (int i = 0; i < AVecPos.Count; i++)
                    {
                        biAVecPos.Add(AVecPos[i]);
                    }
                    for (int i = 0; i < EVecPos.Count; i++)
                    {
                        biEVecPos.Add(EVecPos[i]);
                    }

                    //loop to eliminate any A vectors that have no cross-coupling
                    for (int i = 0; i < biAVecPos.Count; i++)
                    {
                        for (int j = 0; j < biEVecPos.Count; j++)
                        {
                            if (biAVecPos[i] > biEVecPos[j])
                            {
                                if (input.crossTermMatrix[biEVecPos[j], biAVecPos[i]] != 0)
                                {
                                    break;
                                }
                            }
                            else
                            {
                                if (input.crossTermMatrix[biAVecPos[i], biEVecPos[j]] != 0)
                                {
                                    break;
                                }
                            }
                            if (j == biEVecPos.Count - 1)
                            {
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
                                if (input.crossTermMatrix[biEVecPos[i], biAVecPos[j]] != 0)
                                {
                                    break;
                                }
                            }
                            else
                            {
                                if (input.crossTermMatrix[biAVecPos[j], biEVecPos[i]] != 0)
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

                    if (biAVecPos.Count == 0 && biEVecPos.Count == 0)
                    {
                        bilinear = false;
                    }
                }

                //stuff for cross quadratic
                /*
                for (int i = 0; i < EVecPos.Count; i++)
                {
                    for (int j = i + 1; j < EVecPos.Count; j++)
                    {
                        if (input.crossTermMatrix[EVecPos[i], EVecPos[j]] != 0)
                        {
                            break;
                        }
                        if (j == EVecPos.Count - 1)
                        {
                            EVecPos.RemoveAt(i);
                            i--;
                            break;
                        }
                    }
                }
                */
            }//end if CrossTerms == true
            #endregion

            //generates the Hamiltonian Matrix
            //just the HO terms
            if (input.AT == true)
            {
                for (int n = 0; n < matSize; n++)
                {
                    #region HO Terms + Cross
                    //one mode harmonic and anharmonic terms
                    for (int i = 0; i < input.nModes; i++)
                    {
                        int degeneracy = 2;
                        if (basisVectorsByJ[n].modesInVec[i].symmetryIsA)
                        {
                            degeneracy = 1;
                        }
                        double temp = basisVectorsByJ[n].modesInVec[i].modeOmega * ((double)basisVectorsByJ[n].modesInVec[i].v + (double)degeneracy / 2D) - basisVectorsByJ[n].modesInVec[i].anharmonicity * Math.Pow(((double)basisVectorsByJ[n].modesInVec[i].v + (double)degeneracy / 2), 2);//I deleted the (double) from the degeneracy / 2
                        alglib.sparseadd(A, n, n, temp);//changed from (A, n, m, temp)
                        //temp = 0;
                    }

                    //cross anharmonic terms for interaction between two A modes
                    if (input.AT == true)
                    {
                        #region CrossAnharmonic
                        for (int i = 1; i < input.nModes; i++)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (basisVectorsByJ[n].modesInVec[i].symmetryIsA == false || basisVectorsByJ[n].modesInVec[j].symmetryIsA == false)
                                {
                                    continue;
                                }
                                if (input.crossTermMatrix[i, j] != 0D)
                                {
                                    //added the *4 to this
                                    double temp = input.crossTermMatrix[i, j] * 4D * ((double)basisVectorsByJ[n].modesInVec[i].v + (double)1 / 2) * ((double)basisVectorsByJ[n].modesInVec[j].v + (double)1 / 2);//took off the -1
                                    alglib.sparseadd(A, n, n, temp);//changed from (A, n, m, temp)
                                    //temp = 0;
                                }//end if
                            }//end j for loop
                        }//end i for loop
                        #endregion
                    }//end if AT == true
                    continue;
                    #endregion
                }
            }//end for if AT = true

            else
            {
                for (int n = 0; n < matSize; n++)
                {
                    #region HO Terms
                    //one mode harmonic and anharmonic terms
                    for (int i = 0; i < input.nModes; i++)
                    {
                        int degeneracy = 2;
                        if (basisVectorsByJ[n].modesInVec[i].symmetryIsA)
                        {
                            degeneracy = 1;
                        }
                        double temp = basisVectorsByJ[n].modesInVec[i].modeOmega * ((double)basisVectorsByJ[n].modesInVec[i].v + (double)degeneracy / 2D) - basisVectorsByJ[n].modesInVec[i].anharmonicity * Math.Pow(((double)basisVectorsByJ[n].modesInVec[i].v + (double)degeneracy / 2), 2);//I deleted the (double) from the degeneracy / 2
                        alglib.sparseadd(A, n, n, temp);//changed from (A, n, m, temp)
                        //temp = 0;
                    }//end if AT == true
                    continue;
                    #endregion
                }
            }

            //now the JT terms

            var rangePartitioner = Partitioner.Create(0, matSize);
            ParallelOptions parOp = new ParallelOptions();
            parOp.MaxDegreeOfParallelism = par;
            Parallel.ForEach(rangePartitioner, parOp, (range, loopState) =>
            //Parallel.For(0, par, hh =>
                {
                    //for (int n = hh * (matSize / par); n < ((hh + 1) * matSize) / par; n++)
                    //for(int n = 0; n < matSize; n++)
                    for(int n = range.Item1; n < range.Item2; n++)
                    {
                        for (int m = n + 1; m < matSize; m++)//changed from r + 1
                        {
                            double temp;
                            int numberOfChanges = modeChangeCounter(basisVectorsByJ, input.nModes, n, m);
                            //makes it skip any places where matrix elements are 0.
                            if (numberOfChanges > 2)
                            {
                                continue;
                            }

                            //indexes n and m are for the rows and columns of the matrix respectively
                            if (basisVectorsByJ[n].Lambda == basisVectorsByJ[m].Lambda && input.Special == false && input.AT == false)//checks to see if it can have an off diagonal term              
                            {
                                continue;
                            }

                            if (input.Special == true && numberOfChanges <= 2 && basisVectorsByJ[n].Lambda == basisVectorsByJ[m].Lambda)
                            {
                                #region Special
                                if (numberOfChanges == 1)
                                {
                                    if (basisVectorsByJ[n].modesInVec[0].v + 1 == basisVectorsByJ[m].modesInVec[0].v)
                                    {
                                        temp = input.crossTermMatrix[0, 0] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[0].v + 1D) * (2D * (double)basisVectorsByJ[n].modesInVec[1].v + 1D);//mode two stuff
                                        Tuple<int, int, double> tTemp = new Tuple<int,int,double>(n, m, temp);
                                        matPos.Add(tTemp);
                                        continue;
                                    }
                                    /*
                                    if (basisVectorsByJ[n].modesInVec[0].v - 1 == basisVectorsByJ[m].modesInVec[0].v)
                                    {
                                        temp = input.crossTermMatrix[0, 0] * Math.Sqrt(basisVectorsByJ[n].modesInVec[0].modeOmega * (double)basisVectorsByJ[n].modesInVec[0].v) * basisVectorsByJ[n].modesInVec[1].modeOmega * (2D * (double)basisVectorsByJ[n].modesInVec[1].v + 1D);//mode two stuff
                                        alglib.sparseadd(A, n, m, temp);
                                        temp = 0;
                                        continue;
                                    }
                                    */
                                }
                                if (numberOfChanges == 2)
                                {
                                    if (basisVectorsByJ[n].modesInVec[0].v + 1 == basisVectorsByJ[m].modesInVec[0].v)
                                    {
                                        /*
                                        if (basisVectorsByJ[n].modesInVec[1].v - 2 == basisVectorsByJ[m].modesInVec[1].v)
                                        {
                                            temp = input.crossTermMatrix[0, 0] * Math.Sqrt(basisVectorsByJ[n].modesInVec[0].modeOmega * (double)basisVectorsByJ[n].modesInVec[0].v + 1D) * basisVectorsByJ[n].modesInVec[1].modeOmega * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[1].v - 1D) * (double)basisVectorsByJ[n].modesInVec[1].v);
                                            alglib.sparseadd(A, n, m, temp);
                                            temp = 0;
                                            continue;
                                        }
                                        */
                                        if (basisVectorsByJ[n].modesInVec[1].v + 2 == basisVectorsByJ[m].modesInVec[1].v)
                                        {
                                            temp = input.crossTermMatrix[0, 0] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[0].v + 1D) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[1].v + 2D) * ((double)basisVectorsByJ[n].modesInVec[1].v + 1D));
                                            Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                            matPos.Add(tTemp);
                                            continue;
                                        }
                                    }
                                    if (basisVectorsByJ[n].modesInVec[0].v - 1 == basisVectorsByJ[m].modesInVec[0].v)
                                    {
                                        if (basisVectorsByJ[n].modesInVec[1].v + 2 == basisVectorsByJ[m].modesInVec[1].v)
                                        {
                                            temp = input.crossTermMatrix[0, 0] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[0].v) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[1].v + 2D) * ((double)basisVectorsByJ[n].modesInVec[1].v + 1D));
                                            Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                            matPos.Add(tTemp);
                                            continue;
                                        }
                                        /*
                                        if (basisVectorsByJ[n].modesInVec[1].v - 2 == basisVectorsByJ[m].modesInVec[1].v)
                                        {
                                            temp = input.crossTermMatrix[0, 0] * Math.Sqrt(basisVectorsByJ[n].modesInVec[0].modeOmega * (double)basisVectorsByJ[n].modesInVec[0].v) * basisVectorsByJ[n].modesInVec[1].modeOmega * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[1].v - 1D) * (double)basisVectorsByJ[n].modesInVec[1].v);
                                            alglib.sparseadd(A, n, m, temp);
                                            temp = 0;
                                            continue;
                                        }
                                        */
                                    }

                                    continue;
                                }
                                #endregion
                            }

                            if (input.AT == true && basisVectorsByJ[n].Lambda == basisVectorsByJ[m].Lambda && numberOfChanges <= 2)
                            {
                                #region Cross-AT
                                for (int i = 1; i < input.nModes; i++)//rows
                                {
                                    if (basisVectorsByJ[n].modesInVec[i].symmetryIsA == false)
                                    {
                                        continue;
                                    }
                                    for (int j = 0; j < i; j++)//columns
                                    {
                                        if (basisVectorsByJ[n].modesInVec[j].symmetryIsA == false)
                                        {
                                            continue;
                                        }
                                        if (input.crossTermMatrix[i, j] != 0D)
                                        {
                                            if (basisVectorsByJ[n].modesInVec[i].v + 2 == basisVectorsByJ[m].modesInVec[i].v)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[j].v + 2 == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[i].v * (double)basisVectorsByJ[n].modesInVec[i].v + 3D * (double)basisVectorsByJ[n].modesInVec[i].v + 2D) * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[j].v * (double)basisVectorsByJ[n].modesInVec[j].v + 3D * (double)basisVectorsByJ[n].modesInVec[j].v + 2D);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                                if (basisVectorsByJ[n].modesInVec[j].v == basisVectorsByJ[m].modesInVec[j].v && numberOfChanges == 1)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[i].v * (double)basisVectorsByJ[n].modesInVec[i].v + 3D * (double)basisVectorsByJ[n].modesInVec[i].v + 2D) * 2D * ((double)basisVectorsByJ[n].modesInVec[j].v + 0.5);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                                if (basisVectorsByJ[n].modesInVec[j].v - 2 == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[i].v * (double)basisVectorsByJ[n].modesInVec[i].v + 3D * (double)basisVectorsByJ[n].modesInVec[i].v + 2D) * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[j].v * (double)basisVectorsByJ[n].modesInVec[j].v - (double)basisVectorsByJ[n].modesInVec[j].v);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                            }
                                            if (basisVectorsByJ[n].modesInVec[i].v == basisVectorsByJ[m].modesInVec[i].v && numberOfChanges == 1)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[j].v + 2 == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * 2D * ((double)basisVectorsByJ[n].modesInVec[i].v + 0.5) * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[j].v * (double)basisVectorsByJ[n].modesInVec[j].v + 3D * (double)basisVectorsByJ[n].modesInVec[j].v + 2D);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                                /*
                                                if (basisVectorsByJ[n].modesInVec[j].v == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * 4D * ((double)basisVectorsByJ[n].modesInVec[i].v + 0.5) * ((double)basisVectorsByJ[n].modesInVec[j].v + 0.5);
                                                    alglib.sparseadd(A, n, m, temp);
                                                    temp = 0;
                                                }
                                                */
                                                if (basisVectorsByJ[n].modesInVec[j].v - 2 == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * 2D * ((double)basisVectorsByJ[n].modesInVec[i].v + 0.5) * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[j].v * (double)basisVectorsByJ[n].modesInVec[j].v - (double)basisVectorsByJ[n].modesInVec[j].v);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                            }
                                            if (basisVectorsByJ[n].modesInVec[i].v - 2 == basisVectorsByJ[m].modesInVec[i].v)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[j].v + 2 == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[i].v * (double)basisVectorsByJ[n].modesInVec[i].v - (double)basisVectorsByJ[n].modesInVec[i].v) * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[j].v * (double)basisVectorsByJ[n].modesInVec[j].v + 3D * (double)basisVectorsByJ[n].modesInVec[j].v + 2D);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                                if (basisVectorsByJ[n].modesInVec[j].v == basisVectorsByJ[m].modesInVec[j].v && numberOfChanges == 1)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[i].v * (double)basisVectorsByJ[n].modesInVec[i].v - (double)basisVectorsByJ[n].modesInVec[i].v) * 2D * ((double)basisVectorsByJ[n].modesInVec[j].v + 0.5);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                                if (basisVectorsByJ[n].modesInVec[j].v - 2 == basisVectorsByJ[m].modesInVec[j].v)
                                                {
                                                    temp = input.crossTermMatrix[i, j] * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[i].v * (double)basisVectorsByJ[n].modesInVec[i].v - (double)basisVectorsByJ[n].modesInVec[i].v) * Math.Sqrt((double)basisVectorsByJ[n].modesInVec[j].v * (double)basisVectorsByJ[n].modesInVec[j].v - (double)basisVectorsByJ[n].modesInVec[j].v);
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }
                                            }
                                            //go through possible nonzero terms
                                        }
                                    }
                                }
                                #endregion
                            }

                            //calculate JT terms
                            if (numberOfChanges == 1)//changed this to one from zero                    
                            {
                                for (int i = 0; i < input.nModes; i++)
                                {
                                    if (basisVectorsByJ[n].modesInVec[i].v == basisVectorsByJ[m].modesInVec[i].v && basisVectorsByJ[n].modesInVec[i].l == basisVectorsByJ[m].modesInVec[i].l)
                                    {
                                        continue;
                                    }//this just makes it skip the extra conditionals if it's not the mode with the potential JT term                            
                                    #region linear terms
                                    if (basisVectorsByJ[n].modesInVec[i].v + 1 == basisVectorsByJ[m].modesInVec[i].v && basisVectorsByJ[n].modesInVec[i].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[i].l)
                                    {
                                        temp = basisVectorsByJ[n].modesInVec[i].modeOmega * Math.Sqrt(basisVectorsByJ[n].modesInVec[i].DBasis * ((double)basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[i].l + 2D));
                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                        matPos.Add(tTemp);
                                        break;//I think I can do this
                                    }//end first if
                                    if (basisVectorsByJ[n].modesInVec[i].v + 1 == basisVectorsByJ[m].modesInVec[i].v && basisVectorsByJ[n].modesInVec[i].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[i].l)
                                    {
                                        temp = basisVectorsByJ[n].modesInVec[i].modeOmega * Math.Sqrt(basisVectorsByJ[n].modesInVec[i].DBasis * ((double)basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[i].l + 2D));
                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                        matPos.Add(tTemp);
                                        break;
                                    }//end second if

                                    //I don't think this code is ever used since I only have the upper triangle of the matrix being used
                                    /*
                                    if (basisVectorsByJ[n].modesInVec[i].v - 1 == basisVectorsByJ[m].modesInVec[i].v && basisVectorsByJ[n].modesInVec[i].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[i].l)
                                    {
                                        matrix[n, m] += basisVectorsByJ[n].modesInVec[i].modeOmega * Math.Sqrt(basisVectorsByJ[n].modesInVec[i].DBasis * ((double)basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[i].l));//don't think this 2 should be here
                                        break;
                                    }//end third if
                                    if (basisVectorsByJ[n].modesInVec[i].v - 1 == basisVectorsByJ[m].modesInVec[i].v && basisVectorsByJ[n].modesInVec[i].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[i].l)
                                    {
                                        matrix[n, m] += basisVectorsByJ[n].modesInVec[i].modeOmega * Math.Sqrt(basisVectorsByJ[n].modesInVec[i].DBasis * ((double)basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[i].l));//don't think this 2 should be here
                                        break;
                                    }//end fourth if
                                    */
                                    #endregion
                                    if (isQuad == true)
                                    {
                                        #region quadratic terms
                                        if (basisVectorsByJ[n].modesInVec[i].v - 2 == basisVectorsByJ[m].modesInVec[i].v)
                                        {

                                            if (basisVectorsByJ[n].modesInVec[i].l + 2 * Math.Pow(-1D, input.S2) == basisVectorsByJ[m].modesInVec[i].l)
                                            {
                                                temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].KBasis / 4D * Math.Sqrt((basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l) * (basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l - 2)));
                                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                matPos.Add(tTemp);
                                            }//end if for top terms in +/-
                                            if (basisVectorsByJ[n].modesInVec[i].l - 2 * Math.Pow(-1D, input.S2) == basisVectorsByJ[m].modesInVec[i].l)
                                            {
                                                temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].KBasis / 4D * Math.Sqrt((basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l) * (basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l - 2)));
                                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                matPos.Add(tTemp);
                                            }//end if for bottom terms in +/-
                                            break;
                                        }//end if for top quadratic term (1 of 3)


                                        if (basisVectorsByJ[n].modesInVec[i].v == basisVectorsByJ[m].modesInVec[i].v)
                                        {
                                            if (basisVectorsByJ[n].modesInVec[i].l + 2 * Math.Pow(-1D, input.S2) == basisVectorsByJ[m].modesInVec[i].l)
                                            {
                                                temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].KBasis / 2D * Math.Sqrt((basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l + 2) * (basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l)));
                                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                matPos.Add(tTemp);
                                            }//end if for top terms in +/-
                                            if (basisVectorsByJ[n].modesInVec[i].l - 2 * Math.Pow(-1D, input.S2) == basisVectorsByJ[m].modesInVec[i].l)
                                            {
                                                temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].KBasis / 2D * Math.Sqrt((basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l + 2) * (basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l)));
                                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                matPos.Add(tTemp);
                                            }//end if for bottom terms in +/-
                                            break;
                                        }//end if for top quadratic term (2 of 3)

                                        if (basisVectorsByJ[n].modesInVec[i].v + 2 == basisVectorsByJ[m].modesInVec[i].v)
                                        {

                                            if (basisVectorsByJ[n].modesInVec[i].l + 2 * Math.Pow(-1D, input.S2) == basisVectorsByJ[m].modesInVec[i].l)
                                            {
                                                temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].KBasis / 4D * Math.Sqrt((basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l + 4D) * (basisVectorsByJ[n].modesInVec[i].v + Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l + 2)));
                                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                matPos.Add(tTemp);
                                            }//end if for top terms in +/-
                                            if (basisVectorsByJ[n].modesInVec[i].l - 2 * Math.Pow(-1D, input.S2) == basisVectorsByJ[m].modesInVec[i].l)
                                            {
                                                temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].KBasis / 4D * Math.Sqrt((basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l + 4D) * (basisVectorsByJ[n].modesInVec[i].v - Math.Pow(-1, input.S2) * basisVectorsByJ[n].modesInVec[i].l + 2)));
                                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                matPos.Add(tTemp);
                                            }//end if for bottom terms in +/-
                                            break;
                                        }//end if for top quadratic term (3 of 3)
                                        #endregion
                                    }//end if                            
                                }//end for loop
                            }//end of JT terms part

                            #region Cross Terms
                            if (input.includeCrossTerms == true && numberOfChanges == 2)
                            {
                                //A-E interaction terms
                                if (containsAVecs == true)
                                {
                                    if (bilinear == true)
                                    {
                                        #region Bilinear
                                        for (int a = 0; a < biAVecPos.Count; a++)
                                        {
                                            for (int e = 0; e < biEVecPos.Count; e++)
                                            {
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
                                                if (input.crossTermMatrix[row, column] == 0D)
                                                {
                                                    continue;
                                                }
                                                if (basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1 == basisVectorsByJ[m].modesInVec[biAVecPos[a]].v)
                                                {
                                                    ////*
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;//I think I can do this
                                                        //Math.Sqrt(basisVectorsByJ[n].modesInVec[biAVecPos[a]].modeOmega * basisVectorsByJ[n].modesInVec[biEVecPos[e]].modeOmega) * 
                                                        //this is what was removed from these matrix elements
                                                    }//end first if
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;
                                                    }//end second if
                                                    //*/
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;
                                                    }//end third if
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;
                                                    }//end fourth if
                                                }
                                                if (basisVectorsByJ[n].modesInVec[biAVecPos[a]].v - 1 == basisVectorsByJ[m].modesInVec[biAVecPos[a]].v)
                                                {
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;//I think I can do this
                                                    }//end first if
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;
                                                    }//end second if
                                                    ////*
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;
                                                    }//end third if
                                                    if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                    {
                                                        temp = 0.5 * input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                        matPos.Add(tTemp);
                                                        continue;
                                                    }//end fourth if
                                                }
                                            }
                                        }
                                        #endregion
                                    }//end bilinear if                                    
                                }//end A-E interaction

                                if (EVecPos.Count > 1)
                                {
                                    #region Cross-Quadratic
                                    int ps = -1;
                                    if (input.S1 == 0)
                                    {
                                        ps = 1;
                                    }
                                    for (int i = 0; i < EVecPos.Count; i++)
                                    {
                                        for (int j = i + 1; j < EVecPos.Count; j++)
                                        {
                                            if (basisVectorsByJ[n].modesInVec[EVecPos[i]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[i]].v && basisVectorsByJ[n].modesInVec[EVecPos[i]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[i]].l)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end first if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end second if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end third if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end fourth if
                                            }//end first if
                                            if (basisVectorsByJ[n].modesInVec[EVecPos[i]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[i]].v && basisVectorsByJ[n].modesInVec[EVecPos[i]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[i]].l)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end first if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end second if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end third if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l + 2D)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end fourth if
                                            }//end second if
                                            if (basisVectorsByJ[n].modesInVec[EVecPos[i]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[i]].v && basisVectorsByJ[n].modesInVec[EVecPos[i]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[i]].l)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end first if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end second if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end third if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end fourth if
                                            }//end third if
                                            if (basisVectorsByJ[n].modesInVec[EVecPos[i]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[i]].v && basisVectorsByJ[n].modesInVec[EVecPos[i]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[i]].l)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end first if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v + 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end second if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l - ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v + (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end third if
                                                if (basisVectorsByJ[n].modesInVec[EVecPos[j]].v - 1 == basisVectorsByJ[m].modesInVec[EVecPos[j]].v && basisVectorsByJ[n].modesInVec[EVecPos[j]].l + ps == basisVectorsByJ[m].modesInVec[EVecPos[j]].l)
                                                {
                                                    temp = input.crossTermMatrix[EVecPos[i], EVecPos[j]] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[i]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[i]].l)) * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[EVecPos[j]].v - (double)ps * (double)basisVectorsByJ[n].modesInVec[EVecPos[j]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end fourth if
                                            }//end fourth if
                                        }
                                    }
                                    #endregion
                                }//end Cross-quadratic

                            }//end if CrossTerms = true
                            #endregion

                        }//column for loop                
                        r++;
                    }//row for loop
                }
            );//end parallel for

            foreach (Tuple<int, int, double> spot in matPos)
            {
                alglib.sparseadd(A, spot.Item1, spot.Item2, spot.Item3);
                alglib.sparseadd(A, spot.Item2, spot.Item1, spot.Item3);
            }

            //alglib.sparseconverttocrs(A);
            return A;
        }//end method genMatrix


        /// <summary>
        /// Generates an alglib sparsematrix object to be used in Lanczos diagonalization.
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
        /// <param name="Large">
        /// boolean value so that genMatrix method could be overloaded.
        /// </param>
        /// <returns></returns>
        public static alglib.sparsematrix genMatrix2(List<BasisFunction> basisVectorsByJ, bool isQuad, FileInfo input, out int nColumns, bool Large, int par)
        {
            int matSize = basisVectorsByJ.Count;
            nColumns = matSize;

            alglib.sparsematrix A = new alglib.sparsematrix();
            alglib.sparsecreate(matSize, matSize, 10, out A);
            bool containsAVecs = false;
            //bool bilinear = true;
            int nModes = input.nModes;
            List<int> AVecPos = new List<int>();
            List<int> EVecPos = new List<int>();
            List<int> biAVecPos = new List<int>();
            List<int> biEVecPos = new List<int>();
            ConcurrentBag<Tuple<int, int, double>> matPos = new ConcurrentBag<Tuple<int, int, double>>();
            int[] change = new int[3];

            #region CrossTermInitializationStuff
            if (input.includeCrossTerms == true)
            {

                for (int i = 0; i < input.nModes; i++)
                {
                    if (basisVectorsByJ[0].modesInVec[i].symmetryIsA == true)
                    {
                        AVecPos.Add(i);
                        containsAVecs = true;
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
                    //new lists just for bilinear coupling
                    //biAVecPos = AVecPos;
                    //biEVecPos = EVecPos;
                    for (int i = 0; i < AVecPos.Count; i++)
                    {
                        biAVecPos.Add(AVecPos[i]);
                    }
                    for (int i = 0; i < EVecPos.Count; i++)
                    {
                        biEVecPos.Add(EVecPos[i]);
                    }

                    //loop to eliminate any A vectors that have no cross-coupling
                    for (int i = 0; i < biAVecPos.Count; i++)
                    {
                        for (int j = 0; j < biEVecPos.Count; j++)
                        {
                            if (biAVecPos[i] > biEVecPos[j])
                            {
                                if (input.crossTermMatrix[biEVecPos[j], biAVecPos[i]] != 0)
                                {
                                    break;
                                }
                            }
                            else
                            {
                                if (input.crossTermMatrix[biAVecPos[i], biEVecPos[j]] != 0)
                                {
                                    break;
                                }
                            }
                            if (j == biEVecPos.Count - 1)
                            {
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
                                if (input.crossTermMatrix[biEVecPos[i], biAVecPos[j]] != 0)
                                {
                                    break;
                                }
                            }
                            else
                            {
                                if (input.crossTermMatrix[biAVecPos[j], biEVecPos[i]] != 0)
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
                    /*
                    if (biAVecPos.Count == 0 || biEVecPos.Count == 0)
                    {
                        bilinear = false;
                    }
                    */
                }
            }//end if CrossTerms == true
            #endregion

            //generates the Hamiltonian Matrix
            //just the HO terms
                for (int n = 0; n < matSize; n++)
                {
                    #region HO Terms
                    //one mode harmonic and anharmonic terms
                    for (int i = 0; i < input.nModes; i++)
                    {
                        int degeneracy = 2;
                        if (basisVectorsByJ[n].modesInVec[i].symmetryIsA)
                        {
                            degeneracy = 1;
                        }
                        double temp = basisVectorsByJ[n].modesInVec[i].modeOmega * (basisVectorsByJ[n].modesInVec[i].v + (double)degeneracy / 2D) - basisVectorsByJ[n].modesInVec[i].anharmonicity * Math.Pow((basisVectorsByJ[n].modesInVec[i].v + (double)degeneracy / 2), 2);//I deleted the (double) from the degeneracy / 2
                        alglib.sparseadd(A, n, n, temp);
                    }//end if AT == true
                    continue;
                    #endregion                
                }

                int[,] vlLambda = new int[matSize, nModes * 2 + 2];
                vlLambda.AsParallel();
                for (int i = 0; i < matSize; i++)
                {
                    for (int j = 0; j < nModes; j++)
                    {
                        vlLambda[i, j] = basisVectorsByJ[i].modesInVec[j].v;
                        vlLambda[i, j + nModes] = basisVectorsByJ[i].modesInVec[j].l;
                    }
                    vlLambda[i, nModes * 2] = basisVectorsByJ[i].Lambda;
                    vlLambda[i, nModes * 2 + 1] = (int) (basisVectorsByJ[i].J - 0.5M);
                }
                //indexes n and m are for the rows and columns of the matrix respectively

            var rangePartitioner = Partitioner.Create(0, matSize);
            ParallelOptions parOp = new ParallelOptions();
            parOp.MaxDegreeOfParallelism = par;
            Parallel.ForEach(rangePartitioner, parOp, (range, loopState) =>
            {
                int[] vdiff = new int[nModes];
                int[] ldiff = new int[nModes];
                for (int n = range.Item1; n < range.Item2; n++)
                {
                    for (int m = n + 1; m < matSize; m++)//changed from r + 1
                    {
                        double temp;
                        if(vlLambda[n, nModes *2] == vlLambda[m, nModes * 2])//Delta Lambda must be +/- 1
                        {
                            continue;
                        }                        
                        for (int b = 0; b < nModes; b++)
                        {
                            vdiff[b] = vlLambda[n, b] - vlLambda[m, b];
                            ldiff[b] = vlLambda[n, b + nModes] - vlLambda[m, b + nModes];
                        }

                        //Linear JT elements
                        if (vlLambda[n, nModes * 2 + 1] - vlLambda[m, nModes * 2 + 1] == 0)//means Delta J = 0, possible linear term
                        {
                            int[] vabs = new int[nModes];
                            int[] labs = new int[nModes];
                            for (int h = 0; h < nModes; h++)
                            {
                                vabs[h] = Math.Abs(vdiff[h]);
                                labs[h] = Math.Abs(ldiff[h]);
                            }
                            if (vabs.Sum() != 1)//Delta v = +/- 1 in only one mode
                            {
                                continue;
                            }
                            if (labs.Sum() != 1)//Delta l = +/- 1 in only one mode
                            {
                                continue;
                            }
                            int[] vlprod = new int[nModes];
                            for (int h = 0; h < nModes; h++)
                            {
                                vlprod[h] = vabs[h] * labs[h];
                            }
                            if (vlprod.Sum() != 1)
                            {
                                continue;
                            }
                            int lval = ldiff.Sum() * -1;
                            int mode = 0;
                            for (int h = 0; h < nModes; h++)
                            {
                                if (vdiff[h] == 0)
                                {
                                    continue;
                                }
                                else
                                {
                                    mode = h;
                                    break;
                                }
                            }
                            temp = basisVectorsByJ[n].modesInVec[mode].modeOmega * Math.Sqrt(basisVectorsByJ[n].modesInVec[mode].DBasis * ((double)basisVectorsByJ[n].modesInVec[mode].v + lval * (double)basisVectorsByJ[n].modesInVec[mode].l + 2D));
                            Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                            matPos.Add(tTemp);
                            continue;
                        }

                        if (Math.Abs(vlLambda[n, nModes * 2 + 1] - vlLambda[m, nModes * 2 + 1]) == 3)//means Delta J = 3, possible quadratic term
                        {
                            if (Math.Abs(ldiff.Sum()) != 2)//Delta l = 2 or -2
                            {
                                continue;
                            }
                            if (Math.Abs(vdiff.Sum()) != 2 && vdiff.Sum() != 0)//Delta v = 2, -2, or 0
                            {
                                continue;
                            }
                            int count = 0;
                            int pos = 0;
                            int count2 = 0;
                            int pos2 = 0;
                            int sign = 1;
                            for (int h = 0; h < nModes; h++)
                            {
                                if (vdiff[h] != 0)
                                {
                                    count++;
                                    pos = h;
                                }
                                if (ldiff[h] != 0)
                                {
                                    count2++;
                                    pos2 = h;
                                }
                            }
                            if (count > 1 || count2 > 1)//says only one mode changes for each
                            {
                                continue;
                            }
                            if (count == 1)//means Delta v = +/- 2
                            {
                                if (pos != pos2)//same mode has changes in v and l
                                {
                                    continue;
                                }
                                else
                                {
                                    if (vdiff.Sum() == 2)//top matrix elements on Ham page
                                    {
                                        if (ldiff.Sum() == 2)
                                        {
                                            sign = -1;
                                        }
                                        temp = basisVectorsByJ[n].modesInVec[pos].modeOmega * (basisVectorsByJ[n].modesInVec[pos].KBasis / 4D * Math.Sqrt((basisVectorsByJ[n].modesInVec[pos].v - sign * basisVectorsByJ[n].modesInVec[pos].l) * (basisVectorsByJ[n].modesInVec[pos].v - sign * basisVectorsByJ[n].modesInVec[pos].l - 2)));
                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                        matPos.Add(tTemp);
                                        continue;
                                    }
                                    else if (vdiff.Sum() == -2)//bottom matrix elements on Ham page
                                    {
                                        if (ldiff.Sum() == 2)
                                        {
                                            sign = -1;
                                        }
                                        temp = basisVectorsByJ[n].modesInVec[pos].modeOmega * (basisVectorsByJ[n].modesInVec[pos].KBasis / 4D * Math.Sqrt((basisVectorsByJ[n].modesInVec[pos].v + sign * basisVectorsByJ[n].modesInVec[pos].l + 4D) * (basisVectorsByJ[n].modesInVec[pos].v + sign * basisVectorsByJ[n].modesInVec[pos].l + 2)));
                                        Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                        matPos.Add(tTemp);
                                        continue;
                                    }
                                }
                            }
                            else//means Delta v = 0
                            {
                                if (ldiff.Sum() > 0)
                                {
                                    sign = -1;
                                }
                                temp = basisVectorsByJ[n].modesInVec[pos2].modeOmega * (basisVectorsByJ[n].modesInVec[pos2].KBasis / 2D * Math.Sqrt((basisVectorsByJ[n].modesInVec[pos2].v + sign * basisVectorsByJ[n].modesInVec[pos2].l + 2) * (basisVectorsByJ[n].modesInVec[pos2].v - sign * basisVectorsByJ[n].modesInVec[pos2].l)));
                                Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                matPos.Add(tTemp);
                                continue;
                            }
                        }

                        #region Cross Terms
                        /*
                        if (input.includeCrossTerms == true && numberOfChanges == 2)
                        {
                            //A-E interaction terms
                            if (containsAVecs == true)
                            {
                                if (bilinear == true)
                                {
                                    #region Bilinear
                                    //Delta J = 0
                                    //lambdas not equal
                                    //only v + 1 for both
                                    //make loop that checks the column row positions of vlLambda for correct changes in v and l
                                    //check v for v4, l for v4 and v for v1 
                                    for (int a = 0; a < biAVecPos.Count; a++)
                                    {
                                        for (int e = 0; e < biEVecPos.Count; e++)
                                        {
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
                                            if (input.crossTermMatrix[row, column] == 0D)
                                            {
                                                continue;
                                            }
                                            if (basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1 == basisVectorsByJ[m].modesInVec[biAVecPos[a]].v)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;//I think I can do this
                                                    //Math.Sqrt(basisVectorsByJ[n].modesInVec[biAVecPos[a]].modeOmega * basisVectorsByJ[n].modesInVec[biEVecPos[e]].modeOmega) * 
                                                    //this is what was removed from these matrix elements
                                                }//end first if
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end second if
                                                //
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end third if
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v + 1D) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end fourth if
                                            }
                                            if (basisVectorsByJ[n].modesInVec[biAVecPos[a]].v - 1 == basisVectorsByJ[m].modesInVec[biAVecPos[a]].v)
                                            {
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;//I think I can do this
                                                }//end first if
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + 2D));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end second if
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l - (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v + Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end third if
                                                if (basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - 1 == basisVectorsByJ[m].modesInVec[biEVecPos[e]].v && basisVectorsByJ[n].modesInVec[biEVecPos[e]].l + (int)Math.Pow(-1D, (double)input.S1) == basisVectorsByJ[m].modesInVec[biEVecPos[e]].l)
                                                {
                                                    temp = input.crossTermMatrix[row, column] * Math.Sqrt(((double)basisVectorsByJ[n].modesInVec[biAVecPos[a]].v) * ((double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].v - Math.Pow(-1D, (double)input.S1) * (double)basisVectorsByJ[n].modesInVec[biEVecPos[e]].l));
                                                    Tuple<int, int, double> tTemp = new Tuple<int, int, double>(n, m, temp);
                                                    matPos.Add(tTemp);
                                                    continue;
                                                }//end fourth if
                                            }
                                        }
                                    }
                                    #endregion
                                }//end bilinear if                                
                            }//end A-E interaction
                        }//end if CrossTerms = true
                        */
                        #endregion

                    }//column for loop
                }//row for loop
            }
            );//end parallel for

            foreach (Tuple<int, int, double> spot in matPos)
            {
                alglib.sparseadd(A, spot.Item1, spot.Item2, spot.Item3);
                alglib.sparseadd(A, spot.Item2, spot.Item1, spot.Item3);
            }
            return A;
        }//end method genMatrix


        /// <summary>
        /// Counts the number of Basis objects two JBasisVectors that are different from one another
        /// </summary>
        /// <param name="jList">
        /// jList is the List of JBasisVectors used as basis vectors in this J block.
        /// </param>
        /// <param name="nModes">
        /// Number of modes (Basis objects) in each JBasisVector
        /// </param>
        /// <param name="n">
        /// Row of matrix.
        /// </param>
        /// <param name="m">
        /// Column of matrix.
        /// </param>
        /// <returns>
        /// Integer value of number of Basis objects not equal to one another.
        /// </returns>
        private static int modeChangeCounter(List<BasisFunction> jList, int nModes, int n, int m)
        {
            int change = 0;
            for (int i = 0; i < nModes; i++)//first check to see if more than one mode changes, if so then continue outer for loop
            {
                if (jList[n].modesInVec[i].v == jList[m].modesInVec[i].v && jList[n].modesInVec[i].l == jList[m].modesInVec[i].l)// && basisVectorsByJ[n].Lambda == basisVectorsByJ[m].Lambda)
                {
                    continue;
                }
                else
                {                    
                    if ((int)Math.Abs(jList[n].modesInVec[i].v - jList[m].modesInVec[i].v) > 2 || (int)Math.Abs(jList[n].modesInVec[i].l - jList[m].modesInVec[i].l) > 2)// && basisVectorsByJ[n].Lambda == basisVectorsByJ[m].Lambda)
                    {
                        change = 10;
                        break;
                    }
                    change++;
                    if (change > 2)
                    {
                        break;
                    }
                    continue;
                }                
            }
            return change;
        }//end method modeChangeCounter

        private static bool modeChangeCounter(List<BasisFunction> jList, int nModes, int n, int m, out int[] change)
        {
            bool nzero = true;
            change = new int[3];
            for (int i = 0; i < nModes; i++)//first check to see if more than one mode changes, if so then continue outer for loop
            {
                change[1] += Math.Abs(jList[n].modesInVec[i].v - jList[m].modesInVec[i].v);
                change[2] += Math.Abs(jList[n].modesInVec[i].l - jList[m].modesInVec[i].l);
                if (jList[n].modesInVec[i].v != jList[m].modesInVec[i].v || jList[n].modesInVec[i].l != jList[m].modesInVec[i].l)// && basisVectorsByJ[n].Lambda == basisVectorsByJ[m].Lambda)
                {
                    change[0]++;
                }
            }
            if (change[0] > 2 || change[1] > 2 || change[2] > 2)
            {
                nzero = false;
            }
            return nzero;
        }//end method modeChangeCounter

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
        public static List<BasisFunction> sortByJ(List<BasisFunction> jVecs, decimal J)
        {
            List<BasisFunction> outList = new List<BasisFunction>();
            //count = 0;
            foreach (BasisFunction vector in jVecs)
            {
                if (vector.J == J)
                {
                    outList.Add(vector);
                    //count++;
                }//end if
            }//end foreach
            return outList;
        }//end method sortByJ
    }//end class GenHamMat
}
