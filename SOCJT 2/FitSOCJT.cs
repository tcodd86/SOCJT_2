﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace ConsoleApplication1
{
    class FitSOCJT
    {
        static int count = 0;

        public static List<string> fit(List<ModeInfo> Modes, bool isQuad, string[] inputFile, FileInfo input, String filepath)
        {
            //string to return
            List<string> output = new List<string>();

            //reads and parses fit file
            string[] fitF = {};
            try
            {
                using (StreamReader FitIn = new StreamReader(filepath))
                {
                    string lineS;
                    string[] FitNewLine;
                    while ((lineS = FitIn.ReadLine()) != null)
                    {
                        char[] delimiters = new char[] { '\t', '\r', '=', ' ' };
                        FitNewLine = lineS.Split(delimiters, StringSplitOptions.RemoveEmptyEntries);
                        fitF = fitF.Concat(FitNewLine).ToArray();
                    }
                }
            }
            catch(FileNotFoundException)
            {
                throw new FitFileNotFoundException();
            }
            int nToFit = Convert.ToInt16(fitF[0]);
            Eigenvalue[] userInput = new Eigenvalue[nToFit];
            for (int i = 0; i < nToFit; i++)
            {
                userInput[i] = new Eigenvalue(Convert.ToDecimal(fitF[i * 4 + 2]), Convert.ToInt16(fitF[i * 4 + 3]), Convert.ToDecimal(fitF[i * 4 + 4]), Convert.ToDouble(fitF[i * 4 + 1]));
            }
            
            //initializes the X vector, boundary conditions and variable scales
            List<double> xList = new List<double>();
            List<double> bndL = new List<double>();
            List<double> bndU = new List<double>();
            List<double> lScale = new List<double>();

            if (input.fitAzeta == true)
            {
                xList.Add(input.Azeta);
                bndL.Add(double.NegativeInfinity);
                bndU.Add(double.PositiveInfinity);
                lScale.Add(1.0);
            }
            if (input.fitOrigin == true)
            {
                xList.Add(0.0);
                bndL.Add(double.NegativeInfinity);
                bndU.Add(double.PositiveInfinity);
                lScale.Add(200.0);//changed from 75
            }
            for (int i = 0; i < input.nModes; i++)
            {                
                if (Modes[i].fitOmega == true)
                {
                    xList.Add(Modes[i].modeOmega);
                    bndL.Add(0.0);
                    bndU.Add(double.PositiveInfinity);
                    lScale.Add(200.0);//changed from 200
                }
                if (Modes[i].fitD == true)
                {
                    xList.Add(Modes[i].D);
                    bndL.Add(0.0);
                    bndU.Add(double.PositiveInfinity);
                    lScale.Add(10.0);//changed from 1000
                }
                if (Modes[i].fitK == true)
                {
                    xList.Add(Modes[i].K);
                    bndL.Add(double.NegativeInfinity);
                    bndU.Add(double.PositiveInfinity);
                    lScale.Add(1.0);//changed from 100
                }
                if (Modes[i].fitWEXE == true)
                {
                    xList.Add(Modes[i].wExe);
                    bndL.Add(double.NegativeInfinity);
                    bndU.Add(double.PositiveInfinity);
                    lScale.Add(50.0);//change from 250
                }
            }
            if (input.includeCrossTerms == true)
            {
                for (int i = 0; i < input.nModes; i++)
                {
                    for (int j = 0; j < input.nModes; j++)
                    {
                        if (input.crossTermFit[i, j] == true)
                        {
                            xList.Add(input.crossTermMatrix[i, j]);
                            bndL.Add(double.NegativeInfinity);
                            bndU.Add(double.PositiveInfinity);
                            if (j > i)
                            {
                                lScale.Add(100.0);//scale for bilinear coupling + cross quadratic, changed from 1
                            }
                            if (j == i)
                            {
                                lScale.Add(50.0);//change from 250
                            }
                            else
                            {
                                lScale.Add(50.0);//scale for cross anharmonic, change from 250
                            }
                        }
                    }
                }
            }

            //xVec for minLM optimizer
            double[] xVec = xList.ToArray();

            //Upper/lower boundary arrays for minLM optimizer
            double[] lowBound = bndL.ToArray();
            double[] upBound = bndU.ToArray();

            //scaling array for minLM optimizer
            double[] scale = lScale.ToArray();

            //Generate new Master Object
            SOCJT run = new SOCJT();            
            MasterObject Masterly = new MasterObject();            
            Masterly.Initialize(run, input, Modes, inputFile, isQuad, userInput);

            //initialize and run MinLM algorithm
            double epsg = input.gTol;
            double epsf = input.fTol;
            double epsx = input.xTol;
            int maxits = input.maxFev;
            alglib.minlmstate state;
            alglib.minlmreport rep;

            alglib.minlmcreatev(userInput.Length, xVec, input.factor, out state);
            alglib.minlmsetbc(state, lowBound, upBound);
            alglib.minlmsetscale(state, scale);
            alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
            alglib.minlmsetxrep(state, true);
            //acc = 1, meaning it's on by default.  Just leave that.
            //alglib.minlmsetacctype(state, 0);

            alglib.minlmoptimize(state, function, null, Masterly);
            alglib.minlmresults(state, out xVec, out rep);
            int report = rep.terminationtype;
            int iter = rep.iterationscount;

            //this calculates the covariance matrix from the Jacobian using cov[,] = {J^T * J}^-1
            //initialize covariance matrix
            double[,] resultM = new double[state.x.Length, state.x.Length];
            for (int u = 0; u < state.x.Length; u++)
            {
                for (int uu = 0; uu < state.x.Length; uu++)
                {
                    resultM[u, uu] = 0.0;
                }
            }
            //calculate J^T * J and store in resultM
            alglib.rmatrixgemm(state.j.GetLength(1), state.j.GetLength(1), state.j.GetLength(0), 1.0, state.j, 0, 0, 1, state.j, 0, 0, 0, 0.0, ref resultM, 0, 0);
            //take inverse of resultM, replaces resultM
            int invInfo;
            alglib.matinvreport invReport = new alglib.matinvreport();
            alglib.rmatrixinverse(ref resultM, out invInfo, out invReport);
            //now make correlation coefficient matrix
            double[,] corMat = new double[resultM.GetLength(0), resultM.GetLength(0)];
            for (int i = 0; i < resultM.GetLength(0); i++)
            {
                for (int j = 0; j < resultM.GetLength(0); j++)
                {
                    corMat[i, j] = resultM[i, j] / (Math.Sqrt(resultM[i, i] * resultM[j, j]));
                }
            }

            //make output
            output = Masterly.nSoc.outp;
            //add something showing RMS error and parameters for each JT mode
            double[] error = ComparerVec(userInput, Masterly.nSoc.finalList, Masterly.nInput.origin, true);

            StringBuilder file = new StringBuilder();
            file.AppendLine("Fit report below...");
            file.AppendLine("Termination Type: ");
            if (report == -9)
            {
                file.Append("Derivative correctness check failed");
            }
            if (report == 1)
            {
                file.Append("Relative function improvement is no more than EpsF");
            }
            if (report == 2)
            {
                file.Append("Relative step is no more than EpsX");
            }
            if (report == 4)
            {
                file.Append("Gradient is no more than EpsG");
            }
            if (report == 5)
            {
                file.Append("Maximum number of iterations was exceeded");
            }
            if (report == 7)
            {
                file.Append("Stopping conditions are too stringent, further improvement is impossible");
            }

            file.AppendLine(" ");

            file.AppendLine("Number of iterations: " + Convert.ToString(iter));
            file.AppendLine("Number of times the eigenvalues were calculated: " + Convert.ToString(rep.nfunc));

            file.AppendLine(" ");

            file.AppendLine("A * zeta e = " + Convert.ToString(Masterly.nInput.Azeta));
            file.AppendLine("Final Parameters for Each Mode:");
            file.AppendLine("Mode #" + "\t" + "V(min)" + "\t" + "V(max)" + "\t" + "Omega(E)" + "\t" + "wexe" + "\t" + "D" + "\t" + "K" + "\t" + "JTSE" + "\t" + "Omega(A)" + "\t" + "A Type?");
            for (int i = 0; i < Masterly.nInput.nModes; i++)
            {
                double JTSE;
                if (Masterly.nModes[i].IsAType == true)
                {
                    JTSE = 0.0;
                }
                else
                {
                    JTSE = Masterly.nModes[i].D * Masterly.nModes[i].modeOmega * (Masterly.nModes[i].K + 1.0);
                }
                file.AppendLine(Convert.ToString(i + 1) + "\t" + String.Format("{0,4:0}", 0) + "\t" + String.Format("{0,4:0}", Masterly.nModes[i].modeVMax) + "\t" + String.Format("{0,7:0.00}", Masterly.nModes[i].modeOmega) + "\t" + "\t" + String.Format("{0,4:0.00}", Masterly.nModes[i].wExe) + "\t" + String.Format("{0,5:0.0000}", Masterly.nModes[i].D) + "\t" + String.Format("{0,5:0.0000}", Masterly.nModes[i].K) + "\t" + String.Format("{0,4:0.00}", JTSE) + "\t" + String.Format("{0,7:0.00}", Masterly.nModes[i].modeAOmega) + "\t" + "\t" + Convert.ToString(Masterly.nModes[i].IsAType));
            }
            file.AppendLine("  ");
            if (input.fitOrigin == true)
            { 
                file.AppendLine("Origin Shift = " + String.Format("{0,8:0.000}", -1.0 * Masterly.nInput.origin));
            }
            file.AppendLine("  ");
            if (input.includeCrossTerms == true)
            {
                for (int i = 0; i < input.nModes; i++)
                {
                    for (int j = 0; j < input.nModes; j++)
                    {
                        if (input.crossTermMatrix[i, j] != 0.0 || input.crossTermFit[i, j] == true)
                        {
                            if (i < j)
                            {
                                file.AppendLine("Mode " + Convert.ToString(i + 1) + " Mode " + Convert.ToString(j + 1) + " JT cross-term = " + String.Format("{0,10:0.0000}", input.crossTermMatrix[i, j]));
                            }
                            else
                            {
                                file.AppendLine("Mode " + Convert.ToString(j + 1) + " Mode " + Convert.ToString(i + 1) + " AT cross-term = " + String.Format("{0,10:0.00}", input.crossTermMatrix[i, j]));
                            }
                        }
                    }
                }
                file.AppendLine(" ");
            }            

            file.AppendLine("Fitting Results:");
            if (input.fitOrigin == true)
            {
                file.AppendLine("Experimental values are shifted by " + String.Format("{0,8:0.000}", Masterly.nInput.origin) + " wavenumbers.");
            }
            file.AppendLine("FitFile Value" + "\t" + "Calculated Value" + "\t" + "Exp - Calc" + "\t" + "(Exp - Calc)^2");
            for (int i = 0; i < error.Length; i++)
            {
                file.AppendLine(String.Format("{0,13:0.000}", userInput[i].Ev + Masterly.nInput.origin) + "\t" + String.Format("{0,13:0.000}", userInput[i].Ev + Masterly.nInput.origin - error[i]) + "\t" + "\t" + String.Format("{0,9:0.000}", error[i]) + "\t" + String.Format("{0,11:0.000}", Math.Pow(error[i], 2D)));
            }
            file.AppendLine("  ");
            file.AppendLine("RMS Error = " + String.Format("{0,10:0.000}", (Math.Sqrt(FitSOCJT.Comparer(userInput, Masterly.nSoc.finalList, Masterly.nInput.origin) / userInput.Length))));
            file.AppendLine("  ");

            int l = 0;
            if (Masterly.nInput.fitAzeta == true)
            {
                file.AppendLine("Azeta StdDev = " + String.Format("{0,10:0.000000}", Math.Sqrt(resultM[l, l])));
                l++;
            }
            if (Masterly.nInput.fitOrigin == true)
            {
                file.AppendLine("Origin StdDev = " + String.Format("{0,10:0.00}", Math.Sqrt(resultM[l, l])));
                l++;
            }
            for (int i = 0; i < Masterly.nInput.nModes; i++)
            {
                if (Masterly.nModes[i].fitOmega == true)
                {
                    file.AppendLine("Mode " + Convert.ToString(i + 1) + " Omega StdDev = " + String.Format("{0,10:0.00}", Math.Sqrt(resultM[l, l])));
                    l++;
                }
                if (Masterly.nModes[i].fitD == true)
                {
                    file.AppendLine("Mode " + Convert.ToString(i + 1) + " D StdDev = " + String.Format("{0,10:0.0000}", Math.Sqrt(resultM[l, l])));
                    l++;
                }
                if (Masterly.nModes[i].fitK == true)
                {
                    file.AppendLine("Mode " + Convert.ToString(i + 1) + " K StdDev = " + String.Format("{0,10:0.0000}", Math.Sqrt(resultM[l, l])));
                    l++;
                }
                if (Masterly.nModes[i].fitWEXE == true)
                {
                    file.AppendLine("Mode " + Convert.ToString(i + 1) + " wexe StdDev = " + String.Format("{0,10:0.00}", Math.Sqrt(resultM[l, l])));
                    l++;
                }
            }
            if (Masterly.nInput.includeCrossTerms == true)
            {
                for (int i = 0; i < Masterly.nInput.nModes; i++)
                {
                    for (int h = 0; h < Masterly.nInput.nModes; h++)
                    {
                        if (Masterly.nInput.crossTermFit[i, h] == true)
                        {
                            if (i < h)
                            {
                                file.AppendLine("Mode " + Convert.ToString(i + 1) + " Mode " + Convert.ToString(h + 1) + " JT Term StdDev = " + String.Format("{0,10:0.0000}", Math.Sqrt(resultM[l, l])));
                                l++;
                            }
                            else
                            {
                                file.AppendLine("Mode " + Convert.ToString(h + 1) + " Mode " + Convert.ToString(i + 1) + " AT Term StdDev = " + String.Format("{0,10:0.00}", Math.Sqrt(resultM[l, l])));
                                l++;
                            }
                        }
                    }
                }
            }
            file.AppendLine(" ");
            file.AppendLine("Correlation coefficient matrix:");
            for (int i = 0; i < resultM.GetLength(0); i++)
            {
                file.AppendLine("");
                file.Append("(");
                for (int j = 0; j < resultM.GetLength(0); j++)
                {
                    if (j != 0)
                    {
                        file.Append(",  ");
                    }
                    file.Append(String.Format("{0,8:0.000}", corMat[i, j]));
                }
                file.Append(")");
            }
            file.AppendLine(" ");
            file.AppendLine(" ");

            output.Add(file.ToString());

            output.AddRange(OutputFile.inputFileMaker(Masterly.nInput, Masterly.nModes));
            return output;
        }

        /// <summary>
        /// Function called by MinLM algorithm.
        /// </summary>
        /// <param name="x">
        /// Contains values of paramters being fit.  They are assigned within the function.
        /// </param>
        /// <param name="fi">
        /// The vector containing the errors for a given set of parameters.
        /// </param>
        /// <param name="obj">
        /// This passes an instance of MasterObject class.  Object is recast within function.
        /// </param>
        public static void function(double[] x, double[] fi, Object obj)//, MasterObject Master)
        {
            //MasterObject Master = ConsoleApplication1.FitSOCJT.Masterly;
            MasterObject Master = (MasterObject) obj;//(ConsoleApplication1.MasterObject)obj;
            int j = 0;
            double[] temp;
            if (Master.nInput.fitAzeta == true)
            {
                Master.nInput.Azeta = x[j];
                j++;
            }
            if (Master.nInput.fitOrigin == true)
            {
                Master.nInput.origin = x[j];
                j++;
            }
            for (int i = 0; i < Master.nInput.nModes; i++)
            {
                if (Master.nModes[i].fitOmega == true)
                {
                    Master.nModes[i].modeOmega = x[j];
                    j++;
                }
                if (Master.nModes[i].fitD == true)
                {
                    Master.nModes[i].D = x[j];
                    j++;
                }
                if (Master.nModes[i].fitK == true)
                {
                    Master.nModes[i].K = x[j];
                    j++;
                }
                if (Master.nModes[i].fitWEXE == true)
                {
                    Master.nModes[i].wExe = x[j];
                    j++;
                }
            }
            if (Master.nInput.includeCrossTerms == true)
            {
                for (int i = 0; i < Master.nInput.nModes; i++)
                {
                    for (int h = 0; h < Master.nInput.nModes; h++)
                    {
                        if (Master.nInput.crossTermFit[i, h] == true)
                        {
                            Master.nInput.crossTermMatrix[i, h] = x[j];
                            j++;
                        }
                    }
                }
            }
            Master.nSoc.SOCJTroutine(Master.nModes, Master.nIsQuad, Master.nInputFile, Master.nInput);
            temp = FitSOCJT.ComparerVec(Master.nFitFile, Master.nSoc.finalList, Master.nInput.origin, true);
            for (int i = 0; i < temp.Length; i++)
            {
                fi[i] = temp[i];
            }
            double rmsError = 0.0;
            for (int w = 0; w < temp.Length; w++)
            {
                rmsError += Math.Pow(temp[w], 2D);
            }
            rmsError /= temp.Length;

            Console.WriteLine("Eigenvalues have been calculated {0} times.", count++);
            Console.WriteLine("RMS Error is {0}.", rmsError);
        }

        /// <summary>
        /// Returns a vector containing the error between experimental and calculated eigenvalues.
        /// </summary>
        /// <param name="exp">
        /// Experimental eigenvalues.
        /// </param>
        /// <param name="socjtOut">
        /// Calculated eigenvalues.
        /// </param>
        /// <param name="raw">
        /// If true, returns raw (exp - calc) eigenvalues.  If false, returns (exp - calc)^2.
        /// </param>
        /// <returns>
        /// double[] containing error values.
        /// </returns>
        public static double[] ComparerVec(Eigenvalue[] exp, Eigenvalue[] socjtOut, double origin, bool raw)
        {
            List<double> errors = new List<double>();
            double[] errorsvec;
            //double ZPE = socjtOut[0].Ev;
            for (int i = 0; i < exp.Length; i++)
            {
                for (int j = 0; j < socjtOut.Length; j++)
                {
                    if (exp[i].pJ == socjtOut[j].pJ)
                    {
                        if (exp[i].nJ == socjtOut[j].nJ)
                        {
                            if (exp[i].Sig == socjtOut[j].Sig)
                            {
                                if (raw == false)
                                {
                                    errors.Add(Math.Pow((exp[i].Ev + origin) - socjtOut[j].Ev, 2D));//took out ZPE
                                    break;
                                }
                                if (raw == true)
                                {
                                    errors.Add((exp[i].Ev + origin) - socjtOut[j].Ev);//took out ZPE
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            errorsvec = errors.ToArray();
            return errorsvec;
        }

        public static double Comparer(Eigenvalue[] exp, Eigenvalue[] socjtOut, double origin)
        {
            double RMS = 0D;
            //int h = 0;
            for (int i = 0; i < exp.Length; i++)
            {
                for (int j = 0; j < socjtOut.Length; j++)
                {
                    if (exp[i].pJ == socjtOut[j].pJ)
                    {
                        if (exp[i].nJ == socjtOut[j].nJ)
                        {
                            RMS += Math.Pow((exp[i].Ev + origin) - socjtOut[j].Ev, 2D);
                            //h = j;
                            break;
                        }
                    }
                }
            }
            return RMS;
        }
    }
}
