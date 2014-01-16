using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Reflection;

namespace ConsoleApplication1
{
    class MainS
    {
        static void Main(string[] args)
        {
            try
            {
                //variable to store filepath
                string fileDirectory;

                //check to see if running on Mono for linux cluster
                //bool runningOnMono = Type.GetType("Mono.Runtime") != null;
                //commented this out and just made it the default method for Release operation

                //meaning it's running .NET on a Windows system -- UPDATE: meaning it's running in DEBUG mode
                //consider making this the default for Release jobs  -- done
                //if (!runningOnMono)
#if DEBUG
                //prompt user for input directory.  Default value is C:\SOCJT 2
                Console.WriteLine("Enter file directory or press enter to use C:\\SOCJT 2");
                fileDirectory = Console.ReadLine();
                if (fileDirectory == "")
                {
                    fileDirectory = "C:\\SOCJT 2";
                }

                //if entered directory doesn't exist provide option to create it.  If not, throw error and end execution.
                if (Directory.Exists(fileDirectory) == false)
                {
                    throw new DirectoryNotFoundException();
                }
#else

                fileDirectory = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);

#endif

                //set the directory for reading and writing files
                Directory.SetCurrentDirectory(fileDirectory);

                //prompt user to enter input file name
                Console.WriteLine("Enter file name including extension:");
                string inFileName = Console.ReadLine();

                //obtain output file name.  If blank, default value of input + .out is used
                Console.WriteLine("Enter output file name or press enter to use " + inFileName + ".out:");
                string outFile = Console.ReadLine();
                if (outFile == "" || outFile == " ")
                {
                    outFile = string.Concat(inFileName, ".out");
                }

                //start timer for overall program execution
                Stopwatch totalTime = new Stopwatch();
                totalTime.Start();

                //set input, output, and fit file values.
                string filepath = string.Concat(fileDirectory);
                //filepath = string.Concat("\\");
                filepath += "\\";
                string filepathIN = string.Copy(filepath);
                string filepathOUT = string.Copy(filepath);
                string filepathFIT = string.Copy(filepath);
                //filepathIN = string.Concat(inFileName);
                filepathIN += inFileName;
                //filepathOUT = string.Concat(outFile);
                filepathOUT += outFile;
                if (filepathOUT == filepathIN)
                {
                    throw new FileNameError("outFile");
                }

                //read input file
                string[] inputFile = FileInfo.fileRead(filepathIN);

                //create new FileInfo object
                FileInfo input = new FileInfo();

                //set input object data from input file
                input.setFileInfo(inputFile);
                input.filePath = filepath;

                //make the fitfile point to something
                filepathFIT = string.Concat(input.fitFile);

                //see if matFile is true, and if so if the matfile exists or not.
                if (input.useMatFile)
                {
                    //create file pointer for the matrix file
                    input.matFilePath = string.Copy(filepath);
                    input.matFilePath = string.Concat(input.matFile);
                    //check that the mat file has a valid path
                    if (input.matFilePath == filepathIN || input.matFilePath == filepathOUT || input.matFilePath == filepathFIT)
                    {
                        throw new FileNameError("matFile");
                    }
                    //if this file already exists, then use it for the matrix generation
                    if (File.Exists(input.matFilePath))
                    {
                        input.matMade = true;
                    }
                }

                //check that spin is integer or half integer only
                if (input.S % 0.5M != 0M)
                {
                    throw new SpinInvalidException();
                }

                //iterates through the Mode info in the input file and initializes the modes
                bool tempB = false;
                List<ModeInfo> Modes = new List<ModeInfo>();
                for (int i = 0; i < input.nModes; i++)
                {
                    Modes.Add(new ModeInfo(i, inputFile, out tempB));
                    if (tempB)
                    {
                        Modes[i].D = 0.5 * Math.Pow(Modes[i].kappa / Modes[i].modeOmega, 2D);
                        Modes[i].fitD = Modes[i].fitKappa;
                        Modes[i].K = Modes[i].eta / Modes[i].modeOmega;
                        Modes[i].fitK = Modes[i].fitEta;
                        input.useKappaEta = true;
                    }
                }//end for

                //sets bool isQuad = false if K is zero and fitK is false for all modes.
                bool isQuad = false;
                for (int i = 0; i < Modes.Count; i++)
                {
                    if (Modes[i].fitK == false && Modes[i].K == 0)
                    {
                        continue;
                    }//end if
                    else
                    {
                        isQuad = true;
                        break;
                    }
                }//end for

                //this sets isQuad = true if there is a cross quadratic term between an A and E mode
                if (isQuad == false)
                {
                    for (int i = 0; i < input.nModes; i++)
                    {
                        if (input.crossTermMatrix == null)
                        {
                            break;
                        }
                        for (int j = i + 1; j < input.nModes; j++)
                        {
                            if (input.crossTermMatrix[j, i] != 0)
                            {
                                if (Modes[i].IsAType == true)
                                {
                                    if (Modes[j].IsAType == false)
                                    {
                                        isQuad = true;
                                        break;
                                    }
                                }
                                else
                                {
                                    if (Modes[j].IsAType == true)
                                    {
                                        isQuad = true;
                                    }
                                }//end else
                            }//end if
                        }//end inner for
                    }//end outer for
                }//end if isQuad == false

                //Sets the fit boolean value to true if any values are to be fit
                bool fit = false;
                for (int i = 0; i < input.nModes; i++)
                {
                    if (input.fitAzeta == true)
                    {
                        fit = true;
                        break;
                    }
                    if (input.fitOrigin == true)
                    {
                        fit = true;
                        break;
                    }
                    if (Modes[i].fitOmega == true)
                    {
                        fit = true;
                        break;
                    }
                    if (Modes[i].fitD == true)
                    {
                        fit = true;
                        break;
                    }
                    if (Modes[i].fitK == true)
                    {
                        fit = true;
                        break;
                    }
                    if (Modes[i].fitWEXE == true)
                    {
                        fit = true;
                        break;
                    }
                }

                //this checks the cross term fit boolean values
                if (input.includeCrossTerms == true)
                {
                    foreach (bool fitter in input.crossTermFit)
                    {
                        if (fitter == true)
                        {
                            fit = true;
                            break;
                        }
                    }
                }

                //turns off the scan function if any values are being fit
                if (fit == true)
                {
                    input.Scan = false;
                }

                SOCJT runner = new SOCJT();
                FitSOCJT fitt = new FitSOCJT();

                //main subroutine execution when not running a scan
                if (input.Scan == false)
                {
                    //initialize variable to hold the output file and initialize it with input file containing original values
                    List<string> linesToWrite = new List<string>();
                    linesToWrite = OutputFile.inputFileMaker(input, Modes);

                    //if not fitting any variables run SOCJT routine directly
                    if (!fit)
                    {
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input));
                    }
                    else//else run FitSOCJT routine which will run LM optimizer which will call SOCJT
                    {
                        linesToWrite.AddRange(fitt.fit(Modes, isQuad, inputFile, input, filepathFIT));
                    }

                    //end stopwatch for total program time execution and appends the total run time to the output file
                    totalTime.Stop();
                    double TIME = totalTime.ElapsedMilliseconds / 1000D;
                    linesToWrite.Add(" ");
                    linesToWrite.Add("SOCJT 2 has completed. Total time elapsed = " + String.Format("{0,11:0.0000}", TIME) + " seconds.");
#if DEBUG
                    if (input.blockLanczos)
                    {
                        linesToWrite.Add("Orthog took " + Lanczos.reorthogTime / 1000L + " seconds");
                    }
#endif
                    //writes all info to the output file
                    File.WriteAllLines(filepathOUT, linesToWrite);
                    Console.WriteLine("SOCJT Has Completed");
                }//end no scan

                else//means a scan is being run
                {
                    #region SCAN
                    //create a list to store the final output eigenvalues from each iteration of SOCJT
                    List<Eigenvalue[]> scanList = new List<Eigenvalue[]>();

                    //create a List to hold the output file information
                    List<string> finalOut = new List<string>();

                    //for loop to run steps of scan
                    for (int h = 0; h < input.Steps; h++)
                    {
                        //basically run the SOCJT subroutine for each step of the scan and save the output file at the end of the loop.
                        List<string> linesToWrite = new List<string>();
                        //SOCJT runner = new SOCJT();

                        //run loop over each variable to be scanned and increment by step size times loop iteration.
                        for (int n = 0; n < input.scanList.Count; n++)
                        {
                            if (input.scanList[n].varToFit.ToUpper() == "OMEGA")
                            {
                                Modes[input.scanList[n].Mode - 1].modeOmega = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "WEXE")
                            {
                                Modes[input.scanList[n].Mode - 1].wExe = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "D")
                            {
                                Modes[input.scanList[n].Mode - 1].D = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "K")
                            {
                                Modes[input.scanList[n].Mode - 1].K = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "B")
                            {
                                int tRow = input.scanList[n].Mode;
                                int tCol = input.scanList[n].Cross;
                                int temp;
                                if (tRow > tCol)
                                {
                                    temp = tRow;
                                    tRow = tCol;
                                    tCol = temp;
                                }
                                input.crossTermMatrix[tRow - 1, tCol - 1] = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "SPECIAL")
                            {
                                input.crossTermMatrix[0, 0] = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "AT")
                            {
                                int tRow = input.scanList[n].Mode;
                                int tCol = input.scanList[n].Cross;
                                int temp;
                                if (tRow < tCol)
                                {
                                    temp = tRow;
                                    tRow = tCol;
                                    tCol = temp;
                                }
                                input.crossTermMatrix[tRow - 1, tCol - 1] = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                            if (input.scanList[n].varToFit.ToUpper() == "JT")
                            {
                                int tRow = input.scanList[n].Mode;
                                int tCol = input.scanList[n].Cross;
                                int temp;
                                if (tRow > tCol)
                                {
                                    temp = tRow;
                                    tRow = tCol;
                                    tCol = temp;
                                }
                                input.crossTermMatrix[tRow - 1, tCol - 1] = input.scanList[n].Start + input.scanList[n].Step * (double)h;
                            }
                        }//end adjust each variable to fit

                        //now make output file for each step
                        linesToWrite = OutputFile.inputFileMaker(input, Modes);

                        //run SOCJT subroutine on this step's values
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input));

                        //write output file for this step
                        string stepFile = filepathOUT + "_step_" + Convert.ToString(h + 1) + ".out";
                        File.WriteAllLines(stepFile, linesToWrite);

                        //add this steps eigenvalues to the scan output file
                        scanList.Add(runner.finalList);

                        //clear values for the next iteration
                        linesToWrite = null;
                        stepFile = null;
                        runner = null;
                    }//end of steps loop

                    //write scan output file including total elapsed time
                    List<string> scanOut = OutputFile.ScanOutput(scanList);
                    totalTime.Stop();
                    double TIME = totalTime.ElapsedMilliseconds / 1000D;
                    scanOut.Add(" ");
                    scanOut.Add("SOCJT 2 has completed. Total time elapsed = " + String.Format("{0,11:0.0000}", TIME) + " seconds.");
                    File.WriteAllLines(filepathOUT, scanOut);
                    #endregion
                }//end else for scan

                //code to write matrix file to disc here?
                if (input.useMatFile && !input.matMade)
                {
                    OutputFile.writeMatFile(input);
                }//end if to write matrix to file

                //Here, if necessary, the eigenvectors for very large matrices are calculated after the rest of the calculations have been done.
                //This is done because this process may take a large amount of time
                if (runner.lanczosEVectors != null || fitt.lanczosEVectors != null)
                {
                    Console.WriteLine("\n Eigenvectors being calculated...");
                    Console.WriteLine("This may take some time.");
                    if (fit)
                    {
                        eVecGenerator(fitt.lanczosEVectors, filepath, fitt.basisSet, input);
                    }//end if
                    else
                    {
                        eVecGenerator(runner.lanczosEVectors, filepath, runner.basisSet, input);
                    }//end else
                }//end if

            }//end try block

            //Exception handling
            catch (DirectoryNotFoundException)
            {
                Console.WriteLine("The directory does not exist");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (FileNotFoundException a)
            {
                Console.WriteLine(a.Message);
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
#if !DEBUG
            catch (IndexOutOfRangeException)
            {
                Console.WriteLine("An index out of range exception has occurred");
                Console.WriteLine("Check to make sure that NMODES is correct and");
                Console.WriteLine("that any CROSS_TERM arguments are correct.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (BasisSetTooSmallException a)
            {
                Console.Write(a.eMessage);
                Console.ReadLine();
            }
            catch (SpinInvalidException)
            {
                Console.WriteLine("S has an incorrect value.  S must be either integer or half integer.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (AEAnharmonicTermException)
            {
                Console.WriteLine("Cross anharmonic terms can only be included between two A modes.");
                Console.WriteLine("Take out any cross anharmonic terms between A and E modes.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (RepeaterError)
            {
                Console.WriteLine("NaN Error in Naive Lanczos routine.");
                Console.WriteLine("Try decreasing NOITS or using Block Lanczos instead.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (FileNameError err)
            {
                Console.WriteLine(err.eMessage);
                Console.ReadLine();
            }
            catch (MatrixFileError)
            {
                Console.WriteLine("The basis sets of the input file and matrix file are different sizes.");
                Console.WriteLine("Try regenerating the matrix.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }

            catch (Exception ex)
            {
                Console.WriteLine("An exception has occurred: " + ex.Message);
                Console.WriteLine("Check your input file and try again.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
#endif
        }//end Main

        /// <summary>
        /// Reads a text file into memory and parses it into a string array.
        /// </summary>
        /// <param name="filepath">
        /// Filepath of file to be read and parsed.
        /// </param>
        /// <returns>
        /// Array containing all values from text file.
        /// </returns>
        public static double[] vecRead(string filepath, int n)
        {
            List<double> vector = new List<double>();
            //string[] inputFa = { };
            using (StreamReader vecIn = new StreamReader(filepath))
            {
                string lineS;
                //string[] SOCJTNewLine;
                //char[] delimiters = new char[] { '\t', '\r', '=', ' ' };
                while ((lineS = vecIn.ReadLine()) != ("START_VEC" + n))
                {
                    continue;
                }
                while ((lineS = vecIn.ReadLine()) != "END_VEC")
                {
                    vector.Add(FileInfo.parseDouble(lineS));
                }//end while
            }//end StreamReader
            return vector.ToArray();
        }//end method vecRead


        private static void eVecGenerator(List<double[,]> lanczosEVectors, string filepath, List<List<BasisFunction>> basisSet, FileInfo input)
        {
            List<double[,]> eVecs = new List<double[,]>();
            string file;
            double evMin = input.evMin;
            double[] lanczosVector;
            for (int i = 0; i < lanczosEVectors.Count; i++)
            {
                //this uses the same convention as in the Lanczos method which saves the vectors. Not as general as I'd like but probably ok.
                file = filepath + "temp_vecs_" + i + ".tmp";
                lanczosVector = vecRead(file, 0);
                int numberOfEigenvalues = lanczosEVectors[i].GetLength(1);
                int iterations = lanczosEVectors[i].GetLength(0);
                int dimension = lanczosVector.Length;
                eVecs[i] = new double[dimension, numberOfEigenvalues];
                for (int j = 0; j < iterations; j++)
                {
                    lanczosVector = vecRead(file, j);
                    for (int m = 0; m < numberOfEigenvalues; m++)
                    {
                        //read in the right vector to memory from the lanczos vector file                    
                        for (int n = 0; n < lanczosVector.Length; n++)
                        {                            
                            eVecs[i][n, m] += lanczosVector[n] * lanczosEVectors[i][j, m];
                        }//end loop over rows of lanczos Vector
                    }//end loop over columns of lanczosEVectors
                }
                //here delete the temp lanczos vector file
                File.Delete(file);
            }//end loop to generate eigenvectors
            //here use basis sets and eigenvectors and write them to file
            StringBuilder output = new StringBuilder();
            for (int i = 0; i < lanczosEVectors.Count; i++)
            {
                decimal jblock = (decimal)i + 0.5M;
                output.AppendLine("J-Block " + jblock);
                output.AppendLine(" ");
                for (int j = 0; j < lanczosEVectors[i].GetLength(0); j++)
                {
                    output.AppendLine(" " + "\r");
                    output.AppendLine("Eigenvalue" + "\t" + Convert.ToString(j + 1));
                    output.AppendLine(" " + "\r");
                    output.AppendLine("Eigenvector: (Only vectors with coefficients larger than " + Convert.ToString(evMin) + " are shown)");
                    output.AppendLine(" ");

                    bool a1 = SOCJT.isA(basisSet[i], lanczosEVectors[i], j, input);
                    if (a1)
                    {
                        output.AppendLine("Vector is Type 1");
                    }
                    else
                    {
                        output.AppendLine("Vector is Type 2");
                    }
                    output.AppendLine(" ");

                    output.Append("Coefficient" + "\t");
                    for (int h = 0; h < input.nModes; h++)
                    {
                        output.Append("v(" + Convert.ToString(h + 1) + ")" + "\t" + "l(" + Convert.ToString(h + 1) + ")" + "\t");
                    }
                    output.Append("lambda");
                    for (int h = 0; h < basisSet[i].Count; h++)//goes through basis vectors
                    {
                        if (lanczosEVectors[i][h, j] > evMin || lanczosEVectors[i][h, j] < -1.0 * input.evMin)
                        {
                            output.AppendLine("\t");
                            output.Append(String.Format("{0,10:0.000000}", lanczosEVectors[i][h, j]));
                            for (int m = 0; m < input.nModes; m++)//goes through each mode
                            {
                                output.Append("\t" + "  " + Convert.ToString(basisSet[i][h].modesInVec[m].v) + "\t" + String.Format("{0,3}", basisSet[i][h].modesInVec[m].l));//  "  " + Convert.ToString(jBasisVecsByJ[i][h].modesInVec[m].l));
                            }
                            output.Append("\t" + String.Format("{0,4}", basisSet[i][h].Lambda));
                        }
                    }
                    output.AppendLine("\r");
                }//end j for loop
            }//end i loop
            //code here to write to file
            string fileName = filepath + input.title + "_EVecs.out";
            List<string> ou = new List<string>();
            ou.Add(output.ToString());
            File.WriteAllLines(fileName, ou);
        }//end eVecGenerator
    }//end class Program
}//end namespace
