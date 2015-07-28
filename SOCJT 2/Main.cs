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
                string fileDirectory;
#if DEBUG
                //prompt user for input directory.  Default value is C:\SOCJT 2
                Console.WriteLine("Enter file directory or press enter to use C:\\SOCJT 2");
                fileDirectory = Console.ReadLine();
                if (fileDirectory == "")
                {
                    fileDirectory = "C:\\SOCJT 2";
                }

                if (Directory.Exists(fileDirectory) == false)
                {
                    throw new DirectoryNotFoundException();
                }
#else

                fileDirectory = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);

#endif

                Directory.SetCurrentDirectory(fileDirectory);

                /* string inFileName = args[0]; // Used for NFG SOCJT file. I don't know why the commented section below doesn't work.
                string outFile = args[1];

                //if (args == null || args.Length == 0)
                //{
                //    Console.WriteLine("Enter file name including extension:");
                //    inFileName = Console.ReadLine();

                //    //obtain output file name.  If blank, default value of input + .out is used
                //    Console.WriteLine("Enter output file name or press enter to use " + inFileName + ".out:");
                //    outFile = Console.ReadLine();
                //    if (outFile == "" || outFile == " ")
                //    {
                //        outFile = string.Concat(inFileName, ".out");
                //    }
                //}
                */

                //Original Input
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
                filepath += "\\";
                string filepathIN = string.Copy(filepath);
                string filepathOUT = string.Copy(filepath);
                string filepathFIT = string.Copy(filepath);
                filepathIN += inFileName;
                filepathOUT += outFile;
                if (filepathOUT == filepathIN)
                {
                    throw new FileNameError("outFile");
                }

                //read and parse input file, then initialize FileInfo object
                string[] inputFile = FileInfo.FileRead(filepathIN);
                FileInfo input = new FileInfo();
                input.SetFileInfo(inputFile, filepath);
                input.FilePath = filepath;

                //make the fitfile point to something
                filepathFIT = string.Concat(input.FitFile);

                //see if matFile is true, and if so if the matfile exists or not.
                if (input.UseMatrixFile)
                {
                    input.MatrixFilePath = string.Copy(filepath);
                    input.MatrixFilePath += input.MatrixFile;
                    //check that the mat file has a valid path
                    if (input.MatrixFilePath == filepathIN || input.MatrixFilePath == filepathOUT || input.MatrixFilePath == filepathFIT)
                    {
                        throw new FileNameError("matFile");
                    }
                    //if this file already exists, then use it for the matrix generation
                    if (File.Exists(input.MatrixFilePath))
                    {
                        input.MatrixMade = true;
                    }
                }

                //initialize the modes
                List<ModeInfo> Modes = ModeInitialization(inputFile, input);

                //Determines if the quadratic basis set should be used
                bool isQuad = IsQuad(input, Modes);

                //Determines if any values are being fit
                bool fit = IsFit(input, Modes);

                //main subroutine execution when not running a scan
                if (input.Scan == false)
                {
                    List<string> linesToWrite = new List<string>();
                    linesToWrite = OutputFile.inputFileMaker(input, Modes);

                    if (!fit)
                    {
                        SOCJT runner = new SOCJT();
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input, input.useAbsoluteEV));
                        WriteOutputFile(totalTime, filepathOUT, linesToWrite);
                        //Here, if necessary, the eigenvectors for very large matrices are calculated after the rest of the calculations have been done.
                        //This is done because this process may take a large amount of time
                        if (runner.lanczosEVectors != null)
                        {
                            CalcLargeEvecs(filepath, input, runner.lanczosEVectors, runner.basisSet);
                        }
                    }
                    else
                    {
                        FitSOCJT fitt = new FitSOCJT();
                        linesToWrite.AddRange(fitt.fit(Modes, isQuad, inputFile, input, filepathFIT));
                        WriteOutputFile(totalTime, filepathOUT, linesToWrite);
                        if (fitt.lanczosEVectors != null)
                        {
                            CalcLargeEvecs(filepath, input, fitt.lanczosEVectors, fitt.basisSet);
                        }
                    }
                    //WriteOutputFile(totalTime, filepathOUT, linesToWrite);
                }//end no scan

                else//means a scan is being run
                {
                    #region SCAN
                    SOCJT runner = new SOCJT();
                    List<Eigenvalue[]> scanList = new List<Eigenvalue[]>();
                    List<string> finalOut = new List<string>();

                    for (int h = 0; h < input.ScanSteps; h++)
                    {
                        List<string> linesToWrite = new List<string>();

                        //run loop over each variable to be scanned and increment by step size times loop iteration.
                        for (int n = 0; n < input.ScanList.Count; n++)
                        {
                            if (input.ScanList[n].varToFit.ToUpper() == "OMEGA")
                            {
                                Modes[input.ScanList[n].Mode - 1].modeOmega = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "WEXE")
                            {
                                Modes[input.ScanList[n].Mode - 1].wExe = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "D")
                            {
                                Modes[input.ScanList[n].Mode - 1].D = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "K")
                            {
                                Modes[input.ScanList[n].Mode - 1].K = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "B" || input.ScanList[n].varToFit.ToUpper() == "C")
                            {
                                int tRow = input.ScanList[n].Mode;
                                int tCol = input.ScanList[n].Cross;
                                int temp;
                                if (tRow > tCol)
                                {
                                    temp = tRow;
                                    tRow = tCol;
                                    tCol = temp;
                                }
                                input.CrossTermMatrix[tRow - 1, tCol - 1] = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "SPECIAL")
                            {
                                input.CrossTermMatrix[0, 0] = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "AT")
                            {
                                int tRow = input.ScanList[n].Mode;
                                int tCol = input.ScanList[n].Cross;
                                int temp;
                                if (tRow < tCol)
                                {
                                    temp = tRow;
                                    tRow = tCol;
                                    tCol = temp;
                                }
                                input.CrossTermMatrix[tRow - 1, tCol - 1] = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                            if (input.ScanList[n].varToFit.ToUpper() == "JT")
                            {
                                int tRow = input.ScanList[n].Mode;
                                int tCol = input.ScanList[n].Cross;
                                int temp;
                                if (tRow > tCol)
                                {
                                    temp = tRow;
                                    tRow = tCol;
                                    tCol = temp;
                                }
                                input.CrossTermMatrix[tRow - 1, tCol - 1] = input.ScanList[n].Start + input.ScanList[n].Step * (double)h;
                            }
                        }//end adjust each variable to fit

                        //now make output file for each step
                        linesToWrite = OutputFile.inputFileMaker(input, Modes);
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input, input.useAbsoluteEV));

                        string stepFile = filepathOUT + "_step_" + Convert.ToString(h + 1) + ".out";
                        File.WriteAllLines(stepFile, linesToWrite);

                        //add this steps eigenvalues to the scan output file
                        scanList.Add(runner.finalList);

                        linesToWrite = null;
                        stepFile = null;
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

                //Writes the matrix to file if needed
                if (input.UseMatrixFile && !input.MatrixMade)
                {
                    OutputFile.writeMatFile(input);
                }//end if to write matrix to file
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
            catch (InvalidInput a)
            {
                Console.WriteLine(a.EMessage);
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
                Console.Write(a.EMessage);
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
                Console.WriteLine(err.EMessage);
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
        /// Appends timing information to final output file and writes to disk.
        /// </summary>
        /// <param name="totalTime">
        /// Timer tracking total execution time.
        /// </param>
        /// <param name="filepathOUT">
        /// Filepath for output file to be written to.
        /// </param>
        /// <param name="linesToWrite">
        /// List with content already generated by SOCJT and FITSOCJT routines.
        /// </param>
        private static void WriteOutputFile(Stopwatch totalTime, string filepathOUT, List<string> linesToWrite)
        {
            totalTime.Stop();
            double TIME = totalTime.ElapsedMilliseconds / 1000D;
            linesToWrite.Add(" ");
            linesToWrite.Add("SOCJT 2 has completed. Total time elapsed = " + String.Format("{0,11:0.0000}", TIME) + " seconds.");
            File.WriteAllLines(filepathOUT, linesToWrite);
            Console.WriteLine("SOCJT Has Completed");
        }

        /// <summary>
        /// Calculates the eigenvectors for very large calculations after main routines have run.
        /// </summary>
        /// <param name="filepath">
        /// Directory all files are being stored in.
        /// </param>
        /// <param name="input">
        /// FileInfo object
        /// </param>
        /// <param name="lanczosEVectors">
        /// Eigenvectors of the LanczosMatrix to be transformed to the eigenvectors of the Hamiltonian.
        /// </param>
        /// <param name="basisSet">
        /// Basis sets of the Hamiltonian blocks.
        /// </param>
        private static void CalcLargeEvecs(string filepath, FileInfo input, List<double[,]> lanczosEVectors, List<List<BasisFunction>> basisSet)
        {
            Console.WriteLine("\nEigenvectors being calculated...");
            Console.WriteLine("This may take some time.");
            eVecGenerator(lanczosEVectors, filepath, basisSet, input);
        }

        /// <summary>
        /// Checks if any values are being fit.
        /// </summary>
        /// <param name="input">
        /// FileInfo object.
        /// </param>
        /// <param name="Modes">
        /// Modes in this calculations
        /// </param>
        /// <returns>
        /// True if a fit is being run. False if not.
        /// </returns>
        private static bool IsFit(FileInfo input, List<ModeInfo> Modes)
        {
            //Sets the fit boolean value to true if any values are to be fit
            bool fit = false;
            for (int i = 0; i < input.nModes; i++)
            {
                if (input.FitAzeta == true)
                {
                    fit = true;
                    break;
                }
                if (input.FitOrigin == true)
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
            if (input.IncludeCrossTerms == true)
            {
                foreach (bool fitter in input.CrossTermFit)
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
            return fit;
        }

        /// <summary>
        /// Tests if the quadratic basis set needs to be used.
        /// </summary>
        /// <param name="input">
        /// FileInfo object.
        /// </param>
        /// <param name="Modes">
        /// List of Modes.
        /// </param>
        /// <returns>
        /// True if quadratic, false if linear only.
        /// </returns>
        private static bool IsQuad(FileInfo input, List<ModeInfo> Modes)
        {
            //checks to see if quadratic basis set should be used.
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

            //this sets isQuad = true if there is a cross quadratic term between an two E modes
            if (isQuad == false)
            {
                for (int i = 0; i < input.nModes; i++)
                {
                    if (input.CrossTermMatrix == null)
                    {
                        break;
                    }
                    for (int j = i + 1; j < input.nModes; j++)
                    {
                        if (input.CrossTermMatrix[i, j] != 0)
                        {
                            if (!Modes[i].IsAType && !Modes[j].IsAType)
                            {
                                isQuad = true;
                            }
                        }//end if
                    }//end inner for
                }//end outer for
            }//end if isQuad == false
            return isQuad;
        }

        /// <summary>
        /// Creates and initializes all modes.
        /// </summary>
        /// <param name="inputFile">
        /// The parsed input file.
        /// </param>
        /// <param name="input">
        /// FileInfo object
        /// </param>
        /// <returns>
        /// List of initialized modes.
        /// </returns>
        private static List<ModeInfo> ModeInitialization(string[] inputFile, FileInfo input)
        {
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
                    input.UseKappaEta = true;
                }
            }//end for
            return Modes;
        }

        /// <summary>
        /// Reads from a lanczos vector file and returns the specified vector
        /// </summary>
        /// <param name="filepath">
        /// Filepath of file to be read and parsed.
        /// </param>
        /// <param name="n">
        /// Which vector is needed
        /// </param>
        /// <param name="vecIn">
        /// StreamReader object from reading previous vectors
        /// </param>
        /// <returns>
        /// Array containing all values from text file.
        /// </returns>
        public static double[] vecRead(string filepath, int n, ref StreamReader vecIn)
        {
            List<double> vector = new List<double>();
            string lineS;
            while ((lineS = vecIn.ReadLine()) != ("START_VEC" + n))
            {
                continue;
            }
            while ((lineS = vecIn.ReadLine()) != "END_VEC")
            {
                vector.Add(FileInfo.ParseDouble(lineS));
            }//end while
            return vector.ToArray();
        }//end method vecRead

        /// <summary>
        /// Generates the eigenvectors for calculations where the basis set was too large to store all of the lanczos vectors in memory.
        /// In that case they are written to disk in a temp file and here read and used for the proper transformations.
        /// </summary>
        /// <param name="lanczosEVectors">
        /// Eigenvectors of the Lanczos matrix to be transformed.
        /// </param>
        /// <param name="filepath">
        /// Directory to read Lanczos vectors from and to write results to.
        /// </param>
        /// <param name="basisSet">
        /// Basis set of the hamiltonian.
        /// </param>
        /// <param name="input">
        /// FileInfo object
        /// </param>
        private static void eVecGenerator(List<double[,]> lanczosEVectors, string filepath, List<List<BasisFunction>> basisSet, FileInfo input)
        {
            List<double[,]> eVecs = new List<double[,]>();
            string file;
            double evMin = input.EigenvectorCoefficientMinimum;
            double[] lanczosVector;
            for (int i = 0; i < lanczosEVectors.Count; i++)
            {
                //this uses the same convention as in the Lanczos method which saves the vectors. Not as general as I'd like but probably ok.
                file = filepath + "temp_vecs_" + i + ".tmp";
                StreamReader vecIn = new StreamReader(file);
                lanczosVector = vecRead(file, 0, ref vecIn);
                int numberOfEigenvalues = lanczosEVectors[i].GetLength(1);
                int iterations = lanczosEVectors[i].GetLength(0);
                int dimension = lanczosVector.Length;
                vecIn.Close();
                vecIn = new StreamReader(file);
                double[,] temp = new double[dimension, numberOfEigenvalues];
                for (int j = 0; j < iterations; j++)
                {
                    lanczosVector = vecRead(file, j, ref vecIn);
                    for (int m = 0; m < numberOfEigenvalues; m++)
                    {
                        //read in the right vector to memory from the lanczos vector file
                        for (int n = 0; n < lanczosVector.Length; n++)
                        {
                            temp[n, m] += lanczosVector[n] * lanczosEVectors[i][j, m];
                        }//end loop over rows of lanczos Vector
                    }//end loop over columns of lanczosEVectors
                }
                //here delete the temp lanczos vector file
                vecIn.Close();
                File.Delete(file);
                Lanczos.normalize(ref temp);
                eVecs.Add(temp);
            }//end loop to generate eigenvectors
            //here use basis sets and eigenvectors and write them to file
            SOCJT.writeVecFile(input, eVecs, basisSet, evMin, filepath + input.Title + "_EVecs.out");
            if (input.EVectorFile)
            {
                SOCJT.writeVecFile(input, eVecs, basisSet, 0.0);
            }
        }//end eVecGenerator
    }//end class Program
}//end namespace
