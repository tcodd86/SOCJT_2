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
                //{
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
                        //Commenting this out unless I add an input file maker helper method
                        /*
                        Console.WriteLine("The directory does not exist.  Would you like to create it? Y/N");
                        string YN = Console.ReadLine();
                        if (YN.ToUpper() == "Y" || YN.ToUpper() == "YES")
                        {
                            Directory.CreateDirectory(fileDirectory);
                        }
                        else
                        {
                            */
                        throw new DirectoryNotFoundException();
                        //}
                    }
                //}
#else
                //else//meaning it's running mono on the Linux cluster
                //{
                    fileDirectory = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
                //}
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
                filepath = string.Concat("\\");
                string filepathIN = string.Copy(filepath);
                string filepathOUT = string.Copy(filepath);
                string filepathFIT = string.Copy(filepath);
                filepathIN = string.Concat(inFileName);
                filepathOUT = string.Concat(outFile);

                //read input file
                string[] inputFile = FileInfo.fileRead(filepathIN);

                //create new FileInfo object
                FileInfo input = new FileInfo();

                //set input object data from input file
                input.setFileInfo(inputFile);

                //make the fitfile point to something
                filepathFIT = string.Concat(input.fitFile);

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

                //main subroutine execution when not running a scan
                if (input.Scan == false)
                {
                    //initialize variable to hold the output file and initialize it with input file containing original values
                    List<string> linesToWrite = new List<string>();
                    linesToWrite = OutputFile.inputFileMaker(input, Modes);

                    //if not fitting anyvariables run SOCJT routine directly
                    if (fit == false)
                    {
                        SOCJT runner = new SOCJT();
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input));
                    }
                    else//else run FitSOCJT routine which will run LM optimizer which will call SOCJT
                    { 
                        linesToWrite.AddRange(FitSOCJT.fit(Modes, isQuad, inputFile, input, filepathFIT));
                    }

                    //end stopwatch for total program time execution and appends the total run time to the output file
                    totalTime.Stop();
                    double TIME = totalTime.ElapsedMilliseconds / 1000D;
                    linesToWrite.Add(" ");
                    linesToWrite.Add("SOCJT 2 has completed. Total time elapsed = " + String.Format("{0,11:0.0000}", TIME) + " seconds.");
#if DEBUG
                    linesToWrite.Add("Orthog took " + Lanczos.reorthogTime / 1000L + " seconds");
#endif

                    //writes all info to the output file
                    File.WriteAllLines(filepathOUT, linesToWrite);
                }//end no scan

                //means a scan is being run
                else
                {
                    //create a list to store the final output eigenvalues from each iteration of SOCJT
                    List<Eigenvalue[]> scanList = new List<Eigenvalue[]>();

                    //create a List to hold the output file information
                    List<string> finalOut = new List<string>();

                    //for loop to run steps of scan
                    for (int h = 0; h < input.Steps; h++)
                    {
                        //basically run the SOCJT subroutine for each step of the scan and save the output file at the end of the loop.
                        List<string> linesToWrite = new List<string>();
                        SOCJT runner = new SOCJT();

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
                }
            }//end try block

            //Exception handling
            catch (DirectoryNotFoundException)
            {
                Console.WriteLine("The directory does not exist");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (FileNotFoundException)
            {
                Console.WriteLine("The file does not exist");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (IndexOutOfRangeException)
            {
                Console.WriteLine("An index out of range exception has occurred");
                Console.WriteLine("Check to make sure that NMODES is correct and");
                Console.WriteLine("that any CROSS_TERM arguments are correct.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (FitFileNotFoundException)
            {
                Console.WriteLine("The fit file does not exist.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            catch (BasisSetTooSmallException)
            {
                Console.WriteLine("At least one of the j-blocks has a dimension smaller than M.");
                Console.WriteLine("Increase the basis set or decrease M.");
                Console.WriteLine("Press enter to terminate the program.");
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
#if !DEBUG
            catch (Exception ex)
            {
                Console.WriteLine("An exception has occurred: " + ex.Message);
                Console.WriteLine("Check your input file and try again.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
#endif
        }//end Main
    }//end class Program
}//end namespace
