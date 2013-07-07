using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace ConsoleApplication1
{
    class MainS
    {
        static void Main(string[] args)
        {
            try
            {
                Console.WriteLine("Enter file directory or press enter to use C:\\SOCJT 2");
                string fileDirectory = Console.ReadLine();
                if(fileDirectory == "")
                {
                    fileDirectory = "C:\\SOCJT 2";
                }
                if(Directory.Exists(fileDirectory) == false)
                {
                    Console.WriteLine("The directory does not exist.  Would you like to create it? Y/N");
                    string YN = Console.ReadLine();
                    if (YN.ToUpper() == "Y" || YN.ToUpper() == "YES")
                    {
                        Directory.CreateDirectory(fileDirectory);
                    }
                    else
                    {
                        throw new DirectoryNotFoundException();
                    }
                }
                Directory.SetCurrentDirectory(fileDirectory);
                Console.WriteLine("Enter file name including extension:");
                string inFileName = Console.ReadLine();
                Console.WriteLine("Enter output file name or press enter to use " + inFileName + ".out:");
                string outFile = Console.ReadLine();
                if (outFile == "" || outFile == " ")
                {
                    outFile = string.Concat(inFileName, ".out");
                }

                Stopwatch totalTime = new Stopwatch();
                totalTime.Start();
                string filepath = string.Concat(fileDirectory);
                filepath = string.Concat("\\");
                string filepathIN = string.Copy(filepath);
                string filepathOUT = string.Copy(filepath);
                string filepathFIT = string.Copy(filepath);
                filepathIN = string.Concat(inFileName);
                filepathOUT = string.Concat(outFile);
                string[] inputFile = FileInfo.fileRead(filepathIN);
                FileInfo input = new FileInfo();

                input.Scan = false;//added this

                input.setFileInfo(inputFile);
                filepathFIT = string.Concat(input.fitFile);

                if (input.S % 0.5M != 0M)
                {
                    throw new SpinInvalidException();
                }

                //This iterates through the Mode info in the input file and initializes the modes
                List<Mode> Modes = new List<Mode>();
                for (int i = 0; i < input.nModes; i++)
                {
                    Modes.Add(new Mode(i, inputFile));
                    /*
                    Modes.Add(new Mode());
                    Modes[i].fitOmega = false;
                    Modes[i].fitD = false;
                    Modes[i].fitK = false;
                    Modes[i].fitWEXE = false;
                    Modes[i].IsAType = false;
                    Modes[i].setMode(Modes[i], i, inputFile);
                    */
                }//end for

                //This sets bool isQuad = false if K is zero and fitK is false for all modes.
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

                //this sets isQuad = true if there is a cross quartic term between an A and E mode
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
                                }
                            }
                        }
                    }
                }

                //This sets the fit boolean value to true if any values are to be fit
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


                if (input.Scan == false)
                {
                    List<string> linesToWrite = new List<string>();
                    linesToWrite = OutputFile.inputFileMaker(input, Modes);
                    if (fit == false)
                    {
                        SOCJT runner = new SOCJT();
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input));
                    }
                    else
                    {
                        linesToWrite.AddRange(FitSOCJT.fit(Modes, isQuad, inputFile, input, filepathFIT));
                    }
                    totalTime.Stop();
                    double TIME = totalTime.ElapsedMilliseconds / 1000D;
                    linesToWrite.Add(" ");
                    linesToWrite.Add("SOCJT 2 has completed. Total time elapsed = " + String.Format("{0,11:0.0000}", TIME) + " seconds.");
                    File.WriteAllLines(filepathOUT, linesToWrite);
                }//end no scan

                //means a scan is being run
                else
                {
                    List<Eigenvalue[]> scanList = new List<Eigenvalue[]>();
                    List<string> finalOut = new List<string>();
                    for (int h = 0; h < input.Steps; h++)
                    {
                        List<string> linesToWrite = new List<string>();
                        SOCJT runner = new SOCJT();
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

                        linesToWrite = OutputFile.inputFileMaker(input, Modes);
                        linesToWrite.AddRange(runner.SOCJTroutine(Modes, isQuad, inputFile, input));
                        string stepFile = filepathOUT + "_step_" + Convert.ToString(h + 1) + ".out";
                        File.WriteAllLines(stepFile, linesToWrite);
                        scanList.Add(runner.finalList);

                        linesToWrite = null;
                        stepFile = null;
                        runner = null;
                    }//end of steps loop

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
            /*
            catch (Exception ex)
            {
                Console.WriteLine("An exception has occurred: " + ex.Message);
                Console.WriteLine("Check your input file and try again.");
                Console.WriteLine("Press enter to terminate the program.");
                Console.ReadLine();
            }
            */
        }//end Main
    }//end class Program
}//end namespace
