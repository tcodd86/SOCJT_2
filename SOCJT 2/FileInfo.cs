using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;

namespace ConsoleApplication1
{
    /// <summary>
    /// Class for reading input files and storing values from it.
    /// </summary>
    class FileInfo
    {
        #region properties

        /// <summary>
        /// Number of modes in the calculation
        /// </summary>
        private int _nModes;
        public int nModes 
        {
            get { return _nModes; }
            private set
            {
                NonnegativeCheck(value, 0, "NMODES");
                _nModes = value;
            }
        }

        /// <summary>
        /// Title of the calculations
        /// </summary>
        public string title { get; private set;}

        /// <summary>
        /// Value of spin
        /// </summary>
        private decimal _S;
        public decimal S
        {
            get { return _S; } 
            private set
            {
                //check that spin is integer or half integer only
                if (value % 0.5M != 0M)
                {
                    throw new InvalidInput("S");
                }
                _S = value;
            } 
        }

        /// <summary>
        /// Value of a*zeta, the SO coupling constant
        /// </summary>
        public double Azeta { get; set; }

        /// <summary>
        /// True if the SO coupling constant is being fit
        /// </summary>
        public bool fitAzeta { get; private set; }

        /// <summary>
        /// Maximum value of j to use in the basis set
        /// </summary>
        private decimal _maxJ;
        public decimal maxJ
        {
            get { return _maxJ; }
            set 
            {
                NonnegativeCheck(value, 0M, "MAXJ");
                if (value % 0.5M != 0M || value < 0.5M)
                {
                    throw new InvalidInput("MAXJ");
                }
                _maxJ = value;
            }
        }

        /// <summary>
        /// Minimum value of j that should be used in calculations
        /// </summary>
        private decimal _minJ;
        public decimal minJ
        {
            get { return _minJ; }
            set
            {
                NonnegativeCheck(value, 0M, "MINJ");
                if (value % 0.5M != 0M)
                {
                    throw new InvalidInput("MINJ");
                }
                _minJ = value;
            }
        }

        /// <summary>
        /// True if min J is user specified
        /// </summary>
        public bool minJBool { get; private set; }

        /// <summary>
        /// Value of S1
        /// </summary>
        private int _S1;
        public int S1
        {
            get { return _S1; }
            set
            {
                if (value != 0 && value != 1)
                {
                    throw new InvalidInput("S1");
                }
                _S1 = value;
            }
        }

        /// <summary>
        /// Value of S2
        /// </summary>
        private int _S2;
        public int S2
        {
            get { return _S2; }
            set
            {
                if (value != 0 && value != 1)
                {
                    throw new InvalidInput("S2");
                }
                _S2 = value;
            }
        }

        /// <summary>
        /// True if the basis set should be printed in the output file.
        /// </summary>
        public bool printBasis { get; private set; }

        /// <summary>
        /// True if the Hamiltonian should be printed in the output file.
        /// </summary>
        public bool pMatrix { get; private set; }

        /// <summary>
        /// True if the eigenvectors should be printed in the output file.
        /// </summary>
        public bool pVector { get; set; }  

        /// <summary>
        /// True if a file with the eigenvectors should be printed.
        /// </summary>
        public bool vecFile { get; private set; }

        /// <summary>
        /// Value of the origin to add to the eigenvalues
        /// </summary>
        public double origin { get; set; }

        /// <summary>
        /// True if the origin should be fit.
        /// </summary>
        public bool fitOrigin { get; private set; }

        /// <summary>
        /// True if you want to save a file with the basis set in it.
        /// </summary>
        public bool basisFile { get; private set; }

        /// <summary>
        /// Number of eigenvalues/eigenvectors to find
        /// </summary>
        private int _M;
        public int M 
        {
            get { return _M; }
            private set
            {
                NonnegativeCheck(value, 0, "M");
                _M = value;
            }
        }

        /// <summary>
        /// Size of block (number of columns) in Block Lanczos
        /// </summary>
        private int _kFactor;
        public int kFactor 
        {
            get { return _kFactor; }
            private set
            {
                NonnegativeCheck(value, 0, "K_FACTOR");
                _kFactor = value;
            }
        }

        /// <summary>
        /// Max number of iterations to run the block Lanczos or size of Lanczos matrix to be generated in Naive Lanczos
        /// </summary>
        private int _noIts;
        public int noIts 
        {
            get { return _noIts; }
            private set 
            { 
                NonnegativeCheck(value, 0, "NOITS");
                _noIts = value;
            }
        }

        /// <summary>
        /// Tolerance used in Block Lanczos for convergance or in Naive Lanczos for eigenvalue comparison
        /// </summary>
        private double _tol;
        public double tol 
        {
            get { return _tol; }
            private set
            {
                NonnegativeCheck(value, 0.0, "TOL");
                _tol = value;
            }
        }

        /// <summary>
        /// Name of the fit file to be used for a fit.
        /// </summary>
        public string fitFile { get; private set; }

        /// <summary>
        /// F-Tolerance value for LM optimization
        /// </summary>
        private double _fTol;
        public double fTol 
        {
            get { return _fTol; }
            private set
            {
                NonnegativeCheck(value, 0.0, "FTOL");
                _fTol = value;
            }
        }

        /// <summary>
        /// X-Tolerance value for LM optimization
        /// </summary>
        private double _xTol;
        public double xTol
        {
            get { return _xTol; }
            private set
            {
                NonnegativeCheck(value, 0.0, "XTOL");
                _xTol = value;
            }
        }

        /// <summary>
        /// G-Tolerance value for LM optimization
        /// </summary>
        private double _gTol;
        public double gTol
        {
            get { return _gTol; }
            private set
            {
                NonnegativeCheck(value, 0.0, "GTOL");
                _gTol = value;
            }
        }
        
        /// <summary>
        /// Max number of iterations (steps) in the LM optimizer
        /// </summary>
        private int _maxFev;
        public int maxFev
        {
            get { return _maxFev; }
            private set 
            {
                NonnegativeCheck(value, 0, "MAXFEV");
                _maxFev = value;
            }
        }

        /// <summary>
        /// Factor used for step size in LM optimizer
        /// </summary>
        private double _factor;
        public double factor
        {
            get { return _factor; }
            private set
            {
                NonnegativeCheck(value, 0.0, "FACTOR");
                _factor = value;
            }
        }

        /// <summary>
        /// True if there are cross-terms in the Hamiltonian
        /// </summary>
        public bool includeCrossTerms { get; private set; }

        /// <summary>
        /// True if a scan is being run
        /// </summary>
        public bool scan { get; set; }

        /// <summary>
        /// Number of steps to run in a scan
        /// </summary>
        private int _steps;
        public int steps
        {
            get { return _steps; }
            set
            {
                NonnegativeCheck(value, 0, "SCAN STEPS");
                _steps = value;
            }
        }

        /// <summary>
        /// How long (in seconds) the Hamiltonian generation took.
        /// </summary>
        public double matGenTime { get; set; }

        /// <summary>
        /// How long (in seconds) the diagonalization took.
        /// </summary>
        public double diagTime { get; set; }

        /// <summary>
        /// How much to parallelize the Hamiltonian matrix generation
        /// </summary>
        private int _parMat;
        public int parMat
        {
            get { return _parMat; }
            private set
            {
                NonnegativeCheck(value, 0, "PARMAT");
                _parMat = value;
            }
        }

        /// <summary>
        /// How much to parallelize the matrix vector multiplication
        /// </summary>
        private int _parVec;
        public int parVec
        {
            get { return _parVec; }
            private set
            {
                NonnegativeCheck(value, 0, "PARVEC");
                _parVec = value;
            }
        }

        /// <summary>
        /// How many j-blocks should be run in parallel
        /// </summary>
        private int _parJ;
        public int parJ
        {
            get { return _parJ; }
            private set
            {
                NonnegativeCheck(value, 0, "PARJ");
                _parJ = value;
            }
        }

        /// <summary>
        /// Minimum value of coefficients to print in eigenvectors
        /// </summary>
        private double _evMin;
        public double evMin
        {
            get { return _evMin; }
            private set
            {
                NonnegativeCheck(value, 0.0, "EVMIN");
                _evMin = value;
            }
        }

        /// <summary>
        /// True if SO coupling is nonzero
        /// </summary>
        public bool inclSO { get; private set; }

        public List<Tuple<decimal, int, int>> eVecs { get; private set; }
        
        /// <summary>
        /// True if using kappa and eta for linear and quadratic JT coupling instead of D and K
        /// </summary>
        public bool useKappaEta { get; set; }//end useKappaEta

        /// <summary>
        /// List of Scanner objects for scan. Contains, variable to scan, step size, and start value.
        /// </summary>
        public List<Scanner> scanList;

        /// <summary>
        /// Array containing cross-terms
        /// </summary>
        public double[,] crossTermMatrix;

        /// <summary>
        /// Array containing booleans for fitting various cross-terms
        /// </summary>
        public bool[,] crossTermFit;

        /// <summary>
        /// True if using block lanczos instead of naive Lanczos.
        /// </summary>
        public bool blockLanczos { get; private set; }//end blockLanczos

        /// <summary>
        /// Filepath of matrix file to read.
        /// </summary>
        public string matFilePath { get; set; }//end matFilePath

        /// <summary>
        /// String value. File name of file containing the previous matrix or the filename desired to write the matrix to.
        /// </summary>
        public string matFile { get; private set; }//end matFile

        /// <summary>
        /// Boolean indicating that a matrix file should be used or generated.
        /// </summary>
        public bool useMatFile { get; private set; }//end useMatFile

        /// <summary>
        /// Boolean indicating whether the matrix file exists or not.
        /// </summary>
        public bool matMade { get; set; }//end matMade

        /// <summary>
        /// String to store the directory of the input and output files.
        /// </summary>
        public string filePath { get; set; }

        /// <summary>
        /// Boolean indicating whether or not an output file of the eigenvectors using the complete basis set should be printed.
        /// </summary>
        public bool vecFileComplete { get; private set; }

        #endregion properties

        /// <summary>
        /// Constructor for FileInfo object. Sets reasonable values for most properties in case user forgets something in input file.
        /// </summary>
        public FileInfo()
        {            
            fitAzeta = false;
            minJBool = false;
            fitOrigin = false;
            inclSO = false;
            printBasis = false;
            pMatrix = false;
            pVector = false;
            includeCrossTerms = false;
            useKappaEta = false;
            blockLanczos = false;
            useMatFile = false;
            matMade = false;
            vecFile = false;
            vecFileComplete = false;
            
            matFile = "matrix.txt";
            title = "TITLE";
            fitFile = "fit.fit";

            //these are reasonable values of J for a basic quadratic problem
            maxJ = 7.5M;
            minJ = 0.5M;            

            //S1 and S2 values for 3-fold symmetry
            S1 = 0;
            S2 = 1;
            
            origin = 0.0;
            parMat = 1;
            nModes = 1;
            S = 0.5M;
            Azeta = 0.0;
            evMin = 0.2;            
            parVec = 1;
            parMat = 1;
            parJ = 2;
            M = 5;
            kFactor = 2;
            noIts = 2000;
            tol = 0.000001;
            fTol = 0.0;
            xTol = 0.0;
            gTol = 0.0;
            maxFev = 25;
            factor = 0.001;        
        }
        
        /// <summary>
        /// Reads a text file into memory and parses it into a string array.
        /// </summary>
        /// <param name="filepath">
        /// Filepath of file to be read and parsed.
        /// </param>
        /// <returns>
        /// Array containing all values from text file.
        /// </returns>
        public static string[] fileRead(string filepath)
        {
            List<string> inputF = new List<string>();
            string[] inputFa = { };
            using (StreamReader SOCJTin = new StreamReader(filepath))
            {
                string lineS;
                string[] SOCJTNewLine;
                char[] delimiters = new char[] { '\t', '\r', '=', ' ' };
                while ((lineS = SOCJTin.ReadLine()) != null)
                {                    
                    SOCJTNewLine = lineS.Split(delimiters, StringSplitOptions.RemoveEmptyEntries);
                    for (int i = 0; i < SOCJTNewLine.Length; i++)
                    {
                        inputF.Add(SOCJTNewLine[i]);
                    }
                }//end while
            }//end StreamReader
            inputFa = inputF.ToArray();
            return inputFa;
        }//end method fileRead

        public void setFileInfo(string[] inputf, string filepath)
        {
            this.filePath = filePath;
            for (int i = 0; i < inputf.Length; i++)
            {
                if (inputf[i].ToUpper() == "&GENERAL")
                {
                    #region &GENERAL
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "TITLE")
                        {
                            title = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "NMODES")
                        {
                            nModes = Convert.ToInt16(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "S")
                        {
                            S = parseDecimal(inputf[u + 1]);
                            if (S < 0M)
                            {
                                S = S * -1M;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "AZETA")
                        {
                            Azeta = parseDouble(inputf[u + 1]);
                            if (Azeta != 0D)
                            {
                                inclSO = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FIT_AZETA")
                        {
                            fitAzeta = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MAXJ")
                        {
                            maxJ = parseDecimal(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MINJ")
                        {
                            minJ = parseDecimal(inputf[u + 1]);
                            minJBool = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FIT_ORIGIN")
                        {
                            fitOrigin = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "ORIGIN")
                        {
                            origin = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "S1")
                        {
                            S1 = Convert.ToInt16(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "S2")
                        {
                            S2 = Convert.ToInt16(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "USE_KAPPA_ETA")
                        {
                            useKappaEta = TorF(inputf[u + 1]);
                        }
                        if (inputf[u] == "/")
                        {
                            break;
                        }
                    }//end for loop
                    continue;
                    #endregion                    
                }//end GENERAL if

                if (inputf[i].ToUpper() == "&IO_INFO")
                {
                    #region &IO_INFO
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "PRINT_BASIS")
                        {
                            printBasis = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PRINT_MATRIX")
                        {
                            pMatrix = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PRINT_VEC")
                        {
                            pVector = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "USE_MATRIX_FILE")
                        {
                            useMatFile = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "VEC_FILE")
                        {
                            vecFile = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "VEC_FILE_COMPLETE")
                        {
                            vecFileComplete = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MATRIX_FILE")
                        {
                            matFile = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "EV_MIN")
                        {
                            evMin = Convert.ToDouble(inputf[u + 1]);
                            if (evMin < 0.0)
                            {
                                evMin = -1.0 * evMin;
                            }
                            continue;
                        }
                        if (inputf[u] == "/")
                        {
                            break;
                        }
                    }//end for loop
                    continue;
                    #endregion                    
                }//end IO INFO if

                if (inputf[i].ToUpper() == "&SOLVE_INFO")
                {
                    #region &SOLVE_INFO
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "M")
                        {
                            M = Convert.ToInt32(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "K_FACTOR")
                        {
                            kFactor = Convert.ToInt32(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "NOITS")
                        {
                            noIts = Convert.ToInt32(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "TOL")
                        {
                            tol = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PARVEC")
                        {
                            parVec = Convert.ToInt16(inputf[u + 1]);
                            if (parVec < 1)
                            {
                                parVec = 1;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PARMAT")
                        {
                            parMat = Convert.ToInt16(inputf[u + 1]);
                            if (parMat < 1)
                            {
                                parMat = 1;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PARJ")
                        {
                            parJ = Convert.ToInt16(inputf[u + 1]);
                            if (parJ < 1)
                            {
                                parJ = 1;
                            }
                        }
                        if (inputf[u].ToUpper() == "BLOCK_LANCZOS")
                        {
                            blockLanczos = TorF(inputf[u + 1]);
                        }
                        if (inputf[u] == "/")
                        {
                            break;
                        }
                    }//end for loop                    
                    continue;
                    #endregion
                }//end SOLVE INFO if

                if (inputf[i].ToUpper() == "&FIT_INFO")
                {
                    #region &FIT_INFO
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "FITFILE")
                        {
                            fitFile = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FTOL")
                        {
                            fTol = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "XTOL")
                        {
                            xTol = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "GTOL")
                        {
                            gTol = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MAXFEV")
                        {
                            maxFev = Convert.ToInt32(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FACTOR")
                        {
                            factor = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "/")
                        {
                            break;
                        }
                    }//end for
                    #endregion
                }//end FIT INFO if                

                if (inputf[i].ToUpper() == "&CROSS_TERMS")
                {
                    #region &CROSS_TERMS
                    int row;
                    int column;
                    int temp;
                    bool tbool;
                    includeCrossTerms = true;
                    crossTermMatrix = new double[nModes, nModes];
                    crossTermFit = new bool[nModes, nModes];
                    for (int l = 0; l < nModes; l++)
                    {
                        for (int m = 0; m < nModes; m++)
                        {
                            crossTermMatrix[l, m] = 0D;
                            crossTermFit[l, m] = false;
                        }
                    }//end loops to initalize crossTermMatrix to all 0's

                    for (int j = i; ; j++)
                    {
                        if (inputf[j].ToUpper() == "JT")
                        {
                            tbool = false;
                            row = Convert.ToInt16(inputf[j + 2]) - 1;
                            column = Convert.ToInt16(inputf[j + 4]) - 1;
                            if (row > column)
                            {
                                temp = row;
                                row = column;
                                column = temp;
                            }//end if
                            crossTermMatrix[row, column] = parseDouble(inputf[j + 5]);
                            tbool = TorF(inputf[j + 7]);
                            crossTermFit[row, column] = tbool;
                            continue;
                        }
                        if (inputf[j] == "/")
                        {
                            break;
                        }//end end of CROSS_TERM if
                    }//end CROSS_TERM for
                    #endregion
                }//end CROSS_TERM if
                
                if (inputf[i].ToUpper() == "&SCAN")
                {
                    #region &SCAN
                    scanList = new List<Scanner>();
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "STEPS")
                        {
                            steps = Convert.ToInt32(inputf[u + 1]);
                        }
                        if (inputf[u].ToUpper() == "MODE")
                        {
                            Scanner tempScan = new Scanner();
                            tempScan.Mode = Convert.ToInt32(inputf[u + 1]);
                            tempScan.varToFit = inputf[u + 2];
                            tempScan.Start = parseDouble(inputf[u + 3]);
                            tempScan.Step = parseDouble(inputf[u + 4]);
                            scanList.Add(tempScan);
                            tempScan = null;
                            scan = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "CROSS")
                        {
                            Scanner tempScan = new Scanner();
                            tempScan.Mode = Convert.ToInt32(inputf[u + 1]);
                            tempScan.Cross = Convert.ToInt32(inputf[u + 2]);
                            tempScan.varToFit = inputf[u + 3];
                            tempScan.Start = parseDouble(inputf[u + 4]);
                            tempScan.Step = parseDouble(inputf[u + 5]);
                            scanList.Add(tempScan);
                            tempScan = null;
                            scan = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "/")
                        {
                            break;
                        }
                    }
                    #endregion
                }
            }//end for

            //this is to initialize the cross-term matrix no matter what
            if (crossTermMatrix == null)
            {
                crossTermMatrix = new double[nModes, nModes];
                crossTermFit = new bool[nModes, nModes];
            }
        }//end method setFileInfo

        /// <summary>
        /// Checks a string to see if it's T or TRUE, converts to boolean value
        /// </summary>
        /// <param name="val">
        /// String to be checked.
        /// </param>
        /// <returns>
        /// True if the string is T or TRUE, false if not
        /// </returns>
        public static bool TorF(string val)
        {
            bool trfa = false;
            if (val.ToUpper() == "T" || val.ToUpper() == "TRUE")
            {
                trfa = true;
            }
            return trfa;
        }

        /// <summary>
        /// To parse a string containing a double that may or may not have scientific notation in it.
        /// </summary>
        /// <param name="s">
        /// String to be parsed
        /// </param>
        /// <returns>
        /// Double value of parsed string
        /// </returns>
        public static double parseDouble(string s)
        {
            try
            {
                return Double.Parse(s, NumberStyles.AllowExponent | NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint);
            }
            catch
            {
                throw new InvalidInput("A numerical item");
            }
        }

        /// <summary>
        /// To parse a string containing a decimal that may or may not have scientific notation in it.
        /// </summary>
        /// <param name="s">
        /// String to be parsed
        /// </param>
        /// <returns>
        /// Decimal value of parsed string
        /// </returns>
        public static decimal parseDecimal(string s)
        {
            try
            {
                return Decimal.Parse(s, NumberStyles.AllowExponent | NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint);
            }
            catch
            {
                throw new InvalidInput("A numerical item");
            }
        }

        /// <summary>
        /// To save me from writing negative value checks for every single parameter.
        /// </summary>
        /// <typeparam name="T">
        /// Anyhing implementing the IComparable Interface
        /// </typeparam>
        /// <param name="val">
        /// Parameter to be checked.
        /// </param>
        /// <param name="zero">
        /// Value to be checked against. Should be zero of the same type as val
        /// </param>
        /// <param name="s">
        /// String representation of the parameter so that if an error is thrown a description can be given.
        /// </param>
        public static void NonnegativeCheck<T>(T val, T zero, string s) where T : IComparable<T>
        {
            if (val.CompareTo(zero) < 0)
            {
                throw new InvalidInput(s);
            }
        }
    }//class FileInfo
}//end namespace