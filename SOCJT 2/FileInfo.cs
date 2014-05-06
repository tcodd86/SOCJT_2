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
        public string Title { get; private set;}

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
        public bool FitAzeta { get; private set; }

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
        public bool MinJBool { get; private set; }

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
        public bool PrintBasis { get; private set; }

        /// <summary>
        /// True if the Hamiltonian should be printed in the output file.
        /// </summary>
        public bool PrintMatrix { get; private set; }

        /// <summary>
        /// True if the eigenvectors should be printed in the output file.
        /// </summary>
        public bool PrintVector { get; set; }  

        /// <summary>
        /// True if a file with the eigenvectors should be printed.
        /// </summary>
        public bool EVectorFile { get; private set; }

        /// <summary>
        /// Value of the origin to add to the eigenvalues
        /// </summary>
        public double Origin { get; set; }

        /// <summary>
        /// True if the origin should be fit.
        /// </summary>
        public bool FitOrigin { get; private set; }

        /// <summary>
        /// True if you want to save a file with the basis set in it.
        /// </summary>
        public bool BasisSetFile { get; private set; }

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
        public int NumberOfIts 
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
        public double Tolerance 
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
        public string FitFile { get; private set; }

        /// <summary>
        /// F-Tolerance value for LM optimization
        /// </summary>
        private double _fTol;
        public double FTol 
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
        public double XTol
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
        public double GTol
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
        public int MaxOptimizerSteps
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
        public double Factor
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
        public bool IncludeCrossTerms { get; private set; }

        /// <summary>
        /// True if a scan is being run
        /// </summary>
        public bool Scan { get; set; }

        /// <summary>
        /// Number of steps to run in a scan
        /// </summary>
        private int _steps;
        public int ScanSteps
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
        public double MatrixGenerationTime { get; set; }

        /// <summary>
        /// How long (in seconds) the diagonalization took.
        /// </summary>
        public double DiagonalizationTime { get; set; }

        /// <summary>
        /// How much to parallelize the Hamiltonian matrix generation
        /// </summary>
        private int _parMat;
        public int ParMatrix
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
        public int ParVectorMultiplication
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
        public int ParJ
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
        public double EigenvectorCoefficientMinimum
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
        public bool IncludeSO { get; private set; }

        public List<Tuple<decimal, int, int>> eVecs { get; private set; }
        
        /// <summary>
        /// True if using kappa and eta for linear and quadratic JT coupling instead of D and K
        /// </summary>
        public bool UseKappaEta { get; set; }//end useKappaEta

        /// <summary>
        /// List of Scanner objects for scan. Contains, variable to scan, step size, and start value.
        /// </summary>
        public List<Scanner> ScanList;

        /// <summary>
        /// Array containing cross-terms
        /// </summary>
        public double[,] CrossTermMatrix;

        /// <summary>
        /// Array containing booleans for fitting various cross-terms
        /// </summary>
        public bool[,] CrossTermFit;

        /// <summary>
        /// True if using block lanczos instead of naive Lanczos.
        /// </summary>
        public bool BlockLanczos { get; private set; }//end blockLanczos

        /// <summary>
        /// Filepath of matrix file to read.
        /// </summary>
        public string MatrixFilePath { get; set; }//end matFilePath

        /// <summary>
        /// String value. File name of file containing the previous matrix or the filename desired to write the matrix to.
        /// </summary>
        public string MatrixFile { get; private set; }//end matFile

        /// <summary>
        /// Boolean indicating that a matrix file should be used or generated.
        /// </summary>
        public bool UseMatrixFile { get; private set; }//end useMatFile

        /// <summary>
        /// Boolean indicating whether the matrix file exists or not.
        /// </summary>
        public bool MatrixMade { get; set; }//end matMade

        /// <summary>
        /// String to store the directory of the input and output files.
        /// </summary>
        public string FilePath { get; set; }

        /// <summary>
        /// Boolean indicating whether or not an output file of the eigenvectors using the complete basis set should be printed.
        /// </summary>
        public bool VectorFileComplete { get; private set; }

        /// <summary>
        /// Boolean indicating whether or not an eigenvector will be checked for accuracy
        /// </summary>
        public bool CheckEigenvector { get; private set; }
        
        /// <summary>
        /// File name of the eigenvector being checked
        /// </summary>
        public string EigenvectorFileName { get; private set; }

        /// <summary>
        /// Which j block the eigenvector belongs to and which eigenvector in that block is being checked.
        /// </summary>
        public Tuple<decimal, int> JBlockEigenvector { get; private set; }

        /// <summary>
        /// The tolerance that will define if a vector is considered to be an eigenvector or not
        /// </summary>
        public double EigenvectorTolerance { get; set; }

        /// <summary>
        /// Boolean to indicate that a vector from John Stanton is being checked.
        /// </summary>
        public bool JSInten { get; private set; }

        /// <summary>
        /// Name of the vector file to be read.
        /// </summary>
        public string VectorName { get; private set; }

        /// <summary>
        /// Boolean indicating whether intensity is being checked or not.
        /// </summary>
        public bool Intensity { get; private set; }

        /// <summary>
        /// Int indicating which vector from an eigenvector file should be read.
        /// </summary>
        public int VectorIndex { get; private set; }

        /// <summary>
        /// Decimal indicating which jblock the vector to be read is in.
        /// </summary>
        public decimal VectorJBlock { get; private set; }

        #endregion properties

        /// <summary>
        /// Constructor for FileInfo object. Sets reasonable values for most properties in case user forgets something in input file.
        /// </summary>
        public FileInfo()
        {            
            FitAzeta = false;
            MinJBool = false;
            FitOrigin = false;
            IncludeSO = false;
            PrintBasis = false;
            PrintMatrix = false;
            PrintVector = false;
            IncludeCrossTerms = false;
            UseKappaEta = false;
            BlockLanczos = false;
            UseMatrixFile = false;
            MatrixMade = false;
            EVectorFile = false;
            VectorFileComplete = false;
            CheckEigenvector = false;
            JSInten = false;
            Intensity = false;
            
            MatrixFile = "matrix.txt";
            Title = "TITLE";
            FitFile = "fit.fit";
            VectorName = "vector.txt";

            //these are reasonable values of J for a basic quadratic problem
            maxJ = 7.5M;
            minJ = 0.5M;            

            //S1 and S2 values for 3-fold symmetry
            S1 = 0;
            S2 = 1;
            
            Origin = 0.0;
            ParMatrix = 1;
            nModes = 1;
            S = 0.5M;
            Azeta = 0.0;
            EigenvectorCoefficientMinimum = 0.2;            
            ParVectorMultiplication = 1;
            ParMatrix = 1;
            ParJ = 2;
            M = 5;
            kFactor = 2;
            NumberOfIts = 2000;
            Tolerance = 0.000001;
            FTol = 0.0;
            XTol = 0.0;
            GTol = 0.0;
            MaxOptimizerSteps = 25;
            Factor = 0.001;
            EigenvectorTolerance = 0.0001;
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
        public static string[] FileRead(string filepath)
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

        public void SetFileInfo(string[] inputf, string filepath)
        {
            this.FilePath = FilePath;
            for (int i = 0; i < inputf.Length; i++)
            {
                if (inputf[i].ToUpper() == "&GENERAL")
                {
                    #region &GENERAL
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "TITLE")
                        {
                            Title = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "NMODES")
                        {
                            nModes = Convert.ToInt16(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "S")
                        {
                            S = ParseDecimal(inputf[u + 1]);
                            if (S < 0M)
                            {
                                S = S * -1M;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "AZETA")
                        {
                            Azeta = ParseDouble(inputf[u + 1]);
                            if (Azeta != 0D)
                            {
                                IncludeSO = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FIT_AZETA")
                        {
                            FitAzeta = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MAXJ")
                        {
                            maxJ = ParseDecimal(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MINJ")
                        {
                            minJ = ParseDecimal(inputf[u + 1]);
                            MinJBool = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FIT_ORIGIN")
                        {
                            FitOrigin = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "ORIGIN")
                        {
                            Origin = ParseDouble(inputf[u + 1]);
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
                            UseKappaEta = TorF(inputf[u + 1]);
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
                            PrintBasis = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PRINT_MATRIX")
                        {
                            PrintMatrix = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PRINT_VEC")
                        {
                            PrintVector = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "USE_MATRIX_FILE")
                        {
                            UseMatrixFile = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "VEC_FILE")
                        {
                            EVectorFile = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "VEC_FILE_COMPLETE")
                        {
                            VectorFileComplete = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MATRIX_FILE")
                        {
                            MatrixFile = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "EV_MIN")
                        {
                            EigenvectorCoefficientMinimum = Convert.ToDouble(inputf[u + 1]);
                            if (EigenvectorCoefficientMinimum < 0.0)
                            {
                                EigenvectorCoefficientMinimum = -1.0 * EigenvectorCoefficientMinimum;
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
                            NumberOfIts = Convert.ToInt32(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "TOL")
                        {
                            Tolerance = ParseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PARVEC")
                        {
                            ParVectorMultiplication = Convert.ToInt16(inputf[u + 1]);
                            if (ParVectorMultiplication < 1)
                            {
                                ParVectorMultiplication = 1;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PARMAT")
                        {
                            ParMatrix = Convert.ToInt16(inputf[u + 1]);
                            if (ParMatrix < 1)
                            {
                                ParMatrix = 1;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PARJ")
                        {
                            ParJ = Convert.ToInt16(inputf[u + 1]);
                            if (ParJ < 1)
                            {
                                ParJ = 1;
                            }
                        }
                        if (inputf[u].ToUpper() == "BLOCK_LANCZOS")
                        {
                            BlockLanczos = TorF(inputf[u + 1]);
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
                            FitFile = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FTOL")
                        {
                            FTol = ParseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "XTOL")
                        {
                            XTol = ParseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "GTOL")
                        {
                            GTol = ParseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "MAXFEV")
                        {
                            MaxOptimizerSteps = Convert.ToInt32(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "FACTOR")
                        {
                            Factor = ParseDouble(inputf[u + 1]);
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
                    IncludeCrossTerms = true;
                    CrossTermMatrix = new double[nModes, nModes];
                    CrossTermFit = new bool[nModes, nModes];
                    for (int l = 0; l < nModes; l++)
                    {
                        for (int m = 0; m < nModes; m++)
                        {
                            CrossTermMatrix[l, m] = 0D;
                            CrossTermFit[l, m] = false;
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
                            CrossTermMatrix[row, column] = ParseDouble(inputf[j + 5]);
                            tbool = TorF(inputf[j + 7]);
                            CrossTermFit[row, column] = tbool;
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
                    ScanList = new List<Scanner>();
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "STEPS")
                        {
                            ScanSteps = Convert.ToInt32(inputf[u + 1]);
                        }
                        if (inputf[u].ToUpper() == "MODE")
                        {
                            Scanner tempScan = new Scanner();
                            tempScan.Mode = Convert.ToInt32(inputf[u + 1]);
                            tempScan.varToFit = inputf[u + 2];
                            tempScan.Start = ParseDouble(inputf[u + 3]);
                            tempScan.Step = ParseDouble(inputf[u + 4]);
                            ScanList.Add(tempScan);
                            tempScan = null;
                            Scan = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "CROSS")
                        {
                            Scanner tempScan = new Scanner();
                            tempScan.Mode = Convert.ToInt32(inputf[u + 1]);
                            tempScan.Cross = Convert.ToInt32(inputf[u + 2]);
                            tempScan.varToFit = inputf[u + 3];
                            tempScan.Start = ParseDouble(inputf[u + 4]);
                            tempScan.Step = ParseDouble(inputf[u + 5]);
                            ScanList.Add(tempScan);
                            tempScan = null;
                            Scan = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "/")
                        {
                            break;
                        }
                    }
                    #endregion
                }

                if (inputf[i].ToUpper() == "&INTENSITY")
                {
                    #region &INTENSITY
                    Intensity = true;
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "JS")
                        {
                            JSInten = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "VECTOR_FILE")
                        {
                            VectorName = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "INDEX")
                        {
                            VectorIndex = Convert.ToInt32(inputf[u + 1]);
                        }
                        if (inputf[u].ToUpper() == "JBLOCK")
                        {
                            VectorJBlock = ParseDecimal(inputf[u + 1]);
                        }
                        if (inputf[u] == "/")
                        {
                            break;
                        }
                    }//end for
                    #endregion
                }

                if (inputf[i].ToUpper() == "&SPECIAL")
                {
                    #region &SPECIAL
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "CHECK_EIGENVECTOR")
                        {
                            CheckEigenvector = TorF(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "EIGENVECTOR_FILE")
                        {
                            EigenvectorFileName = inputf[u + 1];
                            continue;
                        }
                        if (inputf[u].ToUpper() == "INDICES")
                        {
                            JBlockEigenvector = new Tuple<decimal, int>(ParseDecimal(inputf[u + 1]), Convert.ToInt32(inputf[u + 2]));
                            continue;
                        }
                        if (inputf[u].ToUpper() == "EIGENVECTOR_TOLERANCE")
                        {
                            EigenvectorTolerance = ParseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u] == "/")
                        {
                            break;
                        }
                    }//end for
                    #endregion
                }
            }//end for

            //this is to initialize the cross-term matrix no matter what
            if (CrossTermMatrix == null)
            {
                CrossTermMatrix = new double[nModes, nModes];
                CrossTermFit = new bool[nModes, nModes];
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
        public static double ParseDouble(string s)
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
        public static decimal ParseDecimal(string s)
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