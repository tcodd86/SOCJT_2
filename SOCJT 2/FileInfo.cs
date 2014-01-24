using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;

namespace ConsoleApplication1
{
    class FileInfo
    {
        #region properties

        private int nNModes;
        public int nModes
        { 
            get{ return nNModes; }
            set { nNModes = value; }
        }//end property nNModes

        private string nTitle;
        public string title
        {
            get { return nTitle; }
            set { nTitle = value; }
        }//end property nTitle

        private decimal nS;
        public decimal S
        {
            get { return nS; }
            set { nS = value; }
        }//end property S

        private double nAzeta;
        public double Azeta
        {
            get { return nAzeta; }
            set { nAzeta = value; }
        }//end property Azeta

        private bool nfitAzeta;
        public bool fitAzeta
        {
            get { return nfitAzeta; }
            set { nfitAzeta = value; }
        }//end property fitAzeta

        private decimal nmaxJ;
        public decimal maxJ
        {
            get { return nmaxJ; }
            set { nmaxJ = value; }
        }//end property maxJ

        private decimal nminJ;
        public decimal minJ
        {
            get { return nminJ; }
            set { nminJ = value; }
        }//end property minJ

        private bool nminJBool;
        public bool minJBool
        {
            get { return nminJBool; }
            set { nminJBool = value; }
        }//end property calcDeriv

        private bool nCalcDeriv;
        public bool calcDeriv
        {
            get { return nCalcDeriv; }
            set { nCalcDeriv = value; }
        }//end property calcDeriv

        private decimal nZetaE;
        public decimal zetaE
        {
            get { return nZetaE; }
            set { nZetaE = value; }
        }//end property zetaE

        private int nS1;
        public int S1
        {
            get { return nS1; }
            set { nS1 = value; }
        }//end property S1

        private int nS2;
        public int S2
        {
            get { return nS2; }
            set { nS2 = value; }
        }//end property S2

        private bool nPrintBasis;
        public bool printBasis
        {
            get { return nPrintBasis; }
            set { nPrintBasis = value; }
        }//end property printBasis

        private bool npMatrix;
        public bool pMatrix
        {
            get { return npMatrix; }
            set { npMatrix = value; }
        }//end property pMatrix

        private bool npVec;
        public bool pVector
        {
            get { return npVec; }
            set { npVec = value; }
        }//end property pVector

        private bool npMonit;
        public bool pMonit
        {
            get { return npMonit; }
            set { npMonit = value; }
        }//end property pMonit

        private bool npDer;
        public bool pDerivs
        {
            get { return npDer; }
            set { npDer = value; }
        }//end property pDerivs

        private bool nvecFile;
        public bool vecFile
        {
            get { return nvecFile; }
            set { nvecFile = value; }
        }//end property vecFile

        private double nOrigin;
        public double origin
        {
            get { return nOrigin; }
            set { nOrigin = value; }
        }//end property origin

        private bool nfitOrigin;
        public bool fitOrigin
        {
            get { return nfitOrigin; }
            set { nfitOrigin = value; }
        }//end property fitOrigin

        private bool nbasisFile;
        public bool basisFile
        {
            get { return nbasisFile; }
            set { nbasisFile = value; }
        }//end property basisFile

        private int nM;
        public int M
        {
            get { return nM; }
            set { nM = value; }
        }//end property M

        private int nkFactor;
        public int kFactor
        {
            get { return nkFactor; }
            set { nkFactor = value; }
        }//end property kFactor

        private int nnoIts;
        public int noIts
        {
            get { return nnoIts; }
            set { nnoIts = value; }
        }//end property noIts

        private double ntol;
        public double tol
        {
            get { return ntol; }
            set { ntol = value; }
        }//end property tol

        private bool nguesses;
        public bool guesses
        {
            get { return nguesses; }
            set { nguesses = value; }
        }//end property guesses

        private string nfitFile;
        public string fitFile
        {
            get { return nfitFile; }
            set { nfitFile = value; }
        }//end property fitFile

        private double nfTol;
        public double fTol
        {
            get { return nfTol; }
            set { nfTol = value; }
        }//end property fTol

        private double nxTol;
        public double xTol
        {
            get { return nxTol; }
            set { nxTol = value; }
        }//end property xTol

        private double ngTol;
        public double gTol
        {
            get { return ngTol; }
            set { ngTol = value; }
        }//end property gTol
        
        private int nmaxFev;
        public int maxFev
        {
            get { return nmaxFev; }
            set { nmaxFev = value; }
        }//end property maxFev

        private double nfactor;
        public double factor
        {
            get { return nfactor; }
            set { nfactor = value; }
        }//end property factor

        private int nprint;
        public int print
        {
            get { return nprint; }
            set { nprint = value; }
        }//end property print
        
        private bool nIncludeCrossTerms;
        public bool includeCrossTerms
        {
            get { return nIncludeCrossTerms; }
            set { nIncludeCrossTerms = value; }
        }//end property includeCrossTerms

        private bool nAT;
        public bool AT
        {
            get { return nAT; }
            set { nAT = value; }
        }//end property AT

        private bool nSpecial;
        public bool Special
        {
            get { return nSpecial; }
            set { nSpecial = value; }
        }//end property Special

        private bool nScan;
        public bool Scan
        {
            get { return nScan; }
            set { nScan = value; }
        }//end property Scan

        private int nSteps;
        public int Steps
        {
            get { return nSteps; }
            set { nSteps = value; }
        }//end property Steps

        private double nmatGenTime;
        public double matGenTime
        {
            get { return nmatGenTime; }
            set { nmatGenTime = value; }
        }//end property matGenTime

        private double ndiagTime;
        public double diagTime
        {
            get { return ndiagTime; }
            set { ndiagTime = value; }
        }//end property diagTime

        private int nPar;
        public int parMat
        {
            get { return nPar; }
            set { nPar = value; }
        }//end property par

        private int nParVec;
        public int parVec
        {
            get { return nParVec; }
            set { nParVec = value; }
        }//end property pVec

        private int nParJ;
        public int parJ
        {
            get { return nParJ; }
            set { nParJ = value; }
        }//end property parJ

        private double nEvMin;
        public double evMin
        {
            get { return nEvMin; }
            set { nEvMin = value; }
        }//end property evMin

        private bool nInclSO;
        public bool inclSO
        {
            get { return nInclSO; }
            set { nInclSO = value; }
        }//end property inclSO

        private List<Tuple<decimal, int, int>> neVecs;
        public List<Tuple<decimal, int, int>> eVecs
        {
            get { return neVecs; }
            set { neVecs = value; }
        }//end property eVecs

        private bool nbeVecs;
        public bool beVecs
        {
            get { return nbeVecs; }
            set { nbeVecs = value; }
        }//end property nbeVecs

        public bool useKappaEta { get; set; }//end useKappaEta

        public List<Scanner> scanList;

        public double[,] crossTermMatrix;

        public bool[,] crossTermFit;

        public bool blockLanczos { get; private set; }//end naiveLanczos

        public bool oldRandom { get; set; }//end newRandom

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

        #endregion properties

        /// <summary>
        /// Constructor for FileInfo object. Sets reasonable values for most properties in case user forgets something in input file.
        /// </summary>
        public FileInfo()
        {
            //initialize all booleans to false by default            
            fitAzeta = false;
            calcDeriv = false;
            minJBool = false;
            fitOrigin = false;
            inclSO = false;
            printBasis = false;
            pMatrix = false;
            pVector = false;
            AT = false;
            Special = false;
            includeCrossTerms = false;
            beVecs = false;
            useKappaEta = false;
            blockLanczos = false;
            oldRandom = false;
            matFile = "matrix.txt";
            useMatFile = false;
            matMade = false;

            title = "TITLE";
            origin = 0.0;
            parMat = 1;
            nModes = 1;
            S = 0.5M;
            Azeta = 0.0;
            zetaE = 0.0M;

            //these are reasonable values of J for a basic quadratic problem
            maxJ = 7.5M;
            minJ = 0.5M;            

            //S1 and S2 values for 3-fold symmetry
            S1 = 0;
            S2 = 1;
                        
            evMin = 0.2;            
            parVec = 1;
            parMat = 1;
            parJ = 2;
            M = 5;
            kFactor = 2;
            noIts = 10000;
            tol = 0.0001;
            fitFile = "fit.fit";
            fTol = 0.0;
            xTol = 0.0;
            gTol = 0.0;
            maxFev = 25;
            factor = 0.001;
            nprint = 0;            
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

        public void setFileInfo(string[] inputf)
        {
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
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                fitAzeta = true;
                            }
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
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                fitOrigin = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "ORIGIN")
                        {
                            origin = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "CALC_DERIV")
                        {
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                calcDeriv = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "ZETAE")
                        {
                            zetaE = parseDecimal(inputf[u + 1]);
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
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                useKappaEta = true;
                            }
                            else
                            {
                                useKappaEta = false;
                            }
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
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                printBasis = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PRINT_MATRIX")
                        {
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                pMatrix = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "PRINT_VEC")
                        {
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                pVector = true;
                            }
                            continue;
                        }
                        if (inputf[u].ToUpper() == "USE_MATRIX_FILE")
                        {
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                useMatFile = true;
                            }
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
                            //tol = Convert.ToDouble(inputf[u + 1]);
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
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                blockLanczos = true;
                            }
                        }
                        if (inputf[u].ToUpper() == "OLD_RANDOM")
                        {
                            if (inputf[u + 1].ToUpper() == "T" || inputf[u + 1].ToUpper() == "TRUE")
                            {
                                oldRandom = true;
                            }
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
                            //fTol = Convert.ToDouble(inputf[u + 1]);
                            fTol = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "XTOL")
                        {
                            //xTol = Convert.ToDouble(inputf[u + 1]);
                            xTol = parseDouble(inputf[u + 1]);
                            continue;
                        }
                        if (inputf[u].ToUpper() == "GTOL")
                        {
                            //gTol = Convert.ToDouble(inputf[u + 1]);
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
                            //factor = Convert.ToDouble(inputf[u + 1]);
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
                            //crossTermMatrix[row, column] = Convert.ToDouble(inputf[j + 5]);
                            crossTermMatrix[row, column] = parseDouble(inputf[j + 5]);
                            if (inputf[j + 7].ToUpper() == "T" || inputf[j + 7].ToUpper() == "TRUE")
                            {
                                tbool = true;
                            }//end if
                            crossTermFit[row, column] = tbool;
                            //j += 8;
                            continue;
                        }

                        if (inputf[j].ToUpper() == "AT")
                        {
                            AT = true;
                            tbool = false;
                            row = Convert.ToInt16(inputf[j + 2]) - 1;
                            column = Convert.ToInt16(inputf[j + 4]) - 1;
                            if (row < column)
                            {
                                temp = row;
                                row = column;
                                column = temp;
                            }//end if
                            //crossTermMatrix[row, column] = Convert.ToDouble(inputf[j + 5]);
                            crossTermMatrix[row, column] = parseDouble(inputf[j + 5]);
                            if (inputf[j + 7].ToUpper() == "T" || inputf[j + 7].ToUpper() == "TRUE")
                            {
                                tbool = true;
                            }//end if
                            crossTermFit[row, column] = tbool;
                            //j += 8;
                            continue;
                        }

                        if (inputf[j].ToUpper() == "SPECIAL")
                        {
                            tbool = false;
                            Special = true;
                            //crossTermMatrix[0, 0] = Convert.ToDouble(inputf[j + 1]);
                            crossTermMatrix[0, 0] = parseDouble(inputf[j + 1]);
                            if (inputf[j + 3].ToUpper() == "T" || inputf[j + 3].ToUpper() == "TRUE")
                            {
                                tbool = true;
                            }
                            crossTermFit[0, 0] = tbool;
                        }
                        if (inputf[j] == "/")
                        {
                            break;
                        }//end end of CROSS_TERM if
                    }//end CROSS_TERM for
                    #endregion
                }//end CROSS_TERM if

                if (inputf[i].ToUpper() == "&COMPARE_EIGENVECTORS")
                {
                    #region &COMPARE_EIGENVECTORS
                    beVecs = true;
                    eVecs = new List<Tuple<decimal,int,int>>();
                    for (int u = i + 1; ; u += 3)
                    {
                        eVecs.Add(new Tuple<decimal, int, int>(parseDecimal(inputf[u]), Convert.ToInt16(inputf[u + 1]),  Convert.ToInt16(inputf[u + 2])));
                        if (inputf[u + 3].ToUpper() == "/")
                        {
                            break;
                        }
                    }
                    #endregion
                }//end COMPARE_EIGENVECTORS if

                if (inputf[i].ToUpper() == "&SCAN")
                {
                    #region &SCAN
                    scanList = new List<Scanner>();
                    for (int u = i; ; u++)
                    {
                        if (inputf[u].ToUpper() == "STEPS")
                        {
                            nSteps = Convert.ToInt32(inputf[u + 1]);
                        }
                        if (inputf[u].ToUpper() == "MODE")
                        {
                            Scanner tempScan = new Scanner();
                            tempScan.Mode = Convert.ToInt32(inputf[u + 1]);
                            tempScan.varToFit = inputf[u + 2];
                            //tempScan.Start = Convert.ToDouble(inputf[u + 3]);
                            tempScan.Start = parseDouble(inputf[u + 3]);
                            //tempScan.Step = Convert.ToDouble(inputf[u + 4]);
                            tempScan.Step = parseDouble(inputf[u + 4]);
                            scanList.Add(tempScan);
                            tempScan = null;
                            nScan = true;
                            continue;
                        }
                        if (inputf[u].ToUpper() == "CROSS")
                        {
                            Scanner tempScan = new Scanner();
                            tempScan.Mode = Convert.ToInt32(inputf[u + 1]);
                            tempScan.Cross = Convert.ToInt32(inputf[u + 2]);
                            tempScan.varToFit = inputf[u + 3];
                            //tempScan.Start = Convert.ToDouble(inputf[u + 4]);
                            tempScan.Start = parseDouble(inputf[u + 4]);
                            //tempScan.Step = Convert.ToDouble(inputf[u + 5]);
                            tempScan.Step = parseDouble(inputf[u + 5]);
                            scanList.Add(tempScan);
                            tempScan = null;
                            nScan = true;
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
            return Double.Parse(s, NumberStyles.AllowExponent | NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint);
        }

        public static decimal parseDecimal(string s)
        {
            return Decimal.Parse(s, NumberStyles.AllowExponent | NumberStyles.AllowLeadingSign | NumberStyles.AllowDecimalPoint);
        }
    }//class FileInfo
}//end namespace
