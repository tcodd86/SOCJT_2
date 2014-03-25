using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;

namespace ConsoleApplication1
{
    /// <summary>
    /// Stores the information for a mode in the input file.  All info in &amp;MODE_INFO section stored here.
    /// </summary>
    class ModeInfo
    {
        #region properties
        private double nModeOmega;
        public double modeOmega
        {
            get { return nModeOmega; }
            set 
            {
                if (value < 0.0)
                {
                    throw new InvalidInput("omega");
                }                
                nModeOmega = value; 
            }
        }//end property modeOmega

        private double nD;
        public double D
        {
            get { return nD; }
            set 
            {
                if (value < 0.0)
                {
                    throw new InvalidInput("D (linear JT constant)");
                }
                nD = value; 
            }
        }//end property D

        private double nK;
        public double K
        {
            get { return nK; }
            set { nK = value; }
        }//end property K

        private double nwExe;
        public double wExe
        {
            get { return nwExe; }
            set { nwExe = value; }
        }//end property wExe

        private int nmodeVMax;
        public int modeVMax
        {
            get { return nmodeVMax; }
            private set
            {
                if (value < 0)
                {
                    throw new InvalidInput("MODEVMAX");
                }
                nmodeVMax = value;
            }
        }//end property modeVMax

        private bool nfitOmega;
        public bool fitOmega
        {
            get { return nfitOmega; }
        }//end property nfitOmega

        private bool nfitD;
        public bool fitD
        {
            get { return nfitD; }
            set { nfitD = value; }
        }//end property nfitD

        private bool nfitK;
        public bool fitK
        {
            get { return nfitK; }
            set { nfitK = value; }
        }//end property nfitK

        private bool nfitWEXE;
        public bool fitWEXE
        {
            get { return nfitWEXE; }
        }//end property nfitWEXE

        private bool nIsAType;
        public bool IsAType
        {
            get { return nIsAType; }
        }//end property nIsAType

        private double nmodeAOmega;
        public double modeAOmega
        {
            get { return nmodeAOmega; }
        }//end property modeAOmega

        private double nmodeZeta;
        public double modeZeta
        {
            get { return nmodeZeta; }
        }//end property modeZeta
                
        public double eta { get; private set; }//end property eta

        public bool fitEta { get; private set; }//end fitEta

        public double kappa { get; private set; }//end property kappa

        public bool fitKappa { get; private set; }//end fitKappa

        #endregion properties

        /// <summary>
        /// The constructor for an InputMode object.
        /// </summary>
        /// <param name="modeN">
        /// What number mode is being initialized
        /// </param>
        /// <param name="inputF">
        /// The string array containing the parsed input file
        /// </param>
        /// <param name="tReturn">
        /// Boolean to indicate if Kappa and Eta are being used.
        /// </param>
        public ModeInfo(int modeN, string[] inputF, out bool tReturn)
        {
            int whatMode = -1;
            tReturn = false;
            for (int i = 0; i < inputF.Length; i++)
            {
                if (inputF[i] == "&MODE_INFO")
                {
                    whatMode++;
                }//end check to see if it's the right mode

                //if it's the right mode then set the values for this ModeInfo object.
                if (whatMode == modeN)
                {
                    for (int u = i; ; u++)
                    {
                        if (inputF[u] == "MODEOMEGA")
                        {
                            modeOmega = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODED")
                        {
                            D = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEK")
                        {
                            K = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEWEXE")
                        {
                            wExe = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEVMAX")
                        {
                            modeVMax = Convert.ToInt32(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEA_OMEGA")
                        {
                            nmodeAOmega = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEZETA")
                        {
                            nmodeZeta = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "FIT_OMEGA")
                        {
                            nfitOmega = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "FIT_D")
                        {
                            nfitD = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "FIT_K")
                        {
                            nfitK = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "FIT_WEXE")
                        {
                            nfitWEXE = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u].ToUpper() == "ISATYPE")
                        {
                            nIsAType = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u].ToUpper() == "KAPPA")
                        {
                            tReturn = true;
                            kappa = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u].ToUpper() == "FIT_KAPPA")
                        {
                            fitKappa = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u].ToUpper() == "ETA")
                        {
                            tReturn = true;
                            eta = FileInfo.ParseDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u].ToUpper() == "FIT_ETA")
                        {
                            fitEta = FileInfo.TorF(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "/")
                        {
                            break;
                        }
                    }//end u for loop
                }//end if
            }//end for
        }//end method setMode
    }//end class Mode
}
