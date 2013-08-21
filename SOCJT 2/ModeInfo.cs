using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    /// <summary>
    /// Stores the information for a mode in the input file.  All info in &MODE_INFO section stored here.
    /// </summary>
    class ModeInfo
    {
        #region properties
        private double nModeOmega;
        public double modeOmega
        {
            get { return nModeOmega; }
            set { nModeOmega = value; }
        }//end property modeOmega

        private double nD;
        public double D
        {
            get { return nD; }
            set { nD = value; }
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
            //set { nmodeVMax = value; }
        }//end property modeVMax

        private bool nfitOmega;
        public bool fitOmega
        {
            get { return nfitOmega; }
            //set { nfitOmega = value; }
        }//end property nfitOmega

        private bool nfitD;
        public bool fitD
        {
            get { return nfitD; }
            //set { nfitD = value; }
        }//end property nfitD

        private bool nfitK;
        public bool fitK
        {
            get { return nfitK; }
            //set { nfitK = value; }
        }//end property nfitK

        private bool nfitWEXE;
        public bool fitWEXE
        {
            get { return nfitWEXE; }
            //set { nfitWEXE = value; }
        }//end property nfitWEXE

        private bool nIsAType;
        public bool IsAType
        {
            get { return nIsAType; }
            //set { nIsAType = value; }
        }//end property nIsAType

        private double nmodeAOmega;
        public double modeAOmega
        {
            get { return nmodeAOmega; }
            //set { nmodeAOmega = value; }
        }//end property modeAOmega

        private double nmodeZeta;
        public double modeZeta
        {
            get { return nmodeZeta; }
            //set { nmodeZeta = value; }
        }//end property modeZeta
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
        public ModeInfo(int modeN, string[] inputF)
        {
            int whatMode = -1;
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
                            nModeOmega = Convert.ToDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODED")
                        {
                            nD = Convert.ToDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEK")
                        {
                            nK = Convert.ToDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEWEXE")
                        {
                            nwExe = Convert.ToDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEVMAX")
                        {
                            nmodeVMax = Convert.ToInt32(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEA_OMEGA")
                        {
                            nmodeAOmega = Convert.ToDouble(inputF[u + 1]);
                            continue;
                        }
                        if (inputF[u] == "MODEZETA")
                        {
                            nmodeZeta = Convert.ToDouble(inputF[u + 1]);
                        }
                        if (inputF[u] == "FIT_OMEGA")
                        {
                            if (inputF[u + 1].ToUpper() == "T" || inputF[u + 1].ToUpper() == "TRUE")
                            {
                                nfitOmega = true;
                            }
                            else
                            {
                                nfitOmega = false;
                            }
                        }
                        if (inputF[u] == "FIT_D")
                        {
                            if (inputF[u + 1].ToUpper() == "T" || inputF[u + 1].ToUpper() == "TRUE")
                            {
                                nfitD = true;
                            }
                            else
                            {
                                nfitD = false;
                            }
                        }
                        if (inputF[u] == "FIT_K")
                        {
                            if (inputF[u + 1].ToUpper() == "T" || inputF[u + 1].ToUpper() == "TRUE")
                            {
                                nfitK = true;
                            }
                            else
                            {
                                nfitK = false;
                            }
                            continue;
                        }
                        if (inputF[u] == "FIT_WEXE")
                        {
                            if (inputF[u + 1].ToUpper() == "T" || inputF[u + 1].ToUpper() == "TRUE")
                            {
                                nfitWEXE = true;
                            }
                            else
                            {
                                nfitWEXE = false;
                            }
                            continue;
                        }
                        if (inputF[u].ToUpper() == "ISATYPE")
                        {
                            if (inputF[u + 1].ToUpper() == "T" || inputF[u + 1].ToUpper() == "TRUE")
                            {
                                nIsAType = true;
                            }
                            else
                            {
                                nIsAType = false;
                            }
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
