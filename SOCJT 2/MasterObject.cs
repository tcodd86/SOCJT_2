using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class MasterObject
    {
        public SOCJT nSoc = new SOCJT();
        public FileInfo nInput = new FileInfo();
        public string[] nInputFile;
        public bool nIsQuad;
        public List<Mode> nModes = new List<Mode>();
        public Eigenvalue[] nFitFile;

        public void Initialize(SOCJT nSocjt, FileInfo input, List<Mode> Modes, string[] inputFile, bool isQuad, Eigenvalue[] fitFile)
        {
            nSoc = nSocjt;
            nInput = input;
            nModes = Modes;
            nInputFile = inputFile;
            nIsQuad = isQuad;
            nFitFile = fitFile;
        }
    }
}
