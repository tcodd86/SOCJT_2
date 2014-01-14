using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class BasisSetTooSmallException : ApplicationException
    {
        public string eMessage { get; private set; }

        /// <summary>
        /// Constructs BasisSetTooSmallException
        /// </summary>
        /// <param name="blockLanczos">
        /// True if blockLanczos, false if not.
        /// </param>
        public BasisSetTooSmallException(bool blockLanczos)
        {
            if (blockLanczos)
            {
                eMessage = "At least one of the j-blocks has a dimension smaller than M.";
                eMessage += "\r" + "Increase the basis set or decrease M.";
            }
            else
            {
                eMessage = "M is larger than noits.  Either increase noits or decrease M.";
            }
            eMessage += "\r" + "Press enter to terminate the program.";
        }
    }

    class AEAnharmonicTermException : ApplicationException
    {
        public AEAnharmonicTermException()
        { 
            //empty body
        }
    }

    class SpinInvalidException : ApplicationException
    {
        public SpinInvalidException()
        { 
            //empty body
        }
    }

    class RepeaterError : ApplicationException
    {
        public RepeaterError()
        { 
            //empty body
        }
    }

    class FileNameError : ApplicationException
    {
        public string eMessage { get; private set; }
        public FileNameError(string error)
        {
            if (error == "matFile")
            {
                eMessage = "The matrix file may not have the same name as the input, \n output, or fit files. Please change the file name.";
            }
            else if (error == "outFile")
            {
                eMessage = "The output file may not have the same name as the input file. \n Please change the output file name.";
            }
            else
            {
                eMessage = "There is a file name error. Please check your input file name.";
            }
            eMessage += "\r" + "Press enter to terminate the program.";
        }
    }

    class MatrixFileError : ApplicationException
    {
        public MatrixFileError()
        {
            //empty body
        }
    }
}
