using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ExceptionLibrary
{
    public class BasisSetTooSmallException : ApplicationException
    {
        public string EMessage { get; private set; }

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
                EMessage = "At least one of the j-blocks has a dimension smaller than M.";
                EMessage += "\r" + "Increase the basis set or decrease M.";
            }
            else
            {
                EMessage = "M is larger than noits.  Either increase noits or decrease M.";
            }
            EMessage += "\r" + "Press enter to terminate the program.";
        }
    }

    public class AEAnharmonicTermException : ApplicationException
    {
        public AEAnharmonicTermException()
        { 
            //empty body
        }
    }

    public class InvalidInput : ApplicationException
    {
        public string EMessage { get; private set; }
        public InvalidInput(string mess)
        {
            EMessage = mess + " has an invalid value. Please check your input file.";
        }
    }

    public class RepeaterError : ApplicationException
    {
        public RepeaterError()
        { 
            //empty body
        }
    }

    public class FileNameError : ApplicationException
    {
        public string EMessage { get; private set; }
        public FileNameError(string error)
        {
            if (error == "matFile")
            {
                EMessage = "The matrix file may not have the same name as the input, \n output, or fit files. Please change the file name.";
            }
            else if (error == "outFile")
            {
                EMessage = "The output file may not have the same name as the input file. \n Please change the output file name.";
            }
            else
            {
                EMessage = "There is a file name error. Please check your input file name.";
            }
            EMessage += "\r" + "Press enter to terminate the program.";
        }
    }

    public class MatrixFileError : ApplicationException
    {
        public MatrixFileError()
        {
            //empty body
        }
    }
}
