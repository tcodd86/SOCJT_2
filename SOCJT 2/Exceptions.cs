using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class FitFileNotFoundException : ApplicationException
    {
        public FitFileNotFoundException()
        { 
            //empty body
        }
    }

    class BasisSetTooSmallException : ApplicationException
    {
        public BasisSetTooSmallException()
        { 
            //empty body
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
}
