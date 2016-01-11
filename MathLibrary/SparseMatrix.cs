using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathLibrary
{
    public class SparseMatrix
    {
        private alglib.sparsematrix sparseMatrix;
        public bool IsCRS { get; set; }
        public int Rows
        {
            get
            {
                return sparseMatrix.innerobj.m;
            }
        }
        public int Columns
        {
            get
            {
                return sparseMatrix.innerobj.n;
            }
        }

        private void Initialize(alglib.sparsematrix matrix, bool isCRS)
        {
            sparseMatrix = matrix;
            IsCRS = isCRS;
        }

        private void Initialize()
        {
            sparseMatrix = new alglib.sparsematrix();
            IsCRS = false;
        }

        private SparseMatrix(alglib.sparsematrix matrix, bool isCRS)
        {
            Initialize(matrix, isCRS);
        }

        public SparseMatrix(int rows, int columns)
        {
            Initialize();
            alglib.sparsecreate(rows, columns, out sparseMatrix);            
        }

        public SparseMatrix(int rows, int columns, int numberOfNonzeroElements)
        {
            Initialize();
            alglib.sparsecreate(rows, columns, numberOfNonzeroElements, out sparseMatrix);
        }

        public void AddElement(int row, int column, double value)
        {
            alglib.sparseadd(sparseMatrix, row, column, value);
        }

        public void ConvertToCRS()
        {
            if (!IsCRS)
            {
                alglib.sparseconverttocrs(sparseMatrix);
                IsCRS = true;
            }
        }

        public bool EnumerateElements(ref int counter1, ref int counter2, out int row, out int column, out double value)
        {
            return alglib.sparseenumerate(sparseMatrix, ref counter1, ref counter2, out row, out column, out value);
        }

        /// <summary>
        /// Generates a copy of the current SparseMatrix object.
        /// </summary>
        /// <returns>
        /// new SparseMatrix object containing same internal values.
        /// </returns>
        public SparseMatrix Copy()
        {
            return new SparseMatrix(sparseMatrix, IsCRS);
        }
        
        /// <summary>
        /// Takes a list of upper diagonal SparseMatrix objects and combines them. Assumes element 0 is the diagonal elements.
        /// </summary>
        /// <param name="matrixList">
        /// List of matrices to be combined
        /// </param>
        /// <returns>
        /// SparseMatrix which contains all elements of all matrices in the list.
        /// </returns>
        private static SparseMatrix Aggregator(List<SparseMatrix> matrixList)
        {
            int row;
            int column;
            double value;
            var combinedMatrix = new SparseMatrix(matrixList[0].Rows, matrixList[0].Columns);

            for (int m = 0; m < matrixList.Count; m++)
            {
                int t0 = 0;
                int t1 = 0;
                while (matrixList[m].EnumerateElements(ref t0, ref t1, out row, out column, out value))
                {
                    combinedMatrix.AddElement(row, column, value);
                    //adds lower diagonal matrix elements for all off-diagonal matrices (first element in mat is the diagonal elements so don't do this for it). 
                    if (m != 0)
                    {
                        combinedMatrix.AddElement(column, row, value);
                    }
                }
            }
            return combinedMatrix;
        }

        /// <summary>
        /// Used to multiply all elements in a sparse matrix A by the value val
        /// </summary>
        /// <param name="A">
        /// Sparse matrix to be multiplied.
        /// </param>
        /// <param name="val">
        /// Value to be multiplied
        /// </param>
        /// <returns>
        /// Sparesmatrix containing the original matrix A times the double val.
        /// </returns>
        public static SparseMatrix ConstantTimesSparse(SparseMatrix A, double val)
        {
            int row;
            int column;
            double oldVal;
            int t0 = 0;
            int t1 = 0;
            
            var B = new SparseMatrix(A.Rows, A.Columns);
            while (A.EnumerateElements(ref t0, ref t1, out row, out column, out oldVal))
            {
                B.AddElement(row, column, oldVal * val);
            }

            return B;
        }

    }
}
