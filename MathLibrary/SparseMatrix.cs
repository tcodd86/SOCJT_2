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
        /// Takes a list of alglib.sparsematrix objects and combines them.
        /// </summary>
        /// <param name="mat">
        /// List of alglib.sparsematrix objects to be combined
        /// </param>
        /// <returns>
        /// alglib.sparsematrix object which contains all elements of all alglib.sparsematrix objects in the List mat.
        /// </returns>
        private static SparseMatrix Aggregator(List<SparseMatrix> mat)
        {
            int row;
            int column;
            double oldVal;
            var B = new SparseMatrix(mat[0].Rows, mat[0].Columns);

            for (int m = 0; m < mat.Count; m++)
            {
                int t0 = 0;
                int t1 = 0;
                while (mat[m].EnumerateElements(ref t0, ref t1, out row, out column, out oldVal))
                {
                    B.AddElement(row, column, oldVal);
                    alglib.sparseadd(B, row, column, oldVal);
                    //adds lower diagonal matrix elements for all off-diagonal matrices (first element in mat is the diagonal elements so don't do this for it). 
                    if (m != 0)
                    {
                        alglib.sparseadd(B, column, row, oldVal);
                    }
                }
            }
            return B;
        }

    }
}
