using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathLibrary
{
    public class MatrixFunctions
    {

        /// <summary>
        /// Normalizes the vector X
        /// </summary>
        /// <param name="X">
        /// Vector to be normalized.
        /// </param>
        public static void normalize(double[] X)
        {
            double sum = 0.0;
            for (int i = 0; i < X.Length; i++)
            {
                sum += X[i] * X[i];
            }
            sum = Math.Sqrt(sum);
            for (int i = 0; i < X.Length; i++)
            {
                X[i] /= sum;
            }
        }//end normalize

        /// <summary>
        /// Normalizes a collection of vectors stored as a 2D array.
        /// </summary>
        /// <param name="X">
        /// 2D array containing the vectors to be normalized
        /// </param>
        public static void normalize(ref double[,] X)
        {
            double[] temp = new double[X.GetLength(0)];
            for (int j = 0; j < X.GetLength(1); j++)
            {
                for (int i = 0; i < X.GetLength(0); i++)
                {
                    temp[i] = X[i, j];
                }
                normalize(temp);
                for (int i = 0; i < X.GetLength(0); i++)
                {
                    X[i, j] = temp[i];
                }
            }//end loop over columns
        }//end normalize

    }
}
