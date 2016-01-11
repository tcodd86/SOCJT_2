using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathLibrary
{
    public static class MatrixFunctions
    {
        /// <summary>
        /// Normalizes the vector X
        /// </summary>
        /// <param name="X">
        /// Vector to be normalized.
        /// </param>
        public static void Normalize(double[] X)
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
        }

        /// <summary>
        /// Normalizes a collection of vectors stored as a 2D array.
        /// </summary>
        /// <param name="X">
        /// 2D array containing the vectors to be normalized
        /// </param>
        public static void Normalize(ref double[,] X)
        {
            double[] temp = new double[X.GetLength(0)];
            for (int j = 0; j < X.GetLength(1); j++)
            {
                for (int i = 0; i < X.GetLength(0); i++)
                {
                    temp[i] = X[i, j];
                }
                Normalize(temp);
                for (int i = 0; i < X.GetLength(0); i++)
                {
                    X[i, j] = temp[i];
                }
            }
        }


        /// <summary>
        /// Dot product of two vectors
        /// </summary>
        /// <param name="v">
        /// Vector 1
        /// </param>
        /// <param name="u">
        /// Vector 2
        /// </param>
        /// <returns>
        /// Scalar product of vectors v and u
        /// </returns>
        public static double DotProduct(double[] v, double[] u)
        {
            double product = 0.0;
            for (int i = 0; i < v.Length; i++)
            {
                product += v[i] * u[i];
            }
            return product;
        }

        /// <summary>
        /// Computes the projection of v onto u.
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double[] Projection(double[] u, double[] v)
        {
            double uv = MatrixFunctions.DotProduct(u, v);
            double uu = MatrixFunctions.DotProduct(u, u);
            uv /= uu;
            double[] proj = new double[u.Length];
            for (int i = 0; i < u.Length; i++)
            {
                proj[i] = uv * u[i];
            }
            return proj;
        }

        /// <summary>
        /// Calculates the magnitude of a vector
        /// </summary>
        /// <param name="vec">
        /// Vector whose magnitude is to be found
        /// </param>
        /// <returns>
        /// Magnitude of vec.
        /// </returns>
        public static double Magnitude(double[] vec)
        {
            double sum = 0.0;
            for (int i = 0; i < vec.Length; i++)
            {
                sum += vec[i] * vec[i];
            }
            return Math.Sqrt(sum);
        }
    }
}
