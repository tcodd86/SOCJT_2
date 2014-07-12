using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

namespace ConsoleApplication1
{
    /// <summary>
    /// This class contains the functions and subroutines for the block Lanczos routine.
    /// This code was translated from FORTRAN by Terrance J. Codd in 2012.  The original
    /// FORTRAN routine (which was translated from ALGOL) was taken from the source code 
    /// for SOCJT (Barckholtz, T. Miller, T. Int. Rev. Phys. Chem., 1998, Vol. 17, No. 4, 435-524)
    /// and no authorship information was provided in the source code or documentation.
    /// 
    /// Function summaries were taken directly from the SOCJT source code while comments
    /// throughout the code are my own.  There are two key differences between this code
    /// and the original: 1. arrays are indexed to 0 in C# and 1 in FORTRAN so the loop
    /// bounds are often shifted by one as are some variables; 2. I have removed most of
    /// the goto statements found in the original to increase the readability of the code
    /// and replaced them by using conditional blocks of code.  I have labeled the for loops
    /// using the same numbering found in do loops in SOCJT for easy comparison.
    /// 
    /// Where the FORTRAN used functions TRED and TRED2 I have used functions from ALGLIB
    /// PROJECT by Bochkanov, Sergey.
    /// </summary>
    static class Lanczos
    {
        public static int basisSetLimit = 100000000;//40000 * 2000;

        /// <summary>
        /// This subroutine implements the block lanczos method with reorthogonalization.
        /// BKLANC computes a block tridiagonal matrix MS which it stores in rows and 
        /// columns M to M + P * S of the array C, and an orthonormal matrix XS which it
        /// stores in columns M to M + P * S of the N by Q array X.  MS is a symmetric
        /// matrix with P by P symmetric matrices M(1)....M(S) on its diagonal and P by P
        /// upper triangular matries R(2)......R(S) along its' lower diagonal.  Since MS
        /// is symmetric and banded, only its lower triangle (P + 1 diagonals) is stored
        /// in C.  XS is composed of S N by P orthonormal matrices X(1)......X(S) where
        /// X(1) is given and should be stored in columns M to M + P of X.  Furthermore,
        /// X(1) is assumed to satisfy X(1)*A*X(1) = diag(D(M), D(M + 1)....,D(M + P - 1))
        /// and if M > 0 then X(1) is assumed to be orthogonal to the vectors stored in
        /// columns 0 to M of X.  OP is the name of an external subroutine used to define
        /// the matrix A.  (OP was replaced with inline calls to an ALGLIB routine for 
        /// sparse symmetric matrix vector multiplication.  OP is not needed anymore.)
        /// During the first step, the subroutine ERR s called and the 
        /// quantities EJ are computed where EJ = ||A*X1J-D(M+J)*X1J||, X1J is the J-th
        /// column of X(1), and ||*|| denotes the Euclidean norm.  EJ is stored in the 
        /// E(M + J) J = 1, 2,......, P. U and V are auxilliary vectors used by OP.
        /// 
        /// The last two parameters, alglib.sparsematrix A and int par were added to the 
        /// C# version of the code.  In the original FORTRAN version the sparsematrix was
        /// stored in memory as a global variable while in C# it is passed as a parameter.
        /// The int par is a variable added for cases of chunking the matrix vector
        /// multiplications for parallelization.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="Q"></param>
        /// <param name="M"></param>
        /// <param name="P"></param>
        /// <param name="S"></param>
        /// <param name="D"></param>
        /// <param name="C"></param>
        /// <param name="X"></param>
        /// <param name="E"></param>
        /// <param name="U"></param>
        /// <param name="V"></param>
        /// <param name="A"></param>
        /// <param name="par"></param>
        private static void BKLANC(int N, int Q, int M, int P, int S, double[] D, ref double[,] C, ref double[,] X, ref double[] E, ref double[] U, ref double[] V, alglib.sparsematrix A, int par)
        {
            int LL;
            int LU;
            double T;
            int IT;
            int KP1;
            int IL;
#if DEBUG
            System.Diagnostics.Stopwatch orthogTimer = new System.Diagnostics.Stopwatch();
#endif            
            for (int L = 0; L < S; L++)//do 90
            {
                LL = M + (L) * P;
                LU = M + (L + 1) * P;
                for (int K = LL; K < LU; K++)//do 70
                {
                    for (int I = 0; I < N; I++)//do 10
                    {
                        U[I] = X[I, K];
                    }//do 10

                    OP(A, U, ref V, par);
                    
                    if (L == 0)
                    {
                        for (int I = K; I < LU; I++)//do 12
                        {
                            C[I, K] = 0.0;
                        }//do 12
                        C[K, K] = D[K];
                        for (int I = 0; I < N; I++)//do 14
                        {
                            V[I] -= D[K] * X[I, K];
                        }//do 14
                    }//end if 19 
                    else
                    {
                        for (int I = K; I < LU; I++)//do 30
                        {
                            T = 0.0;
                            for (int J = 0; J < N; J++)//do 20
                            {
                                T += V[J] * X[J, I];
                            }//do 20
                            C[I, K] = T;
                        }//do 30
                        IT = K - P;
                        for (int I = 0; I < N; I++)//do 60
                        {
                            T = 0.0;
                            for (int J = IT; J <= K; J++)//do 40
                            {
                                T += X[I, J] * C[K, J];
                            }//do 40
                            if (K != LU - 1)
                            {
                                KP1 = K + 1;
                                for (int J = KP1; J < LU; J++)//do 50
                                {
                                    T += X[I, J] * C[J, K];
                                }//do 50
                            }
                            V[I] -= T;
                        }//do 60
                    }//end else
                    if (L != S - 1)
                    {
                        for (int I = 0; I < N; I++)//do 63
                        {
                            X[I, K + P] = V[I];
                        }//do 63
                    }
                }//do 70
                if (L == 0)
                {
                    ERR(N, Q, M, P, X, ref E);
                }
                if (L < S - 1)
                {
                    ORTHG(N, Q, LU, P, ref C, ref X);
                    IL = LU;
                    int Ir = LU - 1;
                    for (int J = 0; J < P; J++)
                    {
                        Ir++;
                        for (int I = IL; I <= Ir; I++)
                        {
                            C[I, Ir - P] = C[I, Ir];
                        }
                    }
                }
            }//do 90
        }

        /// <summary>
        /// CNVTST determines which of the P eigenvalues stored in elements M to M + P
        /// of D have converged.  ERRC is a measure of the accumulated error in the M
        /// previously computed eigenvalues and eigenvectors. ERRC is updated if more
        /// approximations have converged.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="Q"></param>
        /// <param name="M"></param>
        /// <param name="P"></param>
        /// <param name="ERRC"></param>
        /// <param name="EPS"></param>
        /// <param name="D"></param>
        /// <param name="E"></param>
        /// <param name="NCONV"></param>
        private static void CNVTST(int N, int Q, int M, int P, ref double ERRC, double EPS, double[] D, double[] E, ref int NCONV)
        {
            //double CHEPS = 2.22e-16;  This is the value in the original FORTRAN
            //This should be machine epsilon.  I've estimated it here but it could just be calculated ahead of time and passed as a parameter.
            double CHEPS = 1.11e-16;
            int K = 0;
            double T;
            for (int I = 0; I < P; I++)
            {
                T = Math.Abs(D[M + I]);
                if (T < 1.0)
                {
                    T = 1.0;
                }
                if (E[M + I] > (T * (EPS + 10.0 * (double)N * CHEPS) + ERRC))
                {
                    break;
                }
                K++;
            }
            NCONV = K;
            if (K != 0)
            {
                T = 0.0;
                for (int I = 0; I < K; I++)
                {
                    T += Math.Pow(E[M + I], (double)2);
                }
                ERRC = Math.Sqrt(ERRC * ERRC + T);
            }
        }//end method CNVTST

        /// <summary>
        /// EIGEN solves the eigenproblem for the symmetric matrix MS of order PS stored
        /// in rows and columns M to M + PS of C.  The eigenvalues of MS are stored in 
        /// elements M to M + PS of D and the eigenvectors are stored in rows and columns
        /// 0 to PS of C possibly overwriting MS.
        /// </summary>
        /// <param name="Q"></param>
        /// <param name="M"></param>
        /// <param name="P"></param>
        /// <param name="PS"></param>
        /// <param name="C"></param>
        /// <param name="D"></param>
        private static void EIGEN(int Q, int M, int P, int PS, ref double[,] C, ref double[] D)
        {
            double[,] tempMat = new double[PS, PS];
            int LIM;

            for (int I = 0; I < PS; I++)
            {
                LIM = I - P;
                int LM1;
                if (I <= P)
                {
                    LIM = 0;
                }
                if (LIM > 0)
                {
                    LM1 = LIM - 1;
                    for (int J = 0; J <= LM1; J++)
                    {
                        C[I, J] = 0.0;
                    }
                }
                for (int J = LIM; J <= I; J++)
                {
                    C[I, J] = C[I + M, J + M];
                }
            }

            for (int i = 0; i < PS; i++)
            {
                for (int j = 0; j < PS; j++)
                {
                    tempMat[i, j] = C[i, j];
                }
            }

            double[] evs = new double[PS];
            double[,] z = new double[PS, PS];

            //solve eigenproblem using alglib library
            alglib.evd.smatrixevd(tempMat, PS, 1, false, ref evs, ref z);         
            //put results back into the expected format for MINVAL routine
            for (int i = 0; i < PS; i++)
            {
                D[i + M] = evs[i];
            }
            for (int i = 0; i < PS; i++)
            {
                for (int j = 0; j < PS; j++)
                {
                    C[i, j] = z[i, j];
                }
            }
        }//end EIGEN

        /// <summary>
        /// ERR computes the Euclidean lengths of the vectors stored in the columns
        /// M + P + 1 through M + P + P of the N by Q array X and stores them in
        /// elements M + 1 through M + P of E.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="Q"></param>
        /// <param name="M"></param>
        /// <param name="P"></param>
        /// <param name="X"></param>
        /// <param name="E"></param>
        private static void ERR(int N, int Q, int M, int P, double[,] X, ref double[] E)
        {
            int MP1 = M + P;
            int MPP = M + P + P;
            double T;
            for (int K = MP1; K < MPP; K++)
            {
                T = 0.0;
                for (int I = 0; I < N; I++)
                {
                    T += X[I, K] * X[I, K];
                }
                E[K - P] = Math.Sqrt(T);
            }
        }
        //ERR done new

        /// <summary>
        /// This is the function called to perform the matrix vector multiplications.
        /// I've here used functions from the ALGLIB library which is available at
        /// ALGLIB.net and has both C# and C++ implementations.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="U"></param>
        /// <param name="V"></param>
        /// <param name="par"></param>
        private static void OP(alglib.sparsematrix A, double[] U, ref double[] V, int par)
        {
            //alglib.sparsesmv(A, true, U, ref V);
            
            //if par == 1 then no parallelization for M*v products, use regular.  Else use parallel version.
            if (par == 1)
            {
                alglib.sparsemv(A, U, ref V);
            }
            else
            {
                alglib.sparsemvTC(A, U, ref V, par);
            }
            
        }


        /// <summary>
        /// ORTHG reorthogonalizes the N by P matrix Z stored in columns F to F + P of the
        /// N by Q array X with respect to the vectors stored in the first F columns of X
        /// and then decomposes the resulting matrix into the product of an N by P
        /// orthonormal matrix XORTH and and P by P upper triangular matrix R.  XORTH is
        /// stored over Z and the upper triangle of R is stored in rows and columns F 
        /// to F + P of the Q by Q array B.  A stable variant of the Gram-Schmidt
        /// orthogonalization method is utilised.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="Q"></param>
        /// <param name="F"></param>
        /// <param name="P"></param>
        /// <param name="B"></param>
        /// <param name="X"></param>
        private static void ORTHG(int N, int Q, int F, int P, ref double[,] B, ref double[,] X)
        {
            bool orig;
            int FP1;
            int FPP;
            int KM1;
            double T;
            double S;

            if (P != 0)
            {
                FP1 = F;
                FPP = F + P;
                for (int K = FP1; K < FPP; K++)//do 50
                {
                    orig = true;
                    KM1 = K;

                    S = 0.0;
                    T = 0.0;
                    while (S <= T / 100.0)
                    {
                        T = 0.0;
                        if (KM1 >= 1)
                        {
                            for (int I = 0; I < KM1; I++)//do 20
                            {
                                S = 0.0;
                                for (int J = 0; J < N; J++)//do 15
                                {
                                    S += X[J, I] * X[J, K];
                                }//do 15
                                if (orig == true & I >= F)
                                {
                                    B[I, K] = S;
                                }
                                T += S * S;
                                for (int J = 0; J < N; J++)//do 20
                                {
                                    X[J, K] -= S * X[J, I];
                                }//do 20
                            }//end do 20
                        }
                        S = 0.0;
                        for (int J = 0; J < N; J++)//do 30                    
                        {
                            S += X[J, K] * X[J, K];
                        }//do 30                    
                        T += S;
                        if (S <= T / 100.0)
                        {
                            orig = false;
                        }
                    }

                    S = Math.Sqrt(S);
                    B[K, K] = S;
                    if (S != 0)
                    {
                        S = 1.0 / S;
                    }
                    for (int J = 0; J < N; J++)//do 50
                    {
                        X[J, K] = S * X[J, K];
                    }//do 50
                }//end do 50                
            }//end if P != 0;            
        }//end ORTHG

        /// <summary>
        /// Based on the values of N, Q, M, R and NCONV, PCH chooses new values for P and S,
        /// the block size and number of steps for the block Lanczos method.  The strategy
        /// used here is to choose P to be the smaller of the two following values:
        /// 1. The previous block size
        /// 2. The number of values left to be computed.
        /// S is chosen as large as possible subject to the constraints imposed by the limits
        /// of storage.  In any event, S is greater than or equal to 2.  N is the order
        /// of the problem and Q is the number of vectors available for storing eigenvectors
        /// and applying the block lanczos method.  M is the number of eigenvalues and 
        /// eigenvectors that have already been computed and R is the required number.
        /// Finally, NCONV is the number of eigenvalues and eigenvectors that have 
        /// converged in the current iteration.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="Q"></param>
        /// <param name="M"></param>
        /// <param name="R"></param>
        /// <param name="NCONV"></param>
        /// <param name="P"></param>
        /// <param name="S"></param>
        private static void PCH(int N, int Q, int M, int R, int NCONV, ref int P, ref int S)
        {
            int PT;
            int ST;
            int MT;
            MT = M + NCONV;
            PT = R - MT;
            if (PT > P)
            {
                PT = P;
            }
            if (PT > 0)
            {
                ST = (Q - MT) / PT;
                if (ST <= 2)
                {
                    ST = 2;
                    PT = (Q - MT) / 2;
                }
                P = PT;
                S = ST;
            }
            else
            {
                P = 0;
            }
        }
        
        /// <summary>
        /// Generates random numbers and puts them in column L of X
        /// </summary>
        /// <param name="N">
        /// Number of rows in X
        /// </param>
        /// <param name="L">
        /// The column of X to be filled with random numbers
        /// </param>
        /// <param name="X">
        /// Matrix to hold the random numbers in column L
        /// </param>
        private static void RANDOM(int N, int L, ref double[,] X)
        {
            double[] T = new double[100];
            int X1;
            int F1 = 71416;
            int F2 = 27183;
            int FT;
            int K;
            int A = 6821;
            int C = 5327;
            int X0 = 5328;
            for (int I = 0; I < 100; I++)//do 100
            {
                X1 = A * X0 + C;
                if (X1 >= 10000)
                {
                    X1 -= 10000;
                }
                T[I] = (double)X1 / 9999.0 - 0.5;
                X0 = X1;
            }//do 100
            for (int I = 0; I < N; I++)//do 200
            {
                FT = F1 + F2;
                if (FT >= 1000000)
                {
                    FT -= 1000000;
                }
                F1 = F2;
                F2 = FT;
                K = FT / 10000;
                X[I, L] = T[K];
                X1 = A * X0 + C;
                if (X1 >= 10000)
                {
                    X1 -= 10000;
                }
                T[K] = (double)X1 / 9999.0 - 0.5;
                X0 = X1;
            }//do 200
            T = null;
        }//end function Random for BlockLanczos Routine

        /// <summary>
        /// Function to generate a normalized vector of length N filled with random numbers
        /// </summary>
        /// <param name="N">
        /// Length of the vector to be returned.
        /// </param>
        /// <returns>
        /// Normalized, random vector of length N.
        /// </returns>
        private static double[] RANDOM(int N)
        {
            var X = new double[N];            
            Random randy = new Random(6821);
            for (int i = 0; i < N; i++)
            {
                X[i] = randy.NextDouble();
            }
            normalize(X);
            return X;
        }

        /// <summary>
        /// Rotate computes the first L columns of the matrix XS*QS where XS is an
        /// N by PS orthonormal matrix stored in columns M + 1 through M + PS of the
        /// N by Q array X and QS is a PS by PS orthonormal matrix stored in rows and
        /// columns 0 to PS of the array C.  The result is stored in columns M + 1
        /// through M + L of X overwriting part of XS.
        /// </summary>
        /// <param name="N">Number of rows in X.</param>
        /// <param name="Q">Number of columns in X.</param>
        /// <param name="M">Where XS is stored in X.</param>
        /// <param name="PS">Dim of array QS.</param>
        /// <param name="L"></param>
        /// <param name="C">Array to contain QS.</param>
        /// <param name="X">Starting array.</param>
        private static void ROTATE(int N, int Q, int M, int PS, int L, double[,] C, ref double[,] X)
        {
            double[] V = new double[Q];

            for (int I = 0; I < N; I++)
            {
                for (int K = 0; K < L; K++)
                {
                    double T = 0.0;
                    for (int J = 0; J < PS; J++)
                    {
                        T += X[I, M + J] * C[J, K];
                    }
                    V[K] = T;
                }
                for (int K = 0; K < L; K++)
                {
                    X[I, M + K] = V[K];
                }
            }
        }//end Rotate

        /// <summary>
        /// SECTN transforms the N by P orthonormal matrix X1, say, stored in columns M 
        /// to M + p of the N by Q array X so that X1'*A*X1 = diag(D1, D2,...,DP), where
        /// ' denotes transpose and A is a symmetric matrix of order N defined by the
        /// subroutine OP.  The values D1,...,DP are stored in elements M to M + p of D.
        /// SECTN forms the matrix X1'*A*X1 = CP, storing CP in the array for C.  The
        /// values D1, D2,...,DP and the eigenvectors QP of CP are computed by EIGEN
        /// and stored in D and C respectively.  ROTATE then carries out the
        /// transformation X1*QP to X1.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="Q"></param>
        /// <param name="M"></param>
        /// <param name="P"></param>
        /// <param name="X"></param>
        /// <param name="C"></param>
        /// <param name="D"></param>
        /// <param name="U"></param>
        /// <param name="V"></param>
        private static void SECTN(int N, int Q, int M, int P, ref double[,] X, ref double[,] C, ref double[] D, ref double[] U, ref double[] V, alglib.sparsematrix A, int par)
        {
            int ICOL1 = M - 1;
            int ICOL2;
            double T;
            for (int J = 0; J < P; J++)//do 300
            {
                ICOL1++;
                for (int I = 0; I < N; I++)//do 100
                {
                    U[I] = X[I, ICOL1];
                }//do 100

                OP(A, U, ref V, par);

                ICOL2 = M - 1;

                for (int I = 0; I <= J; I++)//do 300
                {
                    ICOL2++;
                    T = 0.0;
                    for (int K = 0; K < N; K++)//do 200
                    {
                        T += V[K] * X[K, ICOL2];
                    }//do 200
                    C[ICOL1, ICOL2] = T;
                }//do 300
            }//do 300
            EIGEN(Q, M, P, P, ref C, ref D);
            ROTATE(N, Q, M, P, P, C, ref X);
        }//end method SECTN

        /// <summary>
        /// This subroutine is the main subroutine implementing the iterative block lanczos
        /// method for computing the eigenvalues and eigenvectors of symmetric matrices.
        /// </summary>
        /// <param name="N">Order of the symmetric matrix A whose eigenvalues and 
        /// eigenvectors are being computed.</param>
        /// <param name="Q">The number of vectors of length N contained in the array
        /// X.  The value of Q should be less than or equal to 25, at least one
        /// greater than R and less than or equal to N.</param>
        /// <param name="PINIT">The initial block size to be used in the block 
        /// lanczos method.</param>
        /// <param name="R">The number of eigenvalues and eigenvectors being computed.</param>
        /// <param name="MMAX">The maximum number of matrix vector products A*X allowed
        /// in one call of this subroutine.</param>
        /// <param name="EPS">The relative precision to which MINVAL will attempt to 
        /// compute the eigenvalues and eigenvectors of A.</param>
        /// <param name="M">The number of eigenvalues and eigenvectors already computed.  This
        /// should initially be zero.</param>
        /// <param name="D">D is a vector containing the computed eigenvalues.  Should have
        /// at least Q elements.</param>
        /// <param name="X">X contains the computed eigenvectors.  X should be N x Q
        /// array.  Also used as working storage while computing.</param>
        /// <param name="IECODE">
        /// Status code to return.
        /// </param>
        /// <param name="A">A is the sparse matrix being diagonalized.</param>
        public static int MINVAL(int N, int Q, int PINIT, int R, int MMAX, double EPS, int M, ref double[] D, ref double[,] X, ref int IECODE, alglib.sparsematrix A, int par)
        {
            double[] E = new double[Q];
            double[,] C = new double[Q, Q];
            double[] U = new double[N];
            double[] V = new double[N];
            int PS;
            int ITER = 0;
            int NCONV = 0;
            double ERRC = 0.0;
            int IMM;
            int IMMURN = 0;

            //this checks that the values are in the proper ranges. If they're not, then skip the lanczos routine.
            bool go = true;
            if (N < 2)
            {
                go = false;
                IECODE = 1;
            }
            if (R < 1)
            {
                go = false;
                IECODE = 3;
            }
            if (Q <= R)
            {
                go = false;
                IECODE = 4;
            }
            if (Q > N)
            {
                go = false;
                IECODE = 6;
            }

            if (go)
            {
                //chooses initial values for block size P, the number of steps that the block
                //lanczos method is carried out, and chooses an initial N by P orthonormal
                //matrix X1 used to start the block lanczos method
                int P = PINIT;
                if (P < 0)
                {
                    P = -P;
                }
                int S = (Q - M) / P;

                if (S <= 2)
                {
                    S = 2;
                    P = Q / 2;
                }

                if (PINIT >= 0)
                {
                    for (int K = 0; K < P; K++)//check what RANDOM does to see what ought to go here            
                    {
                        RANDOM(N, K, ref X);
                    }
                }

                if (M == 0)
                {
                    ORTHG(N, Q, M, P, ref C, ref X);
                    SECTN(N, Q, M, P, ref X, ref C, ref D, ref U, ref V, A, par);
                    ERRC = 0.0D;
                }

                IMM = 0;

                //The main body of the subroutine starts here.  IMM counts the number of 
                //matrix-vector products computed which is the number of times the 
                //subroutine named by OP is called.  ERRC measures the accumulated
                //error in the eigenvalues and eigenvectors.

                while (M < R)
                {
                    if (IMM > MMAX)
                    {
                        break;
                    }

                    //add a conditional here to be able to kill it after a certain number of Lanczos iterations

                    ITER++;
                    PS = P * S;

                    //BKLANC carries out the block lanczos method and stores the resulting
                    //block tridiagonal matrix MS in C and the N by PS orthonormal matrix
                    //XS in X.  The initial N by P orthonormal matrix is assumed to be
                    //stored in columns M to M + PS of X.  The residuals for these
                    //vectors and the eigenvalue approximations in D are computed and stored
                    //in E.
                    BKLANC(N, Q, M, P, S, D, ref C, ref X, ref E, ref U, ref V, A, par);

                    //Eigen solves the eigenproblem for MS, storing the eigenvalues in
                    //elements M to M + PS of D and the eigenvectors in the first
                    //P * S rows and columns of C (overwriting MS, possibly)
                    EIGEN(Q, M, P, PS, ref C, ref D);

                    //CNVTST determines how many of the eigenvalues and eigenvectors have
                    //converged using the error estimates stored in E.  The number that have
                    //converged is stoered in NCONV.  If NCONV = 0 then none have converged.
                    CNVTST(N, Q, M, P, ref ERRC, EPS, D, E, ref NCONV);

                    //PCH chooses new values for P and S, the block size and the number of 
                    //steps for the block Lanczos subprogram, respectively.
                    PCH(N, Q, M, R, NCONV, ref P, ref S);

                    //Rotate computes the eigenvectors of the restricted matrix using XS
                    //stored in X and the eigenvectors of MS stored in C.  These vectors
                    //serve both as eigenvector approximations and to form the matrix used
                    //to start the block Lanczos method in the next iteration.
                    ROTATE(N, Q, M, PS, NCONV + P, C, ref X);//see if any of this should be ref.

                    M += NCONV;
                    IMM += P * S;
                }//end while

                if (M >= R)
                {
                    IECODE = 0;
                }
                if (IMM > MMAX)
                {
                    IECODE = 7;
                    PINIT = -P;
                }
            }//end if(go = true)
        
            MMAX = IMMURN;
            return ITER;

        }//end method MinVal

        /// <summary>
        /// Single vector Lanczos routine with no reorthogonalization. 'Ghost' or 'repeat' eigenvalues are removed using the technique of Cullum
        /// </summary>
        /// <param name="evs">
        /// Vector to store eigenvalues
        /// </param>
        /// <param name="z">
        /// Array to store eigenvectors if needed, in the case that NTooBig == true, this array will contain the eigenvectors of the Lanczos matrix
        /// </param>
        /// <param name="A">
        /// Matrix to be diagonalized.
        /// </param>
        /// <param name="its">
        /// Size of the Lanczos matrix to be generated.
        /// </param>
        /// <param name="tol">
        /// Tolerance to be used when using Cullum technique to check if an eigenvalue is real or a 'ghost'
        /// </param>
        /// <param name="evsNeeded">
        /// True if eigenvectors are needed.
        /// </param>
        /// <param name="n">
        /// Which j-block is being diagonalized. Used in case the Lanczos vectors need to be written to file.
        /// </param>
        /// <param name="file">
        /// FileInfo object
        /// </param>
        public static void NaiveLanczos(ref double[] evs, ref double[,] z, alglib.sparsematrix A, int its, double tol, bool evsNeeded, int n, string file)
        {
            int N = A.innerobj.m;
            int M = evs.Length;
            
            var alphas = new double[its];
            var betas = new double[its];
            betas[0] = 0.0;

            var lanczosVecs = new double[0,0];
            bool NTooBig = false;
            if (N * its > basisSetLimit)
            {
                NTooBig = true;
            }
            //This will create the temp files needed to store the lanczos vectors for the eigenvector calculation        
            if (NTooBig && evsNeeded)
            {
                //create file to store the eigenvectors in the directory
                string fileDirectory = file + "temp_vecs_" + n + ".tmp";
                StreamWriter writer = new StreamWriter(fileDirectory);
                writer.WriteLine("Temporary storage of Lanczos Vectors. \n");
                LanczosIterations(A, its, evsNeeded, ref alphas, ref betas, ref lanczosVecs, NTooBig, writer);
                writer.Close();
            }
            else     //means either the eigenvectors are not needed or they are needed and the basis set is sufficiently small that lanczos vectors will be kept in memory
            {
                //initialize array for eigenvectors if necessary
                if (evsNeeded)
                {
                    lanczosVecs = new double[N, its];
                }
                LanczosIterations(A, its, evsNeeded, ref alphas, ref betas, ref lanczosVecs, NTooBig);
            }
                        
            double[] nBetas = new double[its - 1];
            //tAlphas and tBetas are diagonal and off diagonal for matix "T^2"
            double[] tAlphas = new double[its - 1];
            double[] tBetas = new double[its - 2];
            for (int i = 0; i < nBetas.Length; i++)
            {
                nBetas[i] = betas[i + 1];
                tAlphas[i] = alphas[i + 1];
                if (i == nBetas.Length - 1)
                {
                    break;
                }
                tBetas[i] = betas[i + 2];
            }

            //Diagonalize Lanczos Matrix to find M correct eigenvalues
            LanczosMatrixDiagonalization(ref evs, ref z, tol, M, alphas, nBetas, tAlphas, tBetas, evsNeeded);

            //if needed, build array of eigenvectors to return
            if (evsNeeded)
            {       
                /*
                //generate the eigenvectors by matrix multiplication
                //temporary storage for eigenvectors
                var tempEvecs = new double[its, correctEvs.Count];
                for (int i = 0; i < correctEvs.Count; i++)
                {
                    for (int j = 0; j < its; j++)
                    {
                        //pull only eigenvectors which correspond to a true eigenvalue
                        tempEvecs[j, i] = z[j, correctEvs[i].Item1];
                    }
                }
                */
                if (!NTooBig)
                {
                    //do matrix multiplication of tempEvecs and laczosVecs, results stored in transEvecs which are true eigenvectors.
                    double[,] transEvecs = new double[N, evs.Length];
                    alglib.rmatrixgemm(N, evs.Length, its, 1.0, lanczosVecs, 0, 0, 0, z, 0, 0, 0, 0.0, ref transEvecs, 0, 0);//changed tempEvecs to z

                    //then normalize
                    normalize(ref transEvecs);

                    //then set equal to z
                    z = transEvecs;
                }
                    /*
                else
                {
                    //if evectors will be generated after the fact then pass back the eigenvectors of the Lanczos matrix for later use
                    z = tempEvecs;
                }
                */
            }//end if evsNeeded
        }//end NaiveLanczos

        /// <summary>
        /// Diagonalizes the tridiagonal Lanczos matrix T and the matrix T^2 and finds the true or 'non-spuriou' eigenvalues of T
        /// </summary>
        /// <param name="evs">
        /// Will contain correct eigenvalues on return
        /// </param>
        /// <param name="z">
        /// If eigenvectors are requested, will contain eigenvectors on return
        /// </param>
        /// <param name="tol">
        /// Tolerance used for eigenvalue comparison
        /// </param>
        /// <param name="M">
        /// Number of eigenvalues needed
        /// </param>
        /// <param name="alphas">
        /// Diagonal elements of T
        /// </param>
        /// <param name="nBetas">
        /// Off diagonal elements of T
        /// </param>
        /// <param name="tAlphas">
        /// Diagonal elements of T^2
        /// </param>
        /// <param name="tBetas">
        /// Off diagonal elements of T^2
        /// </param>
        /// <param name="zz">
        /// Flag for ALGLIB function call for diagonalizing T to requeste eigenvalues or not
        /// </param>
        /// <param name="correctEvs">
        /// Tuple to store the positions of the correct eigenvalues so that the correct eigenvectors may be extracted.
        /// </param>
        private static void LanczosMatrixDiagonalization(ref double[] evs, ref double[,] z, double tol, int M, double[] alphas, double[] nBetas, double[] tAlphas, double[] tBetas, bool evsNeeded)
        {
            //dummy matrix for diagonalization of T^2
            var ZZ = new double[0, 0];

            int its = alphas.Length;
            int evsRequested = M * 4;
            if (evsRequested >= its)
            {
                evsRequested = its - its / 5;
            }

            //Tuple so that we know which eigenvectors to pull if necessary
            List<Tuple<int, double>> correctEvs = new List<Tuple<int, double>>();

            //since the vectors alphas and tAlphas are overwritten by the diagonaliztion routine but may be needed at a later point they are stored here so that if the diagonalization needs to run
            //for larger M to get all requested eigenvalues then they can be regenerated.
            double[] alphOriginal = new double[alphas.Length];
            double[] tAlphaOriginal = new double[tAlphas.Length];
            for (int i = 0; i < alphas.Length; i++)
            {
                alphOriginal[i] = alphas[i];
                if (i < tAlphas.Length)
                {
                    tAlphaOriginal[i] = tAlphas[i];
                }
            }

            int evalsFound = 0;
            while (evalsFound < M)
            {
                correctEvs = new List<Tuple<int, double>>();
                //set flag for eigenvectors and initialize eigenvector array if necessary
                int zz = 0;
                if (evsNeeded)
                {
                    zz = 2;
                    z = new double[its, evsRequested];
                }

                //if the diagonalization has run then alphas has been replaced by the eigenvalues
                //reset it and tAlphas to their original values
                if (alphas.Length != alphOriginal.Length)
                {
                    alphas = new double[alphOriginal.Length];
                    tAlphas = new double[tAlphaOriginal.Length];
                    for (int i = 0; i < alphOriginal.Length; i++)
                    {
                        alphas[i] = alphOriginal[i];
                        if (i < tAlphaOriginal.Length)
                        {
                            tAlphas[i] = tAlphaOriginal[i];
                        }
                    }
                }

                //this is the diagonalization of the complete Lanczos matrix
                bool test = alglib.smatrixtdevdi(ref alphas, nBetas, alphas.Length, zz, 0, evsRequested, ref z);//changed M - 1 to evsRequested
                evs = alphas;

                //this is the diagonalization of the Lanczos matrix without the first row and column
                bool test2 = alglib.smatrixtdevdi(ref tAlphas, tBetas, tAlphas.Length, 0, 0, evsRequested, ref ZZ);//changed M to evsRequested

                //here I run test from Lanczos book to see if evs are good.  Briefly, 
                //there are 3 cases to consider:
                //1. The ev is in T^2 and is a multiple ev in Tm: accept ev as good.
                //2. The ev is in T^2 and is in Tm but not as a multiple: reject ev.
                //3. The ev is not in T^2 and is not a multiple ev in Tm: accept ev.
                //evs evaluated as correct will be stored in this list

                for (int i = 0; i < tAlphas.Length - 1; i++)
                {
                    //check to see if the value is in T^2 eigenvalues
                    bool temp = checkInTT(alphas[i], tAlphas, tol);

                    //tells how many times alphas[i] is in alphas, always at least 1
                    int repeater = repeat(i, alphas, tol);

                    if (repeater == 0)    //means there is some problem in the alphas vectors (probably NaN error)
                    {
                        throw new RepeaterError();
                    }
                    //this evaluates to true if the ev is not a repeat by evaluating function repeat for alphas[i], this is condition 3.
                    if (repeater == 1)
                    {
                        //loop over elements of tAlphas to see if this ev is also in tAlphas, if not add it to the output list
                        if (!temp)
                        {
                            correctEvs.Add(new Tuple<int, double>(i, alphas[i]));
                            continue;
                        }
                    }
                    else    //means this evalue has a repeat and is a true eigenvalue of A
                    {
                        correctEvs.Add(new Tuple<int, double>(i, alphas[i]));
                        i += repeater - 1;
                    }
                }
                evalsFound = correctEvs.Count();
                if (evalsFound < M)
                {
                    evsRequested += 10;
                }
            }//end while loop

            //build array of final eigenvalues to return
            double[] checkedEvs = new double[M];
            for (int i = 0; i < checkedEvs.Length; i++)
            {
                checkedEvs[i] = correctEvs[i].Item2;
            }
            evs = checkedEvs;

            if (evsNeeded)
            {
                //temporary storage for eigenvectors
                var tempEvecs = new double[its, M];
                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < its; j++)
                    {
                        //pull only eigenvectors which correspond to a true eigenvalue
                        tempEvecs[j, i] = z[j, correctEvs[i].Item1];
                    }
                }
                z = tempEvecs;
            }
        }//end LanczosMatrixDiagonalization

        /// <summary>
        /// Function which generates the Lanczos matrix and optionally stores or writes to file the Lanczos vectors if necessary.
        /// </summary>
        /// <param name="A">
        /// Sparse matrix being diagonalized
        /// </param>
        /// <param name="its">
        /// Number of iterations for the lanczos routine
        /// </param>
        /// <param name="evsNeeded">
        /// True if eigenvectors will be calculated
        /// </param>
        /// <param name="alphas">
        /// Vector for storing the diagonal of the Lanczos matrix.
        /// </param>
        /// <param name="betas">
        /// Vector for storing the off diagonal of the Lanczos matrix.
        /// </param>
        /// <param name="lanczosVecs">
        /// Array to store the lanczosVecs, should be initialized if evsNeeded and NTooBig == false
        /// </param>
        /// <param name="NTooBig">
        /// True if the basis set is too large to store the Lanczos Vectors in memory and they must be written to file.
        /// </param>
        /// <param name="writer">
        /// StreamWriter for writing LanczosVectors to disc if necessary, null by default
        /// </param>
        private static void LanczosIterations(alglib.sparsematrix A, int its, bool evsNeeded, ref double[] alphas, ref double[] betas, ref double[,] lanczosVecs, bool NTooBig, StreamWriter writer = null)
        {
            int N = A.innerobj.m;
            //initialize vectors to store the various vectors used in the Lanczos iterations
            var vi = RANDOM(N);
            var viminusone = new double[N];
            var viplusone = new double[N];
            double[] Axvi = new double[N];
            betas[0] = 0.0;

            //loop to generate the Lanczos matrix
            for (int i = 0; i < its; i++)
            {
                //Axvi will contain the product of A and vi                
                OP(A, vi, ref Axvi, 1);

                //assign alpha value for this iteration
                //***************************************************************
                alphas[i] = Alpha_i(Axvi, viminusone, vi, betas[i]);
                //***************************************************************

                //if calculating the eigenvectors then store the lanczos vectors here
                if (evsNeeded)
                {
                    //store in memory if small enough
                    if (!NTooBig)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            lanczosVecs[j, i] = vi[j];
                        }
                    }
                    //if not, then write them to file
                    else
                    {
                        //write which iteration this is
                        writer.WriteLine(" ");
                        writer.WriteLine("START_VEC" + i);
                        //write the vector to the evFile
                        for (int j = 0; j < N; j++)
                        {
                            writer.WriteLine(vi[j]);
                        }
                        writer.WriteLine("END_VEC");
                        writer.WriteLine(" ");
                    }
                }//end if evsNeeded

                //to keep from calculating meaningless values for beta and vi
                if (i == its - 1)
                {
                    Console.WriteLine("Lanczos iterations completed. \nEntering diagonalization.");
                    break;
                }

                //calculate viplusone and beta i + 1.
                viplusone = betavplusone(Axvi, alphas[i], vi, betas[i], viminusone);
                //***************************************************************
                betas[i + 1] = Magnitude(viplusone);
                //***************************************************************
                for (int j = 0; j < viplusone.Length; j++)
                {
                    viplusone[j] /= betas[i + 1];
                }
                //now reassign vi vectors for next iteration.
                viminusone = vi;
                vi = viplusone;

                //let the user know things are happening
                Console.WriteLine("Lanczos iteration " + (i + 1) + " done.");
            }//end loop to generate Lanczos matrix
        }//end LanczosIterations

        /// <summary>
        /// Calculates alpha according to Cullum 4.3.1
        /// </summary>
        /// <param name="Axvi">
        /// Product of A and vi, vector
        /// </param>
        /// <param name="v_iminusone">
        /// Lanczos Vector from previous iteration
        /// </param>
        /// <param name="v_i">
        /// Lanczos vector from current iteration
        /// </param>
        /// <param name="beta_i">
        /// Beta_i
        /// </param>
        /// <returns>
        /// Alpha_i
        /// </returns>
        private static double Alpha_i(double[] Axvi, double[] v_iminusone, double[] v_i, double beta_i)
        {
            double[] temp = new double[Axvi.Length];
            for (int i = 0; i < temp.Length; i++)
            {
                temp[i] = Axvi[i] - beta_i * v_iminusone[i];
            }
            return vxv(v_i, temp);
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
        private static double vxv(double[] v, double[] u)
        {
            double product = 0.0;
            for (int i = 0; i < v.Length; i++)
            {
                product += v[i] * u[i];
            }
            return product;
        }//end vxv

        private static double[] betavplusone(double[] avi, double alphai, double[] vi, double betai, double[] viminusone)
        {
            var product = new double[avi.Length];
            for (int i = 0; i < product.Length; i++)
            {
                product[i] = avi[i] - alphai * vi[i] - betai * viminusone[i];
            }
            return product;
        }//end betavplusone

        /// <summary>
        /// Computes the projection of v onto u.
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static double[] projection(double[] u, double[] v)
        {
            double uv = vxv(u, v);
            double uu = vxv(u, u);
            uv /= uu;
            double[] proj = new double[u.Length];
            for (int i = 0; i < u.Length; i++)
            {
                proj[i] = uv * u[i];
            }
            return proj;
        }//end projection

        /// <summary>
        /// Checks if alphas is in tAlphas.
        /// </summary>
        /// <param name="alphas">
        /// Value to check against.
        /// </param>
        /// <param name="tAlphas">
        /// Array to look in for value.
        /// </param>
        /// <returns>
        /// bool. True if alphas is in tAlphas, false if not.
        /// </returns>
        private static bool checkInTT(double alphas, double[] tAlphas, double tol)
        {
            bool temp = false;
            for (int j = 0; j < tAlphas.Length; j++)
            {
                if (Math.Abs(alphas - tAlphas[j]) < tol)
                {
                    temp = true;
                    return temp;
                }
            }
            return temp;
        }//end checkInTT

        /// <summary>
        /// Checks a certain entry in an array to see how many times that same value, to within the tolerance, is in the array.
        /// </summary>
        /// <param name="j">
        /// Index of array value to be checked for repeats
        /// </param>
        /// <param name="alphvec">
        /// Vector being checked
        /// </param>
        /// <param name="tol">
        /// Tolerance to define a repeat so that if value1 - value2 is less than tol it is treated as a repeat
        /// </param>
        /// <returns>
        /// Number of times the value in alphvec[j] is repeated.  Always at least 1.
        /// </returns>
        private static int repeat(int j, double[] alphvec, double tol)
        {
            int count = 0;
            for (int i = 0; i < alphvec.Length; i++)
            {
                if (Math.Abs(alphvec[j] - alphvec[i]) < tol)
                {
                    count++;
                }//end if < tol
            }//end for
            return count;
        }//end repeat

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
    }//end class Lanczos
}
