using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
    /// PROJECT by Bochkanov Sergey.
    /// </summary>
    static class Lanczos
    {
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
        private static void BKLANC(int N, int Q, int M, int P, int S, double[] D, ref double[,] C, ref double[,] X, ref double[] E, ref double[] U, ref double[] V, alglib.sparsematrix A, int par)
        {
            int LL;
            int LU;
            double T;
            int IT;
            int KP1;
            int IL;
            //I need to pay special attention to this routines loops and their bounds
            for (int L = 0; L < S; L++)//do 90
            {
                LL = M + (L) * P;//I changed this to L from (L - 1) and took off the + 1
                LU = M + (L + 1) * P;//I changed this to (L + 1) from L
                for (int K = LL; K < LU; K++)//do 70 --must be K < LU--
                {
                    for (int I = 0; I < N; I++)//do 10
                    {
                        U[I] = X[I, K];
                    }//do 10

                    OP(A, U, ref V, par);
                    
                    if (L == 0)//changed this to 0
                    {
                        for (int I = K; I < LU; I++)//do 12 --changed to I < LU from I <= LU--
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
                        for (int I = K; I < LU; I++)//do 30 --must be I < LU--
                        {
                            T = 0.0;
                            for (int J = 0; J < N; J++)//do 20
                            {
                                T += V[J] * X[J, I];
                            }//do 20
                            C[I, K] = T;
                        }//do 30
                        IT = K - P;//does this value go negative?--no
                        for (int I = 0; I < N; I++)//do 60
                        {
                            T = 0.0;
                            for (int J = IT; J <= K; J++)//do 40 --changed back to J < K from J <= K--
                            {
                                T += X[I, J] * C[K, J];
                            }//do 40
                            if (K != LU - 1)//changed this from K < LU - 1
                            {
                                KP1 = K + 1;//took off the + 1
                                for (int J = KP1; J < LU; J++)//do 50 --must be J < LU--
                                {
                                    T += X[I, J] * C[J, K];
                                }//do 50
                            }
                            V[I] -= T;
                        }//do 60
                    }//end else
                    //label 61
                    if (L != S - 1)//changed from L < S - 1
                    {
                        for (int I = 0; I < N; I++)//do 63
                        {
                            X[I, K + P] = V[I];//changed back to original and left - 1 off of K + p
                        }//do 63
                    }
                }//do 70
                if (L == 0)//changed to 0 from 1
                {
                    ERR(N, Q, M, P, X, ref E);
                }
                if (L < S - 1)//changed from L != S - 1 to L < S - 1
                {
                    ORTHG(N, Q, LU, P, ref C, ref X);
                    IL = LU;//took off + 1
                    int Ir = LU - 1;//changed back to original
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
            double CHEPS = 1.11e-16;//SOCJT_2 seems to converge faster no matter what

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
                //K = I + 1;//added the plus 1
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
        //CNVTST done new

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
                if (I <= P)//took off the -1
                {
                    LIM = 0;
                }
                if (LIM > 0)
                {
                    LM1 = LIM - 1;
                    for (int J = 0; J <= LM1; J++)//Changed to <=
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
            alglib.evd.smatrixevd(tempMat, PS, 1, false, ref evs, ref z);//changed to false            
            //put results back into the expected format for MINVAL routine
            for (int i = 0; i < PS; i++)
            {
                D[i + M] = evs[i];
            }
            for (int i = 0; i < PS; i++)
            {
                for (int j = 0; j < PS; j++)
                {
                    C[i, j] = z[i, j];//changed this
                }
            }
        }//end EIGEN
        //EIGEN done new I think

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
            int MP1 = M + P;//took off the + 1
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
                            for (int I = 0; I < KM1; I++)//do 20---changed this to I < KM1 from I <= KM1
                            {
                                S = 0.0;
                                for (int J = 0; J < N; J++)//do 15
                                {
                                    S += X[J, I] * X[J, K];
                                }//do 15
                                if (orig == true & I >= F)//changed to I>=F from I>F CHANGED
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
        //ORTHG done new I think

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
        //PCH done new
        
        private static void RANDOM(int N, int Q, int L, ref double[,] X)
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
                X1 = A * X0 + C;//added + 1
                if (X1 >= 10000)
                {
                    X1 -= 10000;
                }
                T[I] = (double) X1 / 9999.0 - 0.5;
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
                K = FT / 10000;//probably should not have the + 1 but does better with it.
                X[I, L] = T[K];
                X1 = A * X0 + C;//try a + 1 here too
                if (X1 >= 10000)
                {
                    X1 -= 10000;
                }
                T[K] = (double)X1 / 9999.0 - 0.5;
                X0 = X1;
            }//do 200
            T = null;
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
                for (int K = 0; K < L; K++)//changed this to K < L from K <= L.
                {
                    X[I, M + K] = V[K];
                }
            }
        }//end Rotate
        //rotate done new

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
            int ICOL1 = M - 1;//added this so that ICOL1++ for 1st iteration is 0.
            int ICOL2;
            double T;
            for (int J = 0; J < P; J++)//do 300
            {
                ICOL1++;  //PUT IT BACK, -- old comment -- commented this out because not necessary when indexed to 0.
                for (int I = 0; I < N; I++)//do 100
                {
                    U[I] = X[I, ICOL1];
                }//do 100

                OP(A, U, ref V, par);

                ICOL2 = M - 1;//fixed like above

                for (int I = 0; I <= J; I++)//do 300, changed to I <= J from I < J
                {
                    ICOL2++;//same as above, indexed to 0
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
        //SECTN done new

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
        /// <param name="IECODE"></param>
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

            //this checks that the values are in the proper ranges
            //need to remove these last gotos but low priority.
            if (N < 2)
            {
                goto onehundred;
            }
            if (R < 1)
            {
                goto onehundred;
            }
            if (Q <= R)
            {
                goto onehundred;
            }
            if (Q > N)
            {
                goto onehundred;
            }

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

                    RANDOM(N, Q, K, ref X);
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
         
            onehundred:
            if (N < 2)
            {
                IECODE = 1;
            }
            if (R < 1)
            {
                IECODE = 3;
            }
            if (Q <= R)
            {
                IECODE = 4;
            }
            if (Q > N)
            {
                IECODE = 6;
            }
        
            MMAX = IMMURN;
            return ITER;

        }//end method MinVal

    }
}
