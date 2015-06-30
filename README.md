SOCJT_2
=======

A program for calculating the vibronic eigenvalues and eigenvectors in a Jahn-Teller active electronic state. 

SOCJT 2 is a a follow up program to SOCJT, written by Terry A. Miller and Timothy Barckholtz. Compared with SOCJT, SOCJT 2 includes more possible types of coupling and the possibility of including nondegenerate modes in calculations. SOCJT 2 runs between 10 and 100 times faster than SOCJT due to the use of much more efficient numerical routines for the Hamiltonian diagonalization and parallelization throughout.

SOCJT 2 may be used to calculate spectra to see the effect of various types of coupling or experimental spectra may be fit by including experimental frequencies in a .fit file. For fits, user specified parameters will be fit using a Levenberg-Marquardt optimizer to minimize the error in the energies. This optimizer along with other numerical routines are from ALGLIB.net and are used under the educational license

SOCJT 2 includes two diagonalization routines, both of which are Lanczos type algorithms. The first uses the Underwood method which is a block Lanczos algorithm. This method was translated from FORTRAN 77 to C# for this project. It has been well tested in both SOCJT and SOCJT 2. The second is a so called 'naive' Lanczos which is the single vector Lanczos algorithm with no reorthogonalization. Spurious eigenvalues arising from loss of orthogonality are eliminated using the method proposed by Jane Cullum. The benefit of the naive Lanczos is that it is orders of magnitude faster in most cases than the block method when only the eigenvalues are needed (e.g. for fits).

SOCJT 2 has been tested in both Windows and on a Linux cluster runnng Red Hat. An instruction manual is included in this repository.
