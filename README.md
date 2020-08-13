# OTOC
Out-of-time-ordered-correlators for the anharmonic (quartic) quantum oscillator

Create directories "matrices" and "results"

Use "compile.sh" to compile the three executables (gen_mat, evs, and otoc)

Running: first you need to generate the Matrices, syntax

./gen_mat [N]

where N is dimension of matrix you want to generate. If the matrix file already exists in the "matrices" subdirectory, the code will abort unless you force it to overwrite the file with ./gen_mat [N] f

Next, solve the eigensystem using

./evs [N] 

If you want to generate also the matrix elements <x^2>, then you can force it to do that with ./evs [N] f

Finally, run

./otoc [N] [MAX]

where [MAX] is the trunction of the matrix system for calculating the otocs. Generically, you should not need to go beyond MAX=32.

The results folder then will contain the following files:

* En[N].dat for the eigenvalues
* Vn[N].dat for the spectral eigenfunctions
* xn[N].dat for the matrix elements <n|x|m> with n even, m odd
* Ax2[N].dat for the matrix elements <n|x^2|n> (if requested)
* cn[N].dat for the microcanonical OTOCs c_n(t)
* otoc[N]-MAX[MAX].dat for the thermal OTOCs O_T(t)
* SSF[N]-MAX[MAX].dat for the spectral form factor g(beta,t)
* Tcomp[N]-MAX[MAX].dat for partition function Z(beta) and others as a function of temperature
