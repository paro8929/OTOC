# OTOC
Out-of-time-ordered-correlators for the anharmonic (quartic) quantum oscillator

Create directories "matrices" and "results"

Use "compile.sh" to compile the three executables (gen_mat, evs, and otoc)

Running: first you need to generate the Matrices, syntax

./gen_mat [N]

where N is dimension of matrix you want to generate. If the matrix file already exists in the "matrices" subdirectory, the code will abort unless you force it to overwrite the file with ./gen_mat [N] f
