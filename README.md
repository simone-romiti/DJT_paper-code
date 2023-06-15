# su2_discrete-code

DJT algorithm for the lattice $SU(2)$ Hamiltonian. 
See pdf note here for more details: https://github.com/simone-romiti/su2_discrete

This repository contains a python implementation of the gauge degrees of freedom for a $SU(2)$ gauge theory in the magnetic basis.
The main modules needed to use it are `operators.py` and `DJT_matrix.py`.
The files `test-*.py` test the main features of the algorithm.

## Things done

- Implementation of the DJT + test of orthogonality of columns
- pseudo change of basis matrix $V$ + test invertibility
- Canonical momenta built from the similarity transformation with the matrix $V$
- Gauge links: $N_\alpha \times N_\alpha$ components in color space + test of unitarity
- Test of the canonical commutation relations for the $L_a$ up to $N_{q-1/2}$.
- Test of the Lie algebra commutation relations (this should be trivial to check because $V^{-1} V = 1$)

## Things to be done

**Achtung**: I am using numpy.matrix arrays to represent matrices. This means that the `*` operator is the matrix multiplication. One can also use `np.dot()`. If one uses numpy.array the "*" will not do the matrix multiplication, but the one component by component.

- Test of the canonical commutation relations for the $R_a$ up to $N_{q-1/2}$.
- Plaquette operator: Hilbert space tensor product of links at multiple points and usual matrix produc in color space
- Hamiltonian for some given lattice
- Locality of the canonical momenta (heatmap plot?)
- Gauge invariant states and Gauss law

