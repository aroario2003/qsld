module quantum.impure_state.observable;

import std.stdio;
import std.complex;
import std.typecons : Nullable;

import linalg.vector;
import linalg.matrix;

struct Observable {
    Complex!real[] coeffs;
    string[] observable_decomp;
    int num_qubits;

    /**
    * The constructor for the observable type, an observable in quantum mechanics is 
    * something you can measure within a quantum system.
    * 
    * params:
    * observable_decomp = The observable deomposed into its tensor products of pauli terms
    *
    * coeffs = the weight of each pauli term in the observable decomposition
    */
    this(string[] observable_decomp, Complex!real[] coeffs, int num_qubits) {
        this.observable_decomp = observable_decomp;
        this.coeffs = coeffs;
        this.num_qubits = num_qubits;

        assert(this.observable_decomp.length == this.coeffs.length,
            "Array of coefficients and array of pauli terms have to be the same length");
    }

    /**
    * Applies the observable object to the density matrix rho
    *
    * params:
    * rho = The density matrix of the quantum system
    *
    * returns: An array of full pauli matrices applied to rho
    */
    Matrix!(Complex!real)[] apply(Matrix!(Complex!real) rho) {
        Matrix!(Complex!real) pauli_x = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) pauli_y = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, -1)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 1), Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) pauli_z = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(-1, 0)
                    ])
            ]);

        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

        Matrix!(Complex!real)[] full_pauli_ops;
        foreach (j, term; this.observable_decomp) {
            assert(term.length == this.num_qubits,
                "Each term in the observable decomposition must be the same length as the number of qubits");

            Matrix!(Complex!real) pauli;
            Matrix!(Complex!real) full_pauli_term;
            for (int i = 0; i < term.length; i++) {
                if (term[i] == 'I') {
                    pauli = identity;
                } else if (term[i] == 'X') {
                    pauli = pauli_x;
                } else if (term[i] == 'Y') {
                    pauli = pauli_y;
                } else if (term[i] == 'Z') {
                    pauli = pauli_z;
                }

                if (i == 0) {
                    full_pauli_term = pauli;
                } else {
                    full_pauli_term = full_pauli_term.kronecker(pauli);
                }
            }
            full_pauli_ops[full_pauli_ops.length++] = full_pauli_term;
        }

        Matrix!(Complex!real)[] result_mats;
        foreach (i, pauli_mat; full_pauli_ops) {
            result_mats[result_mats.length++] = rho.mult_mat(pauli_mat);
        }

        return result_mats;
    }
}
