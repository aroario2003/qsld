module quantum.pure_state.observable;

import std.stdio;
import std.complex;
import std.format;
import std.conv;

import linalg.vector;

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

        assert(this.observable_decomp.length == this.coeffs.length, "Array of coefficients and array of pauli terms have to be the same length");
    }

    /**
    * Applies the observable object to the state vector psi 
    *
    * params:
    * psi = the state vector of the quantum system
    *
    * returns: A new complex valued vector
    */
    Vector!(Complex!real) apply(Vector!(Complex!real) psi) {
        Vector!(Complex!real) phi = Vector!(Complex!real)(cast(int) psi.length(), new Complex!real[psi.length()]);
        // this is neccessary because the compiler doesnt 0 initialize the vector automatically
        for (int i = 0; i < phi.length(); i++) {
            phi[i] = Complex!real(0, 0);
        }

        foreach (i, term; this.observable_decomp) {
            assert(term.length == this.num_qubits, "Each term in the observable decomposition must be the same length as the number of qubits");

            int j_prime = 0;

            // j is the basis state index
            for (int j = 0; j < psi.length(); j++) {
                Complex!real phase_acc = Complex!real(1, 0);
                string bit_str = format("%0*b", this.num_qubits, j);
                char[] bit_arr = bit_str.dup;
                // q is the qubit index in bit_str
                for (int q = 0; q < this.num_qubits; q++) {
                    if (term[q] == 'I') {
                        continue;
                    } else if (term[q] == 'X') {
                        // convert the char in the bit_str which is '0' or '1' to an int
                        int bit_q = bit_arr[q] - '0';

                        // flip that bit since we are applying pauli-x
                        char bit_prime = cast(char)((bit_q ^ 1) + '0');
                        bit_arr[q] = bit_prime;

                    } else if (term[q] == 'Y') {
                        // convert the char in the bit_str which is '0' or '1' to an int
                        int bit_q = bit_arr[q] - '0';

                        // flip that bit since we are applying pauli-y
                        char bit_prime = cast(char)((bit_q ^ 1) + '0');
                        bit_arr[q] = bit_prime;

                        phase_acc *= Complex!real(0, 1);
                    } else if (term[q] == 'Z') {
                        if (bit_arr[q] == '1') {
                            phase_acc *= Complex!real(-1, 0);
                        }
                    }
                }
                j_prime = to!int(bit_arr, 2);
                phi[j_prime] = phi[j_prime] + this.coeffs[i] * phase_acc * psi[j];
            }
        }
        return phi;
    }
}
