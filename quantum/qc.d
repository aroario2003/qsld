module quantum.qc;

import std.stdio;
import std.complex;
import std.math;
import std.typecons;
import std.format;
import std.algorithm;
import std.random;

import linalg.matrix;
import linalg.vector;

struct QuantumCircuit {
    int num_qubits;
    Vector!(Complex!real) state;

    this(int num_qubits) {
        this.num_qubits = num_qubits;

        int num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];
        //initialize state vector to all 0+0i amplitudes
        state_arr[] = Complex!real(0, 0);
        // start with a valid classical state by setting one of the amplitudes probabilities to 100%
        state_arr[0] = Complex!real(1, 0);
        this.state = Vector!(Complex!real)(num_probabilities, state_arr);
    }

    /// The hadamard quantum gate puts the state into superposition with equal probabilities for each state in 
    /// superposition if applied to all qubits in the system. Otherwise, Some states will have different probability
    /// amplitudes then others.
    void hadamard(int qubit_idx) {
        // make sure that the 1/sqrt(2) is scalar multiplied by the hadamard matrix
        Matrix!(Complex!real) hadamard = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(-1, 0)
                    ])
            ]).mult_scalar(Complex!real(1 / sqrt(2.0), 0));

        auto pairs = new Vector!int[(this.state.length() / 2)];
        int pairs_idx = 0;

        // Find bit flipped unique pairs to apply hadamard without having to construct full multi-qubit hadamard matrix.
        // Doing so is O(n^2^n) whereas this method is O(n^2n).
        for (int i = 0; i < this.state.length(); i++) {
            bool qubit_is_one = (i & (1 << qubit_idx)) != 0;
            if (!qubit_is_one) {
                int j = i ^ (1 << qubit_idx);
                pairs[pairs_idx] = Vector!int(2, [
                        i, j
                    ]);
                pairs_idx++;
            }
        }

        foreach (vec; pairs) {
            Vector!(Complex!real) amplitudes = Vector!(Complex!real)(2, [
                    this.state[vec[0]], this.state[vec[1]]
                ]);
            Vector!(Complex!real) updated_amplitudes = hadamard.mult_vec(amplitudes);

            this.state[cast(ulong) vec[0]] = updated_amplitudes[0];
            this.state[cast(ulong) vec[1]] = updated_amplitudes[1];
        }
    }

    /// The pauli-x gate or NOT gate, negates the current state of the qubit so |0> -> |1> and |1> -> |0>.
    /// More concisely it performs a bit flip. 
    void pauli_x(int qubit_idx) {
        Matrix!(Complex!real) pauli_x = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ])
            ]);

        auto pairs = new Vector!int[(this.state.length() / 2)];
        int pairs_idx = 0;

        for (int i = 0; i < this.state.length(); i++) {
            bool qubit_is_zero = (i & (1 << qubit_idx)) == 0;
            if (qubit_is_zero) {
                int j = i ^ (1 << qubit_idx);
                pairs[pairs_idx] = Vector!int(2, [
                        i, j
                    ]);
                pairs_idx++;
            }
        }

        foreach (vec; pairs) {
            Complex!real temp = this.state[vec[0]];
            this.state[vec[0]] = this.state[vec[1]];
            this.state[vec[1]] = temp;
        }
    }

    /// The pauli-y gate applies an imaginary relative phase to a state when flipping the state, for |1> -> |0> multiply by i.
    /// And for |0> -> |1> multiply by -i.
    void pauli_y(int qubit_idx) {
        for (int i = 0; i < this.state.length(); i++) {
            int j = i ^ (1 << qubit_idx);
            if (i < j) {
                Complex!real temp = this.state[i];
                this.state[i] = this.state[j] * Complex!real(0, 1);
                this.state[j] = temp * Complex!real(0, -1);
            }
        }
    }

    /// The pauli-z gate puts a relative phase on the |1> state and leaves |0> untouched
    void pauli_z(int qubit_idx) {
        for (int i = 0; i < this.state.length(); i++) {
            if ((i & (1 << qubit_idx)) != 0) {
                this.state[i] = this.state[i] * Complex!real(-1, 0);
            }
        }
    }

    string measure() {
        Vector!float probs = Vector!float(cast(int) this.state.length(), new float[this.state.length()]);
        // Take the magnitude of each complex probability amplitude
        foreach (i, c; this.state.elems) {
            float magnitude = sqrt(pow(c.re, 2) + pow(c.im, 2));
            float prob = pow(magnitude, 2);
            probs.elems[i] = prob;
        }

        // Perform inverse transform sampling on probabilities since measurement is non-algorithmic
        auto rng = Random(unpredictableSeed);
        auto r = uniform(0.0, 1.0f, rng);

        float sum = 0;
        string binary_result;
        foreach (i, elem; probs.elems) {
            sum += elem;
            if (r < sum) {
                binary_result = format("%0*b", this.num_qubits, i);
                break;
            }
        }
        return binary_result;
    }
}
