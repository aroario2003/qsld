module algos.grovers;

import std.complex;
import std.math;
import std.format;

import quantum.pure_state.qc;

struct Grovers {
    int num_qubits;
    QuantumCircuit qc;

    /**
    * The constructor for the grovers algorithm object
    *
    * params:
    * num_qubits = The number of qubits to use in the algorithms internal circuit
    */
    this(int num_qubits) {
        this.num_qubits = num_qubits;
        this.qc = QuantumCircuit(this.num_qubits);
    }

    // The grovers algorithm oracle to be used
    private void oracle(int function(string) f) {
        for (int i = 0; i < this.qc.state.elems.length; i++) {
            if (f(format("%0*b", this.num_qubits, i)) == 1) {
                this.qc.state.elems[i] = this.qc.state.elems[i] * Complex!real(-1, 0);
            }
        }
    }

    // The diffusion operator to be used 
    private void diffusion() {
        Complex!real sum = Complex!real(0, 0);
        for (int i = 0; i < this.qc.state.elems.length; i++) {
            sum = sum + this.qc.state.elems[i];
        }

        Complex!real mean = sum / Complex!real(this.qc.state.elems.length, 0);

        for (int i = 0; i < this.qc.state.elems.length; i++) {
            this.qc.state.elems[i] = 2 * mean - this.qc.state.elems[i];
        }
    }

    /**
    * The main grovers algorithm which searches for a solution in an unsorted 
    * search space.
    *
    * params:
    * f = The function which encodes the solution to the search
    *
    * returns: An associative array of basis state to number of times measured
    */
    int[string] grovers(int function(string) f) {
        for (int i = 0; i < this.num_qubits; i++) {
            this.qc.hadamard(i);
        }

        real pi_over_four = PI / 4;
        int num_iterations = cast(int)(pi_over_four * sqrt(cast(real) pow(2, this.num_qubits)));

        for (int i = 0; i < num_iterations; i++) {
            oracle(f);
            diffusion();
        }

        return this.qc.measure_all(2000);
    }
}
