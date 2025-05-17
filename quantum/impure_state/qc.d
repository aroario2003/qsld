// NOTE: This file contains the code for the impure subsystem of QSLD. this means that 
// the modules within this folder (quantum/impure_state/) work more correctly but not
// neccessarily efficiently for mixed or impure quantum states. The lack of efficiency
// is due to the fact that in order to maintain correctness when decohering an entangled
// state, I must construct a density matrix and full multi-qubit gates. This leads to
// inefficiency and high memory usage. When you are using this module it is important
// to keep this in mind. It is recommended that you do not use more than 5 or 6 qubits 
// in a single quantum circuit because otherwise you might kill your computer due to 
// high memory usage. If you would like to have efficiency and correctness for pure states
// with the ability to use more qubits, you should use the pure subsystem (quantum/pure_state/).

module quantum.impure_state.qc;

// standard library modules
import std.stdio;
import std.complex;
import std.math;
import std.typecons;

// linear algebra modules
import linalg.vector;
import linalg.matrix;

struct QuantumCircuit {
    // These are for the circuit itself
    int num_qubits;
    Matrix!(Complex!real) density_mat;
    int num_probabilities;

    // These are for circuit visualization
    int timestep;
    Tuple!(string, int[], int)[] visualization_arr;

    /**
    * Constructor for quantum circuit object (impure subsystem) 
    * 
    * params: 
    * num_qubits = The number of qubits for the circuit to have
    */
    this(int num_qubits) {
        this.num_qubits = num_qubits;

        this.num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];

        state_arr[] = Complex!real(0, 0);
        state_arr[0] = Complex!real(1, 0);

        Vector!(Complex!real) ket_psi = Vector!(Complex!real)(num_probabilities, state_arr);
        Matrix!(Complex!real) bra_psi = ket_psi.dagger();
        this.density_mat = ket_psi.outer_prod(bra_psi);
    }

    /**
    * Overload of the constructor for the quantum circuit object (impure subsystem)
    *
    * params:
    * num_qubits = The number of qubits for the circuit to have
    *
    * starting_state_idx = The index in the state vector of the amplitude to have 100% 
    *                      probability when starting out
    */
    this(int num_qubits, int starting_state_idx) {
        this.num_qubits = num_qubits;

        this.num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];

        state_arr[] = Complex!real(0, 0);
        state_arr[starting_state_idx] = Complex!real(1, 0);

        Vector!(Complex!real) ket_psi = Vector!(Complex!real)(num_probabilities, state_arr);
        Matrix!(Complex!real) bra_psi = ket_psi.dagger();
        this.density_mat = ket_psi.outer_prod(bra_psi);
    }

    // Updates the visualization internal representation 
    private void update_visualization_arr(string gate_name, int[] qubit_idxs) {
        this.visualization_arr[this.visualization_arr.length++] = tuple(gate_name, qubit_idxs, this
                .timestep);
        this.timestep += 1;
    }

    /**
    * The hadamard quantum gate puts the state into superposition with equal probabilities for each state in 
    * superposition if applied to all qubits in the system. Otherwise, Some states will have different probability
    * amplitudes then others. (impure subsystem)
    * 
    * params:
    * qubit_idx = the index of the qubit to affect
    */
    void hadamard(int qubit_idx) {
        update_visualization_arr("H", [qubit_idx]);

        Matrix!(Complex!real) hadamard = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(-1, 0)
                    ])
            ]).mult_scalar(Complex!real(1 / sqrt(2.0), 0));

        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);
        Matrix!(Complex!real)[] kronecker_chain;

        for (int i = 0; i < this.num_qubits; i++) {
            if (i == qubit_idx) {
                kronecker_chain[kronecker_chain.length++] = hadamard;
            } else {
                kronecker_chain[kronecker_chain.length++] = identity;
            }
        }

        Matrix!(Complex!real) result = kronecker_chain[0];
        for (int i = 1; i < kronecker_chain.length; i++) {
            result = result.kronecker(kronecker_chain[i]);
        }

        this.density_mat = result.mult_mat(this.density_mat).mult_mat(result.dagger());
    }
}
