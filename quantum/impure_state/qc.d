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
    int initial_state_idx;

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
        this.initial_state_idx = 0;

        this.num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];

        state_arr[] = Complex!real(0, 0);
        state_arr[this.initial_state_idx] = Complex!real(1, 0);

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
        this.initial_state_idx = starting_state_idx;

        this.num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];

        state_arr[] = Complex!real(0, 0);
        state_arr[this.initial_state_idx] = Complex!real(1, 0);

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

    // Builds the full matrix operator for a given controlled gate
    private Matrix!(Complex!real) build_full_controlled_gate(Matrix!(
            Complex!real) gate_matrix, int control_qubit_idx, int target_qubit_idx) {

        Matrix!(Complex!real) projection_0 = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) projection_1 = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ])
            ]);

        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

        // Make a kronecker product chain for when the control qubit is in 
        // the |0> state
        Matrix!(Complex!real)[] kronecker_chain_p0;
        for (int i = this.num_qubits - 1; i >= 0; i--) {
            if (i == control_qubit_idx) {
                kronecker_chain_p0[kronecker_chain_p0.length++] = projection_0;
            } else {
                kronecker_chain_p0[kronecker_chain_p0.length++] = identity;
            }
        }

        // Build a new unitary operator from the kronecker chain for control
        // in |0> state
        Matrix!(Complex!real) control_0_op = kronecker_chain_p0[0];
        for (int i = 1; i < kronecker_chain_p0.length; i++) {
            control_0_op = control_0_op.kronecker(kronecker_chain_p0[i]);
        }

        // Make a kronecker product chain for when the control qubit is in the 
        // |1> state
        Matrix!(Complex!real)[] kronecker_chain_p1;
        for (int i = this.num_qubits - 1; i >= 0; i--) {
            if (i == control_qubit_idx) {
                kronecker_chain_p1[kronecker_chain_p1.length++] = projection_1;
            } else if (i == target_qubit_idx) {
                kronecker_chain_p1[kronecker_chain_p1.length++] = gate_matrix;
            } else {
                kronecker_chain_p1[kronecker_chain_p1.length++] = identity;
            }
        }

        // Build a new unitary operator from the kronecker chain for control
        // in |1> state
        Matrix!(Complex!real) control_1_op = kronecker_chain_p1[0];
        for (int i = 1; i < kronecker_chain_p1.length; i++) {
            control_1_op = control_1_op.kronecker(kronecker_chain_p1[i]);
        }

        return control_0_op + control_1_op;
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

        for (int i = this.num_qubits - 1; i >= 0; i--) {
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

    /**
    * The controlled hadamard gate applies a hadamard transformation to the target qubit when the 
    * control qubit is in the state |1>
    *
    * params:
    * control_qubit_idx = the index of the qubit which determines if the other qubit is affected or not
    *
    * target_qubit_idx = the index of the qubit which is affected by the control 
    */
    void ch(int control_qubit_idx, int target_qubit_idx) {
        assert(this.num_qubits >= 2, "The number of qubits must be greater than or equal to two in order to use controlled gates");

        update_visualization_arr("CH", [control_qubit_idx, target_qubit_idx]);

        Matrix!(Complex!real) hadamard = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(-1, 0)
                    ])
            ]).mult_scalar(Complex!real(1 / sqrt(2.0), 0));

        Matrix!(Complex!real) ch_op = build_full_controlled_gate(hadamard, control_qubit_idx, target_qubit_idx);
        this.density_mat = ch_op.mult_mat(this.density_mat).mult_mat(ch_op.dagger());
    }

    /**
    * The controlled NOT gate checks if the control qubit is |1> if so it flips the target qubit.
    * 
    * params:
    * control_qubit_idx = the index of the qubit which determines if the target will be affected
    *
    * target_qubit_idx = the index of the qubit which is affected based on the state of the control 
    */
    void cnot(int control_qubit_idx, int target_qubit_idx) {
        assert(this.num_qubits >= 2, "The number of qubits must be greater than or equal to two in order to use controlled gates");

        update_visualization_arr("CX", [control_qubit_idx, target_qubit_idx]);

        Matrix!(Complex!real) pauli_x = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) cnot_op = build_full_controlled_gate(pauli_x, control_qubit_idx, target_qubit_idx);
        this.density_mat = cnot_op.mult_mat(this.density_mat).mult_mat(cnot_op.dagger());
    }
}
