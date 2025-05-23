// NOTE: This file contains the code for the pure subsystem of QSLD, this means that
// the modules within this folder (quantum/pure_state/) work more correctly and efficiently 
// with pure quantum states. A pure state is one that does not suffer from decoherence 
// or noise of any sort. You can use decoherence with this module, however, you should 
// make sure that your quantum state isn't heavily entangled, otherwise, the result of 
// decoherence will be incorrect to varying degrees depending on the level of entanglement.
// The reason for this has to do with the code simulating state vector evolution and not
// density matrix evolution and therefore the state vector does not store and maintian 
// entanglement information as well as the density matrix. If you would like accurate
// results with an entangled state, you should use the impure subsystem (quantum/impure_state/).

module quantum.pure_state.qc;

// standard library modules
import std.stdio;
import std.complex;
import std.math;
import std.typecons;
import std.format;
import std.random;
import std.typecons;

// linear algebra modules
import linalg.matrix;
import linalg.vector;

// quantum related modules
import quantum.pure_state.observable;
import quantum.pure_state.decoherence;

// visualization related modules
import viz.visualization;

struct QuantumCircuit {
    // These are for the circuit itself
    int num_qubits;
    Vector!(Complex!real) state;
    int initial_state_idx;

    // These are for circuit visualization
    int timestep;
    Tuple!(string, int[], int)[] visualization_arr;

    // This is for decoherence with T1/T2 decay
    DecoherenceConfig decoherence_conf;

    /**
    * Constructor for quantum circuit object (pure subsystem)
    * 
    * params: 
    * num_qubits = The number of qubits for the circuit to have
    */
    this(int num_qubits) {
        this.num_qubits = num_qubits;
        this.initial_state_idx = 0;

        int num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];
        //initialize state vector to all 0+0i amplitudes
        state_arr[] = Complex!real(0, 0);
        // start with a valid classical state by setting one of the amplitudes probabilities to 100%
        state_arr[0] = Complex!real(1, 0);

        this.state = Vector!(Complex!real)(num_probabilities, state_arr);
        this.timestep = 0;

        this.decoherence_conf = DecoherenceConfig(Nullable!T1Decay.init, Nullable!T2Decay.init, "none");
    }

    /**
    * Overload of the constructor for the quantum circuit object (pure subsystem)
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

        int num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];
        //initialize state vector to all 0+0i amplitudes
        state_arr[] = Complex!real(0, 0);
        // start with a valid classical state by setting one of the amplitudes probabilities to 100%
        state_arr[starting_state_idx] = Complex!real(1, 0);

        this.state = Vector!(Complex!real)(num_probabilities, state_arr);
        this.timestep = 0;

        this.decoherence_conf = DecoherenceConfig(Nullable!T1Decay.init, Nullable!T2Decay.init, "none");
    }

    /**
    * Overload of the constructor for the quantum circuit object
    *
    * params:
    * num_qubits = The number of qubits for the circuit to have
    *
    * decoherence_conf = The config for whether you want to use T1 or T2 decay or both
    *                    and whether you want it be automatic, manual or not happen at all 
    */
    this(int num_qubits, DecoherenceConfig decoherence_conf) {
        this.num_qubits = num_qubits;
        this.initial_state_idx = 0;

        int num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];
        //initialize state vector to all 0+0i amplitudes
        state_arr[] = Complex!real(0, 0);
        // start with a valid classical state by setting one of the amplitudes probabilities to 100%
        state_arr[0] = Complex!real(1, 0);

        this.state = Vector!(Complex!real)(num_probabilities, state_arr);
        this.timestep = 0;

        this.decoherence_conf = decoherence_conf;
    }

    /**
    * Overload of the constructor for the quantum circuit object
    *
    * params:
    * num_qubits = The number of qubits for the circuit to have
    *
    * starting_state_idx = The index in the state vector of the amplitude to have 100% 
    *                      probability when starting out
    *
    * decoherence_conf = The config for whether you want to use T1 or T2 decay or both
    *                    and whether you want it be automatic, manual or not happen at all 
    */
    this(int num_qubits, int starting_state_idx, DecoherenceConfig decoherence_conf) {
        this.num_qubits = num_qubits;
        this.initial_state_idx = starting_state_idx;

        int num_probabilities = pow(2, this.num_qubits);
        Complex!real[] state_arr = new Complex!real[num_probabilities];
        //initialize state vector to all 0+0i amplitudes
        state_arr[] = Complex!real(0, 0);
        // start with a valid classical state by setting one of the amplitudes probabilities to 100%
        state_arr[starting_state_idx] = Complex!real(1, 0);

        this.state = Vector!(Complex!real)(num_probabilities, state_arr);
        this.timestep = 0;

        this.decoherence_conf = decoherence_conf;
    }

    // Updates the visualization internal representation 
    private void update_visualization_arr(string gate_name, int[] qubit_idxs) {
        this.visualization_arr[this.visualization_arr.length++] = tuple(gate_name, qubit_idxs, this
                .timestep);
        this.timestep += 1;
    }

    // Apply T1, T2 decay or both depending on the specification of 
    // DecoherenceConfig
    private void apply_decoherence(int qubit_idx, int gate_duration) {
        if (this.decoherence_conf.decoherence_mode == "automatic") {
            if (!this.decoherence_conf.t1.isNull()) {
                T1Decay t1 = this.decoherence_conf.t1.get();
                this.state = t1.apply(qubit_idx, gate_duration, this.state);
            }

            if (!this.decoherence_conf.t2.isNull()) {
                T2Decay t2 = this.decoherence_conf.t2.get();
                this.state = t2.apply(qubit_idx, gate_duration, this.state);
            }
        }
    }

    /**
    * The hadamard quantum gate puts the state into superposition with equal probabilities for each state in 
    * superposition if applied to all qubits in the system. Otherwise, Some states will have different probability
    * amplitudes then others. (pure subsystem)
    * 
    * params:
    * qubit_idx = the index of the qubit to affect
    */
    void hadamard(int qubit_idx) {
        update_visualization_arr("H", [qubit_idx]);

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

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 20);
    }

    /**
    * Overload for the hadamard gate to apply it to multiple qubits
    *
    * params:
    * qubit_idxs = the qubit indices to apply the hadamard gate to
    */
    void hadamard(int[] qubit_idxs) {
        foreach (idx; qubit_idxs) {
            assert(idx < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.hadamard(idx);
        }
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

        for (int i = 0; i < this.state.length(); i++) {
            bool cntl_qubit_is_one = (i & (1 << control_qubit_idx)) != 0;
            if (cntl_qubit_is_one) {
                int j = i ^ (1 << target_qubit_idx);
                if (i < j) {
                    Complex!real temp_i = this.state[i];
                    Complex!real temp_j = this.state[j];
                    this.state[i] = (temp_i + temp_j) / sqrt(2.0);
                    this.state[j] = (temp_i - temp_j) / sqrt(2.0);
                }
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(control_qubit_idx, 40);
        apply_decoherence(target_qubit_idx, 40);
    }

    /**
    * Overload for the controlled hadamard gate to apply it to multiple pairs of qubits at once
    *
    * params:
    * qubit_idxs = A tuple array of qubit indices with (int, int) pairs where index 0 is control and index 1 is target
    */
    void ch(Tuple!(int, int)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits && idx_tuple[1] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount specified for the system");
            this.ch(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The pauli-x gate or NOT gate, negates the current state of the qubit so |0> -> |1> and |1> -> |0>.
    * More concisely it performs a bit flip. 
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    */
    void pauli_x(int qubit_idx) {
        update_visualization_arr("X", [qubit_idx]);

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

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 20);
    }

    /**
    * Overload for the pauli_x gate to apply it to multiple qubits at once
    *
    * params:
    * qubit_idxs = Array of qubit indices to apply the pauli-x gate to
    */
    void pauli_x(int[] qubit_idxs) {
        foreach (idx; qubit_idxs) {
            assert(idx < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.pauli_x(idx);
        }
    }

    /**
    * The pauli-y gate applies an imaginary relative phase to a state when 
    * flipping the state, for |1> -> |0> multiply by i. And for |0> -> |1> 
    * multiply by -i.
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    */
    void pauli_y(int qubit_idx) {
        update_visualization_arr("Y", [qubit_idx]);

        for (int i = 0; i < this.state.length(); i++) {
            int j = i ^ (1 << qubit_idx);
            if (i < j) {
                Complex!real temp = this.state[i];
                this.state[i] = this.state[j] * Complex!real(0, 1);
                this.state[j] = temp * Complex!real(0, -1);
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 20);
    }

    /**
    * Overload for the pauli_y gate to apply it to multiple qubits at a time
    *
    * params:
    * qubit_idxs = Array of qubit indices to apply the pauli-y gate to
    */
    void pauli_y(int[] qubit_idxs) {
        foreach (idx; qubit_idxs) {
            assert(idx < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.pauli_y(idx);
        }
    }

    /**
    * The pauli-z gate puts a relative phase on the |1> state and leaves |0> untouched
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    */
    void pauli_z(int qubit_idx) {
        update_visualization_arr("Z", [qubit_idx]);

        for (int i = 0; i < this.state.length(); i++) {
            if ((i & (1 << qubit_idx)) != 0) {
                this.state[i] = this.state[i] * Complex!real(-1, 0);
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 0);
    }

    /**
    * Overload for the pauli_z gate to apply to multiple qubits at once
    *
    * params: 
    * qubit_idxs = An array of qubit indices to apply the pauli-z gate to
    */
    void pauli_z(int[] qubit_idxs) {
        foreach (idx; qubit_idxs) {
            assert(idx < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.pauli_z(idx);
        }
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

        for (int i = 0; i < this.state.length(); i++) {
            bool control_is_one = (i & (1 << control_qubit_idx)) != 0;

            if (control_is_one) {
                int j = i ^ (1 << target_qubit_idx);
                if (i < j) {
                    Complex!real temp = this.state[i];
                    this.state[i] = this.state[j];
                    this.state[j] = temp;
                }
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(control_qubit_idx, 100);
        apply_decoherence(target_qubit_idx, 100);
    }

    /**
    * Overload for the cnot gate to apply it multiple qubit pairs at once
    * 
    * params:
    * qubit_idxs = An array of qubit indices as tuples of (int, int) where 
    *              index 0 is control and index 1 is target
    */
    void cnot(Tuple!(int, int)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits && idx_tuple[1] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount specified for the system");
            this.cnot(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The S phase shift gate or PI/4 gate applies a phase shift of PI/4 to the state |1>
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    */
    void s(int qubit_idx) {
        update_visualization_arr("S", [qubit_idx]);

        for (int i = 0; i < this.state.length(); i++) {
            bool qubit_is_one = (i & (1 << qubit_idx)) != 0;
            if (qubit_is_one) {
                this.state[i] = this.state[i] * Complex!real(0, 1);
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 20);
    }

    /**
    * Overload of the s gate to apply it to multiple qubits at once
    *
    * params:
    * qubit_idxs = An array of qubit indices to apply the gate to
    */
    void s(int[] qubit_idxs) {
        foreach (idx; qubit_idxs) {
            assert(idx < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.s(idx);
        }
    }

    /**
    * The T phase shift gate or PI/8 gate applies a phase shift of PI/8 to the state |1>
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    */
    void t(int qubit_idx) {
        update_visualization_arr("T", [qubit_idx]);

        for (int i = 0; i < this.state.length(); i++) {
            bool qubit_is_one = (i & (1 << qubit_idx)) != 0;
            if (qubit_is_one) {
                this.state[i] = this.state[i] * expi(PI / 4);
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 20);
    }

    /**
    * Overload of the t gate to apply it to multiple qubits at once
    * 
    * params:
    * qubit_idxs = An array of qubit indices to apply the gate to
    */
    void t(int[] qubit_idxs) {
        foreach (idx; qubit_idxs) {
            assert(idx < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.t(idx);
        }
    }

    /**
    * The controlled z gate applies a phase flip to the target qubit if both the 
    * control and target are in the state |1>
    *
    * params:
    * control_qubit_idx = the index of the qubit which determines if the target is affected
    *
    * target_qubit_idx = the index of the qubit which is affected
    */
    void cz(int control_qubit_idx, int target_qubit_idx) {
        assert(this.num_qubits >= 2, "The number of qubits must be greater than or equal to two in order to use controlled gates");

        update_visualization_arr("CZ", [control_qubit_idx, target_qubit_idx]);

        for (int i = 0; i < this.state.length(); i++) {
            bool cntl_qubit_is_one = (i & (1 << control_qubit_idx)) != 0;
            bool tgt_qubit_is_one = (i & (1 << target_qubit_idx)) != 0;
            if (cntl_qubit_is_one && tgt_qubit_is_one) {
                this.state[i] = this.state[i] * Complex!real(-1, 0);
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(control_qubit_idx, 40);
        apply_decoherence(target_qubit_idx, 40);
    }

    /**
    * Overload of the controlled z gate to apply to multiple qubit pairs at once
    *
    * params:
    * qubit_idxs = An array of tuples of qubit indices with (int, int) pairs where
    *              index 0 is control and index 1 is target
    */
    void cz(Tuple!(int, int)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits && idx_tuple[1] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount specified for the system");
            this.cz(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The SWAP gate takes two qubits and if their states are different at index i it calculates a
    * new position j to swap the amplitudes of two states.
    *
    * params:
    * qubit1 = the first qubit to be swapped by the gate
    *
    * qubit2 = the second qubit to be swapped by the gate
    */
    void swap(int qubit1, int qubit2) {
        assert(this.num_qubits >= 2, "The number of qubits must be greater than or equal to two in order to use the swap gates");

        update_visualization_arr("SWAP", [qubit1, qubit2]);

        for (int i = 0; i < this.state.length(); i++) {
            int qubit1_val = (i >> qubit1) & 1;
            int qubit2_val = (i >> qubit2) & 1;
            if (qubit1_val != qubit2_val) {
                int j = i ^ ((1 << qubit1) | (1 << qubit2));
                if (i < j) {
                    Complex!real temp = this.state[i];
                    this.state[i] = this.state[j];
                    this.state[j] = temp;
                }
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit1, 50);
        apply_decoherence(qubit2, 50);
    }

    /**
    * Overload of the swap gate to apply it to multiple qubit pairs at once
    *
    * params:
    * qubit_idxs = An array of qubit indices as tuples of (int, int) pairs
    */
    void swap(Tuple!(int, int)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits && idx_tuple[1] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount specified for the system");
            this.swap(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The iSWAP gate does the same thing as the SWAP gate but also multiplies the amplitudes
    * of the states at index i and j by 0+1i
    *
    * params:
    * qubit1 = the first qubit to be swapped by the gate
    *
    * qubit2 = the second qubit to be swapped by the gate
    */
    void iswap(int qubit1, int qubit2) {
        assert(this.num_qubits >= 2, "The number of qubits must be greater than or equal to two in order to use the swap gates");

        update_visualization_arr("iSWAP", [qubit1, qubit2]);

        for (int i = 0; i < this.state.length(); i++) {
            int qubit1_val = (i >> qubit1) & 1;
            int qubit2_val = (i >> qubit2) & 1;
            if (qubit1_val != qubit2_val) {
                int j = i ^ ((1 << qubit1) | (1 << qubit2));
                if (i < j) {
                    Complex!real temp = this.state[i];
                    this.state[i] = this.state[j];
                    this.state[j] = temp;
                    this.state[i] = this.state[i] * Complex!real(0, 1);
                    this.state[j] = this.state[j] * Complex!real(0, 1);

                }
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit1, 50);
        apply_decoherence(qubit2, 50);
    }

    /**
    * Overload of the iswap gate to apply it to multiple qubit pairs at once
    *
    * params:
    * qubit_idxs = An array of qubit indices as tuples of (int, int) pairs
    */
    void iswap(Tuple!(int, int)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits && idx_tuple[1] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount specified for the system");
            this.iswap(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The Rx gate rotates the state vector around the x axis of the bloch sphere by some angle
    * theta in radians. If theta is PI then it rotates the qubit 180 degrees, essentially flipping
    * it like a pauli-x gate. If theta is PI/2 then it creates an equal superposition over the |0>
    * and |1> states but with a specific phase shift applied. If theta is 0 then the qubit does not 
    * change at all.
    *
    * params:
    * qubit_idx = the index of the qubit to be affected by the gate
    *
    * theta = the angle to rotate the qubit by in radians
    */
    void rx(int qubit_idx, real theta) {
        update_visualization_arr("RX", [qubit_idx]);

        Complex!real c = Complex!real(cos(theta / 2), 0);
        Complex!real s = Complex!real(0, -1) * Complex!real(sin(theta / 2), 0);
        Vector!(Complex!real) psi = Vector!(Complex!real)(cast(int) this.state.length(), new Complex!real[this
                .state.length()]);

        // The .init value of psi without this loop will be nan+nani for all elements
        for (int i = 0; i < psi.length(); i++) {
            psi[i] = Complex!real(0, 0);
        }

        for (int i = 0; i < this.state.length(); i++) {
            int j = i ^ (1 << qubit_idx);

            if (i < j) {
                Complex!real a = this.state[i];
                Complex!real b = this.state[j];
                psi[i] = c * a + s * b;
                psi[j] = s * a + c * b;
            }
        }
        this.state = psi;

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 30);
    }

    /**
    * Overload of the Rx gate to apply it to multiple qubits at once with different values of theta
    *
    * params:
    * qubit_idxs = An array of qubit indices and theta values in tuples of (int, real) pairs
    */
    void rx(Tuple!(int, real)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.rx(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The Ry gate rotates the state vector by an angle theta in radiansaround the y-axis in the bloch sphere. 
    * The main difference between the Rx gate and this one is that this one does not introduce any imaginary 
    * values into the amplitudes.
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    * 
    * theta = the angle in radians to rotate the qubit around the y-axis
    */
    void ry(int qubit_idx, real theta) {
        update_visualization_arr("RY", [qubit_idx]);

        Complex!real c = Complex!real(cos(theta / 2), 0);
        Complex!real s = Complex!real(sin(theta / 2), 0);

        Vector!(Complex!real) psi = Vector!(Complex!real)(cast(int) this.state.length(), new Complex!real[this
                .state.length()]);

        // The .init value of psi without this loop will be nan+nani for all elements
        for (int i = 0; i < psi.length(); i++) {
            psi[i] = Complex!real(0, 0);
        }

        for (int i = 0; i < this.state.length(); i++) {
            int j = i ^ (1 << qubit_idx);
            if (i < j) {
                Complex!real a = this.state[i];
                Complex!real b = this.state[j];

                psi[i] = c * a - s * b;
                psi[j] = s * a + c * b;
            }
        }
        this.state = psi;

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 30);
    }

    /**
    * Overload of the Ry gate to apply it to multiple qubits at once with different values of theta
    *
    * params:
    * qubit_idxs = An array of qubit indices and theta values in tuples of (int, real) pairs
    */
    void ry(Tuple!(int, real)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.ry(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The Rz gate applies a phase shift to the target qubit based on its state. If the target qubit is 
    * in the state |0> then it applies a phase shift of e^-i(theta/2). If the qubit is in the state |1>
    * then it applies a phase shift of e^i(theta/2).
    *
    * params:
    * qubit_idx = the index of the qubit to affect
    *
    * theta = the angle in radians to apply to the phase shift exponential
    */
    void rz(int qubit_idx, real theta) {
        update_visualization_arr("RZ", [qubit_idx]);

        Complex!real z0 = exp(Complex!real(0, -1) * Complex!real(theta / 2, 0));
        Complex!real z1 = exp(Complex!real(0, 1) * Complex!real(theta / 2, 0));

        for (int i = 0; i < this.state.length(); i++) {
            int qubit_value = (i >> qubit_idx) & 1;
            if (qubit_value == 0) {
                this.state[i] = this.state[i] * z0;
            } else if (qubit_value == 1) {
                this.state[i] = this.state[i] * z1;
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(qubit_idx, 30);
    }

    /**
    * Overload of the Rz gate to apply it to multiple qubits at once with different values of theta
    *
    * params:
    * qubit_idxs = An array of qubit indices and theta values in tuples of (int, real) pairs
    */
    void rz(Tuple!(int, real)[] qubit_idxs) {
        foreach (idx_tuple; qubit_idxs) {
            assert(idx_tuple[0] < this.num_qubits,
                "One or more of the qubit indices is beyond the amount you specified for the system");
            this.rz(idx_tuple[0], idx_tuple[1]);
        }
    }

    /**
    * The CR_k gate or controlled rotation of order k gate, rotate the phase by e^2 * PI / 2^k
    * if and only if the control and target qubits are 1
    *
    * params:
    * control_qubit_idx = the index of the qubit which determines if the target is affected
    *
    * target_qubit_idx = the index of the qubit which is affected by the control's state
    *
    * k = the exponent k to apply in the phase factor
    *
    * inverse = whether or not to invert the gate, this gate is not hermittian so it is not it's
    *           own inverse
    */
    void cr(int control_qubit_idx, int target_qubit_idx, int k, bool inverse = false) {
        assert(this.num_qubits >= 2, "The number of qubits must be greater than or equal to two in order to use the controlled rotation gate");

        update_visualization_arr("CR", [control_qubit_idx, target_qubit_idx]);

        for (int i = 0; i < this.state.length(); i++) {
            int cntl_qubit_val = (i >> control_qubit_idx) & 1;
            int tgt_qubit_val = (i >> target_qubit_idx) & 1;

            if (cntl_qubit_val == 1 && tgt_qubit_val == 1) {
                if (!inverse) {
                    this.state[i] = this.state[i] * expi(2 * PI / pow(2.0, k));
                } else {
                    this.state[i] = this.state[i] * expi(-2 * PI / pow(2.0, k));
                }
            }
        }

        // This will only happen if DecoherenceConfig.decoherence_mode is
        // set to automatic
        apply_decoherence(control_qubit_idx, 70);
        apply_decoherence(target_qubit_idx, 70);
    }

    /**
    * Computes the expectation value of an observable on the current quantum state of the system
    * 
    * params:
    * observable = The observable affecting the quantum system as a linear combination of weighted 
    *              pauli operator kronecker products
    *
    * returns: A real value, the average measurement value or expectation value
    */
    real expectation_value(Observable observable) {
        Matrix!(Complex!real) psi_dagger = this.state.dagger();
        Vector!(Complex!real) phi = observable.apply(this.state);

        real result = psi_dagger.inner_product(phi);

        return result;
    }

    // measurement of a single qubit internal logic, this function
    // exists solely to prevent code duplication
    private string measure_internal(int qubit_idx) {
        real probability_0 = 0;
        real probability_1 = 0;

        // sum probabilities of the qubit in each state
        for (int i = 0; i < this.state.length(); i++) {
            bool qubit_is_zero = (i & (1 << qubit_idx)) == 0;
            bool qubit_is_one = (i & (1 << qubit_idx)) != 0;

            if (qubit_is_zero) {
                real state_prob = norm(this.state[i]);
                probability_0 += state_prob;
            } else if (qubit_is_one) {
                real state_prob = norm(this.state[i]);
                probability_1 += state_prob;
            }
        }

        // get a random number over a uniform distribution
        auto rng = Random(unpredictableSeed);
        auto r = uniform(0.0, 1.0f, rng);

        int result;

        // determine measurement result
        if (r < probability_0) {
            result = 0;
        } else if (r >= probability_0) {
            result = 1;
        }

        return format("%d", result);
    }

    /**
    * Measure the state of one qubit
    *
    * params:
    * qubit_idx = The index of the qubit to measure
    *
    * returns: A string representing the measured state of the qubit
    */
    string measure(int qubit_idx) {
        string result = measure_internal(qubit_idx);
        return result;
    }

    /**
    * Overload of the measure function to measure the qubit
    * many times to see the probabilistic outcomes
    *
    * params:
    * qubit_idx = The index of the qubit to measure
    * 
    * shots = The amount of times to measure the qubit
    *
    * returns: A string to int map, representing the state 
    *          measured and the amount of times it was measured
    */
    int[string] measure(int qubit_idx, int shots) {
        assert(shots >= 2,
            "using this overload of the measure function requires shots to be greater than or equal to 2, it is recommended to use over a 1000");
        int[string] counts;
        for (int i = 0; i < shots; i++) {
            string result = measure_internal(qubit_idx);
            counts[result] += 1;
        }
        return counts;
    }

    // measurement for the entire system internal logic, this function exists solely 
    // to prevent code duplication
    private string measure_all_internal() {
        Vector!float probs = Vector!float(cast(int) this.state.length(), new float[this
                .state.length()]); // Take the magnitude of each complex probability amplitude
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

    /**
    * Collapses the possible superposition of basis states into one classical state 
    * based on inverse transform sampling (https://en.wikipedia.org/wiki/Inverse_transform_sampling)
    *
    * returns: the bitstring of the state which was measured probabilistically
    */
    string measure_all() {
        string binary_result = measure_all_internal();
        return binary_result;
    }

    /**
    * Overload of the measure() function with shots parameter, to be able to see
    * statistical variation in measurement results
    *
    * params:
    * shots = number of times measurement should be preformed
    *
    * returns: An associative array of bitstring to amount of times it was measured
    */
    int[string] measure_all(int shots) {
        assert(shots >= 2,
            "using this overload of the measure function requires shots to be greater than or equal to 2, it is recommended to use over a 1000");
        int[string] counts;
        for (int i = 0; i < shots; i++) {
            string binary_result = measure_all_internal();
            counts[binary_result] += 1;
        }
        return counts;
    }

    /**
    * Draws the circuit which the user created with latex
    *
    * params:
    * compiler = The name of the latex compiler to use (default: pdflatex)
    *
    * filename = The name of the file to write the latex to and to compile (default: circuit.tex)
    */
    void draw(string compiler = "pdflatex", string filename = "circuit.tex") {
        Visualization vis = Visualization(this.visualization_arr, this.num_qubits, this
                .initial_state_idx);
        vis.parse_and_write_vis_arr(filename);
        vis.compile_tex_and_cleanup(compiler, filename);
    }
}
