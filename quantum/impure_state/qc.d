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
import std.random;
import std.format;
import std.range;

// linear algebra modules
import linalg.vector;
import linalg.matrix;

// quantum related modules
import quantum.impure_state.observable;

//visualization modules
import viz.visualization;

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
        this.initial_state_idx = starting_state_idx;

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

    private Matrix!(Complex!real) build_full_gate(Matrix!(Complex!real) gate_matrix, int qubit_idx) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);
        Matrix!(Complex!real)[] kronecker_chain;

        for (int i = this.num_qubits - 1; i >= 0; i--) {
            if (i == qubit_idx) {
                kronecker_chain[kronecker_chain.length++] = gate_matrix;
            } else {
                kronecker_chain[kronecker_chain.length++] = identity;
            }
        }

        Matrix!(Complex!real) result = kronecker_chain[0];
        for (int i = 1; i < kronecker_chain.length; i++) {
            result = result.kronecker(kronecker_chain[i]);
        }

        return result;
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

    private Matrix!(Complex!real) build_full_pauli_op(Matrix!(Complex!real) pauli_op, int qubit1, int qubit2) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);
        Matrix!(Complex!real) full_pauli_op;

        for (int i = 0; i < this.num_qubits; i++) {
            if (i == 0 && i != qubit1 && i != qubit2) {
                full_pauli_op = identity;
            } else if (i == 0 && i == qubit1 && i != qubit2) {
                full_pauli_op = pauli_op;
            } else if (i == 0 && i != qubit1 && i == qubit2) {
                full_pauli_op = pauli_op;
            } else if (i != 0 && (i == qubit1 || i == qubit2)) {
                full_pauli_op = full_pauli_op.kronecker(pauli_op);
            } else {
                full_pauli_op = full_pauli_op.kronecker(identity);
            }
        }

        return full_pauli_op;
    }

    private Matrix!(Complex!real) build_full_iswap(Matrix!(Complex!real) iswap_mat, int qubit1, int qubit2) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);
        Matrix!(Complex!real) full_iswap;

        for (int i = 0; i < this.num_qubits; i++) {
            if (i == qubit1) {
                if (i == 0) {
                    full_iswap = iswap_mat;
                } else {
                    full_iswap = full_iswap.kronecker(iswap_mat);
                }
            } else if (i == qubit2) {
                continue;
            } else {
                if (i == 0) {
                    full_iswap = identity;
                } else {
                    full_iswap = full_iswap.kronecker(identity);
                }
            }
        }

        return full_iswap;
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

        Matrix!(Complex!real) result = build_full_gate(hadamard, qubit_idx);
        this.density_mat = result.mult_mat(this.density_mat).mult_mat(result.dagger());
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

        Matrix!(Complex!real) pauli_x = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) pauli_x_op = build_full_gate(pauli_x, qubit_idx);
        this.density_mat = pauli_x_op.mult_mat(this.density_mat).mult_mat(pauli_x_op.dagger());
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

        Matrix!(Complex!real) pauli_y = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, -1)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 1), Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) pauli_y_op = build_full_gate(pauli_y, qubit_idx);
        this.density_mat = pauli_y_op.mult_mat(this.density_mat).mult_mat(pauli_y_op.dagger());
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

        Matrix!(Complex!real) pauli_z = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(-1, 0)
                    ])
            ]);

        Matrix!(Complex!real) pauli_z_op = build_full_gate(pauli_z, qubit_idx);
        this.density_mat = pauli_z_op.mult_mat(this.density_mat).mult_mat(pauli_z_op.dagger());
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

    /**
    * The S phase shift gate or PI/4 gate applies a phase shift of PI/4 to the state |1>
    *
    * params:
    * qubit_idx = the index of the qubit to be affected
    */
    void s(int qubit_idx) {
        update_visualization_arr("S", [qubit_idx]);

        Matrix!(Complex!real) s = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, 1)
                    ])
            ]);

        Matrix!(Complex!real) s_op = build_full_gate(s, qubit_idx);
        this.density_mat = s_op.mult_mat(this.density_mat).mult_mat(s_op.dagger());
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

        Matrix!(Complex!real) t = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), expi(PI / 4)
                    ])
            ]);

        Matrix!(Complex!real) t_op = build_full_gate(t, qubit_idx);
        this.density_mat = t_op.mult_mat(this.density_mat).mult_mat(t_op.dagger());
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
        assert(this.num_qubits >= 2,
            "The number of qubits must be greater than or equal to two in order to use controlled gates");

        update_visualization_arr("CZ", [control_qubit_idx, target_qubit_idx]);

        Matrix!(Complex!real) pauli_z = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(-1, 0)
                    ])
            ]);

        Matrix!(Complex!real) cz_op = build_full_controlled_gate(pauli_z, control_qubit_idx, target_qubit_idx);
        this.density_mat = cz_op.mult_mat(this.density_mat).mult_mat(cz_op.dagger());
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
    void swap(int qubit1, int qubit2, bool visualize = true) {
        assert(this.num_qubits >= 2,
            "The number of qubits must be greater than or equal to two in order to use the swap gates");

        if (visualize) {
            update_visualization_arr("SWAP", [qubit1, qubit2]);
        }

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

        Matrix!(Complex!real) full_identity_op = build_full_pauli_op(identity, qubit1, qubit2);
        Matrix!(Complex!real) full_pauli_x_op = build_full_pauli_op(pauli_x, qubit1, qubit2);
        Matrix!(Complex!real) full_pauli_y_op = build_full_pauli_op(pauli_y, qubit1, qubit2);
        Matrix!(Complex!real) full_pauli_z_op = build_full_pauli_op(pauli_z, qubit1, qubit2);

        Matrix!(Complex!real) full_swap = full_identity_op
            .add_mat(full_pauli_x_op)
            .add_mat(full_pauli_y_op)
            .add_mat(full_pauli_z_op);

        full_swap = full_swap.mult_scalar(Complex!real(0.5, 0));

        this.density_mat = full_swap.mult_mat(this.density_mat).mult_mat(full_swap.dagger());
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
        assert(this.num_qubits >= 2,
            "The number of qubits must be greater than or equal to two in order to use the swap gates");

        update_visualization_arr("iSWAP", [qubit1, qubit2]);

        Matrix!(Complex!real) iswap = Matrix!(Complex!real)(4, 4, [
                Vector!(Complex!real)(4, [
                        Complex!real(1, 0), Complex!real(0, 0),
                        Complex!real(0, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(4, [
                        Complex!real(0, 0), Complex!real(0, 0),
                        Complex!real(0, 1), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(4, [
                        Complex!real(0, 0), Complex!real(0, 1),
                        Complex!real(0, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(4, [
                        Complex!real(0, 0), Complex!real(0, 0),
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
            ]);

        Tuple!(int, int)[] swap_seq;

        if (qubit1 > qubit2) {
            int tmp = qubit1;
            qubit1 = qubit2;
            qubit2 = tmp;
        }

        if (qubit1 + 1 < qubit2) {
            while (qubit1 + 1 < qubit2) {
                this.swap(qubit1, qubit1 + 1, false);
                swap_seq[swap_seq.length++] = tuple(qubit1, qubit1 + 1);
                qubit1 += 1;
            }
        }

        Matrix!(Complex!real) full_iswap = build_full_iswap(iswap, qubit1, qubit2);

        this.density_mat = full_iswap.mult_mat(this.density_mat).mult_mat(full_iswap.dagger());

        for (int i = cast(int) swap_seq.length - 1; i >= 0; i--) {
            Tuple!(int, int) item = swap_seq[i];
            this.swap(item[0], item[1], false);
        }
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

        Matrix!(Complex!real) rx = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(cos(theta / 2), 0),
                        Complex!real(0, -1) * sin(theta / 2)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, -1) * sin(theta / 2),
                        Complex!real(cos(theta / 2), 0)
                    ])
            ]);

        Matrix!(Complex!real) rx_op = build_full_gate(rx, qubit_idx);
        this.density_mat = rx_op.mult_mat(this.density_mat).mult_mat(rx_op.dagger());
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

        Matrix!(Complex!real) ry = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(cos(theta / 2), 0),
                        Complex!real(-sin(theta / 2), 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(sin(theta / 2), 0),
                        Complex!real(cos(theta / 2), 0)
                    ])
            ]);

        Matrix!(Complex!real) ry_op = build_full_gate(ry, qubit_idx);
        this.density_mat = ry_op.mult_mat(this.density_mat).mult_mat(ry_op.dagger());
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

        Matrix!(Complex!real) rz = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        exp(Complex!real(0, -1) * theta), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), exp(Complex!real(0, 1) * theta)
                    ])
            ]);

        Matrix!(Complex!real) rz_op = build_full_gate(rz, qubit_idx);
        this.density_mat = rz_op.mult_mat(this.density_mat).mult_mat(rz_op.dagger());
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
        update_visualization_arr("CR", [control_qubit_idx, target_qubit_idx]);

        Complex!real exponential = Complex!real(0, 0);
        if (inverse) {
            exponential = expi(-2 * PI / pow(2.0, k));
        } else {
            exponential = expi(2 * PI / pow(2.0, k));

        }

        Matrix!(Complex!real) cr = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [Complex!real(0, 0), exponential])
            ]);

        Matrix!(Complex!real) cr_op = build_full_controlled_gate(cr, control_qubit_idx, target_qubit_idx);
        this.density_mat = cr_op.mult_mat(this.density_mat).mult_mat(cr_op.dagger());
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
        Matrix!(Complex!real)[] full_pauli_ops = observable.apply(this.density_mat);
        real sum = 0;

        foreach (i, pauli_op; full_pauli_ops) {
            Complex!real trace = Complex!real(pauli_op.trace(), 0);
            Complex!real result = observable.coeffs[i] * trace;
            sum = sum + result.re;
        }

        return sum.re;
    }

    // measurement of a single qubit internal logic, this function exists
    // solely to prevent code duplication
    private string measure_internal(int qubit_idx) {

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

        Matrix!(Complex!real)[] kronecker_chain_p0;
        for (int i = this.num_qubits - 1; i >= 0; i--) {
            if (i == qubit_idx) {
                kronecker_chain_p0[kronecker_chain_p0.length++] = projection_0;
            } else {
                kronecker_chain_p0[kronecker_chain_p0.length++] = identity;
            }
        }

        Matrix!(Complex!real)[] kronecker_chain_p1;
        for (int i = this.num_qubits - 1; i >= 0; i--) {
            if (i == qubit_idx) {
                kronecker_chain_p1[kronecker_chain_p1.length++] = projection_1;
            } else {
                kronecker_chain_p1[kronecker_chain_p1.length++] = identity;
            }
        }

        Matrix!(Complex!real) full_operator_p0 = kronecker_chain_p0[0];
        for (int i = 1; i < kronecker_chain_p0.length; i++) {
            full_operator_p0 = full_operator_p0.kronecker(kronecker_chain_p0[i]);
        }

        Matrix!(Complex!real) full_operator_p1 = kronecker_chain_p1[0];
        for (int i = 1; i < kronecker_chain_p1.length; i++) {
            full_operator_p1 = full_operator_p1.kronecker(kronecker_chain_p1[i]);
        }

        real probability_0 = full_operator_p0.mult_mat(this.density_mat).trace();

        auto rng = Random(unpredictableSeed);
        auto r = uniform(0.0, 1.0f, rng);

        int result;

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
    * returns: A string representing the state of the qubit measured
    */
    string measure(int qubit_idx) {
        update_visualization_arr("M", [qubit_idx]);

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

        update_visualization_arr("M", [qubit_idx]);

        int[string] counts;
        for (int i = 0; i < shots; i++) {
            string result = measure_internal(qubit_idx);
            counts[result] += 1;
        }
        return counts;
    }

    // measurement for the entire system internal logic, this function
    // exists solely to prevent code duplication
    private string measure_all_internal() {
        Vector!(Complex!real) probs = this.density_mat.get_diagonal();

        auto rng = Random(unpredictableSeed);
        auto r = uniform(0.0f, 1.0f, rng);

        float sum = 0;
        string binary_result;
        foreach (i, prob; probs.elems) {
            sum += prob.re;
            if (r < sum) {
                binary_result = format("%0*b", this.num_qubits, i);
                break;
            }
        }

        return binary_result;
    }

    /**
    * Measure the entire system
    *
    * returns: A bitstring representing the final state of the system
    */
    string measure_all() {
        update_visualization_arr("MA", iota(0, this.num_qubits).array);

        string binary_result = measure_all_internal();
        return binary_result;
    }

    /**
    * Overload of the measure_all function to measure the 
    * entire system many times
    *
    * params:
    * shots = The amount of times to measure the entire system
    *
    * returns: An associative array of bitstring to amount of times it was measured
    */
    int[string] measure_all(int shots) {
        assert(shots >= 2,
            "using this overload of the measure function requires shots to be greater than or equal to 2, it is recommended to use over a 1000");

        update_visualization_arr("MA", iota(0, this.num_qubits).array);

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
