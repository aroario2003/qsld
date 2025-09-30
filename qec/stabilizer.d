// NOTE: This module only works with the pure subsystem as doing quantum error correction
// of this kind with the impure subsystem is not possible and would be inaccurate. So if 
// using this subsystem please only use the pure subsystem of QSLD.
// WARNING: This is currently a work in progress, do not use

module qec.stabilizer;

import quantum.pure_state.qc;

import linalg.vector;
import linalg.matrix;

struct QecConfig {
    string qec_mode;

    /**
     * The constructor for the Quantum Error Correction config structure 
     * 
     * params:
     * mode = The way to apply the error correction, "manual" means you do it yourself,
     *        "automatic" means the simulator will do it for you, "none" means no error 
     *        correction.
     */
    this(string mode) {
        assert(mode == "manual" || mode == "automatic" || mode == "none",
            "The mode specified for qec is invalid");

        this.qec_mode = mode;
    }
}

struct Tableau {
    int num_qubits;
    Matrix!int tableau_internal;
    Vector!int error;

    /**
    * The constructor for the tableau object which is part of the stabilizer formalism 
    * simulation. This object keeps track of the current state in a binary matrix and 
    * diagnoses the syndrome based on the binary values.
    *
    * params:
    * num_qubits = The number of qubits in the quantum system.
    */
    this(int num_qubits) {
        this.num_qubits = num_qubits;

        // the number of rows is always 2n, where n is the number of qubits
        int num_rows = 2 * this.num_qubits;
        // the number of cols is always 2n+1 to account for the phase bit
        int num_cols = 2 * this.num_qubits + 1;

        // initialize the matrix which will represent the tableau internally
        Matrix!int tableau_mat = Matrix!int(num_rows, num_cols, []);

        int one_col_idx = 0;
        // build the initial tableau
        for (int i = 0; i < num_rows; i++) {
            Vector!int tableau_row = Vector!int(num_cols, []);
            for (int j = 0; j < num_cols - 1; j++) {
                if (j == one_col_idx) {
                    tableau_row.append(1);
                } else {
                    tableau_row.append(0);
                }
            }
            tableau_row.append(0);
            tableau_mat.append(tableau_row);
            one_col_idx++;
        }

        this.tableau_internal = tableau_mat;
    }

    // flips the phase bit in a given row of the tableau
    private int flip_phase_bit(Vector!int row) {
        int phase_bit = row[row.length() - 1];
        if (phase_bit == 0) {
            phase_bit = 1;
        } else if (phase_bit == 1) {
            phase_bit = 0;
        }
        return phase_bit;
    }

    /** 
     * Update the tableau to represent the state when the hadamard
     * gate is applied
     *
     * params: 
     * qubit_idx = The index of the qubit to update
    */
    void update_hadamard(int qubit_idx) {
        // iterate rows in the tableau
        foreach (row; this.tableau_internal.rows) {
            // get the value of the x and z bits
            int x_bit = row[qubit_idx];
            int z_bit = row[qubit_idx + this.num_qubits];

            // flip the phase bit if both x and z bits are 1
            if (x_bit == 1 && z_bit == 1) {
                int phase_bit = flip_phase_bit(row);
                row[row.length() - 1] = phase_bit;
            }

            // swap all x and z bits
            int temp = x_bit;
            x_bit = z_bit;
            z_bit = temp;

            // this is now the value of the z bit
            row[qubit_idx] = x_bit;
            // this is now the value of the x bit
            row[qubit_idx + this.num_qubits] = z_bit;

        }
    }

    /**
     * Updates the tableau when a pauli-x gate is applied
     * by flipping the phase bit if certain conditions are true
     *
     * params: 
     * qubit_idx = The index of the qubit to update
    */
    void update_pauli_x(int qubit_idx) {
        foreach (row; this.tableau_internal.rows) {
            // get the value of the x and z bits
            int x_bit = row[qubit_idx];
            int z_bit = row[qubit_idx + this.num_qubits];

            if (x_bit == 0 && z_bit == 1) {
                int phase_bit = flip_phase_bit(row);
                row[row.length() - 1] = phase_bit;
            }
        }
    }

    /**
     * Updates the tableau when a pauli-y gate is applied 
     * by flipping the phase bit if certain conditions are true
     * 
     * params;
     * qubit_idx = The index of the qubit to update
    */
    void update_pauli_y(int qubit_idx) {
        foreach (row; this.tableau_internal.rows) {
            // get the value of the x and z bits
            int x_bit = row[qubit_idx];
            int z_bit = row[qubit_idx + this.num_qubits];

            if (x_bit == 1 && z_bit == 0 || x_bit == 0 && z_bit == 1) {
                int phase_bit = flip_phase_bit(row);
                row[row.length() - 1] = phase_bit;
            }
        }
    }

    /**
     * Updates the tableau when a pauli-z gate is applied 
     * by flipping the phase bit if certain conditions are true
     *
     * params:
     * qubit_idx = The index of the qubit to update
    */
    void update_pauli_z(int qubit_idx) {
        foreach (row; this.tableau_internal.rows) {
            // get the value of the x and z bits
            int x_bit = row[qubit_idx];
            int z_bit = row[qubit_idx + this.num_qubits];

            if (x_bit == 1 && z_bit == 0 || x_bit == 1 && z_bit == 1) {
                int phase_bit = flip_phase_bit(row);
                row[row.length() - 1] = phase_bit;
            }
        }
    }

    /**
     * Updates the tableau when the s gate is applied
     * by flipping certain bits in the tableau based on
     * certain conditions
     *
     * params;
     * qubit_idx = The index of the qubit to update
    */
    void update_s(int qubit_idx) {
        foreach (row; this.tableau_internal.rows) {
            // get the value of the x and z bits
            int x_bit = row[qubit_idx];
            int z_bit = row[qubit_idx + this.num_qubits];

            if (x_bit == 1 && z_bit == 0) {
                // set z bit to 1
                row[qubit_idx + this.num_qubits] = 1;
            } else if (x_bit == 1 && z_bit == 1) {
                // set z bit to 0
                row[qubit_idx + this.num_qubits] = 0;
                // flip the phase bit
                int phase_bit = flip_phase_bit(row);
                row[row.length() - 1] = phase_bit;
            }
        }
    }

    /**
     * Updates the talbeau when the cnot gate is applied by flipping
     * certain bits in the tableau based on certain conditions
     *
     * params:
     * control_qubit_idx = The index of the control qubit to update
     * target_qubit_idx = The index of the target qubit to update
    */
    void update_cnot(int control_qubit_idx, int target_qubit_idx) {
        foreach (row; this.tableau_internal.rows) {
            int x_bit_control = row[control_qubit_idx];
            int z_bit_target = row[target_qubit_idx + this.num_qubits];

            if (x_bit_control == 1) {
                row[target_qubit_idx] = 1;
            }

            if (z_bit_target == 1) {
                row[control_qubit_idx + this.num_qubits] = 1;
            }

            if (x_bit_control == 1 && z_bit_target == 1) {
                int phase_bit = flip_phase_bit(row);
                row[row.length() - 1] = phase_bit;
            }
        }
    }

    /**
     * Updates the tableay when the cz gate is applied by xoring 
     * the z control and target bits together to get new values if 
     * applicable
     *
     * params:
     * control_qubit_idx = The index of the control qubit to update
     * target_qubit_idx = The index of the target qubit to update
    */
    void update_cz(int control_qubit_idx, int target_qubit_idx) {
        foreach (row; this.tableau_internal.rows) {
            row[this.num_qubits + control_qubit_idx] = row[this.num_qubits + control_qubit_idx] ^ row[target_qubit_idx];
            row[this.num_qubits + target_qubit_idx] = row[this.num_qubits + target_qubit_idx] ^ row[control_qubit_idx];
        }
    }

    // Used to put the tableau matrix into upper echelon form
    // and find independent stabilizers used to correct the state
    private Matrix!int gaussian_elimination() {
        ulong row_len = this.tableau_internal.rows[0].length - (cast(ulong) 2);
        ulong r = 0;

        bool pivot_found = false;

        Vector!int pivot_row;
        Matrix!int independent_rows;

        for (ulong col_idx = 0; col_idx <= row_len; col_idx++) {
            for (ulong i = this.num_qubits; i < 2 * this.num_qubits; i++) {
                Vector!int row = this.tableau_internal.rows[i];

                if (row[col_idx] == 1) {
                    pivot_found = true;
                    pivot_row = row;
                    independent_rows.append(pivot_row);

                    Vector!int temp = this.tableau_internal.rows[this.num_qubits + r];
                    this.tableau_internal.rows[this.num_qubits + r] = this.tableau_internal.rows[i];
                    this.tableau_internal.rows[i] = temp;
                    break;
                }
            }

            if (pivot_found) {
                for (int i = this.num_qubits; i < 2 * this.num_qubits; i++) {
                    Vector!int row = this.tableau_internal.rows[i];
                    if (i != (this.num_qubits + r)) {
                        if (row[col_idx] == 1) {
                            for (int x = 0; x <= row_len; x++) {
                                this.tableau_internal.rows[i].elems[x] = row[x] ^ pivot_row[x];
                            }
                        }
                    }
                }
                r++;
                pivot_found = false;
            }
        }

        return independent_rows;
    }

    /**
     * Extracts the syndrome from the current state of the tableau
     *
     * params:
     * error = The vector representing the error E applied to the quantum state
     * 
     * returns: The vector of syndrome bits which represent the error
     */
    Vector!int measure() {
        Matrix!int stabilizers = gaussian_elimination();
        Vector!int syndrome_bits;
        foreach (stabilizer; stabilizers.rows) {
            int[] x_bits_s = stabilizer.elems[0 .. this.num_qubits];
            int[] z_bits_s = stabilizer.elems[this.num_qubits .. 2 * this.num_qubits];

            int[] x_bits_e = this.error.elems[0 .. this.num_qubits];
            int[] z_bits_e = this.error.elems[this.num_qubits .. 2 * this.num_qubits];

            int sum = 0;
            for (int bit_idx = 0; bit_idx < this.num_qubits; bit_idx++) {
                int result1 = x_bits_s[bit_idx] * z_bits_e[bit_idx];
                int result2 = x_bits_e[bit_idx] * z_bits_s[bit_idx];

                sum += result1 + result2;
            }

            sum %= 2;
            syndrome_bits.append(sum);
        }

        return syndrome_bits;
    }

    /*
     *--------------------------------------------------------------------------------
     * NOTE: The functions below propogate the errors applied to qubits when certain gates
     * are applied to the quantum state. These functions update the error vector not the
     * tableau object. The pauli operators are not present due to the fact that they commute
     * with themselves.
     *--------------------------------------------------------------------------------
     */

    /**
     * Propogate the error on a qubit if the qubit at the specified index (qubit_idx) 
     * already has an error applied to it. Used for the hadamard gate.
     *
     * params:
     * qubit_idx = The index of the qubit to propogate the error of, if it exists
     */
    void propogate_hadamard(int qubit_idx) {
        int temp = this.error[qubit_idx];
        this.error[qubit_idx] = this.error[this.num_qubits + qubit_idx];
        this.error[this.num_qubits + qubit_idx] = temp;
    }

    /**
     * Propogate the error on a qubit if the qubit at the specified index (qubit_idx)
     * already has an error applied to it. Used for the s gate.
     *
     * params:
     * qubit_idx = The index of the qubit to propogate the error of, if it exists
     */
    void propogate_s(int qubit_idx) {
        this.error[this.num_qubits + qubit_idx] = this.error[this.num_qubits + qubit_idx] ^ this
            .error[qubit_idx];
    }

    /**
     * Propogate the error on a qubit if the qubit at the specified index (control_qubit_idx or target_qubit_idx)
     * already has an error applied to it. Used for the controlled not gate.
     *
     * params:
     * control_qubit_idx = The index of the control qubit to propogate the error of, if it exists
     * target_qubit_idx = The index of the target qubit to propogate the error of, if it exists
     */
    void propogate_cnot(int control_qubit_idx, int target_qubit_idx) {
        this.error[target_qubit_idx] = this.error[target_qubit_idx] ^ this
            .error[control_qubit_idx];

        this.error[this.num_qubits + control_qubit_idx] = this.error[this.num_qubits + control_qubit_idx] ^ this
            .error[this.num_qubits + target_qubit_idx];
    }

    /**
     * Propogate the error on a qubit if the qubit at the specified index (control_qubit_idx or target_qubit_idx)
     * already has an error applied to it. Used for the controlled z gate.
     *
     * params:
     * control_qubit_idx = The index of the control qubit to propogate the error of, if it exists
     * target_qubit_idx = The index of the target qubit to propogate the error of, if it exists
     */
    void propogate_cz(int control_qubit_idx, int target_qubit_idx) {
        this.error[this.num_qubits + control_qubit_idx] = this
            .error[this.num_qubits = control_qubit_idx] ^ this.error[target_qubit_idx];

        this.error[this.num_qubits + target_qubit_idx] = this
            .error[this.num_qubits + target_qubit_idx] ^ this.error[control_qubit_idx];
    }
}
