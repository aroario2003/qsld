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

module quantum.impure_state.gate_noise;

import std.random;
import std.math;
import std.conv;
import std.algorithm.iteration;
import std.complex;
import std.uni;
import std.array;

import std.algorithm : canFind, countUntil;
import std.range : repeat;

import linalg.vector;
import linalg.matrix;

import quantum.impure_state.qc;

struct GateNoise {
    QuantumCircuit* qc;
    Random rng;

    /**
    * The constructor for the GateNoise object
    *
    * params:
    * qc = The QuantumCircuit object being used for simulation so that gate noise can
    *      affect it.
    */
    this(QuantumCircuit* qc) {
        this.qc = qc;
        this.rng = Random(unpredictableSeed);
    }

    // For a given symbol specifically one of the characters representing
    // a pauli gate, return the corresponding 2x2 gate matrix
    private Matrix!(Complex!real) get_gate_matrix(char gate_sym) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

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

        Matrix!(Complex!real) result;
        switch (gate_sym) {
        case 'I':
            result = identity;
            break;
        case 'X':
            result = pauli_x;
            break;
        case 'Y':
            result = pauli_y;
            break;
        case 'Z':
            result = pauli_z;
            break;
        default:
            assert(false, "Invalid gate symbol provided");
        }

        return result;
    }

    // For a given array of qubit indices and pauli combinations with probabilities create the full
    // kraus operators with the probabilities encoded for each pauli combination and associated 
    // probability
    private Matrix!(Complex!real)[] build_full_krauss(int[] qubit_idxs, float[string] probability_map) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

        Matrix!(Complex!real)[] kronecker_chain = new Matrix!(Complex!real)[this.qc.num_qubits];
        Matrix!(Complex!real)[] krauss;

        foreach (combo, prob; probability_map) {
            char[] pauli_ops = combo[].dup;
            for (int i = this.qc.num_qubits - 1; i >= 0; i--) {
                if (qubit_idxs.canFind(i)) {
                    ulong idx = qubit_idxs.countUntil!(idx => idx == i);
                    kronecker_chain[i] = get_gate_matrix(pauli_ops[idx]);
                } else {
                    kronecker_chain[i] = identity;
                }
            }

            Matrix!(Complex!real) result = kronecker_chain[0];
            for (int i = 1; i < kronecker_chain.length; i++) {
                result = result.kronecker(kronecker_chain[i]);
            }

            result = result.mult_scalar(Complex!real(sqrt(prob), 0));

            krauss ~= result;
            kronecker_chain = new Matrix!(Complex!real)[this.qc.num_qubits];
        }

        return krauss;
    }

    // Embed the pauli operator into the full hilbert space affecting multiple qubits
    // with the same operator, this is used for bit and phase flip noise
    private Matrix!(Complex!real) build_full_pauli_op(Matrix!(Complex!real) gate, int[] qubit_idxs, float probability) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);
        Matrix!(Complex!real)[] kronecker_chain = new Matrix!(Complex!real)[this.qc.num_qubits];
        for (int i = this.qc.num_qubits - 1; i >= 0; i--) {
            if (qubit_idxs.canFind(i)) {
                ulong idx = qubit_idxs.countUntil!(idx => idx == i);
                kronecker_chain[i] = gate;
            } else {
                kronecker_chain[i] = identity;
            }
        }

        Matrix!(Complex!real) result = kronecker_chain[0];
        for (int i = 1; i < kronecker_chain.length; i++) {
            result = result.kronecker(kronecker_chain[i]);
        }

        result = result.mult_scalar(Complex!real(sqrt(probability), 0));

        return result;
    }

    /**
    * Generates all possible combinations of pauli operators for a number of qubits
    *
    * params:
    * qubit_num = The number of qubits to generate combinations for 
    *
    * returns: The list of combinations 
    */
    string[] generate_pauli_combos(int qubit_num) {
        assert(qubit_num > 0, "The amount of qubits should be greater than 0");

        string[] pauli_ops = ["I", "X", "Y", "Z"];
        string[] pauli_combos;

        for (int i = 0; i < pow(4, qubit_num); i++) {
            string combo = "";
            int temp = i;
            for (int j = 0; j < qubit_num; j++) {
                combo ~= pauli_ops[temp % 4];
                temp /= 4;
            }
            pauli_combos ~= combo;
        }
        string item = "I".repeat(qubit_num).to!string;
        pauli_combos.filter!(combo => combo != item);
        return pauli_combos;
    }

    /**
    * Apply depolarizing noise to a qubit. The qubit should be the one that the gate
    * acted on. Applying it to an irrelevant qubit will lead to inaccurate results.
    *
    * params:
    * qubit_idx = The qubit acted on by the gate as an index
    *
    * probability = The probability that depolarizing noise is applied to the qubit
    */
    void depolarizing_noise(int qubit_idx, float probability) {
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

        // make the full pauli operator matrices
        Matrix!(Complex!real) full_x = this.qc.build_full_gate(pauli_x, qubit_idx);
        Matrix!(Complex!real) full_y = this.qc.build_full_gate(pauli_y, qubit_idx);
        Matrix!(Complex!real) full_z = this.qc.build_full_gate(pauli_z, qubit_idx);

        // apply the full pauli operator matrices to the density matrix
        Matrix!(Complex!real) applied_x = full_x.mult_mat(this.qc.density_mat)
            .mult_mat(full_x.dagger());
        Matrix!(Complex!real) applied_y = full_y.mult_mat(this.qc.density_mat)
            .mult_mat(full_y.dagger());
        Matrix!(Complex!real) applied_z = full_z.mult_mat(this.qc.density_mat)
            .mult_mat(full_z.dagger());

        // apply the the applied pauli operators to themselves
        Matrix!(Complex!real) full_applied_pauli_mat = applied_x.add_mat(applied_y)
            .add_mat(applied_z);

        // encode probabilities into the matrices
        Matrix!(Complex!real) matrix_one = this.qc.density_mat.mult_scalar(
            Complex!real(1 - probability, 0));
        Matrix!(Complex!real) matrix_two = full_applied_pauli_mat.mult_scalar(
            Complex!real(probability / 3, 0));

        // calculate the result of applying the depolarizing noise
        Matrix!(Complex!real) result = matrix_one.add_mat(matrix_two);

        this.qc.density_mat = result;
    }

    /**
    * Apply depolarizing noise to multiple qubits acted on by a multi-qubit gate. The qubits should be the ones that the gate
    * acted on. Applying it to irrelevant qubits will lead to inaccurate results. Models independent noise on multiple qubits.
    *
    * params:
    * qubit_idxs = The qubits acted on by the gate as indices
    *
    * probabilities = The probabilities that independent depolarizing noise is applied to the qubits
    */
    void depolarizing_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the depolarizing_noise function for single qubits");

        assert(qubit_idxs.length == probabilities.length,
            "The probabilities array should be the same length as the qubit indices");

        foreach (i, qubit_idx; qubit_idxs) {
            depolarizing_noise(qubit_idx, probabilities[i]);
        }
    }

    /**
    * Apply depolarizing noise to multiple qubits acted on by a multi-qubit gate. The qubits should be the ones that the gate
    * acted on. Applying it to irrelevant qubits will lead to inaccurate results. Models correlated noise on multiple qubits.
    *
    * params:
    * qubit_idxs = The qubits acted on by the gate as indices
    *
    * probability_map = The probabilities that correlated depolarizing noise is applied to the qubits
    *                   using the kraus operators generated from pauli combinations embedded into the 
    *                   full hilbert space
    */
    void depolarizing_noise(int[] qubit_idxs, float[string] probability_map) {
        assert(isClose(sum(probability_map.values), 1.0f),
            "The sum of probabilities is not equal to one");

        foreach (key; probability_map.keys) {
            assert(key == toUpper(key), "The key values are not capitalized, they must be");
        }

        float expected_val = 1.0f / probability_map.length;
        foreach (prob; probability_map.values) {
            assert(isClose(prob, expected_val),
                "The probability distribution over pauli terms is not uniform, it should be for depolarizing noise");
        }

        Matrix!(Complex!real)[] krauss = build_full_krauss(qubit_idxs, probability_map);
        Matrix!(Complex!real) result = zeros(pow(2, this.qc.num_qubits));

        foreach (kraus; krauss) {
            result = result.add_mat(kraus.mult_mat(this.qc.density_mat)
                    .mult_mat(kraus.dagger()));
        }

        this.qc.density_mat = result;
    }

    /**
    * Apply pauli noise to a qubit acted on by a gate
    *
    * params:
    * qubit_idx = The index of the qubit acted on by the gate
    *
    * probability_map = An associatrive array of operator name to probability of each 
    *                   pauli operator being applied to the qubit. Operator names should 
    *                   be capitalized and probabilities should add up to 1.
    */
    void pauli_noise(int qubit_idx, float[string] probability_map) {
        assert(probability_map.length == 4,
            "The length of the probability array must be 4");

        assert(isClose(sum(probability_map.values), 1.0f),
            "The sum of probabilities is not equal to one");

        foreach (key; probability_map.keys) {
            assert(key == toUpper(key), "The key values are not capitalized, they must be");
        }

        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

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

        // encode the probabilities into the pauli matrices
        identity = identity.mult_scalar(Complex!real(sqrt(probability_map["I"]), 0));
        pauli_x = pauli_x.mult_scalar(Complex!real(sqrt(probability_map["X"]), 0));
        pauli_y = pauli_y.mult_scalar(Complex!real(sqrt(probability_map["Y"]), 0));
        pauli_z = pauli_z.mult_scalar(Complex!real(sqrt(probability_map["Z"]), 0));

        // embed the pauli matrices with probabilities encoded into the full
        // hilbert space. This is now the full kraus operator for the noise.
        identity = this.qc.build_full_gate(identity, qubit_idx);
        pauli_x = this.qc.build_full_gate(pauli_x, qubit_idx);
        pauli_y = this.qc.build_full_gate(pauli_y, qubit_idx);
        pauli_z = this.qc.build_full_gate(pauli_z, qubit_idx);

        // apply the kraus operator representing the noise channel to the density matrix
        Matrix!(Complex!real) applied_identity = identity.mult_mat(this.qc.density_mat)
            .mult_mat(identity.dagger());
        Matrix!(Complex!real) applied_x = pauli_x.mult_mat(this.qc.density_mat)
            .mult_mat(pauli_x.dagger());
        Matrix!(Complex!real) applied_y = pauli_y.mult_mat(this.qc.density_mat)
            .mult_mat(pauli_y.dagger());
        Matrix!(Complex!real) applied_z = pauli_z.mult_mat(this.qc.density_mat)
            .mult_mat(pauli_z.dagger());

        // replace the density matrix with the new one with noise
        // probablistically applied to it 
        this.qc.density_mat = applied_identity.add_mat(applied_x)
            .add_mat(applied_y).add_mat(applied_z);
    }

    /**
    * Apply pauli noise to qubits acted on by a multi-qubit gate. Models indpendent pauli noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits acted on by the gate
    *
    * probability_maps = The probabilities of each pauli operator being applied to the qubits. 
    *                    Each qubit has it's own set of probabilities for each operator.
    *                    Should add up to 1
    */
    void pauli_noise(int[] qubit_idxs, float[string][] probability_maps) {
        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the pauli_noise function for single qubits");

        assert(probability_maps.length == qubit_idxs.length,
            "The length of the probability map for combinations of pauli operator combinations should be the same as 4 to the number of qubits affected");

        foreach (i, qubit_idx; qubit_idxs) {
            pauli_noise(qubit_idx, probability_maps[i]);
        }
    }

    /**
    * Apply pauli noise to qubits acted on by a multi-qubit gate. Models correlated pauli noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits acted on by the gate
    *
    * probability_map = The probabilities of each pauli operator combination being applied to the qubits
    *                   as an associative array. Should add up to 1
    */
    void pauli_noise(int[] qubit_idxs, float[string] probability_map) {
        assert(isClose(sum(probability_map.values), 1.0f),
            "The sum of probabilities is not equal to one");

        foreach (key; probability_map.keys) {
            assert(key == toUpper(key), "The key values are not capitalized, they must be");
        }

        Matrix!(Complex!real)[] krauss = build_full_krauss(qubit_idxs, probability_map);
        Matrix!(Complex!real) result = zeros(pow(2, this.qc.num_qubits));

        foreach (kraus; krauss) {
            result = result.add_mat(kraus.mult_mat(this.qc.density_mat)
                    .mult_mat(kraus.dagger()));
        }

        this.qc.density_mat = result;
    }

    /**
    * Apply bit flip noise to a qubit acted on by a gate with some probability.
    *
    * params:
    * qubit_idx = The index of the qubit acted on by the gate to apply noise to
    * 
    * probability = The probability that the qubit is flipped
    */
    void bit_flip_noise(int qubit_idx, float probability) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

        Matrix!(Complex!real) pauli_x = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ])
            ]);

        // embed the operators into the full hilbert space
        identity = this.qc.build_full_gate(identity, qubit_idx);
        pauli_x = this.qc.build_full_gate(pauli_x, qubit_idx);

        // encode the probabilites into the matrices
        identity = identity.mult_scalar(Complex!real(sqrt(1 - probability), 0));
        pauli_x = pauli_x.mult_scalar(Complex!real(sqrt(probability), 0));

        // apply the matrices with encoded probabilities to the density matrix
        identity = identity.mult_mat(this.qc.density_mat).mult_mat(identity.dagger());
        pauli_x = pauli_x.mult_mat(this.qc.density_mat).mult_mat(pauli_x.dagger());

        // assign the new density matrix to the addition of the applied matrices
        this.qc.density_mat = identity.add_mat(pauli_x);
    }

    /**
    * Apply bit flip noise qubits acted on by a multi-qubit gate with some probabilities.
    * Models independent bit flip noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits acted on by the gate to apply noise to
    * 
    * probabilities = The probabilities that the qubits are flipped. These do not need to sum to 1
    */
    void bit_flip_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the bit_flip_noise function for single qubits");

        assert(probabilities.length == qubit_idxs.length,
            "The length of the qubit indices array and the probabilities array should be equal, it is not");

        foreach (i, qubit_idx; qubit_idxs) {
            bit_flip_noise(qubit_idx, probabilities[i]);
        }
    }

    /**
    * Apply bit flip noise qubits acted on by a multi-qubit gate with some probabilities. 
    * Models correlated bit flip noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits acted on by the gate to apply noise to
    * 
    * probability = The probability that the qubits are flipped.
    */
    void bit_flip_noise(int[] qubit_idxs, float probability) {
        Matrix!(Complex!real) kraus_identity = Matrix!(Complex!real)(pow(2, this.qc.num_qubits), pow(2, this
                .qc.num_qubits), []).identity(pow(2, this.qc.num_qubits));

        Matrix!(Complex!real) pauli_x = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(1, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ])
            ]);

        kraus_identity = kraus_identity.mult_scalar(Complex!real(sqrt(1 - probability), 0));
        Matrix!(Complex!real) bit_flip_kraus = build_full_pauli_op(pauli_x, qubit_idxs, probability);
        Matrix!(Complex!real) result = zeros(pow(2, this.qc.num_qubits));

        result = result.add_mat(kraus_identity.mult_mat(this.qc.density_mat)
                .mult_mat(kraus_identity.dagger()));

        result = result.add_mat(bit_flip_kraus.mult_mat(this.qc.density_mat)
                .mult_mat(bit_flip_kraus.dagger()));

        this.qc.density_mat = result;
    }

    /**
    * Applies phase flip noise to a qubit with some probability 
    *
    * params:
    * qubit_idx = The index of the qubit to apply the noise to
    *
    * probability = The probability that the phase of the qubit is flipped
    */
    void phase_flip_noise(int qubit_idx, float probability) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

        Matrix!(Complex!real) pauli_z = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(-1, 0)
                    ])
            ]);

        // embed the operators into the full hilbert space
        identity = this.qc.build_full_gate(identity, qubit_idx);
        pauli_z = this.qc.build_full_gate(pauli_z, qubit_idx);

        // encode the probabilites into the matrices
        identity = identity.mult_scalar(Complex!real(sqrt(1 - probability), 0));
        pauli_z = pauli_z.mult_scalar(Complex!real(sqrt(probability), 0));

        // apply the matrices with encoded probabilities to the density matrix
        identity = identity.mult_mat(this.qc.density_mat).mult_mat(identity.dagger());
        pauli_z = pauli_z.mult_mat(this.qc.density_mat).mult_mat(pauli_z.dagger());

        // assign the new density matrix to the addition of the applied matrices
        this.qc.density_mat = identity.add_mat(pauli_z);
    }

    /**
    * Applies phase flip noise to multiple qubits with some probabilities 
    * Models independent phase flip noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits to apply the noise to
    *
    * probabilities = The probabilities that the phase of the qubits is flipped. 
    *                 These do not need to sum to 1.
    */
    void phase_flip_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the bit_flip_noise function for single qubits");

        assert(probabilities.length == qubit_idxs.length,
            "The length of the qubit indices array and the probabilities array should be equal, it is not");

        foreach (i, qubit_idx; qubit_idxs) {
            phase_flip_noise(qubit_idx, probabilities[i]);
        }
    }

    /**
    * Applies phase flip noise to multiple qubits with some probabilities.
    * Models correlated phase flip noise. 
    *
    * params:
    * qubit_idxs = The indices of the qubits to apply the noise to
    *
    * probability = The probability that the phase of the qubits is flipped. 
    *                 These do not need to sum to 1.
    */
    void phase_flip_noise(int[] qubit_idxs, float probability) {
        int size = pow(2, this.qc.num_qubits);
        Matrix!(Complex!real) kraus_identity = Matrix!(Complex!real)(size, size, [
            ]).identity(size);

        Matrix!(Complex!real) pauli_z = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(-1, 0)
                    ])
            ]);

        kraus_identity = kraus_identity.mult_scalar(Complex!real(sqrt(1 - probability), 0));
        Matrix!(Complex!real) phase_flip_kraus = build_full_pauli_op(pauli_z, qubit_idxs, probability);
        Matrix!(Complex!real) result = zeros(pow(2, this.qc.num_qubits));

        result = result.add_mat(kraus_identity.mult_mat(this.qc.density_mat)
                .mult_mat(kraus_identity.dagger()));

        result = result.add_mat(phase_flip_kraus.mult_mat(this.qc.density_mat)
                .mult_mat(phase_flip_kraus.dagger()));

        this.qc.density_mat = result;
    }

    /**
    * Applies bit-phase flip noise to a qubit with some probability
    *
    * params:
    * qubit_idx = The index of the qubit to apply the noise to
    *
    * probability = The probability that the qubit state and phase are flipped
    */
    void bit_phase_flip_noise(int qubit_idx, float probability) {
        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, [
            ]).identity(2);
        Matrix!(Complex!real) pauli_y = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, -1)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 1), Complex!real(0, 0)
                    ])
            ]); // embed the operators into the full hilbert space
        identity = this.qc.build_full_gate(identity, qubit_idx);
        pauli_y = this.qc.build_full_gate(
            pauli_y, qubit_idx); // encode the probabilites into the matrices
        identity = identity.mult_scalar(Complex!real(sqrt(1 - probability), 0));
        pauli_y = pauli_y.mult_scalar(Complex!real(sqrt(probability), 0));

        // apply the matrices with encoded probabilities to the density matrix
        identity = identity.mult_mat(this.qc.density_mat)
            .mult_mat(identity.dagger());
        pauli_y = pauli_y.mult_mat(
            this.qc.density_mat).mult_mat(pauli_y.dagger()); // assign the new density matrix to the addition of the applied matrices
        this.qc.density_mat = identity.add_mat(pauli_y);
    }

    /**
    * Applies bit phase flip noise to multiple qubits with some probabilities.
    * Models independent bit phase flip noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits to apply the noise to
    *
    * probabilities = The probabilities that the qubits states and phases are flipped.
    *                 These do not need to sum to 1.
    */
    void bit_phase_flip_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the bit_flip_noise function for single qubits");

        assert(probabilities.length == qubit_idxs.length,
            "The length of the qubit indices array and the probabilities array should be equal, it is not");

        foreach (i, qubit_idx; qubit_idxs) {
            bit_phase_flip_noise(qubit_idx, probabilities[i]);
        }
    }

    /**
    * Applies bit phase flip noise to multiple qubits with some probabilities.
    * Models correlated bit phase flip noise.
    *
    * params:
    * qubit_idxs = The indices of the qubits to apply the noise to
    *
    * probability = The probability that the qubits states and phases are flipped.
    */
    void bit_phase_flip_noise(int[] qubit_idxs, float probability) {
        int size = pow(2, this.qc.num_qubits);
        Matrix!(Complex!real) kraus_identity = Matrix!(Complex!real)(size, size, [
            ]).identity(size);
        Matrix!(Complex!real) pauli_y = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(0, -1)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 1), Complex!real(0, 0)
                    ])
            ]);

        kraus_identity = kraus_identity.mult_scalar(Complex!real(sqrt(1 - probability), 0));
        Matrix!(Complex!real) phase_flip_kraus = build_full_pauli_op(pauli_y, qubit_idxs, probability);
        Matrix!(Complex!real) result = zeros(pow(2, this.qc.num_qubits));

        result = result.add_mat(kraus_identity.mult_mat(this.qc.density_mat)
                .mult_mat(kraus_identity.dagger()));

        result = result.add_mat(phase_flip_kraus.mult_mat(this.qc.density_mat)
                .mult_mat(phase_flip_kraus.dagger()));

        this.qc.density_mat = result;
    }
}
