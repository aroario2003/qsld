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

module quantum.pure_state.gate_noise;

import std.random;
import std.math;
import std.conv;
import std.algorithm.iteration;
import std.uni;

import std.range : repeat;

import quantum.pure_state.qc;

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
        auto r = uniform(0.0, 1.0f, this.rng);

        float[string] gates = ["X": 1.0f / 3, "Y": 1.0f / 3, "Z": 1.0f / 3];

        if (r < probability) {
            float sum = 0;
            auto r2 = uniform(0.0, 1.0f, this.rng);

            foreach (gate, prob; gates) {
                sum += prob;
                if (r2 < sum) {
                    switch (gate) {
                    case "X":
                        this.qc.pauli_x(qubit_idx);
                        this.qc.tableau.error.elems[qubit_idx] ^= 1;
                        break;
                    case "Y":
                        this.qc.pauli_y(qubit_idx);
                        this.qc.tableau.error.elems[qubit_idx] ^= 1;
                        this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idx] ^= 1;
                        break;
                    case "Z":
                        this.qc.pauli_z(qubit_idx);
                        this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idx] ^= 1;
                        break;
                    default:
                        assert(false, "Invalid gate chosen for depolarizing noise");
                        break;
                    }
                    break;
                }
            }
        }
    }

    /**
    * Apply depolarizing noise to multiple qubits acted on by a multi-qubit gate. The qubits should be the ones that the gate
    * acted on. Applying it to irrelevant qubits will lead to inaccurate results.
    *
    * params:
    * qubit_idxs = The qubits acted on by the gate as indices
    *
    * probability = The probability that depolarizing noise is applied to the qubits
    */
    void depolarizing_noise(int[] qubit_idxs, float probability) {
        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the depolarizing_noise function for single qubits");

        string[] pauli_combos = generate_pauli_combos(cast(int) qubit_idxs.length);
        auto r = uniform(0.0, 1.0f, this.rng);

        if (r < probability) {
            float sum = 0;
            auto r2 = uniform(0.0, 1.0f, rng);
            float pauli_prob = 1.0f / pauli_combos.length;
            foreach (combo; pauli_combos) {
                sum += pauli_prob;
                if (r2 < sum) {
                    char[] pauli_ops = combo[].dup;
                    foreach (i, pauli; pauli_ops) {
                        switch (pauli) {
                        case 'I':
                            continue;
                        case 'X':
                            this.qc.pauli_x(qubit_idxs[i]);
                            this.qc.tableau.error.elems[qubit_idxs[i]] ^= 1;
                            break;
                        case 'Y':
                            this.qc.pauli_y(qubit_idxs[i]);
                            this.qc.tableau.error.elems[qubit_idxs[i]] ^= 1;
                            this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idxs[i]] ^= 1;
                            break;
                        case 'Z':
                            this.qc.pauli_z(qubit_idxs[i]);
                            this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idxs[i]] ^= 1;
                            break;
                        default:
                            assert(false, "Invalid gate encoutered while parsing pauli operator combinations");
                        }
                    }
                    break;
                }
            }
        }
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

        auto r = uniform(0.0, 1.0f, this.rng);

        float sum = 0;
        foreach (op, prob; probability_map) {
            sum += prob;
            if (r < sum) {
                switch (op) {
                case "I":
                    break;
                case "X":
                    this.qc.pauli_x(qubit_idx);
                    this.qc.tableau.error.elems[qubit_idx] ^= 1;
                    break;
                case "Y":
                    this.qc.pauli_y(qubit_idx);
                    this.qc.tableau.error.elems[qubit_idx] ^= 1;
                    this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idx] ^= 1;
                    break;
                case "Z":
                    this.qc.pauli_z(qubit_idx);
                    this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idx] ^= 1;
                    break;
                default:
                    assert(false, "Invalid gate for pauli error noise");
                }

                break;
            }
        }
    }

    /**
    * Apply pauli noise to qubits acted on by a multi-qubit gate
    *
    * params:
    * qubit_idxs = The indices of the qubits acted on by the gate
    *
    * probability_map = The probabilities of each pauli operator combination being applied to the qubits 
    *                   as an associative array, Should add up to 1
    */
    void pauli_noise(int[] qubit_idxs, float[string] probability_map) {
        assert(probability_map.length == pow(4, qubit_idxs.length),
            "The length of the probability map for combinations of pauli operator combinations should be the same as 4 to the number of qubits affected");

        assert(isClose(sum(probability_map.values), 1.0f),
            "The sum of the probabilities in the probability map of pauli operator combinations should sum to 1, it does not");

        assert(qubit_idxs.length > 1,
            "Please use the single qubit version of the pauli_noise function for single qubits");

        foreach (key; probability_map.keys) {
            assert(key == toUpper(key), "The key values are not capitalized, they must be");
        }

        auto r = uniform(0.0, 1.0f, this.rng);

        float sum = 0;
        foreach (combo, prob; probability_map) {
            sum += prob;
            if (r < sum) {
                char[] pauli_ops = combo[].dup;
                foreach (i, pauli; pauli_ops) {
                    switch (pauli) {
                    case 'I':
                        break;
                    case 'X':
                        this.qc.pauli_x(qubit_idxs[i]);
                        this.qc.tableau.error.elems[qubit_idxs[i]] ^= 1;
                        break;
                    case 'Y':
                        this.qc.pauli_y(qubit_idxs[i]);
                        this.qc.tableau.error.elems[qubit_idxs[i]] ^= 1;
                        this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idxs[i]] ^= 1;
                        break;
                    case 'Z':
                        this.qc.pauli_z(qubit_idxs[i]);
                        this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idxs[i]] ^= 1;
                        break;
                    default:
                        assert(false, "Invalid gate encoutered when parsing pauli operator combination");
                    }
                }
                break;
            }
        }
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
        auto r = uniform(0.0, 1.0f, this.rng);

        if (r < probability) {
            this.qc.pauli_x(qubit_idx);
            this.qc.tableau.error.elems[qubit_idx] ^= 1;
        }
    }

    /**
    * Apply bit flip noise qubits acted on by a multi-qubit gate with some probabilities.
    *
    * params:
    * qubit_idxs = The indices of the qubits acted on by the gate to apply noise to
    * 
    * probabilities = The probabilities that the qubits are flipped. These do not need to sum to 1
    */
    void bit_flip_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1 && probabilities.length > 1,
            "Please use the single qubit version of the bit_flip_noise function for single qubits");

        assert(qubit_idxs.length == probabilities.length,
            "The two arrays passed must have equal length");

        foreach (i, qubit_idx; qubit_idxs) {
            auto r = uniform(0.0, 1.0f, this.rng);
            if (r < probabilities[i]) {
                this.qc.pauli_x(qubit_idx);
                this.qc.tableau.error.elems[qubit_idxs[i]] ^= 1;
            }
        }
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
        auto r = uniform(0.0, 1.0f, this.rng);

        if (r < probability) {
            this.qc.pauli_z(qubit_idx);
            this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idx] ^= 1;
        }
    }

    /**
    * Applies phase flip noise to multiple qubits with some probability 
    *
    * params:
    * qubit_idxs = The indices of the qubits to apply the noise to
    *
    * probabilities = The probabilities that the phase of the qubits is flipped. 
    *                 These do not need to sum to 1.
    */
    void phase_flip_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1 && probabilities.length > 1,
            "Please use the single qubit version of the phase_flip_noise function for single qubits");

        assert(qubit_idxs.length == probabilities.length,
            "The two arrays passed must have equal length");

        foreach (i, qubit_idx; qubit_idxs) {
            auto r = uniform(0.0, 1.0f, this.rng);
            if (r < probabilities[i]) {
                this.qc.pauli_z(qubit_idx);
                this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idxs[i]] ^= 1;
            }
        }
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
        auto r = uniform(0.0, 1.0f, this.rng);

        if (r < probability) {
            this.qc.pauli_y(qubit_idx);
            this.qc.tableau.error.elems[qubit_idx] ^= 1;
            this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idx] ^= 1;
        }
    }

    /**
    * Applies bit-phase flip noise to multiple qubits with some probability
    *
    * params:
    * qubit_idxs = The indices of the qubits to apply the noise to
    *
    * probabilities = The probabilities that the qubits states and phases are flipped.
    *                 These do not need to sum to 1.
    */
    void bit_phase_flip_noise(int[] qubit_idxs, float[] probabilities) {
        assert(qubit_idxs.length > 1 && probabilities.length > 1,
            "Please use the single qubit version of the bit_phase_flip_noise function for single qubits");

        assert(qubit_idxs.length == probabilities.length,
            "The two arrays passed must have equal length");

        foreach (i, qubit_idx; qubit_idxs) {
            auto r = uniform(0.0, 1.0f, this.rng);
            if (r < probabilities[i]) {
                this.qc.pauli_y(qubit_idx);
                this.qc.tableau.error.elems[qubit_idxs[i]] ^= 1;
                this.qc.tableau.error.elems[this.qc.num_qubits + qubit_idxs[i]] ^= 1;
            }
        }
    }

    /**
    * Generate a noisy angle for unitary rotation gates, such that it causes
    * a coherent overrotation of the qubit.
    *
    * params:
    * theta = The ideal angle to apply the noise to
    *
    * epsilon = The amount of overrotation to apply to the angle theta
    *
    * returns: The new overrotated angle
    */
    real noisy_angle(real theta, real epsilon) {
        return theta + epsilon;
    }

    /**
    * Generate a noisy angle for multi-qubit unitary rotation gates, such that it
    * causes each qubit to be overrotated by some angle epsilon.
    *
    * params:
    * thetas = The array of ideal angles to rotate the qubits by
    *
    * epsilons = The angles by which to overrotate the angle theta
    *
    * returns: An array of the new overrotated angles for each qubit
    */
    real[] noisy_angle(real[] thetas, real[] epsilons) {
        real[] noisy_thetas;

        foreach (i, epsilon; epsilons) {
            noisy_thetas ~= thetas[i] + epsilon;
        }

        return noisy_thetas;
    }
}
