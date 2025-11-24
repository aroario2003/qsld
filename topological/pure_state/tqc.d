module topological.pure_state.tqc;

import std.complex;
import std.math;
import std.random;

import std.algorithm.mutation : reverse;
import std.algorithm : canFind;
import std.typecons : Tuple, tuple;

import linalg.vector;
import topological.lib;

struct TQC {
    // Numbers which determine the size of the hilbert space
    int num_anyons;
    int num_basis_states;

    // The mapping from valid fusion channels to their indices
    int[string] basis_to_index_map;
    // The inverse of the previous mapping, used for measurement
    string[] index_to_basis_map;

    // The state vector
    Vector!(Complex!real) state;

    /**
     * The constructor for a topological quantum circuit using fibbonacci anyons.
     *
     * params:
     * num_anyons = The number of anyons in the circuit.
     */
    this(int num_anyons) {
        this.num_anyons = num_anyons;
        this.basis_to_index_map = generate_basis_to_index_mapping(this.num_anyons);
        this.index_to_basis_map = generate_index_to_basis_mapping(this.num_anyons);
        this.num_basis_states = cast(int) basis_to_index_map.length;

        Complex!real[] state_arr = new Complex!real[this.num_basis_states];
        state_arr[] = Complex!real(0, 0);
        state_arr[0] = Complex!real(1, 0);

        this.state = Vector!(Complex!real)(this.num_basis_states, state_arr);
    }

    /**
     * The constructor for a topological quantum circuit using fibbonacci anyons.
     * But specifies a different amplitude in the state vector to initialize to 1+0i.
     *
     * params:
     * num_anyons = The number of anyons in the circuit.
     *
     * initial_basis_label_idx = The index of the amplitude in the state vector to initialze to 
     *                           1+0i
     */
    this(int num_anyons, int initial_basis_label_idx) {
        this.num_anyons = num_anyons;
        this.basis_to_index_map = generate_basis_to_index_mapping(this.num_anyons);
        this.index_to_basis_map = generate_index_to_basis_mapping(this.num_anyons);
        this.num_basis_states = cast(int) basis_to_index_map.length;

        Complex!real[] state_arr = new Complex!real[this.num_basis_states];
        state_arr[] = Complex!real(0, 0);
        state_arr[initial_basis_label_idx] = Complex!real(1, 0);

        this.state = Vector!(Complex!real)(this.num_basis_states, state_arr);
    }

    // Compute the valid index and partner index pairs to apply amplitude
    // changes to
    private Tuple!(int, int)[] compute_valid_f_idx_pairs(int bit_idx) {
        char[] partner_label;
        Tuple!(int, int)[] valid_indices;
        int partner_idx;

        foreach (basis_label, idx; this.basis_to_index_map) {
            partner_label = basis_label[].dup;
            partner_label[bit_idx] = '1';

            if (partner_label in this.basis_to_index_map) {
                partner_idx = this.basis_to_index_map[partner_label];

                if (idx < partner_idx) {
                    valid_indices ~= tuple(idx, partner_idx);
                }
            }
        }
        return valid_indices;
    }

    // Applies the F or change of basis matrix
    private void apply_f(int bit_idx) {
        Vector!(Complex!real) psi_prime = Vector!(Complex!real)(
            cast(int) this.state.elems.length, this.state.elems.dup);

        float phi = (1.0f + sqrt(5.0f)) / 2.0f;
        Tuple!(int, int)[] valid_index_pairs = compute_valid_f_idx_pairs(bit_idx);
        foreach (idx_pair; valid_index_pairs) {
            psi_prime[idx_pair[0]] = (1.0f / phi) * this.state[idx_pair[0]] + (
                1.0f / sqrt(phi)) * this.state[idx_pair[1]];

            psi_prime[idx_pair[1]] = (1.0f / sqrt(phi)) * this.state[idx_pair[0]] - (
                1.0f / phi) * this.state[idx_pair[1]];
        }

        this.state = psi_prime;
    }

    // Applies the R matrix, which applies phase to the 
    // quantum state
    private void apply_r(int smaller_anyon_idx, bool inverse) {
        if (!inverse) {
            foreach (basis_label, idx; this.basis_to_index_map) {
                if (basis_label[smaller_anyon_idx] == '0') {
                    this.state[idx] = expi((-4 * PI) / 5) * this.state[idx];
                } else if (basis_label[smaller_anyon_idx] == '1') {
                    this.state[idx] = expi((3 * PI) / 5) * this.state[idx];
                }
            }
        } else {
            foreach (basis_label, idx; this.basis_to_index_map) {
                if (basis_label[smaller_anyon_idx] == '0') {
                    this.state[idx] = conj(expi((-4 * PI) / 5)) * this.state[idx];
                } else if (basis_label[smaller_anyon_idx] == '1') {
                    this.state[idx] = conj(expi((3 * PI) / 5)) * this.state[idx];
                }
            }
        }
    }

    /**
     * Braid two adjacent or non-adjacent anyons
     * 
     * params: 
     * first_anyon_idx = The index of one of the anyons to braid
     * 
     * second_anyon_idx = The index of the second anyon to braid
     *
     * inverse = Whether or not to invert the R matrix when braiding the anyons
    */
    void braid(int first_anyon_idx, int second_anyon_idx, bool inverse = false) {
        int[] f_chain = compute_f_chain(first_anyon_idx, second_anyon_idx);

        foreach (idx; f_chain) {
            apply_f(idx);
        }

        apply_r(f_chain[0], inverse);

        foreach (idx; f_chain.reverse) {
            apply_f(idx);
        }
    }

    /**
     * Braid anyons with a specific sequence
     *
     * params:
     * braid_seq = A 2d array of two indices representing the 
     *             first and second anyon indices
     */
    void braid_sequence(int[][] braid_seq) {
        foreach (anyon_idxs; braid_seq) {
            this.braid(anyon_idxs[0], anyon_idxs[1]);
        }
    }

    /**
    * Collapse the fusion basis state vector to a single fusion channel
    *
    * returns: The bit string corresponding to the resulting fusion channel
    */
    string measure() {
        auto rng = Random(unpredictableSeed);
        auto r = uniform(0.0f, 1.0f, rng);

        real prob_accum = 0.0;
        string result;

        for (int i = 0; i < this.state.length(); i++) {
            prob_accum += norm(this.state[i]);
            if (r < prob_accum) {
                result = this.index_to_basis_map[i];
                break;
            }
        }

        return result;
    }

    int[string] measure(int shots) {
        int[string] counts;
        for (int i = 0; i < shots; i++) {
            string result = this.measure();
            counts[result] += 1;
        }

        return counts;
    }
}
