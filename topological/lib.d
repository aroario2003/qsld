module topological.lib;

import std.stdio;
import std.format;

package bool has_consecutive_zeros(string bit_str) {
    char[] bits = bit_str[].dup;

    for (int i = 0; i < bits.length; i++) {
        if (i == bits.length - (cast(ulong) 1)) {
            break;
        }

        if (bits[0] == '0') {
            return true;
        } else if (bits[i] == '0' && bits[i + 1] == '0') {
            return true;
        }
    }

    return false;
}

package int[string] generate_basis_to_index_mapping(int num_anyons) {
    assert(num_anyons >= 2, "The number of anyons should be greater than 1");

    int basis_state_idx = 0;
    int[string] basis_to_index_map;

    for (int i = 0; i < (1 << (num_anyons - 1)); i++) {
        string bit_str = format("%0*b", num_anyons - 1, i);

        if (!has_consecutive_zeros(bit_str)) {
            basis_to_index_map[bit_str] = basis_state_idx;
            basis_state_idx++;
        }
    }

    return basis_to_index_map;
}

package string[] generate_index_to_basis_mapping(int num_anyons) {
    assert(num_anyons >= 2, "The number of anyons should be greater than 1");

    string[] index_to_basis_map;

    for (int i = 0; i < (1 << (num_anyons - 1)); i++) {
        string bit_str = format("%0*b", num_anyons - 1, i);

        if (!has_consecutive_zeros(bit_str)) {
            index_to_basis_map ~= bit_str;
        }
    }

    return index_to_basis_map;
}

// Compute the chain of F matrices which need to be applied for
// a given braid to be possible
package int[] compute_f_chain(int first_idx, int second_idx) {
    int[] f_chain;
    int larger_idx;

    if (first_idx < second_idx) {
        larger_idx = second_idx;
    } else {
        larger_idx = first_idx;
    }

    for (int k = larger_idx - 1; k >= first_idx; k--) {
        f_chain ~= k;
    }

    return f_chain;
}
