module qec.lib;

import std.algorithm : canFind;

import linalg.vector;
import linalg.matrix;

import qec.stabilizer;

// Used to put the tableau matrix into upper echelon form
// and find independent stabilizers used to correct the state
package Matrix!int gaussian_elimination(Tableau tableau) {
    ulong row_len = tableau.tableau_internal.rows[0].length - (cast(ulong) 1);
    ulong r = 0;

    bool pivot_found = false;

    Vector!int pivot_row = Vector!int(2 * tableau.num_qubits, new int[2 * tableau.num_qubits]);
    Matrix!int independent_rows;
    ulong[] pivoted_cols;
    ulong[] pivoted_rows;

    for (ulong col_idx = 0; col_idx < row_len; col_idx++) {
        if (pivoted_cols.canFind(col_idx))
            continue;

        for (int i = tableau.num_qubits - 1; i < 2 * tableau.num_qubits; i++) {
            if (pivoted_rows.canFind(i))
                continue;

            Vector!int row = tableau.tableau_internal.rows[i];

            if (row[col_idx] == 1) {
                pivot_found = true;

                Vector!int temp = tableau.tableau_internal.rows[tableau.num_qubits + r];
                tableau.tableau_internal.rows[tableau.num_qubits + r] = tableau
                    .tableau_internal.rows[i];
                tableau.tableau_internal.rows[i] = temp;

                pivot_row = Vector!int((2 * tableau.num_qubits), tableau
                        .tableau_internal.rows[tableau.num_qubits + r].elems[0 .. row_len].dup);

                pivoted_rows ~= i;
                pivoted_cols ~= col_idx;
                break;
            }
        }

        if (pivot_found) {
            for (int i = tableau.num_qubits - 1; i < 2 * tableau.num_qubits; i++) {
                Vector!int row = tableau.tableau_internal.rows[i];
                if (i != (tableau.num_qubits + r)) {
                    if (row[col_idx] == 1) {
                        for (int x = 0; x < row_len; x++) {
                            tableau.tableau_internal.rows[i].elems[x] = row[x] ^ pivot_row[x];
                        }
                    }
                }
            }

            r++;
            independent_rows.append(pivot_row);
            pivot_found = false;

            if (r >= tableau.num_qubits) {
                break;
            }
        }
    }

    return independent_rows;
}
