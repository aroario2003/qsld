module viz.visualization;

import std.stdio;
import std.typecons;
import std.algorithm.searching;
import std.process;
import std.format;
import std.array;
import std.file;

import core.stdc.stdlib : exit;
import std.algorithm : canFind;

struct Visualization {
    Tuple!(string, int[], int)[] vis_arr;
    int num_qubits;
    int initial_state_idx;

    /**
    * The constructor for the type that allows for drawing the circuit
    * 
    * params:
    * vis_arr = The array of gates that the user calls as functions in the main program.
    *           This is generated internally by the quantum.qc module.
    * 
    * num_qubits = The number of qubits in the system 
    */
    this(Tuple!(string, int[], int)[] vis_arr, int num_qubits, int initial_state_idx) {
        this.vis_arr = vis_arr;
        this.num_qubits = num_qubits;
        this.initial_state_idx = initial_state_idx;
    }

    /**
    * Parses the entire vis_arr and writes the latex format to a file
    *
    * params:
    * filename = The name of the file to write the latex output to
    */
    void parse_and_write_vis_arr(string filename) {
        // write beginning of file boilerplate
        append(filename, "\\documentclass{standalone}\n");
        append(filename, "\\usepackage{quantikz}\n");
        append(filename, "\\begin{document}\n");
        append(filename, "\\scalebox{1.8} {%\n");
        append(filename, "\\begin{quantikz}\n");

        string[][] lines = [];
        for (int i = 0; i < this.num_qubits; i++) {
            int qubit_val = this.initial_state_idx & (1 << i);
            lines[lines.length++] = [
                format("\\lstick{\\ket{%d}} &", (qubit_val >> i))
            ];
        }

        foreach (i, item; this.vis_arr) {
            string gate_name = item[0];
            int[] qubit_idxs = item[1];
            int timestep = item[2];

            if (!gate_name.startsWith("C") && gate_name != "SWAP" && gate_name != "iSWAP" && gate_name != "TF") {
                if (gate_name != "M" && gate_name != "MA") {
                    lines[qubit_idxs[0]][lines[qubit_idxs[0]].length++] = format(" \\gate{%s} &", gate_name);
                } else {
                    if (gate_name == "M") {
                        lines[qubit_idxs[0]][lines[qubit_idxs[0]].length++] = " \\meter{} &";
                    } else if (gate_name == "MA") {
                        foreach (idx; qubit_idxs) {
                            lines[qubit_idxs[idx]][lines[qubit_idxs[idx]].length++] = " \\meter{} &";
                        }
                    }
                }
            } else {
                for (int j = 0; j < lines.length; j++) {
                    for (int k = 1; k < lines.length; k++) {
                        if (j == k) {
                            break;
                        }

                        if (lines[j].length - 1 < lines[k].length - 1) {
                            while (lines[j].length - 1 < lines[k].length - 1) {
                                lines[j][lines[j].length++] = " \\qw &";
                            }
                        } else if (lines[k].length - 1 < lines[j].length - 1) {
                            while (lines[k].length - 1 < lines[j].length - 1) {
                                lines[k][lines[k].length++] = " \\qw &";
                            }
                        }
                    }
                }
                switch (gate_name) {
                case "CX":
                    foreach (k; 0 .. this.num_qubits) {
                        if (!qubit_idxs.canFind(k)) {
                            lines[k][lines[k].length++] = " \\qw &";
                        }
                    }
                    lines[qubit_idxs[0]][lines[qubit_idxs[0]].length++] = format(" \\ctrl{%d} &", qubit_idxs[1] - qubit_idxs[0]);
                    lines[qubit_idxs[1]][lines[qubit_idxs[1]].length++] = " \\targ{} &";
                    break;
                case "TF":
                    int target_qubit = qubit_idxs[qubit_idxs.length - 1];

                    for (int j = 0; j < qubit_idxs.length - 1; j++) {
                        lines[qubit_idxs[j]][lines[qubit_idxs[j]].length++] = format(" \\ctrl{%d} &", target_qubit - qubit_idxs[j]);
                    }

                    foreach (k; 0 .. this.num_qubits) {
                        if (!qubit_idxs.canFind(k)) {
                            lines[k][lines[k].length++] = " \\qw &";
                        }
                    }

                    lines[target_qubit][lines[target_qubit].length++] = " \\targ{} &";
                    break;
                case "SWAP":
                    foreach (k; 0 .. this.num_qubits) {
                        if (!qubit_idxs.canFind(k)) {
                            lines[k][lines[k].length++] = " \\qw &";
                        }
                    }
                    lines[qubit_idxs[0]][lines[qubit_idxs[0]].length++] = format(" \\swap{%d} &", qubit_idxs[1] - qubit_idxs[0]);
                    lines[qubit_idxs[1]][lines[qubit_idxs[1]].length++] = " \\targX{} &";
                    break;
                default:
                    foreach (k; 0 .. this.num_qubits) {
                        if (!qubit_idxs.canFind(k)) {
                            lines[k][lines[k].length++] = " \\qw &";
                        }
                    }
                    lines[qubit_idxs[0]][lines[qubit_idxs[0]].length++] = format(" \\ctrl{%d} &", qubit_idxs[1] - qubit_idxs[0]);
                    lines[qubit_idxs[1]][lines[qubit_idxs[1]].length++] = format(" \\gate{%s} &", gate_name);
                }
            }
        }

        ulong max_len = 0;
        foreach (line; lines) {
            if (line.length > max_len)
                max_len = line.length;
        }

        foreach (line; lines) {
            while (line.length < max_len) {
                line[line.length++] = " \\qw &";
            }

            line[line.length++] = " \\qw \\\\\n";
            string full_line = line.join();
            append(filename, full_line);
        }

        //write end of file boilerplate
        append(filename, "\\end{quantikz}\n");
        append(filename, "}\n");
        append(filename, "\\end{document}\n");
    }

    /**
    * Compiles the latex file output by the parse_and_write_vis_arr() function
    *
    * params:
    * filename = The name of the file to compile ending with an extension .tex
    */
    void compile_tex_and_cleanup(string compiler, string filename) {
        auto tex_compilation_pid = spawnProcess([compiler, filename]);
        if (wait(tex_compilation_pid) != 0) {
            writeln("The compilation of the the latex file with name ", filename, " failed");
            exit(1);
        }

        string[] filename_split = filename.split(".");
        string filename_no_ext = filename_split[0];

        remove(filename);
        remove(format("./%s.aux", filename_no_ext));
        remove(format("./%s.log", filename_no_ext));
    }
}
