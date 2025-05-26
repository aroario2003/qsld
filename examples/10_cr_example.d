module examples.cr;

import std.stdio;

import quantum.pure_state.qc;

void main() {
    // The CR_k gate or controlled rotation of order k gate, rotates the phase by e^2 * PI / 2^k
    // if and only if the control and target qubits are 1. K determines how much the gate will 
    // rotate the phase

    // Example 1:

    // Initialize a quauntum circuit with 2 qubits in the initial state |00>
    QuantumCircuit qc = QuantumCircuit(2);

    // Apply the CR gate with a k of 2
    qc.cr(0, 1, 2);

    writeln("The state vector after applying a CR gate to a 2 qubit system in the state |00> with a k of 2: ", qc
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 2 qubits in the initial state |01>
    QuantumCircuit qc2 = QuantumCircuit(2, 1);

    // Apply the CR gate with a k of 2
    qc2.cr(0, 1, 2);

    writeln("The state vector after applying a CR gate to a 2 qubit system in the state |01> with a k of 2: ", qc2
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------

    // Example 3: 

    // Initialize a quantum circuit with 2 qubits in the initial state |11>
    QuantumCircuit qc3 = QuantumCircuit(2, 3);

    // Apply a CR gate with a k of 2
    qc3.cr(0, 1, 2);

    writeln("The state vector after applying a CR gate to a 2 qubit system in the state |11> with a k of 2: ", qc3
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------

    // Example 4:

    // Initialize a quantum circuit with 2 qubits in the intial state |11>
    QuantumCircuit qc4 = QuantumCircuit(2, 3);

    // Apply a CR gate with a k of 3
    qc4.cr(0, 1, 3);

    writeln("The state vector after applying a CR gate to a 2 qubit system in the state |11> with a k of 3: ", qc4
            .state.elems);

    // Tip: Try increasing k in different basis states and see how the phase rotation changes ;).
}
