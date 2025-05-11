module examples.pauli_z;

import std.stdio;

import quantum.qc;

void main() {
    // The pauli-z gate is somewhat like the pauli-y gate from the previous example, 
    //the difference is that it does not flip the qubit, instead it only flips the phase 
    // of the probability amplitude if the qubit is in the state |1>

    // Example 1:

    // Initialize a quantum circuii with one qubit in the state |00> 
    QuantumCircuit qc = QuantumCircuit(1);

    // apply pauli-z to the only qubit in the system
    qc.pauli_z(0);

    writeln("Example 1: The state vector after applying pauli-z to a one qubit system in the initial state |0>: ", qc
            .state.elems);

    //----------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 1 qubit in the state |1>
    QuantumCircuit qc2 = QuantumCircuit(1, 1);

    // apply pauli-z once again
    qc2.pauli_z(0);

    writeln("Example 2: The state vector after applying pauli-z to a one qubit system in the initial state |1>: ", qc2
            .state.elems);

    //----------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with 2 qubit and an initial state |10>
    QuantumCircuit qc3 = QuantumCircuit(2, 2);

    // apply pauli-z to the first qubit, think about what will happen based on the initial state
    // of the system ;).
    qc3.pauli_z(0);

    writeln("Example 3: The state vector after applying pauli-z to a two qubit system in the initial state |10>: ", qc3
            .state.elems);
}
