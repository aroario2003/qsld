module examples.pauli_y;

import std.stdio;

import quantum.qc;

void main() {
    // The pauli-y gate is a lot like pauli-x in that it flips the qubit state from either |1> -> |0> or |0> -> |1>.
    // However, it differs in that depending on the initial state of the qubit it applies a phase shift of + or - i.
    // If the state is |0> than it applies a phase shift of +i and if it is |1> it applies a phase shift of -i.

    // Example 1:

    // Initialize a quantum circuit with 1 qubit
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply a pauli-y gate to the only qubit
    qc.pauli_y(0);

    writeln("Example 1: The state vector after applying a pauli-y gate to a one qubit system: ", qc
            .state.elems);

    //----------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 2 qubits
    QuantumCircuit qc2 = QuantumCircuit(2);

    // Apply a pauli-y to the first qubit and not the second
    qc2.pauli_y(0);

    writeln("Example 2: The state vector after applying the pauli-y gate to the first qubit in a 2 qubit system: ", qc2
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with 2 qubits
    QuantumCircuit qc3 = QuantumCircuit(2);

    // Apply a pauli-y to the second qubit and not the first
    qc3.pauli_y(1);

    writeln("Example 3: The state vector after applying the pauli-y gate to the second qubit in a 2 qubit system: ", qc3
            .state.elems);

    //-----------------------------------------------------------------------------------------------------------------

    // Example 4:

    // Initialize a quantum circuit with 2 qubits
    QuantumCircuit qc4 = QuantumCircuit(2);

    // Apply a pauli-y to both qubits
    qc4.pauli_y(0);
    qc4.pauli_y(1);

    writeln("Example 4: The state vector after applying the pauli-y gate to both qubits in a 2 qubit system: ", qc4
            .state.elems);

    // Tip: As always try to increase the number of qubits and apply pauli-y to specific ones to see what happens. 
    //      Also try applying hadamard and than pauli-x and pauli-y and see what happens. Start your code below.

}
