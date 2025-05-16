module examples.pauli_x;

import std.stdio;

import quantum.pure_state.qc;

void main() {
    // Example 1:

    // The pauli-x gate is the quantum equivalent of the classical not gate, that is 
    // if a qubit is in the state |1> it will change to |0> when the gate is applied
    // and vice versa.

    // Iniitalize a quantum circuit with 1 qubit
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply the pauli-x gate to the only qubit in the system
    qc.pauli_x(0);

    writeln("Example 1: The state vector after pauli-x on a one qubit system: ", qc.state.elems);

    //---------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 2 qubits
    QuantumCircuit qc2 = QuantumCircuit(2);

    // Apply a pauli-x gate to qubit at index 1
    qc2.pauli_x(1);

    writeln("Example 2: The state vector after applying pauli-x to one qubit on a 2 qubit system: ", qc2
            .state.elems);

    //---------------------------------------------------------------------------------------------------

    // Example 3:

    //Initialize a quantum circuit with 2 qubits
    QuantumCircuit qc3 = QuantumCircuit(2);

    // Apply pauli-x to both qubits to see the effect
    qc3.pauli_x(0);
    qc3.pauli_x(1);

    writeln("Example 3: The state vector after applying pauli-x to both qubits in a 2 qubit system: ", qc3
            .state.elems);

    // Tip: Try increasing the qubits and applying pauli-x to different combinations. 
    //      Maybe even apply a hadamard to some of them first. Start below this comment ;). 
}
