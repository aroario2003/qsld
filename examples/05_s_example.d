module examples.s;

import std.stdio;

import quantum.pure_state.qc;

void main() {
    // The S gate or the PI/4 phase shift gate, adds a phase of PI/4 to the complex amplitude
    // of a basis state if the state of the qubit is |1>, otherwise nothing happens. Generally,
    // a phase shift can be thought of as a rotation in radians around a unit circle.

    // Example 1:

    // Initialize a quantum ciruit with 1 qubit and in the intial state |0>
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply the S gate to the only qubit in the system
    qc.s(0);

    writeln(
        "Example 1: The state vector after applying the S gate to a one qubit system with an initial state |0>: ", qc
            .state.elems);

    //--------------------------------------------------------------------------------------------------------------------

    //Example 2:

    //Initialize a quantum circuit with 1 qubit in the initial state |1>
    QuantumCircuit qc2 = QuantumCircuit(1, 1);

    // Apply the S gate to the qubit in state |1>
    qc2.s(0);

    writeln(
        "Example 2: The state vector after applying the S gate to a one qubit system with an initial state of |1>: ", qc2
            .state.elems);

    //--------------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with 2 qubits in the initial state |10>
    QuantumCircuit qc3 = QuantumCircuit(2, 2);

    // Apply the S gate to the first qubit
    qc3.s(0);

    writeln(
        "Example 3: The state vector after applying the S gate to a two qubit system on the first qubit with an initial state of |10>: ", qc3
            .state.elems);

    //---------------------------------------------------------------------------------------------------------------------

    // Example 4:

    //Initialize a quantum circuit with 2 qubits in the intial state |01>
    QuantumCircuit qc4 = QuantumCircuit(2, 1);

    // Apply the S gate to the first qubit
    qc4.s(0);

    writeln(
        "Example 4: The state vector after applying the S gate to a two qubit system on the first qubit with an initial state |01>: ", qc4
            .state.elems);

    // Tip: Try applying the S gate to both qubits in a 2 qubit system and try adding some hadamard gates as well.
    //      Maybe even try throwing in some of the pauli x, y and z gates ;). Start your code below.
}
