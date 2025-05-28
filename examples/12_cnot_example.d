module examples.cnot;

import std.stdio;

import quantum.pure_state.qc;

void main() {
    // The controlled not gate flips the target qubit state if and only if the control qubit is in the state |1>

    // Example 1:

    // Initialize a quantum circuit with 2 qubits in the initial state |00>
    QuantumCircuit qc = QuantumCircuit(2);

    // Apply a cnot gate with qubit 0 as control and qubit 1 as the target
    qc.cnot(0, 1);

    writeln(
        "The state vector after applying a cnot gate to qubit 1 with qubit 0 as control with the system in the initial state |00>: ", qc
            .state.elems);

    //-----------------------------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 2 qubits in the initial state |10>
    QuantumCircuit qc2 = QuantumCircuit(2, 2);

    // Apply a cnot gate with qubit 0 as control and qubit 1 as target
    qc2.cnot(0, 1);

    writeln(
        "The state vector after applying a cnot gate to qubit 1 with qubit 0 as control with the system in the initial state |10>: ", qc2
            .state.elems);

    //------------------------------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with 2 qubits in the initial state |01>
    QuantumCircuit qc3 = QuantumCircuit(2, 1);

    // Apply a cnot gate to qubit 1 with qubit 0 as control
    qc3.cnot(0, 1);

    writeln(
        "The state vector after applying a cnot gate to qubit 1 with qubit 0 as control with the system in the initial state |01>: ", qc3
            .state.elems);

    // Tip: Try reversing the order of the qubit indices in the call to the function and see what happens. 
    //      Try to make a maximally entangled bell state, if you know what that is ;).
}
