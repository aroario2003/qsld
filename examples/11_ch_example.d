module examples.ch;

import std.stdio;

import quantum.pure_state.qc;

void main() {
    // The controlled hadamard gate puts specific qubit into an equal superposition of the states
    // |0> and |1> if and only if the control qubit is 0.

    // Example 1:

    // Initialize a quantum circuit with 2 qubits in the initial state |00>
    QuantumCircuit qc = QuantumCircuit(2);

    // Apply a controlled hadamard gate with qubit 0 as the control and qubit 1 as the target
    qc.ch(0, 1);

    writeln("The state vector after applying a controlled hadamard gate to the system in an initial state of |00>: ", qc
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 2 qubits in the initial state |01>
    QuantumCircuit qc2 = QuantumCircuit(2, 1);

    // Apply a controlled hadamard gate with qubit 0 as control and qubit 1 as target
    qc2.ch(0, 1);

    writeln("The state vector after applying a controlled hadamard gate to the system in an initial state of |01>: ", qc2
            .state.elems);

    //--------------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with 2 qubits in the initial state |10>
    QuantumCircuit qc3 = QuantumCircuit(2, 2);

    // Apply a controlled hadamard gate with qubit 1 control and qubit 0 as target
    qc3.ch(1, 0);

    writeln(
        "The state vector after applying a controlled hadamard gate to the system in an initial state of |10> with qubit 1 as control: ", qc3
            .state.elems);

    // Tip; Try another basis state and see what happens, increase the number of qubits and try different basis states with different permutations
    // of qubit indices ;).
}
