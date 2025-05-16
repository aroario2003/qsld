module examples.t;

import std.stdio;

import quantum.pure_state.qc;

void main() {
    // The T phase shift gate or PI/8 gate, adds a phase of PI/8 to the complex probability amplitude
    // of some basis state if the qubit at some index i in that basis state is |1>.

    // Example 1:

    //Initialize a quantum circuit with 1 qubit with initial state |0>
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply T gate to that one qubit
    qc.t(0);

    writeln("Example 1: The state vector after applying the T gate to a one qubit system with initial state |0>: ", qc
            .state.elems);

    //-----------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 1 qubit with initial state |1>
    QuantumCircuit qc2 = QuantumCircuit(1, 1);

    // Apply a T gate to the one qubit
    qc2.t(0);

    writeln("Example 2: The state vector after applying the T gate to a one qubit system with initial state |1>: ", qc2
            .state.elems);

    //------------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with 2 qubits with initial state |01>
    QuantumCircuit qc3 = QuantumCircuit(2, 1);

    // Apply a T gate to the second qubit
    qc3.t(1);

    writeln("Example 3: The state vector after applying the T gate to a two qubit system with initial state |01>: ", qc3
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------

    // Example 4:

    //Initialize a quantum circuit with 2 qubits in the initial state |10>
    QuantumCircuit qc4 = QuantumCircuit(2, 2);

    // Apply a T gate to the second qubit
    qc4.t(1);

    writeln("Example 3: The state vector after applying the T gate to a two qubit system with initial state |10>: ", qc4
            .state.elems);

    // Tip: Try applying a T gate to both qubits and then try adding an S gate. Then maybe add some gates from previous examples.
    //      See what happens and what you'll learn. Start your code below ;).

}
