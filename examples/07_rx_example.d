module examples.rx;

import std.stdio;
import std.math;

import quantum.pure_state.qc;

void main() {
    // The Rx gate rotates the state vector around the x-axis of the bloch sphere by some angle theta
    // in radians. If theta is PI then it rotates the qubit 180 degrees, essentially flipping
    // it like a pauli-x gate. If theta is PI/2 then it creates an equal superposition over the |0>
    // and |1> states but with a specific phase shift applied. If theta is 0 then the qubit does not 
    // change at all.

    // Example 1:

    // Initialize a quantum circuit with 1 qubit
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply a rotation of PI to the qubit
    qc.rx(0, PI);

    writeln(
        "Example 1: The state vector after applying a rotation of PI to a 1 qubit system with initial state |0>: ", qc
            .state.elems);

    //---------------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 1 qubit
    QuantumCircuit qc2 = QuantumCircuit(1);

    // Apply a rotation of PI/2 to the qubit
    qc2.rx(0, PI / 2);

    writeln(
        "Example 2: The state vector after applying a rotation of PI/2 to a 1 qubit system with initial state |0>: ", qc2
            .state.elems);

    //--------------------------------------------------------------------------------------------------------------------

    // Example 3:

    //Initialize a quantum circuit with 1 qubit 
    QuantumCircuit qc3 = QuantumCircuit(1);

    // Apply a rotation of 0 to the qubit
    qc3.rx(0, 0);

    writeln(
        "Example 3: The state vector after applying a rotation of 0 to a 1 qubit system with initial state |0>: ", qc3
            .state.elems);

    // Tip: Try applying other degrees of rotation and also apply a hadamard to the system and see what changes.
    //      Maybe experiment with comparing this gate to pauli-x? ;) Start your code below.
}
