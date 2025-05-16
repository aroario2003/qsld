module examples.ry;

import std.stdio;
import std.math;

import quantum.pure_state.qc;

void main() {
    // The Ry gate applies rotates the state vector by an angle theta in radians around the y-axis of the
    // bloch sphere. The main difference between the Rx gate and this one is that this one does not introduce 
    // any imaginary values into the amplitudes.

    // Example 1:

    // Initialize a quantum circuit with 1 qubit 
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply a rotation of PI to the qubit
    qc.ry(0, PI);

    writeln(
        "Example 1: The state vector after applying the Ry gate with theta = PI on a 1 qubit system with initial state |0>: ", qc
            .state.elems);

    //----------------------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with 1 qubit
    QuantumCircuit qc2 = QuantumCircuit(1);

    // Apply a rotation of PI/2 to the qubit
    qc2.ry(0, PI / 2);

    writeln(
        "Example 2: The state vector after applying the Ry gate with theta = PI/2 on a 1 qubit system with initial state |0>: ", qc2
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------------------

    // Example 3:

    //Initialize a quantum circuit with 1 qubit
    QuantumCircuit qc3 = QuantumCircuit(1);

    // Apply a rotation of 0 to the qubit
    qc3.ry(0, 0);

    writeln(
        "Example 3: The state vector after applying the Ry gate with theta = 0 on a 1 qubit system with initial state |0>: ", qc3
            .state.elems);

    // Tip: Try testing Ry with a different initial state and more qubits. Also, try combining with Rx and 
    //      hadamard gates. Maybe even add some of the pauli gates and phase shift gates to see what happens ;).
    //      State your code below.
}
