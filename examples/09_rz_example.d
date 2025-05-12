module examples.rz;

import std.stdio;
import std.math;

import quantum.qc;

void main() {
    // The Rz gate rotates the state vector around the z-axis of the bloch sphere and applies only a phase shift to 
    // the qubit with the angle of rotation theta in the exponential term. The gate will apply a different phase shift 
    // based on the qubits state, if the qubit is in the state |0> then the phase shift is e^-i(theta/2) and if the qubit 
    // is in the state |1> than the phase shift is e^i(theta/2)

    // Example 1:

    // Initialize a quantum circuit with 1 qubit in the intial state |0>
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply a rotation of PI to the qubit
    qc.rz(0, PI);

    writeln(
        "Example 1: The state vector after applying the Rz gate with theta = PI on a 1 qubit system with initial state |0>: ", qc
            .state.elems);

    //----------------------------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with one qubit in the initial state |1>
    QuantumCircuit qc2 = QuantumCircuit(1, 1);

    // Apply a rotation of PI to the qubit
    qc2.rz(0, PI);

    writeln(
        "Example 2: The state vector after applying the Rz gate with theta = PI on a 1 qubit system with initial state |1>: ", qc2
            .state.elems);

    //-----------------------------------------------------------------------------------------------------------------------------

    // Example 3:

    // Initialize a quantum circuit with one qubit in the initial state |0>
    QuantumCircuit qc3 = QuantumCircuit(1);

    // Apply a rotation of PI/2 to the qubit
    qc3.rz(0, PI / 2);

    writeln(
        "Example 3: The state vector after applying the Rz gate with theta = PI/2 on a 1 qubit system with initial state |0>: ", qc3
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------------------

    // Example 4:

    // Initialize a quantum circuit with one qubit in the initial state |1>
    QuantumCircuit qc4 = QuantumCircuit(1, 1);

    // Apply a rotation of PI/2 to the qubit
    qc4.rz(0, PI / 2);

    writeln(
        "Example 4: The state vector after applying the Rz gate with theta = PI/2 on a 1 qubit system with initial state |1>: ", qc4
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------------------

    // Example 5:

    // Initialize a quantum circuit with one qubit in the initial state |0>
    QuantumCircuit qc5 = QuantumCircuit(1);

    // Apply a rotation of 2*PI to the qubit
    qc5.rz(0, 2 * PI);

    writeln(
        "Example 5: The state vector after applying the Rz gate with theta = 2*PI on a 1 qubit system with initial state |0>: ", qc5
            .state.elems);

    //-------------------------------------------------------------------------------------------------------------------------------

    // Example 6:

    // Initialize a quantum circuit with one qubit in the initial state |1>
    QuantumCircuit qc6 = QuantumCircuit(1, 1);

    // Apply a rotation of 2*PI to the qubit
    qc6.rz(0, 2 * PI);

    writeln(
        "Example 6: The state vector after applying the Rz gate with theta = 2*PI on a 1 qubit system with initial state |1>: ", qc6
            .state.elems);

    // Tip: Try add Rx and Ry gates and see what happens. Also, try combining all the other gates you've learned until now into a 
    //      circuit and see what happens, maybe even use different orders to see the difference ;). Start your code below.
}
