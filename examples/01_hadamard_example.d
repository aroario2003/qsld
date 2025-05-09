module examples.hadamard;

import std.stdio;

import quantum.qc;

void main() {
    // Example 1:

    // The hadamard gate is the most important gate in quantum computing, if applied to all qubits
    // in a system it will put them into equal superposition. This is neccesary because the idea
    // is to go from an equal probability distribution to a state where one state has a greater probability
    // than the others.

    // Initialize a quantum circuit with 1 qubit 
    QuantumCircuit qc = QuantumCircuit(1);

    // Apply the hadamard gate to the qubit at the 0th index, the
    // only qubit
    qc.hadamard(0);

    writeln("Example 1; The state vector after hadamard on 1 qubit: ", qc.state.elems);

    //----------------------------------------------------------------------------------------------------------

    // Example 2:

    // Initialize a quantum circuit with two qubits
    QuantumCircuit qc2 = QuantumCircuit(2);

    // Apply hadamard to only one qubit in the system to see its effect
    // on the superposition.
    qc2.hadamard(0);

    writeln("Example 2: The state vector after hadamard on 1 qubit in a 2 qubit system: ", qc2
            .state.elems);

    //-----------------------------------------------------------------------------------------------------------

    // Example 3:

    //Initialize a quantum circuit with 2 qubits
    QuantumCircuit qc3 = QuantumCircuit(2);

    // Apply hadamard to all qubits in the system to see the effect on the superposition
    qc3.hadamard(0);
    qc3.hadamard(1);

    writeln("Example 3: The state vector after apply hadamard to all qubits in a 2 qubit system: ", qc3
            .state.elems);

    // Tip: try increasing the number qubits and applying hadamards to different indices below to see what happens
}
