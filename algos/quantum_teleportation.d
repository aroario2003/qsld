module algos.quantum_teleportation;

import std.stdio;
import std.format;
import std.conv;

import quantum.pure_state.qc;

struct QuantumTeleportation {
    void quantum_teleportation(void function(QuantumCircuit* qc) randomize_q0) {
        // Initialize a quantum circuit with 3 qubits to represent 
        // Alice and Bob's EPR (Einstein Podolsky Rosen) pair.
        QuantumCircuit qc = QuantumCircuit(3);

        // Randomize the state of qubit 0 before the algorithm starts
        randomize_q0(&qc);

        // Alice puts qubit 1 in superposition
        qc.hadamard(1);

        // Alice entangles qubits 1 and 2
        qc.cnot(1, 2);

        // Entangle qubits 0 and 1 to prepare for telepotation
        qc.cnot(0, 1);

        // Put qubit 0 into superposition
        qc.hadamard(0);

        // Measure qubits 0 and 1 to get what Bob should do with 
        // his qubit
        string q0_measured = qc.measure(0);
        string q1_measured = qc.measure(1);

        int q0_state = to!int(q0_measured);
        int q1_state = to!int(q1_measured);

        // Check which combination of values Bob measured in order
        // to determine what Bob should do with his qubit.
        if (q0_state == 0 && q1_state == 1) {
            qc.pauli_x(2);
        } else if (q0_state == 1 && q1_state == 0) {
            qc.pauli_z(2);
        } else if (q0_state == 1 && q1_state == 1) {
            qc.pauli_x(2);
            qc.pauli_z(2);
        }

        writeln(format("Bob measured: %d%d", q0_state, q1_state));

        // Declare bob's final qubit as a cicuit
        QuantumCircuit bob_qc = QuantumCircuit(1);
        // Give the measurement as a single number between 0 and 3
        // which can be rperesented as a binary number
        int alice_measurement = (q1_state << 1) | q0_state;

        // Loop over each amplitude index in the original state vector
        for (int i = 0; i < 8; i++) {
            // Mask the index to get Alices bits
            int alice_bits = i & 0b11;
            // Compare the bits to the measured value
            if (alice_bits == alice_measurement) {
                // Get the value of Bob's qubit
                int bob_bit = (i >> 2) & 0b1;
                // Put the amplitudes corresponding to Bob's bit into 
                // Bob's state vector for his qubit
                bob_qc.state.elems[bob_bit] = qc.state.elems[i];
            }
        }

        writeln("The full state vector is: ", qc.state.elems);
        writeln("Bob's final state vector is: ", bob_qc.state.elems);
    }
}
