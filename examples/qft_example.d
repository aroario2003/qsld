module examples.qft;

import std.stdio;

import algos.qft;

void main() {
    // Example 1:

    // initialize the QFT quantum circuit with 2 qubits and an
    // initial state of |00> which is default
    QFT qft = QFT(2);

    // apply the QFT to the initial state
    qft.qft();

    // write the state vector to stdout for demonstation of what happens
    writeln("Example 1: state vector after QFT: ", qft.qc.state.elems);

    // invert the QFT (will not be super percise numerically, this is expected)
    qft.qft_inverse();

    // write the state vector to stdout again for demonstation of what happens
    writeln("Example 1: state vector after inverse QFT: ", qft.qc.state.elems);

    //-----------------------------------------------------------------------------

    //Example 2:

    // Initialize the QFT quantum circuit with 2 qubits and an initial state
    // of |01>. In general the first argument is the number of qubits and 
    // the second number is the basis state index which goes up to (2^n)-1 where
    // n is the number of qubits.
    QFT qft2 = QFT(2, 1);

    // apply the QFT to the initial state
    qft2.qft();

    // write the state vector to stdout as before
    writeln("Example 2: state vector after QFT: ", qft2.qc.state.elems);

    //invert the QFT 
    qft2.qft_inverse();

    // write the state vector to stdout as before
    writeln("Example 2: state vector after inverse QFT: ", qft2.qc.state.elems);
}
