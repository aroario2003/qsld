// NOTE: This file contains the code for the impure subsystem of QSLD. this means that 
// the modules within this folder (quantum/impure_state/) work more correctly but not
// neccessarily efficiently for mixed or impure quantum states. The lack of efficiency
// is due to the fact that in order to maintain correctness when decohering an entangled
// state, I must construct a density matrix and full multi-qubit gates. This leads to
// inefficiency and high memory usage. When you are using this module it is important
// to keep this in mind. It is recommended that you do not use more than 5 or 6 qubits 
// in a single quantum circuit because otherwise you might kill your computer due to 
// high memory usage. If you would like to have efficiency and correctness for pure states
// with the ability to use more qubits, you should use the pure subsystem (quantum/pure_state/).

module quantum.impure_state.decoherence;

import std.complex;
import std.math;

import std.typecons : Nullable;

import linalg.vector;
import linalg.matrix;

private Matrix!(Complex!real) build_full_kraus(int qubit_idx, int num_qubits, Matrix!(Complex!real) operator) {
    Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);
    Matrix!(Complex!real) full_kraus;

    for (int i = 0; i < num_qubits; i++) {
        if (i == qubit_idx) {
            if (i == 0) {
                full_kraus = operator;
            } else {
                full_kraus = full_kraus.kronecker(operator);
            }
        } else {
            if (i == 0) {
                full_kraus = identity;
            } else {
                full_kraus = full_kraus.kronecker(identity);
            }
        }
    }

    return full_kraus;
}

struct T1Decay {
    real t1 = 50;

    /**
    * Applies T1 decay or amplitude damping to the density matrix.
    *
    * params:
    * qubit_idx = The index of the qubit which will be affected by the decay
    *
    * gate_duration = The duration of the gate being applied to the qubit
    * 
    * rho = The density matrix for the decay to act on
    */
    Matrix!(Complex!real) apply(int qubit_idx, int num_qubits, int gate_duration, Matrix!(
            Complex!real) rho) {

        int t = gate_duration;
        Matrix!(Complex!real) e0 = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0),
                        Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0),
                        Complex!real(sqrt(exp(-t / this.t1)), 0)
                    ])
            ]);

        Matrix!(Complex!real) e1 = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0),
                        Complex!real(sqrt(1 - exp(-t / this.t1)), 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0),
                        Complex!real(0, 0)
                    ])
            ]);

        Matrix!(Complex!real) full_e0 = build_full_kraus(qubit_idx, num_qubits, e0);
        Matrix!(Complex!real) full_e1 = build_full_kraus(qubit_idx, num_qubits, e1);

        Matrix!(Complex!real) e0_applied = full_e0.mult_mat(rho).mult_mat(full_e0.dagger());
        Matrix!(Complex!real) e1_applied = full_e1.mult_mat(rho).mult_mat(full_e1.dagger());

        return e0_applied.add_mat(e1_applied);
    }
}

struct T2Decay {
    real t2 = 70;

    /**
    * Applies T2 decay or dephasing to the density matrix.
    *
    * params:
    * qubit_idx = The index of the qubit which will be affected by the decay
    *
    * gate_duration = The duration of the gate being applied to the qubit
    * 
    * rho = The density matrix for the decay to act on
    */
    Matrix!(Complex!real) apply(int qubit_idx, int num_qubits, int gate_duration, Matrix!(
            Complex!real) rho) {

        real t = gate_duration;
        real gamma = 1 - exp(-t / this.t2);
        Complex!real kraus_num_1 = sqrt((1 + gamma) / 2);
        Complex!real kraus_num_2 = sqrt((1 - gamma) / 2);

        Matrix!(Complex!real) identity = Matrix!(Complex!real)(2, 2, []).identity(2);

        Matrix!(Complex!real) pauli_z = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0), Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0), Complex!real(-1, 0)
                    ])
            ]);

        Matrix!(Complex!real) e0 = identity.mult_scalar(kraus_num_1);
        Matrix!(Complex!real) e1 = pauli_z.mult_scalar(kraus_num_2);

        Matrix!(Complex!real) full_e0 = build_full_kraus(qubit_idx, num_qubits, e0);
        Matrix!(Complex!real) full_e1 = build_full_kraus(qubit_idx, num_qubits, e1);

        Matrix!(Complex!real) e0_applied = full_e0.mult_mat(rho).mult_mat(full_e0.dagger());
        Matrix!(Complex!real) e1_applied = full_e1.mult_mat(rho).mult_mat(full_e1.dagger());

        return e0_applied.add_mat(e1_applied);
    }
}

struct DecoherenceConfig {
    Nullable!T1Decay t1;
    Nullable!T2Decay t2;
    string decoherence_mode;

    /**
    * The constuctor for the configuration struct for decoherence and decay
    * 
    * params: 
    * t1 = The T1Decay struct, this can be null if you do not want to apply it
    *
    * t2 = The T2Decay struct, this can be null if you do not want to apply it
    *
    * decoherence_mode = The way in which you would like the decay to be applied
    *                    This can be automatic, manual or none. Manual means you 
    *                    have to apply it yourself after each gate or at random points.
    */
    this(Nullable!T1Decay t1, Nullable!T2Decay t2, string decoherence_mode) {
        this.t1 = t1;
        this.t2 = t2;
        this.decoherence_mode = decoherence_mode;

        assert(this.decoherence_mode == "manual" || this.decoherence_mode == "automatic" || this.decoherence_mode == "none",
            "The decoherence mode must be either manual, automatic or none, nothing else is supported");
    }
}
