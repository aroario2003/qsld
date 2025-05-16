module quantum.pure_state.decoherence;

import std.complex;
import std.math;
import std.stdio;

import std.typecons : Nullable;

import linalg.vector;
import linalg.matrix;

// This computes the probability of one of the kraus operators acting on a given 
// state vector
private real[] compute_probs(int qubit_idx, Matrix!(Complex!real) operator1, Matrix!(
        Complex!real) operator2, Vector!(Complex!real) psi) {

    real probability1 = 0;
    real probability2 = 0;

    for (int i = 0; i < psi.length(); i++) {
        int j = i ^ (1 << qubit_idx);

        if (i < j) {
            Complex!real amplitude1 = psi[i];
            Complex!real amplitude2 = psi[j];

            Vector!(Complex!real) amplitudes = Vector!(Complex!real)(2, [
                    amplitude1, amplitude2
                ]);

            Vector!(Complex!real) psi_prime1 = operator1.mult_vec(amplitudes);
            Vector!(Complex!real) psi_prime2 = operator2.mult_vec(amplitudes);

            probability1 += pow(psi_prime1.mag(), 2);
            probability2 += pow(psi_prime2.mag(), 2);
        }
    }
    return [probability1, probability2];
}

struct T1Decay {
    real t1 = 50;

    /**
    * Applies T1 decay or amplitude damping to the state vector.
    *
    * params:
    * qubit_idx = The index of the qubit which will be affected by the decay
    *
    * gate_duration = The duration of the gate being applied to the qubit
    * 
    * psi = The state vector for the decay to act on
    */
    Vector!(Complex!real) apply(int qubit_idx, int gate_duration, Vector!(Complex!real) psi) {
        real t = gate_duration;

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

        real[] probabilities = compute_probs(qubit_idx, e0, e1, psi);
        Matrix!(Complex!real) operator;
        if (probabilities[0] > probabilities[1]) {
            operator = e0;
        } else {
            operator = e1;
        }

        for (int i = 0; i < psi.length(); i++) {
            int j = i ^ (1 << qubit_idx);

            if (i < j) {
                Complex!real amplitude1 = psi[i];
                Complex!real amplitude2 = psi[j];

                Vector!(Complex!real) amplitudes = Vector!(Complex!real)(2, [
                        amplitude1, amplitude2
                    ]);

                Vector!(Complex!real) result = operator.mult_vec(amplitudes);
                psi[i] = result[0];
                psi[j] = result[1];
            }
        }

        real psi_mag = psi.mag();
        for (int i = 0; i < psi.length(); i++) {
            psi[i] = psi[i] / psi_mag;
        }

        return psi;
    }
}

struct T2Decay {
    real t2 = 70;

    /**
    * Applies T2 decay or dephasing to the state vector.
    *
    * params:
    * qubit_idx = The index of the qubit which will be affected by the decay
    *
    * gate_duration = The duration of the gate being applied to the qubit
    * 
    * psi = The state vector for the decay to act on
    */
    Vector!(Complex!real) apply(int qubit_idx, int gate_duration, Vector!(Complex!real) psi) {
        real t = gate_duration;
        real gamma = 1 - exp(-t / this.t2);

        Matrix!(Complex!real) e0 = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(1, 0),
                        Complex!real(0, 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0),
                        Complex!real(sqrt(1 - gamma), 0)
                    ])
            ]);

        Matrix!(Complex!real) e1 = Matrix!(Complex!real)(2, 2, [
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0),
                        Complex!real(sqrt(gamma), 0)
                    ]),
                Vector!(Complex!real)(2, [
                        Complex!real(0, 0),
                        Complex!real(0, 0)
                    ])
            ]);

        real[] probabilities = compute_probs(qubit_idx, e0, e1, psi);
        Matrix!(Complex!real) operator;
        if (probabilities[0] > probabilities[1]) {
            operator = e0;
        } else {
            operator = e1;
        }

        for (int i = 0; i < psi.length(); i++) {
            int j = i ^ (1 << qubit_idx);

            if (i < j) {
                Complex!real amplitude1 = psi[i];
                Complex!real amplitude2 = psi[j];

                Vector!(Complex!real) amplitudes = Vector!(Complex!real)(2, [
                        amplitude1, amplitude2
                    ]);

                Vector!(Complex!real) result = operator.mult_vec(amplitudes);
                psi[i] = result[0];
                psi[j] = result[1];
            }
        }

        real psi_mag = psi.mag();
        for (int i = 0; i < psi.length(); i++) {
            psi[i] = psi[i] / psi_mag;
        }

        return psi;
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
