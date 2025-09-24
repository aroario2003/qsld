import std.complex;
import std.math;
import std.algorithm.iteration;
import std.array;

import qml.qnn;
import quantum.pure_state.qc;
import quantum.pure_state.observable;

DataElement[] train_data = [
    DataElement([-3.0], -1.0),
    DataElement([-0.7], -1.0),
    DataElement([-1.5], -1.0),
    DataElement([-7.0], -1.0),
    DataElement([0.2], 1.0),
    DataElement([1.4], 1.0),
    DataElement([3.0], 1.0),
    DataElement([5.0], 1.0),
    DataElement([2.0], 1.0),
    DataElement([-40.0], -1.0),
    DataElement([-10000.0], -1.0),
    DataElement([-0.3], -1.0),
    DataElement([-0.023], -1.0),
    DataElement([-0.45], -1.0),
    DataElement([-0.27], -1.0),
    DataElement([0.3], 1.0),
    DataElement([0.023], 1.0),
    DataElement([0.45], 1.0),
    DataElement([1000.0], 1.0),
    DataElement([0.27], 1.0),
    DataElement([0.0], 1.0),
    DataElement([1.0], 1.0),
    DataElement([-1.0], -1.0),
    DataElement([-0.4], -1.0),
    DataElement([-0.93], -1.0),
    DataElement([-0.54], -1.0),
    DataElement([-0.32], -1.0),
    DataElement([0.67], 1.0),
    DataElement([0.78], 1.0),
    DataElement([0.93], 1.0),
    DataElement([10000.0], 1.0),
    DataElement([2.32], 1.0),
    DataElement([3.56], 1.0),
    DataElement([-1000.0], -1.0),
    DataElement([-3.56], -1.0),
    DataElement([-5.75], -1.0),

];

real[] norm(real[] input) {
    return input.map!(v => atan(v) / PI + 0.5).array;
}

void encode(QuantumCircuit qc, real[] normalized_angles) {
    assert(qc.num_qubits == 1, "The number of qubits should be 1, it is not");
    foreach (i, angle; normalized_angles) {
        qc.rz(0, angle);
        qc.rx(0, angle);
        qc.ry(0, angle);
    }
}

real vqc(QuantumCircuit qc, real[] trainable_params) {
    Observable obs = Observable(["Z"], [Complex!real(1, 0)], 1);
    qc.rx(0, trainable_params[0]);
    qc.ry(0, trainable_params[1]);
    return qc.expectation_value(obs);
}

void main() {
    QnnConfig qnn_conf = QnnConfig(5, 1, [0.2, -0.2], 0.05, &norm, &encode, &vqc);
    Qnn qnn = Qnn(qnn_conf);
    qnn.train(train_data);
}
