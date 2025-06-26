import std.math;
import std.complex;
import std.array;
import std.algorithm.iteration;
import std.stdio;

import qml.qnn;
import quantum.pure_state.qc;
import quantum.pure_state.observable;

struct Mod3Qnn {
    Qnn qnn;

    static real[] norm(real[] input) {
        return input.map!(v => v * ((2 * PI) / 3)).array;
    }

    static void encode(QuantumCircuit qc, real[] normalized_angles) {
        real total_angle = ((normalized_angles[0] + normalized_angles[1] + normalized_angles[2]) % 3) * (
            2 * PI / 3);

        qc.rz(0, total_angle);
        qc.rz(0, total_angle);
        qc.ry(0, total_angle);

    }

    static real vqc(QuantumCircuit qc, real[] train_params) {
        Observable obs = Observable(["Z"], [Complex!real(1, 0)], 1);

        qc.rz(0, train_params[0]);
        qc.rz(0, train_params[1]);

        qc.rz(0, train_params[2]);
        qc.rz(0, train_params[3]);

        qc.rz(0, train_params[4]);

        return qc.expectation_value(obs);
    }

    void train(DataElement[] train_data) {
        QnnConfig qnn_conf = QnnConfig(15, 1, [
                0.2, -0.2, 0.2, -0.2, 0.2
            ], 0.05, &norm, &encode, &vqc);
        this.qnn = Qnn(qnn_conf);
        this.qnn.train(train_data);
    }

    real predict(DataElement elem) {
        return this.qnn.predict(elem);
    }
}

void main() {
    DataElement[] class_zero_train_data = [
        DataElement([0, 0, 0], 1.0),
        DataElement([0, 0, 2], 1.0),
        DataElement([0, 1, 0], -1.0),
        DataElement([1, 0, 1], -1.0),
        DataElement([2, 1, 0], 1.0),
        DataElement([1, 1, 0], -1.0),
        DataElement([1, 1, 1], 1.0),
        DataElement([2, 0, 0], -1.0),
        DataElement([1, 1, 0], 1.0),
        DataElement([0, 2, 0], -1.0),
        DataElement([1, 0, 0], -1.0),
        DataElement([2, 1, 0], 1.0),
        DataElement([1, 2, 0], 1.0),
        DataElement([0, 1, 1], -1.0),
        DataElement([1, 0, 2], 1.0),
        DataElement([0, 0, 1], -1.0),
        DataElement([0, 0, 2], -1.0),
        DataElement([0, 1, 2], 1.0),
        DataElement([1, 2, 2], -1.0),
        DataElement([1, 2, 0], 1.0),
        DataElement([2, 0, 1], 1.0),
        DataElement([2, 2, 0], -1.0),
        DataElement([2, 2, 2], 1.0),
        DataElement([2, 1, 1], -1.0),
        DataElement([2, 2, 2], 1.0),
        DataElement([1, 2, 1], -1.0),
        DataElement([2, 0, 1], 1.0),
        DataElement([0, 1, 1], 1.0),
        DataElement([1, 1, 2], -1.0),
        DataElement([0, 2, 1], 1.0),
        DataElement([2, 0, 2], -1.0),
        DataElement([2, 2, 1], -1.0),
        DataElement([0, 0, 1], 1.0),
        DataElement([0, 2, 2], -1.0),
        DataElement([2, 1, 2], -1.0),
        DataElement([0, 2, 1], 1.0),
        DataElement([1, 0, 2], 1.0),
    ];

    DataElement[] class_one_train_data = [
        DataElement([0, 0, 1], 1.0),
        DataElement([0, 0, 0], -1.0),
        DataElement([1, 0, 0], 1.0),
        DataElement([2, 2, 0], 1.0),
        DataElement([1, 1, 1], -1.0),
        DataElement([1, 2, 1], 1.0),
        DataElement([1, 0, 1], -1.0),
        DataElement([2, 1, 1], 1.0),
        DataElement([1, 2, 2], -1.0),
        DataElement([1, 1, 2], 1.0),
        DataElement([2, 0, 1], -1.0),
        DataElement([2, 0, 2], 1.0),
        DataElement([2, 2, 2], -1.0),
        DataElement([0, 2, 2], 1.0),
        DataElement([0, 1, 2], -1.0),
        DataElement([2, 1, 2], -1.0),
    ];

    DataElement[] class_two_train_data = [
        DataElement([0, 0, 2], 1.0),
        DataElement([0, 0, 0], -1.0),
        DataElement([1, 1, 0], 1.0),
        DataElement([2, 2, 1], 1.0),
        DataElement([2, 2, 2], -1.0),
        DataElement([2, 1, 2], 1.0),
        DataElement([1, 1, 1], -1.0),
        DataElement([1, 2, 2], 1.0),
        DataElement([1, 1, 2], -1.0),
        DataElement([1, 2, 1], -1.0),
        DataElement([2, 2, 2], 1.0),
        DataElement([0, 2, 0], 1.0),
        DataElement([2, 0, 1], -1.0),
        DataElement([1, 0, 1], 1.0),
        DataElement([0, 1, 1], -1.0),
        DataElement([2, 1, 0], -1.0),
    ];

    Mod3Qnn qnn0 = Mod3Qnn();
    Mod3Qnn qnn1 = Mod3Qnn();
    Mod3Qnn qnn2 = Mod3Qnn();

    writeln("class zero qnn training:");
    qnn0.train(class_zero_train_data);
    writeln("class one qnn training:");
    qnn1.train(class_one_train_data);
    writeln("class two qnn training:");
    qnn2.train(class_two_train_data);

    DataElement elem = DataElement([1, 2, 0], 0.0);
    DataElement elem1 = DataElement([1, 2, 2], 0.0);

    real class_zero_pred = qnn0.predict(elem);
    real class_one_pred = qnn1.predict(elem);
    real class_two_pred = qnn2.predict(elem);

    writeln("Prediction: ", [class_zero_pred, class_one_pred, class_two_pred]);

    real class_zero_pred1 = qnn0.predict(elem1);
    real class_one_pred1 = qnn1.predict(elem1);
    real class_two_pred1 = qnn2.predict(elem1);

    writeln("Prediction: ", [class_zero_pred1, class_one_pred1, class_two_pred1]);
}
