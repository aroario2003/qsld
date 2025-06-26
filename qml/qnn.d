// NOTE: This is a basic API for building Quantum Neural Networks. This is quite limited,
// this really only supports making binary classifiers, however, it is possible to make 
// a multi-class classifier with teh one vs rest approach with an average of ~82% accuracy.
// I will not be adding multi-class classifiers as they do not do well in general and mostly
// get a range of ~75%-105% loss and ~44% accuracy.

import std.math;
import std.stdio;

import std.typecons : Tuple;

import quantum.pure_state.qc;

// The structure used to configure a custom QNN model
struct QnnConfig {
    int epochs; /// the amount of iterations to train for
    int num_qubits; /// the number of qubits in the quantum circuit which should be equal to the number of features
    real[] trainable_params; // the parameters to be modified by the QNN and used for training
    real learning_rate; // the rate at which the QNN will learn the problem, affects the effect of gradient descent

    real[]function(real[]) norm_func; /// normalizes the input
    void function(QuantumCircuit, real[]) encode; /// encodes the normalized angle or angles into the quantum state
    real function(QuantumCircuit, real[]) vqc; /// the variational quantum circuit applied and used for training
}

// The structure used to represent a singular piece of training data
struct DataElement {
    real[] input;
    real label;
}

struct Qnn {
    QnnConfig conf;
    QuantumCircuit qc;

    /**
    * Constructor for the QNN
    *
    * params:
    * conf = The QnnConfig to specify parameters to use for the QNN
    */
    this(QnnConfig conf) {
        this.conf = conf;
    }

    // converts the normalized inputs into normalized angles
    private real[] convert_norms_to_angles(real[] norms) {
        real[] norm_angles;
        foreach (norm; norms) {
            norm_angles ~= PI * norm;
        }

        return norm_angles;
    }

    // Loss function for raw comparison of binary classification results
    private real mean_squared_error(real result, real label) {
        real difference = result - label;
        real difference_sq = difference ^^ 2;
        return difference_sq * 0.5;
    }

    // Used to compute the gradient of the QNN
    private real parameter_shift(ulong i, real[] trainable_params, real[] normalized) {
        // make two new quantum circuits so that it doesnt affect the original
        QuantumCircuit qc_minus = QuantumCircuit(this.conf.num_qubits);
        QuantumCircuit qc_plus = QuantumCircuit(this.conf.num_qubits);

        // make two copies of the trainable parameters so that
        // so that they can be modified without modifying the 
        // original
        real[] train_params_minus = trainable_params.dup;
        real[] train_params_plus = trainable_params.dup;

        // add and subtract PI/2 from a given angle in each copy
        train_params_minus[i] -= PI / 2;
        train_params_plus[i] += PI / 2;

        // encode the normalized angles into the new quantum circuits
        this.conf.encode(qc_minus, normalized);
        this.conf.encode(qc_plus, normalized);

        // get the result of running the variational quantum circuit on the shifted parameters
        real minus_result = this.conf.vqc(qc_minus, train_params_minus);
        real plus_result = this.conf.vqc(qc_plus, train_params_plus);

        // get the difference between the results 
        real result_diff = plus_result - minus_result;
        result_diff = result_diff * 0.5;
        return result_diff;
    }

    // normalizes, encodes and runs the inputs through the VQC
    private Tuple!(real[], real) forward(real[] input) {
        real[] normalized = this.conf.norm_func(input);
        normalized = convert_norms_to_angles(normalized);
        this.conf.encode(this.qc, normalized);
        real result = this.conf.vqc(this.qc, this.conf.trainable_params);

        return Tuple!(real[], real)(normalized, result);
    }

    /**
    * Train the QNN on the dataset provided as an array
    *
    * params:
    * train_data = The data to train the QNN on
    */
    void train(DataElement[] train_data) {
        foreach (epoch; 1 .. this.conf.epochs) {
            real correct = 0;
            foreach (sample; train_data) {
                this.qc = QuantumCircuit(this.conf.num_qubits);

                Tuple!(real[], real) data = forward(sample.input);
                real[] normalized = data[0];
                real result = data[1];

                real prediction = tanh(result);
                // TODO; Allow for use of binary cross-entropy loss
                real loss = mean_squared_error(prediction, sample.label);

                real pred = 0;
                if (prediction > 0) {
                    pred = 1.0;
                } else {
                    pred = -1.0;
                }

                if (pred == sample.label) {
                    correct += 1;
                }

                writefln("Epoch %s | Label: %.2f | Pred: %.2f | Loss: %.4f", epoch, sample.label, prediction, loss);

                real loss_grad = 2 * (prediction - sample.label) * (1 - prediction * prediction);

                foreach (i, param; this.conf.trainable_params) {
                    real raw_gradient = parameter_shift(i, this.conf.trainable_params, normalized);
                    real gradient = loss_grad * raw_gradient;
                    this.conf.trainable_params[i] = param - this.conf.learning_rate * gradient;
                }
            }

            real accuracy = correct / cast(real) train_data.length;
            writefln("Epoch %s | Accuracy: %.4f", epoch, accuracy);
        }
    }

    /**
    * Give a prediction for a specific piece of data
    *
    * params:
    * elem = The piece of data for which to generate the prediction
    *
    * returns: A real number representing whether or not the data is in a 
    * specific class or not. In general, the interpretation of the results
    * of this function will vary based on the problem.
    */
    real predict(DataElement elem) {
        Tuple!(real[], real) data = forward(elem.input);
        real result = data[1];
        real prediction = tanh(result);
        return prediction;
    }
}
