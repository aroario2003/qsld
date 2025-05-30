module algos.shors;

import std.complex;
import std.math;
import std.random;
import std.stdio;

import core.stdc.stdlib : exit;

import algos.qft;

import linalg.vector;

import quantum.pure_state.qc;

struct Shors {
    QuantumCircuit qc;
    int num_qubits;
    int qubits_first_register;
    int qubits_second_register;
    int a;
    int n;

    /**
    * The constructor for shors algorithm
    *
    * params:
    * a = A number less than n but should also be coprime to n
    *
    * n = The number to find the non-trivial factors of, this should not be 1 or prime
    */
    this(int a, int n) {
        this.qubits_first_register = cast(int) floor(2 * log2(cast(float) n));
        this.qubits_second_register = cast(int) floor(log2(cast(float) n));
        this.num_qubits = this.qubits_first_register + this.qubits_second_register;
        this.qc = QuantumCircuit(this.num_qubits);
        this.a = a;
        this.n = n;
    }

    // gets the greatest common divisor of 2 numbers
    private int gcd(int a, int b) {
        int result = 0;
        if (a < b) {
            result = a;
        } else {
            result = b;
        }

        while (result > 0) {
            if (a % result == 0 && b % result == 0) {
                break;
            }
            result--;
        }

        return result;
    }

    // puts the values calculated in register 1 into register 2
    private void oracle() {
        Vector!(Complex!real) psi_prime = Vector!(Complex!real)(pow(2, this.num_qubits), new Complex!real[pow(2, this
                    .num_qubits)]);

        for (int i = 0; i < pow(2, this.num_qubits); i++) {
            psi_prime[i] = Complex!real(0, 0);
        }

        for (int i = 0; i < this.qc.state.elems.length; i++) {
            if (this.qc.state.elems[i] != Complex!real(0, 0)) {
                int x = i >> this.qubits_second_register;
                int y = i & ((1 << this.qubits_second_register) - 1);
                int y_prime = (y * pow(this.a, x)) % this.n;
                int i_prime = (x << this.qubits_second_register) + y_prime;

                psi_prime[i_prime] = psi_prime[i_prime] + this.qc.state.elems[i];
            }
        }

        this.qc.state = psi_prime;
    }

    // measures the first qubit register
    private int measure_reg_1() {
        float[int] prob_dist;
        for (int i = 0; i < this.qc.state.elems.length; i++) {
            int x = i >> this.qubits_second_register;
            prob_dist[x] += norm(this.qc.state.elems[i]);
        }

        auto rng = Random(unpredictableSeed);
        auto r = uniform(0.0, 1.0f, rng);

        real sum = 0;
        int result = 0;
        foreach (k, v; prob_dist) {
            sum += v;
            if (r < sum) {
                result = k;
                break;
            }
        }
        return result;
    }

    // Finds the approximation of the period r of f(x) = a^x mod N
    private int find_r(int x, int q) {
        real[] frac_coeff_list;

        real r = cast(real) x / cast(real) q;
        frac_coeff_list ~= floor(r);

        int[] numerator_list = [1, cast(int) frac_coeff_list[0]];
        int[] denominator_list = [0, 1];

        int iteration = 2;
        while (true) {
            real frac_part = r - floor(r);
            r = 1.0 / frac_part;
            frac_coeff_list ~= floor(r);

            real numerator = frac_coeff_list[iteration - 1] * numerator_list[iteration - 1] + numerator_list[iteration - 2];
            real denominator = frac_coeff_list[iteration - 1] * denominator_list[iteration - 1] + denominator_list[iteration - 2];
            numerator_list ~= cast(int) numerator;
            denominator_list ~= cast(int) denominator;

            if (denominator_list[$ - 1] > this.n || iteration == 25) {
                break;
            }

            iteration++;
        }

        int result = 0;
        foreach (qr; denominator_list) {
            if ((qr % 2) == 0) {
                if ((pow(this.a, qr) % this.n) == 1) {
                    int val = pow(this.a, qr / 2) % this.n;
                    if (val != this.n - 1 && val != 1) {
                        result = qr;
                        break;
                    }
                }
            }
        }

        if (result == 0) {
            writeln("No valid period r was found, try again with a different value of a");
            exit(0);
        }

        return result;
    }

    /**
    * The main shors algorithm which given an a and an n where a is coprime to n
    * will find non-trivial factors of n probabilisitcally
    */
    void shors() {
        assert(gcd(this.a, this.n) == 1, "The gcd of the inputs is not 1, these inputs can be factored classically");

        // put register 1 into superposition
        for (int i = 0; i < this.qubits_first_register; i++) {
            this.qc.hadamard(i);
        }

        // apply oracle
        oracle();

        // apply inverse qft to entire quantum register
        QFT qft = QFT(this.qc);
        qft.qft_inverse();

        // measure register 1
        int x = measure_reg_1();
        int q = pow(2, this.qubits_first_register);

        // find the continued fraction and approximate r
        int result = find_r(x, q);

        // get non-trivial factors of n based on approximated r
        int factor1 = gcd((pow(this.a, result / 2) % this.n) - 1, this.n);
        int factor2 = gcd((pow(this.a, result / 2) % this.n) + 1, this.n);

        writeln("factors of n: ", factor1, " ", factor2);
    }
}
