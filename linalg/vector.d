module linalg.vector;

import std.stdio;
import std.math;

struct Vector(T) {
    int dim;
    T[] elems;

    this(int dimension, T[] args) {
        this.dim = dimension;
        this.elems = args;
    }

    Vector add(Vector v2) {
        assert(this.dim == v2.dim, "Vectors must be equal in dimension to add");

        T[] elems = new T[this.dim];
        Vector!T v3 = Vector!T(this.dim, elems);

        foreach (i, elem; this.elems) {
            T result = elem + v2.elems[i];
            v3.elems[i] = result;
        }

        return v3;
    }

    Vector opBinary(string s : "+")(Vector rhs) {
        return this.add(rhs);
    }

    Vector sub(Vector v2) {
        assert(this.dim == v2.dim, "Vectors must be equal in dimension to subtract");

        T[] elems = new T[this.dim];
        Vector!T v3 = Vector!T(this.dim, elems);

        foreach (i, elem; this.elems) {
            T result = elem - v2.elems[i];
            v3.elems[i] = result;
        }

        return v3;
    }

    Vector opBinary(string s : "-")(Vector rhs) {
        return this.sub(rhs);
    }

    T dot(Vector v2) {
        assert(this.dim == v2.dim, "Vectors must be equal in dimension to take the dot product");

        T sum = 0;

        foreach (i, elem; this.elems) {
            T result = elem * v2.elems[i];
            sum = sum + result;
        }

        return sum;
    }

    T opBinary(string s : "*")(Vector rhs) {
        return this.dot(rhs);
    }

    Vector mult(T scalar) {
        T[] elems = new T[this.dim];
        Vector!T v = Vector!T(this.dim, elems);

        foreach (i, elem; this.elems) {
            v.elems[i] = elem * scalar;
        }

        return v;
    }

    Vector opBinaryRight(string s : "*")(T lhs) {
        return this.mult(lhs);
    }

    Vector opBinary(string s : "*")(T rhs) {
        return this.mult(rhs);
    }

    void append(T elem) {
        this.elems[elems.length++] = elem;
    }

    void clear() {
        this.elems = [];
        this.elems.length = 0;
    }

    // Array operator overloading begin
    T opIndex(size_t i) const {
        return this.elems[i];
    }

    void opIndexAssign(T value, size_t i) {
        this.elems[i] = value;
    }

    size_t length() const {
        return this.dim;
    }
    // Array operator overloading end
}
