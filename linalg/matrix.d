module linalg.matrix;

import std.stdio;
import std.complex;

import linalg.vector;

struct Matrix(T) {
    Vector!T[] rows;
    int row_num;
    int col_num;

    /**
    * Constructor for the Matrix object
    *
    * params:
    * row_num = the number of rows to be in the matrix
    * col_num = the number of columns to be in the matrix
    * rows = The actual rows of the matrix as an array of Vector objects
    */
    this(int row_num, int col_num, Vector!T[] rows) {
        this.row_num = row_num;
        this.col_num = col_num;
        this.rows = rows;
    }

    /**
    * Get the columns of the matrix
    * 
    * returns: A vector of type T array 
    */
    Vector!T[] get_cols() {
        Vector!T[] cols = new Vector!T[this.col_num];
        foreach (i; 0 .. col_num) {
            T[] col = new T[this.rows.length];
            foreach (j, row; this.rows) {
                col[j] = row[i];
            }
            cols[i] = Vector!T(this.row_num, col);
        }
        return cols;
    }

    /**
    * Append a row as a Vector to the matrix
    *
    * params: 
    * row = the row to append to the matrix
    */
    void append(Vector!T row) {
        this.rows[this.rows.length++] = row;
        this.row_num = cast(int) this.rows.length;
    }

    /**
    * Multiplies a matrix a matrix by a matrix
    *
    * params:
    * mat = The matrix to multiply by
    *
    * returns: A new matrix
    */
    Matrix mult_mat(Matrix mat) {
        assert(this.col_num == mat.row_num, "Cannot multiply due to size incompatibility");

        Vector!T[] cols = mat.get_cols();
        Matrix!T result_mat = Matrix!T(this.row_num, mat.col_num, []);
        Vector!T result_row;

        foreach (row; this.rows) {
            foreach (col; cols) {
                T result = row * col;
                result_row.append(result);
            }
            result_mat.append(result_row);
            result_row.clear();
        }

        return result_mat;
    }

    /**
    * Multiplies a matrix by a vector 
    *
    * params:
    * vec = the vector to multiply by
    *
    * returns: A new vector
    */
    Vector!T mult_vec(Vector!T vec) {
        assert(this.col_num == vec.length(),
            "The length of the vector must be equal to length of a row in the matrix");

        Vector!T result = Vector!T(this.col_num, []);

        foreach (row; this.rows) {
            result.append(row.dot(vec));
        }

        return result;
    }

    /**
    * Multiplies a matrix by a scalar value
    *
    * params:
    * scalar = A scalar value of type T
    *
    * returns: A new matrix
    */
    Matrix mult_scalar(T scalar) {
        Matrix!T result_mat = Matrix!T(this.row_num, this.col_num, []);
        Vector!T result_vec = Vector!T(this.col_num, []);

        foreach (row; this.rows) {
            foreach (elem; row.elems) {
                T result = scalar * elem;
                result_vec.append(result);
            }
            result_mat.append(result_vec);
            result_vec.clear();
        }
        return result_mat;
    }

    /**
    * Kronecker or tensor product two matrices
    *
    * params:
    * target = the matrix to operate with
    *
    * returns: A new matrix
    */
    Matrix kronecker(Matrix target) {
        Matrix!T result_mat = Matrix!T(target.row_num * this.row_num, target.col_num * this.col_num, [
            ]);
        Vector!T result;

        foreach (row; this.rows) {
            foreach (target_row; target.rows) {
                foreach (elem; row.elems) {
                    foreach (target_elem; target_row.elems) {
                        result.append(elem * target_elem);
                    }
                }
                result_mat.append(result);
                result.clear();
            }
        }
        return result_mat;
    }

    /**
    * Takes the quantum inner product of a complex conjugated and transposed row vector with a 
    * complex valued row vector. This is used for expectation value calculaiion.
    *
    * params:
    * target = The vector to multiply by
    * 
    * returns: A real number, in the expectation value case, the expectation value
    */
    real inner_product(Vector!T target) {
        assert(this.col_num == 1, "The matrix for the inner product operation must have column dimension of one");

        Complex!real sum = Complex!real(0, 0);

        foreach (i, vec; this.rows) {
            sum += vec[0] * target[i];
        }

        return sum.re;
    }

    /**
    * Generates a square identity matrix 
    * 
    * params:
    * dim = The dimensions of the matrix
    *
    * reuturns: A square identity matrix
    */
    Matrix identity(int dim) {
        Matrix!T result = Matrix!T(dim, dim, []);
        Vector!T row = Vector!T(dim, []);

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (j == i) {
                    row.append(cast(T) 1);
                } else {
                    row.append(cast(T) 0);
                }
            }
            result.append(row);
            row.clear();
        }
        return result;
    }

    /**
    * Take the complex conjugate transpose of a matrix, this only works
    * if the template T is Complex!real
    */
    Matrix dagger()() if (is(T == Complex!real)) {
        foreach (row; this.rows) {
            foreach (i, elem; row.elems) {
                row[i] = conj(elem);
            }
        }
        return this.transpose();
    }

    //operator overload for matrix multiplication
    Matrix opBinary(string s : "*")(Matrix rhs) {
        return this.mult_mat(rhs);
    }

    //operator overloading for matrix-vector multiplication
    Vector!T opBinary(string s : "*")(Vector!T rhs) {
        return this.mult_vec(rhs);
    }

    /**
    * Transpose a matrix
    *
    * returns: A new matrix with the original rows as columns and original columns
    *          as rows
    */
    Matrix transpose() {
        Matrix!T transpose_mat = Matrix!T(this.col_num, this.row_num, this.get_cols());
        return transpose_mat;
    }

}
