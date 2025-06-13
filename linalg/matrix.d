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
    * Add two matrices together
    *
    * params:
    * mat = The matrix to add 
    *
    * returns: A new matrix
    */
    Matrix add_mat(Matrix mat) {
        Matrix!T result_mat = Matrix!T(this.row_num, this.col_num, []);
        Vector!T result_row = Vector!T(this.col_num, []);
        foreach (i, row; this.rows) {
            foreach (j, elem; row.elems) {
                result_row.append(elem + mat.rows[i].elems[j]);
            }
            result_mat.append(result_row);
            result_row.clear();
        }
        return result_mat;
    }

    // operator overload for matrix addition
    Matrix opBinary(string s : "+")(Matrix rhs) {
        return this.add_mat(rhs);
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
    *
    * returns: A transposed matrix with its complex elements cojugated
    */
    Matrix dagger()() if (is(T == Complex!real)) {
        foreach (row; this.rows) {
            foreach (i, elem; row.elems) {
                row[i] = conj(elem);
            }
        }
        return this.transpose();
    }

    /**
    * Get the trace of a square matrix where a result of 1 means the quantum state is pure
    * and less than 1 means the quantum state is mixed.
    * 
    * returns: A real number which is the result of the trace
    */
    real trace()() if (is(T == Complex!real)) {
        assert(this.row_num == this.col_num, "The matrix must be square to take the trace");

        Complex!real trace_sum = Complex!real(0, 0);
        Vector!T diagonal = this.get_diagonal();

        foreach (elem; diagonal.elems) {
            trace_sum = trace_sum + elem;
        }

        return trace_sum.re;
    }

    Vector!T get_diagonal() {
        Vector!T diagonal = Vector!T(this.row_num, []);

        int diagonal_idx = 0;
        foreach (row; this.rows) {
            diagonal.append(row[diagonal_idx]);
            diagonal_idx++;
        }
        return diagonal;
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

// NOTE: It doesn't make sense to have to declare a Matrix instance when
// calling this function:

/**
 * Makes a square matrix with all zeros as elements
 * 
 * params:
 * size = The size of the square matrix
 *
 * returns: The matrix with all zeros
*/
Matrix!(Complex!real) zeros(int size) {
    Matrix!(Complex!real) result = Matrix!(Complex!real)(size, size, []);
    for (int i = 0; i < size; i++) {
        Vector!(Complex!real) row = Vector!(Complex!real)(size, []);
        for (int j = 0; j < size; j++) {
            row.append(Complex!real(0, 0));
        }
        result.append(row);
        row.clear();
    }
    return result;
}
