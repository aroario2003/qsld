module linalg.matrix;

import std.stdio;

import linalg.vector;

struct Matrix(T) {
    Vector!T[] rows;
    int row_num;
    int col_num;

    this(int row_num, int col_num, Vector!T[] rows) {
        this.row_num = row_num;
        this.col_num = col_num;
        this.rows = rows;
    }

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

    void append(Vector!T row) {
        this.rows[this.rows.length++] = row;
    }

    Matrix mult_mat(Matrix mat) {
        assert(this.col_num == mat.row_num, "Cannot multiply due to size incompatibility");

        Vector!T[] cols = mat.get_cols();
        Matrix!T result_mat;
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

    Vector!T mult_vec(Vector!T vec) {
        assert(this.col_num == vec.length(),
            "The length of the vector must be equal to length of a row in the matrix");

        Vector!T result = Vector!T(this.col_num, []);

        foreach (row; this.rows) {
            result.append(row.dot(vec));
        }

        return result;
    }

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

    Matrix kronecker(Matrix target) {
        Matrix!T result_mat;
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

    //operator overload for matrix multiplication
    Matrix opBinary(string s : "*")(Matrix rhs) {
        return this.mult_mat(rhs);
    }

    //operator overloading for matrix-vector multiplication
    Vector!T opBinary(string s : "*")(Vector!T rhs) {
        return this.mult_vec(rhs);
    }

    Matrix transpose() {
        Matrix!T transpose_mat = Matrix!T(this.col_num, this.row_num, this.get_cols());
        return transpose_mat;
    }

}
