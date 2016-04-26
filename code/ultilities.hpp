#ifndef ULTILITIES_HPP_QBDNOQTX
#define ULTILITIES_HPP_QBDNOQTX


double Dot(double* dx, double *dy, int dim);
double Distance(double *point1, double *point2, int dim);
double Mean(int N, double* x);
double Variance(int N, double* x);
void MatrixMulVector(double** matrix,
		double* vector, int dim);


static inline int Delta(int i, int j) {
    return i == j ? 1 : 0;
}


// allocate memory for 2d array
template <typename T>
T** create2DArray(unsigned nrows, unsigned ncols) {
    T** ptr = new T*[nrows];  // allocate pointers
    T* pool = new T[nrows*ncols];  // allocate pool
    for (unsigned i = 0; i < nrows; ++i, pool += ncols )
        ptr[i] = pool;
    return ptr;
}

// free memory for 2d array
template <typename T>
void delete2DArray(T** arr) {
    delete [] arr[0];  // remove the pool
    delete [] arr;     // remove the pointers
}

#endif /* end of include guard: ULTILITIES_HPP_QBDNOQTX */
