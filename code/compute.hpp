#ifndef COMPUTE_HPP_6KHWKDTL
#define COMPUTE_HPP_6KHWKDTL


class Compute
{
public:
    Compute ();
    virtual ~Compute ();

    double gyrationRadius(int N, double **r);
    double gyrationRadius(int N, double* r);
    double theta(double *r1, double *r2);
};


#endif /* end of include guard: COMPUTE_HPP_6KHWKDTL */
