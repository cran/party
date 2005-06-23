
void C_split(const double *x, int p,
             const double *y, int q,
             const double *weights, int n,
             const int *orderx, const double *score_y,
             const int ORDERED, SEXP splitctrl, SEXP linexpcov2sample, 
             SEXP expcovinf, double *cutpoint, double *maxstat, 
             double *statistics);

void C_splitcategorical(const int *codingx, int p,
                        const double *y, int q,
                        const double *weights, int n,
                        const double *score_y,
                        const int ORDERED, double *standstat,
                        SEXP splitctrl, SEXP linexpcov2sample, 
                        SEXP expcovinf, double *cutpoint, int *levelset, 
                        double *maxstat, double *statistics);
