#define PI 3.141592653589793

#define delta(i, j) ((i)==(j) ? 1 : 0) /* \delta_{ij} */
#define neg1toi(i) ((i)%2 ? -1 : 1) /* (-1)^i */
#define cosipiby2(i) ((i)%2 ? 0 : ((i)%4==0 ? 1 : -1)) /* \cos(i\pi/2) */

coeff *coeff_alloc(int z, int n, int L);
void coeff_set(coeff *c, int z, int n, int L, int m, int i, double x);
double coeff_get(coeff *c, int z, int n, int L, int m, int i);
void print_coeff(coeff *c);

bound_coeff *bound_coeff_alloc(int z, int L);
void bound_coeff_set(bound_coeff *b, int z, int L, int m, int i, double x);
double bound_coeff_get(bound_coeff *b, int z, int L, int m, int i);
void print_bound_coeff(bound_coeff *b);
void print_vector(gsl_vector *v);

void legendre_zero_weight(gsl_vector *zeros, gsl_vector *weights);
double p_lm(int l, int m, double x);
double ytheta_lm(int l, int m, double x);
void decompose_boundary(bound_coeff *b_zlm, 
						int nphi,
						int nx,
						double (*boundary)(int, double, double));
