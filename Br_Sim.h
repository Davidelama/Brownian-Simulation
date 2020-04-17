int Ntrial, N, Nsim, Neq, Ncell, Nmax, dim, Nneigh, Nmeas, Nhist, Nlevel, Nblock;
double L, sigma, rho, rc, T, dt, mu;
double *x, *y, *z, *vx, *vy, *vz, *fx, *fy, *fz;
int *cell_list_index, *hist;
int **cell_list, **neighbor;

class RanMars *random_mars;

void rand_init_pos();

void init_vel();

void init_output(FILE* &fout,FILE* &gout,FILE* &sfout, FILE* &zout, FILE* &drout);  //here I pass by reference (&) a FILE* type (i.e. pointer to FILE), that allows me to treat the FILE* variable as a global one 

void rescale_vel();

double dist(double x1, double y1, double z1, double x2, double y2, double z2);

void print_pos(int step);

void print_g(FILE* &gout);

void print_sf(FILE* &sfout);

void print_single_g(int step);

void print_cell(int step);

void print_block(double** &x_blk, double** &y_blk, double** &z_blk);

void latt_init_pos();

double calculate_force();

double calculate_kin();

double calculate_pot();

void print_SF();

double calculate_pressure();

void md_step();

void bd_step();

void init_cell_list();

void reinit_cell_list();

void test_cell(int step);

void close_output(FILE* &fout,FILE* &gout,FILE* &sfout, FILE* &zout, FILE* &drout);

void init_block(double** &x_blk, double** &y_blk, double** &z_blk, double* &Zx, double* &Zy, double* &Zz);

void block_corr(int step, double** &x_blk, double** &y_blk, double** &z_blk, double* &Zx, double* &Zy, double* &Zz, FILE* &zout);

void block_intg(int step, double** &x_blk, double** &y_blk, double** &z_blk, double* &Zx, double* &Zy, double* &Zz, FILE* &drout);

//remember damn ; at the end