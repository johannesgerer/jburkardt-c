void problem_size ( int *chain_num, int *cr_num, int *gen_num, int *pair_num, 
  int *par_num );
void problem_value ( char **chain_filename, char **gr_filename, 
  double *gr_threshold, int *jumpstep, double limits[], int par_num, 
  int *printstep, char **restart_read_filename, char **restart_write_filename );
double prior_density ( int par_num, double zp[] );
double *prior_sample ( int par_num );
double sample_likelihood ( int par_num, double zp[] );

