# include <complex.h>
# include <math.h>
# include <stdbool.h>
# include <stdio.h>
# include <stdlib.h>
# include <strlen.h>

# include "r8lib.h"

int main ( );

void i4int_to_r8int_test ( );

void perm0_check_test ( );
void perm0_uniform_test ( );

void perm1_check_test ( );
void perm1_uniform_test ( );

void r8_abs_test ( );
void r8_acos_test ( );
void r8_acosh_test ( );
void r8_asinh_test ( );
void r8_atan_test ( );
void r8_atanh_test ( );
void r8_big_test ( );
void r8_cas_test ( );
void r8_ceiling_test ( );
void r8_choose_test ( );
void r8_cosd_test ( );
void r8_cotd_test ( );
void r8_cscd_test ( );
void r8_cube_root_test ( );
void r8_diff_test ( );
void r8_digit_test ( );
void r8_e_test ( );
void r8_epsilon_test ( );
void r8_epsilon_compute_test ( );
void r8_factorial_test ( );
void r8_factorial2_test ( );
void r8_fall_test ( );
void r8_fractional_test ( );
void r8_gamma_test ( );
void r8_gamma_log_test ( );
void r8_huge_test ( );
void r8_log_2_test ( );
void r8_log_b_test ( );
void r8_mant_test ( );
void r8_max_test ( );
void r8_min_test ( );
void r8_mod_test ( );
void r8_modp_test ( );
void r8_mop_test ( );
void r8_nint_test ( );
void r8_normal_01_test ( );
void r8_pi_test ( );
void r8_power_test ( );
void r8_power_fast_test ( );
void r8_rise_test ( );
void r8_round2_test ( );
void r8_roundb_test ( );
void r8_roundx_test ( );
void r8_secd_test ( );
void r8_sign_test ( );
void r8_sign3_test ( );
void r8_sind_test ( );
void r8_swap_test ( );
void r8_swap3_test ( );
void r8_tand_test ( );
void r8_to_r8_discrete_test ( );
void r8_to_i4_test ( );
void r8_uniform_01_test ( );
void r8_uniform_ab_test ( );
void r8_walsh_1d_test ( );
void r8_wrap_test ( );

void r82col_print_part_test ( );

void r82poly2_type_test ( );

void r82row_order_type_test ( );
void r82row_part_quick_a_test ( );
void r82row_print_part_test ( );
void r82row_sort_heap_index_a_test ( );
void r82row_sort_quick_a_test ( );

void r83col_print_part_test ( );

void r83row_print_part_test ( );

void r8block_expand_linear_test ( );
void r8block_new_test ( );
void r8block_print_test ( );

void r8cmat_to_r8mat_new_test ( );

void r8col_find_test ( );
void r8col_insert_test ( );
void r8col_sort_heap_a_test ( );
void r8col_sort_heap_index_a_test ( );
void r8col_sort_quick_a_test ( );
void r8col_sorted_tol_unique_test ( );
void r8col_sorted_unique_count_test ( );
void r8col_sorted_tol_undex_test ( );
void r8col_max_test ( );
void r8col_mean_test ( );
void r8col_min_test ( );
void r8col_permute_test ( );
void r8col_sortr_a_test ( );
void r8col_sum_test ( );
void r8col_swap_test ( );
void r8col_to_r8vec_test ( );
void r8col_tol_undex_test ( );
void r8col_undex_test ( );
void r8col_unique_count_test ( );
void r8col_variance_test ( );

void r8int_to_i4int_test ( );

void r8mat_cholesky_inverse_test ( );
void r8mat_cholesky_solve_test ( );
void r8mat_cholesky_solve_upper_test ( );
void r8mat_det_2d_test ( );
void r8mat_det_3d_test ( );
void r8mat_det_4d_test ( );
void r8mat_det_5d_test ( );
void r8mat_expand_linear_test ( );
void r8mat_expand_linear2_test ( );
void r8mat_fs_new_test ( );
void r8mat_fss_new_test ( );
void r8mat_givens_post_test ( );
void r8mat_givens_pre_test ( );
void r8mat_hess_test ( );
double r8mat_hess_f ( int n, double x[] );
double *r8mat_hess_exact ( int n, double x[] );
void r8mat_house_axh_test ( );
void r8mat_house_form_test ( );
void r8mat_house_post_test ( );
void r8mat_house_pre_test ( );
void r8mat_indicator_new_test ( );
void r8mat_inverse_2d_test ( );
void r8mat_inverse_3d_test ( );
void r8mat_inverse_4d_test ( );
void r8mat_jac_test ( );
double *r8mat_jac_f ( int m, int n, double x[] );
double *r8mat_jac_exact ( int m, int n, double x[] );
void r8mat_kronecker_test ( );
void r8mat_l_inverse_test ( );
void r8mat_l_print_test ( );
void r8mat_l1_inverse_test ( );
void r8mat_lu_test ( );
void r8mat_max_test ( );
void r8mat_max_index_test ( );
void r8mat_maxcol_minrow_test ( );
void r8mat_maxrow_mincol_test ( );
void r8mat_min_test ( );
void r8mat_min_index_test ( );
void r8mat_mincol_maxrow_test ( );
void r8mat_minrow_maxcol_test ( );
void r8mat_mm_test ( );
void r8mat_mm_new_test ( );
void r8mat_mv_test ( );
void r8mat_mv_new_test ( );
void r8mat_mtv_new_test ( );
void r8mat_mtv_test ( );
void r8mat_nint_test ( );
void r8mat_nonzeros_test ( );
void r8mat_norm_fro_test ( );
void r8mat_norm_l1_test ( );
void r8mat_nullspace_test ( );
void r8mat_nullspace_size_test ( );
void r8mat_orth_uniform_new_test ( );
void r8mat_plot_test ( );
void r8mat_power_method_test ( );
void r8mat_print_test ( );
void r8mat_print_some_test ( );
void r8mat_ref_test ( );
void r8mat_rref_test ( );
void r8mat_solve_test ( );
void r8mat_solve_2d_test ( );
void r8mat_solve_3d_test ( );
void r8mat_solve2_test ( );
void r8mat_sub_new_test ( );
void r8mat_symm_jacobi_test ( );
void r8mat_to_r8cmat_new_test ( );
void r8mat_to_r8plu_test ( );
void r8mat_to_r8rmat_test ( );
void r8mat_trace_test ( );
void r8mat_transpose_new_test ( );
void r8mat_transpose_print_test ( );
void r8mat_u_inverse_test ( );
void r8mat_u1_inverse_test ( );
void r8mat_uniform_ab_new_test ( );

void r8plu_det_test ( );
void r8plu_inverse_test ( );
void r8plu_mul_test ( );
void r8plu_sol_test ( );
void r8plu_to_r8mat_test ( );

void r8poly_degree_test ( );
void r8poly_deriv_test ( );
void r8poly_lagrange_coef_test ( );
void r8poly_lagrange_0_test ( );
void r8poly_lagrange_1_test ( );
void r8poly_lagrange_2_test ( );
void r8poly_lagrange_factor_test ( );
void r8poly_lagrange_val_test ( );
void r8poly_print_test ( );
void r8poly_value_horner_test ( );
void r8poly_values_horner_test ( );

void r8poly2_ex_test ( );
void r8poly2_ex2_test ( );
void r8poly2_val_test ( );
void r8poly2_f ( double x, double *y, double *yp, double *ypp );
void r8poly2_val2_test ( );

void r8r8vec_index_insert_unique_test ( );

void r8r8r8vec_index_insert_unique_test ( );

void r8rmat_new_test ( );
void r8rmat_to_r8mat_test ( );

void r8row_max_test ( );
void r8row_mean_test ( );
void r8row_min_test ( );
void r8row_sum_test ( );
void r8row_swap_test ( );
void r8row_to_r8vec_test ( );
void r8row_variance_test ( );

void r8slmat_print_test ( );

void r8vec_amax_test ( );
void r8vec_amin_test ( );
void r8vec_bracket_test ( );
void r8vec_bracket2_test ( );
void r8vec_bracket3_test ( );
void r8vec_bracket5_test ( );
void r8vec_chebyspace_new_test ( );
void r8vec_concatenate_new_test ( );
void r8vec_convolution_test ( );
void r8vec_convolution_circ_test ( );
void r8vec_dif_test ( );
double r8vec_dif_f ( double x );
void r8vec_direct_product_test ( );
void r8vec_direct_product2_test ( );
void r8vec_even_test ( );
void r8vec_even2_test ( );
void r8vec_expand_linear_test ( );
void r8vec_frac_test ( );
void r8vec_histogram_test ( );
void r8vec_house_column_test ( );
void r8vec_index_delete_all_test ( );
void r8vec_index_delete_dupes_test ( );
void r8vec_index_delete_one_test ( );
void r8vec_index_insert_test ( );
void r8vec_index_insert_unique_test ( );
void r8vec_index_order_test ( );
void r8vec_index_search_test ( );
void r8vec_index_sorted_range_test ( );
void r8vec_indexed_heap_d_test ( );
void r8vec_indexed_heap_d_extract_test ( );
void r8vec_indexed_heap_d_insert_test ( );
void r8vec_indexed_heap_d_max_test ( );
void r8vec_indicator0_new_test ( );
void r8vec_legendre_test ( );
void r8vec_linspace_new_test ( );
void r8vec_max_test ( );
void r8vec_max_index_test ( );
void r8vec_mean_test ( );
void r8vec_median_test ( );
void r8vec_midspace_new_test ( );
void r8vec_min_test ( );
void r8vec_min_index_test ( );
void r8vec_nint_test ( );
void r8vec_norm_l0_test ( );
void r8vec_norm_l1_test ( );
void r8vec_norm_l2_test ( );
void r8vec_norm_li_test ( );
void r8vec_normal_01_new_test ( );
void r8vec_normalize_l1_test ( );
void r8vec_order_type_test ( );
void r8vec_permute_test ( );
void r8vec_permute_uniform_test ( );
void r8vec_polarize_test ( );
void r8vec_print_test ( );
void r8vec_rotate_test ( );
void r8vec_reverse_test ( );
void r8vec_search_binary_a_test ( );
void r8vec_sort_bubble_a_test ( );
void r8vec_sort_heap_a_test ( );
void r8vec_sort_heap_d_test ( );
void r8vec_sort_heap_index_a_new_test ( );
void r8vec_sort_heap_index_d_new_test ( );
void r8vec_sort_heap_mask_a_test ( );
void r8vec_sort_insert_a_test ( );
void r8vec_sort_insert_index_a_test ( );
void r8vec_sort_quick_a_test ( );
void r8vec_sorted_merge_a_test ( );
void r8vec_sorted_nearest_test ( );
void r8vec_sorted_range_test ( );
void r8vec_sorted_split_test ( );
void r8vec_sorted_undex_test ( );
void r8vec_sorted_unique_test ( );
void r8vec_sorted_unique_count_test ( );
void r8vec_sorted_unique_hist_test ( );
void r8vec_split_test ( );
void r8vec_transpose_print_test ( );
void r8vec_undex_test ( );
void r8vec_uniform_01_new_test ( );
void r8vec_uniform_ab_new_test ( );
void r8vec_variance_test ( );

void r8vec2_sort_a_test ( );
void r8vec2_sort_d_test ( );
void r8vec2_sort_heap_index_a_test ( );
void r8vec2_sorted_unique_test ( );
void r8vec2_sorted_unique_index_test ( );
void r8vec2_sum_max_index_test ( );

void roots_to_r8poly_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for R8LIB_PRB.

  Discussion:

    R8LIB_PRB tests the R8LIB library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 May 2015

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "R8LIB_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the R8LIB library.\n" );

  i4int_to_r8int_test ( );

  perm0_check_test ( );
  perm0_uniform_test ( );

  perm1_check_test ( );
  perm1_uniform_test ( );

  r8_abs_test ( );
  r8_acos_test ( );
  r8_acosh_test ( );
  r8_asinh_test ( );
  r8_atan_test ( );
  r8_atanh_test ( );
  r8_big_test ( );
  r8_cas_test ( );
  r8_ceiling_test ( );
  r8_choose_test ( );
  r8_cosd_test ( );
  r8_cotd_test ( );
  r8_cscd_test ( );
  r8_cube_root_test ( );
  r8_diff_test ( );
  r8_digit_test ( );
  r8_e_test ( );
  r8_epsilon_test ( );
  r8_epsilon_compute_test ( );
  r8_factorial_test ( );
  r8_factorial2_test ( );
  r8_fall_test ( );
  r8_fractional_test ( );
  r8_gamma_test ( );
  r8_gamma_log_test ( );
  r8_huge_test ( );
  r8_log_2_test ( );
  r8_log_b_test ( );
  r8_mant_test ( );
  r8_max_test ( );
  r8_min_test ( );
  r8_mod_test ( );
  r8_modp_test ( );
  r8_mop_test ( );
  r8_nint_test ( );
  r8_normal_01_test ( );
  r8_pi_test ( );
  r8_power_test ( );
  r8_power_fast_test ( );
  r8_rise_test ( );
  r8_round2_test ( );
  r8_roundb_test ( );
  r8_roundx_test ( );
  r8_secd_test ( );
  r8_sign_test ( );
  r8_sign3_test ( );
  r8_sind_test ( );
  r8_swap_test ( );
  r8_swap3_test ( );
  r8_tand_test ( );
  r8_to_r8_discrete_test ( );
  r8_to_i4_test ( );
  r8_uniform_01_test ( );
  r8_uniform_ab_test ( );
  r8_walsh_1d_test ( );
  r8_wrap_test ( );

  r82col_print_part_test ( );

  r82poly2_type_test ( );

  r82row_order_type_test ( );
  r82row_part_quick_a_test ( );
  r82row_print_part_test ( );
  r82row_sort_heap_index_a_test ( );
  r82row_sort_quick_a_test ( );

  r83col_print_part_test ( );

  r83row_print_part_test ( );

  r8block_expand_linear_test ( );
  r8block_new_test ( );
  r8block_print_test ( );

  r8cmat_to_r8mat_new_test ( );

  r8col_find_test ( );
  r8col_insert_test ( );
  r8col_sort_heap_a_test ( );
  r8col_sort_heap_index_a_test ( );
  r8col_sort_quick_a_test ( );
  r8col_sorted_tol_unique_test ( );
  r8col_sorted_unique_count_test ( );
  r8col_sorted_tol_undex_test ( );
  r8col_max_test ( );
  r8col_mean_test ( );
  r8col_min_test ( );
  r8col_permute_test ( );
  r8col_sortr_a_test ( );
  r8col_sum_test ( );
  r8col_swap_test ( );
  r8col_to_r8vec_test ( );
  r8col_tol_undex_test ( );
  r8col_undex_test ( );
  r8col_unique_count_test ( );
  r8col_variance_test ( );

  r8r8vec_index_insert_unique_test ( );

  r8r8r8vec_index_insert_unique_test ( );

  r8int_to_i4int_test ( );

  r8mat_cholesky_inverse_test ( );
  r8mat_cholesky_solve_test ( );
  r8mat_cholesky_solve_upper_test ( );
  r8mat_det_2d_test ( );
  r8mat_det_3d_test ( );
  r8mat_det_4d_test ( );
  r8mat_det_5d_test ( );
  r8mat_expand_linear_test ( );
  r8mat_expand_linear2_test ( );
  r8mat_fs_new_test ( );
  r8mat_fss_new_test ( );
  r8mat_givens_post_test ( );
  r8mat_givens_pre_test ( );
  r8mat_hess_test ( );
  r8mat_house_axh_test ( );
  r8mat_house_form_test ( );
  r8mat_house_post_test ( );
  r8mat_house_pre_test ( );
  r8mat_indicator_new_test ( );
  r8mat_inverse_2d_test ( );
  r8mat_inverse_3d_test ( );
  r8mat_inverse_4d_test ( );
  r8mat_jac_test ( );
  r8mat_kronecker_test ( );
  r8mat_l_inverse_test ( );
  r8mat_l_print_test ( );
  r8mat_l1_inverse_test ( );
  r8mat_lu_test ( );
  r8mat_max_test ( );
  r8mat_max_index_test ( );
  r8mat_maxcol_minrow_test ( );
  r8mat_maxrow_mincol_test ( );
  r8mat_min_test ( );
  r8mat_min_index_test ( );
  r8mat_mincol_maxrow_test ( );
  r8mat_minrow_maxcol_test ( );
  r8mat_mm_test ( );
  r8mat_mm_new_test ( );
  r8mat_mv_test ( );
  r8mat_mv_new_test ( );
  r8mat_mtv_new_test ( );
  r8mat_mtv_test ( );
  r8mat_nint_test ( );
  r8mat_nonzeros_test ( );
  r8mat_norm_fro_test ( );
  r8mat_norm_l1_test ( );
  r8mat_nullspace_test ( );
  r8mat_nullspace_size_test ( );
  r8mat_orth_uniform_new_test ( );
  r8mat_plot_test ( );
  r8mat_power_method_test ( );
  r8mat_print_test ( );
  r8mat_print_some_test ( );
  r8mat_ref_test ( );
  r8mat_rref_test ( );
  r8mat_solve_test ( );
  r8mat_solve_2d_test ( );
  r8mat_solve_3d_test ( );
  r8mat_solve2_test ( );
  r8mat_sub_new_test ( );
  r8mat_symm_jacobi_test ( );
  r8mat_to_r8cmat_new_test ( );
  r8mat_to_r8plu_test ( );
  r8mat_to_r8rmat_test ( );
  r8mat_trace_test ( );
  r8mat_transpose_new_test ( );
  r8mat_transpose_print_test ( );
  r8mat_u_inverse_test ( );
  r8mat_u1_inverse_test ( );
  r8mat_uniform_ab_new_test ( );

  r8plu_det_test ( );
  r8plu_inverse_test ( );
  r8plu_mul_test ( );
  r8plu_sol_test ( );
  r8plu_to_r8mat_test ( );

  r8poly_degree_test ( );
  r8poly_deriv_test ( );
  r8poly_lagrange_coef_test ( );
  r8poly_lagrange_0_test ( );
  r8poly_lagrange_1_test ( );
  r8poly_lagrange_2_test ( );
  r8poly_lagrange_factor_test ( );
  r8poly_lagrange_val_test ( );
  r8poly_print_test ( );
  r8poly_value_horner_test ( );
  r8poly_values_horner_test ( );

  r8poly2_ex_test ( );
  r8poly2_ex2_test ( );
  r8poly2_val_test ( );
  r8poly2_val2_test ( );

  r8rmat_new_test ( );
  r8rmat_to_r8mat_test ( );

  r8row_max_test ( );
  r8row_mean_test ( );
  r8row_min_test ( );
  r8row_sum_test ( );
  r8row_swap_test ( );
  r8row_to_r8vec_test ( );
  r8row_variance_test ( );

  r8slmat_print_test ( );

  r8vec_amax_test ( );
  r8vec_amin_test ( );
  r8vec_bracket_test ( );
  r8vec_bracket2_test ( );
  r8vec_bracket3_test ( );
  r8vec_bracket5_test ( );
  r8vec_chebyspace_new_test ( );
  r8vec_concatenate_new_test ( );
  r8vec_convolution_test ( );
  r8vec_convolution_circ_test ( );
  r8vec_dif_test ( );
  r8vec_direct_product_test ( );
  r8vec_direct_product2_test ( );
  r8vec_even_test ( );
  r8vec_even2_test ( );
  r8vec_expand_linear_test ( );
  r8vec_frac_test ( );
  r8vec_histogram_test ( );
  r8vec_house_column_test ( );
  r8vec_index_delete_all_test ( );
  r8vec_index_delete_dupes_test ( );
  r8vec_index_delete_one_test ( );
  r8vec_index_insert_test ( );
  r8vec_index_insert_unique_test ( );
  r8vec_index_order_test ( );
  r8vec_index_search_test ( );
  r8vec_index_sorted_range_test ( );
  r8vec_indexed_heap_d_test ( );
  r8vec_indexed_heap_d_extract_test ( );
  r8vec_indexed_heap_d_insert_test ( );
  r8vec_indexed_heap_d_max_test ( );
  r8vec_indicator0_new_test ( );
  r8vec_legendre_test ( );
  r8vec_linspace_new_test ( );
  r8vec_max_test ( );
  r8vec_max_index_test ( );
  r8vec_mean_test ( );
  r8vec_median_test ( );
  r8vec_midspace_new_test ( );
  r8vec_min_test ( );
  r8vec_min_index_test ( );
  r8vec_nint_test ( );
  r8vec_norm_l0_test ( );
  r8vec_norm_l1_test ( );
  r8vec_norm_l2_test ( );
  r8vec_norm_li_test ( );
  r8vec_normal_01_new_test ( );
  r8vec_normalize_l1_test ( );
  r8vec_order_type_test ( );
  r8vec_permute_test ( );
  r8vec_permute_uniform_test ( );
  r8vec_polarize_test ( );
  r8vec_print_test ( );
  r8vec_rotate_test ( );
  r8vec_reverse_test ( );
  r8vec_search_binary_a_test ( );
  r8vec_sort_bubble_a_test ( );
  r8vec_sort_heap_a_test ( );
  r8vec_sort_heap_d_test ( );
  r8vec_sort_heap_index_a_new_test ( );
  r8vec_sort_heap_index_d_new_test ( );
  r8vec_sort_heap_mask_a_test ( );
  r8vec_sort_insert_a_test ( );
  r8vec_sort_insert_index_a_test ( );
  r8vec_sort_quick_a_test ( );
  r8vec_sorted_merge_a_test ( );
  r8vec_sorted_nearest_test ( );
  r8vec_sorted_range_test ( );
  r8vec_sorted_split_test ( );
  r8vec_sorted_undex_test ( );
  r8vec_sorted_unique_test ( );
  r8vec_sorted_unique_count_test ( );
  r8vec_sorted_unique_hist_test ( );
  r8vec_split_test ( );
  r8vec_transpose_print_test ( );
  r8vec_undex_test ( );
  r8vec_uniform_01_new_test ( );
  r8vec_uniform_ab_new_test ( );
  r8vec_variance_test ( );

  r8vec2_sort_a_test ( );
  r8vec2_sort_d_test ( );
  r8vec2_sort_heap_index_a_test ( );
  r8vec2_sorted_unique_test ( );
  r8vec2_sorted_unique_index_test ( );
  r8vec2_sum_max_index_test ( );

  roots_to_r8poly_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "R8LIB_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void i4int_to_r8int_test ( )

/******************************************************************************/
/*
  Purpose:

    I4INT_TO_R8INT_TEST tests I4INT_TO_R8INT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2012

  Author:

    John Burkardt
*/
{
  int i;
  int ihi = 11;
  int ilo = 1;
  int ir;
  double r;
  double r2;
  double rhi = 200.0;
  double rhi2;
  double rlo = 100.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "I4INT_TO_R8INT_TEST\n" );
  printf ( "  For data in an interval,\n" );
  printf ( "  I4INT_TO_R8INT converts an integer to a real;\n" );
  printf ( "\n" );
  printf ( "  Integer interval: [%d, %d]\n", ilo, ihi );
  printf ( "  Real interval:    [%g, %g]\n", rlo, rhi );
  printf ( "\n" );
  printf ( "         R          I(R)       R(I(R))\n" );
  printf ( "\n" );

  seed = 123456789;

  rlo2 = rlo - 15.0;
  rhi2 = rhi + 15.0;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, &seed );
    ir = r8int_to_i4int ( rlo, rhi, r, ilo, ihi );
    r2 = i4int_to_r8int ( ilo, ihi, ir, rlo, rhi );
    printf ( "  %12g  %6d  %12g\n", r, ir, r2 );
  }

  return;
}
/******************************************************************************/

void perm0_check_test ( )

/******************************************************************************/
/*
  Purpose:

    PERM0_CHECK_TEST tests PERM0_CHECK.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 May 2015

  Author:

    John Burkardt
*/
{
  int ierror;
  int n = 5;
  int p1[5] = { 5, 2, 3, 4, 1 };
  int p2[5] = { 4, 1, 3, 0, 2 };
  int p3[5] = { 0, 2, 1, 3, 2 };

  printf ( "\n" );
  printf ( "PERM0_CHECK_TEST\n" );
  printf ( "  PERM0_CHECK checks a permutation of 0, ..., N-1.\n" );
  printf ( "\n" );

  i4vec_transpose_print ( n, p1, "  Permutation 1:" );
  ierror = perm0_check( n, p1 );

  i4vec_transpose_print ( n, p2, "  Permutation 2:" );
  ierror = perm0_check( n, p2 );

  i4vec_transpose_print ( n, p3, "  Permutation 3:" );
  ierror = perm0_check( n, p3 );

  return;
}
/******************************************************************************/

void perm0_uniform_test ( )

/******************************************************************************/
/*
  Purpose:

    PERM0_UNIFORM_TEST tests PERM0_UNIFORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 May 2015

  Author:

    John Burkardt
*/
{
  int i;
  int n = 10;
  int *p;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "PERM0_UNIFORM_TEST\n" );
  printf ( "  PERM0_UNIFORM randomly selects a permutation of 0,...,N-1.\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    p = perm0_uniform_new ( n, &seed );
    printf ( "  " );
    for ( i = 0; i < n; i++ )
    {
      printf ( "%4d", p[i] );
    }
    printf ( "\n" );
    free ( p );
  }
  return;
}
/******************************************************************************/

void perm1_check_test ( )

/******************************************************************************/
/*
  Purpose:

    PERM1_CHECK_TEST tests PERM1_CHECK.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 May 2015

  Author:

    John Burkardt
*/
{
  int ierror;
  int n = 5;
  int p1[5] = { 5, 2, 3, 4, 1 };
  int p2[5] = { 4, 1, 3, 0, 2 };
  int p3[5] = { 0, 2, 1, 3, 2 };

  printf ( "\n" );
  printf ( "PERM1_CHECK_TEST\n" );
  printf ( "  PERM1_CHECK checks a permutation of 1, ..., N.\n" );
  printf ( "\n" );

  i4vec_transpose_print ( n, p1, "  Permutation 1:" );
  ierror = perm1_check( n, p1 );

  i4vec_transpose_print ( n, p2, "  Permutation 2:" );
  ierror = perm1_check( n, p2 );

  i4vec_transpose_print ( n, p3, "  Permutation 3:" );
  ierror = perm1_check( n, p3 );

  return;
}
/******************************************************************************/

void perm1_uniform_test ( )

/******************************************************************************/
/*
  Purpose:

    PERM1_UNIFORM_TEST tests PERM1_UNIFORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 May 2015

  Author:

    John Burkardt
*/
{
  int i;
  int n = 10;
  int *p;
  int seed;
  int test;

  printf ( "\n" );
  printf ( "PERM1_UNIFORM_TEST\n" );
  printf ( "  PERM1_UNIFORM randomly selects a permutation of 1,...,N.\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    p = perm1_uniform_new ( n, &seed );
    printf ( "  " );
    for ( i = 0; i < n; i++ )
    {
      printf ( "%4d", p[i] );
    }
    printf ( "\n" );
    free ( p );
  }
  return;
}
/******************************************************************************/

void r8_abs_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ABS_TEST tests R8_ABS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  double r8;
  double r8_absolute;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed;
  int test;
  int test_num = 10;

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8_ABS_TEST\n" );
  printf ( "  R8_ABS returns the absolute value of an R8.\n" );
  printf ( "\n" );
  printf ( "      X         R8_ABS(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    r8_absolute = r8_abs ( r8 );
    printf ( "  %10.6f  %10.6f\n", r8, r8_absolute );
  }

  return;
}
/******************************************************************************/

void r8_acos_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ACOS_TEST tests R8_ACOS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2014

  Author:

    John Burkardt
*/
{
  double c;
  int test;

  printf ( "\n" );
  printf ( "R8_ACOS_TEST\n" );
  printf ( "  R8_ACOS computes the arc-cosine of an angle.\n" ); 
  printf ( "\n" );
  printf ( "          C            R8_ACOS(C)        ACOS(C)\n" );
  printf ( "\n" );

  for ( test = -1; test <= 13; test++ )
  {
    c = ( double ) ( test - 6 ) / ( double ) ( 6 );

    if ( -1.0 <= c && c <= 1.0 )
    {
      printf ( "  %14.6g  %14.6g  %14.6g\n",
        c, r8_acos ( c ), acos ( c ) );
    }
    else
    {
      printf ( "  %14.6g  %14.6g\n", c, r8_acos ( c ) );
    }
  }

  return;
}
/******************************************************************************/

void r8_acosh_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ACOSH_TEST tests R8_ACOSH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  int test;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "R8_ACOSH_TEST\n" );
  printf ( "  R8_ACOSH computes the arc-hyperbolic-cosine of an angle.\n" );
  printf ( "\n" );
  printf ( "       X            A=R8_ACOSH(X)    COSH(A)\n" );
  printf ( "\n" );

  for ( test = 0; test <= 8; test++ )
  {
    x = 1.0 + ( double ) ( test ) / 2.0;
    a = r8_acosh ( x );
    x2 = cosh ( a );
    printf ( "  %14.6g  %14.6g  %14.6g\n", x, a, x2 );
  }
  return;
}
/******************************************************************************/

void r8_asinh_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ASINH_TEST tests R8_ASINH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 May 2012

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "R8_ASINH_TEST\n" );
  printf ( "  R8_ASINH computes the inverse hyperbolic sine\n" );
  printf ( "  of a given value.\n" );
  printf ( "\n" );
  printf ( "         X     R8_ASINH(X)     SINH(R8_ASINH(X))\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = r8_asinh ( x );
    x2 = sinh ( a );

    printf ( "  %10f  %10f  %10f\n", x, a, x2 );
  }

  return;
}
/******************************************************************************/

void r8_atan_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ATAN_TEST tests R8_ATAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 8

  int test;
  double x;
  double xtest[TEST_NUM] = {
     1.0,  1.0,  0.0, -1.0,
    -1.0, -1.0,  0.0,  1.0 };
  double y;
  double ytest[TEST_NUM] = {
     0.0,  1.0,  1.0,  1.0,
     0.0, -1.0, -1.0, -1.0 };

  printf ( "\n" );
  printf ( "R8_ATAN_TEST\n" );
  printf ( "  R8_ATAN computes the arc-tangent given Y and X;\n" );
  printf ( "  ATAN2 is the system version of this routine.\n" );
  printf ( "\n" );
  printf ( "          X              Y          ATAN2(Y,X)    R8_ATAN(Y,X)\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = xtest[test];
    y = ytest[test];
    printf ( "  %14f  %14f  %14f  %14f\n", 
      x, y, atan2 ( y, x ), r8_atan ( y, x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_atanh_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ATANH_TEST tests R8_ATANH.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 April 2012

  Author:

    John Burkardt
*/
{
  double a;
  int i;
  double x;
  double x2;

  printf ( "\n" );
  printf ( "R8_ATANH_TEST\n" );
  printf ( "  R8_ATANH computes the inverse hyperbolic tangent\n" );
  printf ( "  of a given value.\n" );
  printf ( "\n" );
  printf ( "         X     R8_ATANH(X)     TANH(R8_ATANH(X))\n" );
  printf ( "\n" );

  for ( i = -2; i <= 9; i++ )
  {
    x = ( ( double ) i ) / 10.0;
    a = r8_atanh ( x );
    x2 = tanh ( a );

    printf ( "  %14g  %14g  %14g\n", x, a, x2 ); 
  }

  return;
}
/******************************************************************************/

void r8_big_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_BIG_TEST tests R8_BIG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2014

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "R8_BIG_TEST\n" );
  printf ( "  R8_BIG returns a 'big' R8 value;\n" );
  printf ( "\n" );
  printf ( "  R8_BIG = %g\n", r8_big ( ) );

  return;
}
/******************************************************************************/

void r8_cas_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CAS_TEST tests R8_CAS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 12

  int test;
  double x;

  printf ( "\n" );
  printf ( "R8_CAS_TEST\n" );
  printf ( "  R8_CAS evaluates the casine of a number.\n" );
  printf ( "\n" );
  printf ( "          X             R8_CAS ( X )\n" );
  printf ( "\n" );

  for ( test = 0; test <= TEST_NUM; test++ )
  {
    x = r8_pi ( ) * ( double ) ( test ) / ( double ) ( TEST_NUM );
    printf ( "  %14f  %14f\n", x, r8_cas ( x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_ceiling_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CEILING_TEST tests R8_CEILING.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2008

  Author:

    John Burkardt
*/
{
  int i;
  double rval;
  double rval_rounded;

  printf ( "\n" );
  printf ( "R8_CEILING_TEST\n" );
  printf ( "  R8_CEILING rounds a value up.\n" );
  printf ( "\n" );
  printf ( "         X           R8_CEILING(X)\n" );
  printf ( "\n" );
  for ( i = -6; i <= 6; i++ )
  {
    rval = ( double ) ( i ) / 5.0;
    rval_rounded = r8_ceiling ( rval );
    printf ( "  %14f  %14f\n", rval, rval_rounded );
  }

  return;
}
/******************************************************************************/

void r8_choose_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CHOOSE_TEST tests R8_CHOOSE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double cnk;
  int k;
  int n;

  printf ( "\n" );
  printf ( "R8_CHOOSE_TEST\n" );
  printf ( "  R8_CHOOSE evaluates C(N,K).\n" );
  printf ( "\n" );
  printf ( "         N         K       CNK\n" );
 
  for ( n = 0; n <= 5; n++ )
  {
    printf ( "\n" );
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );
      printf ( "  %8d  %8d  %14.6g\n", n, k, cnk );
    }
  }
 
  return;
}
/******************************************************************************/

void r8_cosd_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_COSD_TEST tests R8_COSD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 November 2013

  Author:

    John Burkardt
*/
{
  double angle;
  int i;

  printf ( "\n" );
  printf ( "R8_COSD_TEST\n" );
  printf ( "  R8_COSD computes the cosine of an angle\n" );
  printf ( "  given in degrees.\n" );
  printf ( "\n" );
  printf ( "  ANGLE    R8_COSD(ANGLE)\n" );
  printf ( "\n" );
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    printf ( "  %8.2f  %14.6g\n", angle, r8_cosd ( angle ) );
  }
 
  return;
}
/******************************************************************************/

void r8_cotd_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_COTD_TEST tests R8_COTD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2015

  Author:

    John Burkardt
*/
{
  double angle;
  int i;

  printf ( "\n" );
  printf ( "R8_COTD_TEST\n" );
  printf ( "  R8_COTD computes the cotangent of an angle\n" );
  printf ( "  given in degrees.\n" );
  printf ( "\n" );
  printf ( "  ANGLE    R8_COTD(ANGLE)\n" );
  printf ( "\n" );
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( i % 180 == 0 )
    {
      printf ( "  %8.2f    Undefined\n", angle );
    }
    else
    {
      printf ( "  %8.2f  %14.6g\n", angle, r8_cotd ( angle ) );
    }
  }
 
  return;
}
/******************************************************************************/

void r8_cscd_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CSCD_TEST tests R8_CSCD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2015

  Author:

    John Burkardt
*/
{
  double angle;
  int i;

  printf ( "\n" );
  printf ( "R8_CSCD_TEST\n" );
  printf ( "  R8_CSCD computes the cosecant of an angle\n" );
  printf ( "  given in degrees.\n" );
  printf ( "\n" );
  printf ( "  ANGLE    R8_CSCD(ANGLE)\n" );
  printf ( "\n" );
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( i % 180 == 0 )
    {
      printf ( "  %8.2f    Undefined\n", angle );
    }
    else
    {
      printf ( "  %8.2f  %14.6g\n", angle, r8_cscd ( angle ) );
    }
  }
 
  return;
}
/******************************************************************************/

void r8_cube_root_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_CUBE_ROOT_TEST tests R8_CUBE_ROOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int seed;
  double x1;
  double y;
  double x2;

  printf ( "\n" );
  printf ( "R8_CUBE_ROOT_TEST\n" );
  printf ( "  R8_CUBE_ROOT computes the cube root of an R8.\n" );
  printf ( "\n" );
  printf ( "       X               Y               Y^3\n" );
  printf ( "\n" );

  a = -10.0;
  b = +10.0;
  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x1 = r8_uniform_ab ( a, b, &seed );
    y = r8_cube_root ( x1 );
    x2 = pow ( y, 3 );
    printf ( "  %14.6g  %14.6g  %14.6g\n", x1, y, x2 );
  }

  return;
}
/******************************************************************************/

void r8_diff_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_DIFF_TEST tests R8_DIFF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 15

  int ndig = 3;
  int test;
  double x = 1.0;
  double y;
  double y_test[TEST_NUM] = {
    0.0625, 0.125, 0.25, 0.50,  0.874,
    0.876,  0.90,  0.95, 0.99,  1.0,
    1.01,   1.05,  1.10, 3.0,  10.0 };

  printf ( "\n" );
  printf ( "R8_DIFF_TEST\n" );
  printf ( "  R8_DIFF computes a difference X-Y to a given\n" );
  printf ( "    number of binary places.\n" );
  printf ( "\n" );
  printf ( "  For this test, we use %d binary places.\n", ndig );
  printf ( "\n" );
  printf ( "      X           Y       X-Y     R8_DIFF(X,Y)\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    y = y_test[test];
    printf ( "  %10f  %10f  %10f  %10f\n", 
      x, y, x - y, r8_diff ( x, y, ndig ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_digit_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_DIGIT_TEST tests R8_DIGIT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
# define MAXDIG 20

  int idigit;
  double x;

  x = r8_pi ( );

  printf ( "\n" );
  printf ( "R8_DIGIT_TEST\n" );
  printf ( "  R8_DIGIT extracts decimal digits.\n" );
  printf ( "\n" );
  printf ( "  Here, we get digits of %24.16f\n", x );
  printf ( "\n" );

  printf ( "  " );
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    printf ( " %2d", idigit );
  }
  printf ( "\n" );

  printf ( "  " );
  for ( idigit = -2; idigit <= MAXDIG; idigit++ )
  {
    printf ( "  %1d", r8_digit ( x, idigit ) );
  }
  printf ( "\n" );

  return;
# undef MAXDIG
}
/******************************************************************************/

void r8_e_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_E_TEST tests R8_E.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2015

  Author:

    John Burkardt
*/
{
  int i;
  int n;
  double value1;
  double value2;

  printf ( "\n" );
  printf ( "R8_E_TEST\n" );
  printf ( "  R8_E returns the value of E.\n" );
  printf ( "  Compare E to (1+1/n)^n\n" );
  value1 = r8_e ( );
  printf ( "  R8_E =      %g\n", value1 );
  printf ( "\n" );
  printf ( "         N     Estimate      Error\n" );
  printf ( "\n" );

  n = 1;
  for ( i = 0; i <= 20; i++ )
  {
    value2 = pow ( ( double ) ( n + 1 ) / ( double ) ( n ), n );
    printf ( "  %8d  %14.6g  %14.6g\n", n, value2, fabs ( value1 - value2 ) );
    n = n * 2;
  }

  return;
}
/******************************************************************************/

void r8_epsilon_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON_TEST tests R8_EPSILON.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt
*/
{
  double r;
  double s;
  double t;

  printf ( "\n" );
  printf ( "R8_EPSILON_TEST\n" );
  printf ( "  R8_EPSILON produces the R8 roundoff unit.\n" );
  printf ( "\n" );

  r = r8_epsilon ( );
  printf ( "  R = R8_EPSILON()         = %24.16e\n", r );

  s = 1.0 + r;
  t = s - 1.0;
  printf ( "  ( 1 + R ) - 1            = %24.16e\n", t );

  s = 1.0 + ( r / 2.0 );
  t = s - 1.0;
  printf ( "  ( 1 + (R/2) ) - 1        = %24.16e\n", t );

  return;
}
/******************************************************************************/

void r8_epsilon_compute_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON_COMPUTE_TEST tests R8_EPSILON_COMPUTE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt
*/
{
  double r;
  double s;
  double t;

  printf ( "\n" );
  printf ( "R8_EPSILON_COMPUTE_TEST\n" );
  printf ( "  R8_EPSILON_COMPUTE computes the R8 roundoff unit.\n" );
  printf ( "\n" );

  r = r8_epsilon_compute ( );
  printf ( "  R = R8_EPSILON_COMPUTE() = %24.16e\n", r );

  s = 1.0 + r;
  t = s - 1.0;
  printf ( "  ( 1 + R ) - 1            = %24.16e\n", t );

  s = 1.0 + ( r / 2.0 );
  t = s - 1.0;
  printf ( "  ( 1 + (R/2) ) - 1        = %24.16e\n", t );

  return;
}
/******************************************************************************/

void r8_factorial_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL_TEST tests R8_FACTORIAL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 December 2014

  Author:

    John Burkardt
*/
{
  double f1;
  double f2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "R8_FACTORIAL_TEST\n" );
  printf ( "  R8_FACTORIAL evaluates the factorial.\n" );
  printf ( "\n" );
  printf ( "    N                Exact                  Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_factorial ( n );

    printf ( "  %4d  %24.16g  %24.16g\n", n, f1, f2 );

  }
 
  return;
}
/******************************************************************************/

void r8_factorial2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL2_TEST tests R8_FACTORIAL2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2015

  Author:

    John Burkardt
*/
{
  double f1;
  double f2;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "R8_FACTORIAL2_TEST\n" );
  printf ( "  R8_FACTORIAL2 evaluates the double factorial.\n" );
  printf ( "\n" );
  printf ( "    N                Exact                  Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial2_values ( &n_data, &n, &f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_factorial2 ( n );

    printf ( "  %4d  %24.16g  %24.16g\n", n, f1, f2 );

  }
 
  return;
}
/******************************************************************************/

void r8_fall_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FALL_TEST tests R8_FALL.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 December 2014

  Author:

    John Burkardt
*/
{
  double f1;
  double f2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_FALL_TEST\n" );
  printf ( "  R8_FALL evaluates the falling factorial Fall(X,N).\n" );
  printf ( "\n" );
  printf ( "    X          N                Exact                  Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_fall_values ( &n_data, &x, &n, &f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_fall ( x, n );

    printf ( "  %8f  %4d  %24.16g  %24.16g\n", x, n, f1, f2 );

  }
 
  return;
}
/******************************************************************************/

void r8_fractional_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_FRACTIONAL_TEST tests R8_FRACTIONAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt
*/
{
  double fractional;
  double r8;
  double r8_hi = 5.0;
  double r8_lo = -3.0;
  int seed = 123456789;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "R8_FRACTIONAL_TEST\n" );
  printf ( "  R8_FRACTIONAL returns the fractional part of an R8.\n" );
  printf ( "\n" );
  printf ( "      X           R8_FRACTIONAL(X)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    r8 = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    fractional = r8_fractional ( r8 );
    printf ( "  %10f  %10f\n", r8, fractional );
  }

  return;
}
/******************************************************************************/

void r8_gamma_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA_TEST tests R8_GAMMA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2013

  Author:

    John Burkardt
*/
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " R8_GAMMA_TEST:\n" );
  printf ( "   R8_GAMMA evaluates the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X                GAMMA(X)         R8_GAMMA(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = r8_gamma ( x );

    printf ( "  %24.16e  %24.16e  %24.16e\n", x, fx1, fx2 );
  }
  return;
}
/******************************************************************************/

void r8_gamma_log_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA_LOG_TEST tests R8_GAMMA_LOG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2013

  Author:

    John Burkardt
*/
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " R8_GAMMA_LOG_TEST:\n" );
  printf ( "   R8_GAMMA_LOG_VALUES evaluates\n" );
  printf ( "   the logarithm of the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X            GAMMA_LOG(X)      R8_GAMMA_LOG(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = r8_gamma_log ( x );

    printf ( "  %24.16e  %24.16e  %24.16e\n", x, fx1, fx2 );
  }
  return;
}
/******************************************************************************/

void r8_huge_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE_TEST tests R8_HUGE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
  printf ( "\n" );
  printf ( "R8_HUGE_TEST\n" );
  printf ( "  R8_HUGE returns a large R8 value;\n" );
  printf ( "\n" );
  printf ( "  R8_HUGE =   %e\n", r8_huge ( ) );

  return;
}
/******************************************************************************/

void r8_log_2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_LOG_2_TEST tests R8_LOG_2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 18

  int test;
  double x;
  double x_test[TEST_NUM] = {
    0.0,  1.0,  2.0,   3.0,  9.0,
   10.0, 11.0, 99.0, 101.0, -1.0,
   -2.0, -3.0, -9.0,   0.5,  0.33,
    0.25, 0.20, 0.01 };

  printf ( "\n" );
  printf ( "R8_LOG_2_TEST\n" );
  printf ( "  R8_LOG_2: computes the logarithm base 2.\n" );
  printf ( "\n" );
  printf ( "        X       R8_LOG_2\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = x_test[test];
    printf ( "  %12f  %12g\n", x, r8_log_2 ( x ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_log_b_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_LOG_B_TEST tests R8_LOG_B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 10

  double b;
  double b_test[TEST_NUM] = {
    2.0, 3.0, 4.0, 5.0, 6.0,
    7.0, 8.0, 16.0, 32.0, 256.0 };
  int test;
  double x;

  x = 16.0;

  printf ( "\n" );
  printf ( "R8_LOG_B_TEST\n" );
  printf ( "  R8_LOG_B computes the logarithm base B.\n" );
  printf ( "\n" );
  printf ( "        X     B     R8_LOG_B\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    b = b_test[test];

    printf ( "  %12f  %12f  %12f\n", x, b, r8_log_b ( x, b ) );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_mant_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MANT_TEST tests R8_MANT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  int is;
  int l;
  double r;
  double x;

  x = -314.159;

  printf ( "\n" );
  printf ( "R8_MANT_TEST\n" );
  printf ( "  R8_MANT decomposes a value.\n" );
  printf ( "\n" );
  printf ( "  Number to be decomposed: X = %f\n", x );

  r8_mant ( x, &is, &r, &l );

  printf ( "\n" );
  printf ( "  X = %d * %f * 2 ^ %d\n", is, r, l );

  return;
}
/******************************************************************************/

void r8_max_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MAX_TEST tests R8_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  int i;
  double r8_hi;
  double r8_lo;
  int seed;

  printf ( "\n" );
  printf ( "R8_MAX_TEST\n" );
  printf ( "  R8_MAX returns the maximum of two R8's.\n" );
  printf ( "\n" );
  printf ( "       A       B      C=R8_MAX(A,B)\n" );
  printf ( "\n" );

  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;

  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    b = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    c = r8_max ( a, b );
    printf ( "  %8.4f  %8.4f  %8.4f\n", a, b, c );
  }

  return;
# undef N
}
/******************************************************************************/

void r8_min_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MIN_TEST tests R8_MIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2015

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  int i;
  double r8_hi;
  double r8_lo;
  int seed;

  printf ( "\n" );
  printf ( "R8_MIN_TEST\n" );
  printf ( "  R8_MIN returns the minimum of two R8's.\n" );
  printf ( "\n" );
  printf ( "       A       B      C=R8_MIN(A,B)\n" );
  printf ( "\n" );

  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;

  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    b = r8_uniform_ab ( r8_lo, r8_hi, &seed );
    c = r8_min ( a, b );
    printf ( "  %8.4f  %8.4f  %8.4f\n", a, b, c );
  }

  return;
# undef N
}
/******************************************************************************/

void r8_mod_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MOD_TEST tests R8_MOD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;
  double x_hi = 10.0;
  double x_lo = -10.0;
  double y;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "R8_MOD_TEST\n" );
  printf ( "  R8_MOD returns the remainder after division.\n" );
  printf ( "  R8_MOD ( X, Y ) has the same sign as X.\n" );
  printf ( "\n" );
  printf ( "      X         Y    FMOD(X,Y)    R8_MOD(X,Y)\n" );
  printf ( "\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( x_lo, x_hi, &seed );
    y = r8_uniform_ab ( x_lo, x_hi, &seed );

    z1 =   fmod ( x, y );
    z2 = r8_mod ( x, y );

    printf ( "  %12f  %12f  %12f  %12f\n", x, y, z1, z2 );
  }

  return;
}
/******************************************************************************/

void r8_modp_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MODP_TEST tests R8_MODP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
# define TEST_NUM 20

  int seed = 123456789;
  int test;
  double x;
  double z1;
  double z2;

  printf ( "\n" );
  printf ( "R8_MODP_TEST\n" );
  printf ( "  R8_MODP returns the remainder after division.\n" );
  printf ( "  R8_MODP ( X, Y ) is positive if Y is.\n" );
  printf ( "\n" );
  printf ( "      X             FMOD(x,1.0)  R8_MODP(x,1.0)\n" );
  printf ( "\n" );

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    x = r8_uniform_01 ( &seed );
    x = 20.0 * ( x - 0.24 );
    z1 = fmod ( x, 1.0 );
    z2 = r8_modp ( x, 1.0 );
    printf ( "  %12f  %12f  %12f\n", x, z1, z2 );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_mop_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_MOP_TEST tests R8_MOP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2014

  Author:

    John Burkardt
*/
{
  int i4;
  int i4_max;
  int i4_min;
  double r8;
  int seed = 123456789;
  int test;

  printf ( "\n" );
  printf ( "R8_MOP_TEST\n" );
  printf ( "  R8_MOP evaluates (-1.0)^I4 as an R8.\n" );
  printf ( "\n" );
  printf ( "    I4  R8_MOP(I4)\n" );
  printf ( "\n" );

  i4_min = -100;
  i4_max = +100;

  for ( test = 1; test <= 10; test++ )
  {
    i4 = i4_uniform_ab ( i4_min, i4_max, &seed );
    r8 = r8_mop ( i4 );
    printf ( "  %4d  %4.1f\n", i4, r8 );
  }

  return;
}
/******************************************************************************/

void r8_nint_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_NINT_TEST tests R8_NINT

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double b;
  double c;
  int seed = 123456789;
  int test;
  int test_num = 10;
  double x;

  printf ( "\n" );
  printf ( "R8_NINT_TEST\n" );
  printf ( "  R8_NINT produces the nearest integer.\n" );
  printf ( "\n" );
  printf ( "       X        R8_NINT(X)\n" );
  printf ( "\n" );

  b = -10.0;
  c = +10.0;

  for ( test = 1; test <= test_num; test++ )
  {
    x = r8_uniform_ab ( b, c, &seed );
    printf ( "   %10f  %6d\n", x, r8_nint ( x ) );
  }

  return;
}
/******************************************************************************/

void r8_normal_01_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_NORMAL_01_TEST tests R8_NORMAL_01.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2010

  Author:

    John Burkardt
*/
{
# define TEST_NUM 20

  int seed = 123456789;
  int test;
  double x;

  printf ( "\n" );
  printf ( "R8_NORMAL_01_TEST\n" );
  printf ( "  R8_NORMAL_01 generates normally distributed random values.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    x = r8_normal_01 ( &seed );
    printf ( "  %10f\n", x );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8_pi_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_PI_TEST tests R8_PI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2010

  Author:

    John Burkardt
*/
{
  double four;
  double one;
  double v1;
  double v2;

  four = ( double ) ( 4 );
  one = ( double ) ( 1 );

  printf ( "\n" );
  printf ( "R8_PI_TEST\n" );
  printf ( "  R8_PI returns the value of PI.\n" );
  printf ( "\n" );
  v1 = r8_pi ( );
  printf ( "  R8_PI =     %24.16f\n", v1 );
  v2 = four * atan ( one );
  printf ( "  4*atan(1) = %24.16f\n", v2 );

  return;
}
/******************************************************************************/

void r8_power_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_POWER_TEST tests R8_POWER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2010

  Author:

    John Burkardt
*/
{
  int p;
  double r;
  double value;

  printf ( "\n" );
  printf ( "R8_POWER_TEST\n" );
  printf ( "  R8_POWER computes R^P\n" );
  printf ( "\n" );
  printf ( "      R          P       R^P\n" );
  printf ( "\n" );

  for ( p = -5; p <= 5; p++ )
  {
    r = 2.0;
    value = r8_power ( r, p );
    printf ( "  %12f  %6d  %12f\n", r, p, value );
  }

  return;
}
/******************************************************************************/

void r8_power_fast_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_POWER_FAST_TEST tests R8_POWER_FAST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt
*/
{
  int i;
  int mults;
  int p;
  double r;
  double rp;

  printf ( "\n" );
  printf ( "R8_POWER_FAST_TEST\n" );
  printf ( "  R8_POWER_FAST computes R^P, economizing on\n" );
  printf ( "    multiplications.\n" );
  printf ( "\n" );
  printf ( "      R          P       R^P        Mults\n" );
  printf ( "\n" );

  for ( i = -10; i <= 40; i++ )
  {
    r = 2.0;
    p = i;
    rp = r8_power_fast ( r, p, &mults );
    printf ( "  %12f  %6d  %12g  %6d\n", r, p, rp, mults );
  }

  return;
}
/******************************************************************************/

void r8_rise_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_RISE_TEST tests R8_RISE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 December 2014

  Author:

    John Burkardt
*/
{
  double f1;
  double f2;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "R8_RISE_TEST\n" );
  printf ( "  R8_RISE evaluates the rising factorial Fall(X,N).\n" );
  printf ( "\n" );
  printf ( "    X          N                Exact                  Computed\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_rise_values ( &n_data, &x, &n, &f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_rise ( x, n );

    printf ( "  %8f  %4d  %24.16g  %24.16g\n", x, n, f1, f2 );

  }
 
  return;
}
/******************************************************************************/

void r8_round2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ROUND2_TEST tests R8_ROUND2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt
*/
{
  int i;
  int nplace;
  double x;
  double xround;

  x = r8_pi ( );

  printf ( "\n" );
  printf ( "R8_ROUND2_TEST\n" );
  printf ( "  R8_ROUND2 rounds a number to a\n" );
  printf ( "  specified number of base 2 digits.\n" );
  printf ( "\n" );
  printf ( "  Test effect on PI:\n" );
  printf ( "  X = %24.16f\n", x );
  printf ( "\n" );
  printf ( "  NPLACE      XROUND\n" );
  printf ( "\n" );

  for ( i = 0; i <= 20; i++ )
  {
    nplace = i;
    xround = r8_round2 ( nplace, x );
    printf ( "  %8d  %24.16f\n", i, xround );
  }

  return;
}
/******************************************************************************/

void r8_roundb_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ROUNDB_TEST tests R8_ROUNDB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt
*/
{
  int base;
  int i;
  int nplace;
  double x;
  double xround;

  base = 3;
  x = r8_pi ( );

  printf ( "\n" );
  printf ( "R8_ROUNDB_TEST\n" );
  printf ( "  R8_ROUNDB rounds a number to a \n" );
  printf ( "  specified number of base BASE digits.\n" );
  printf ( "\n" );
  printf ( "  Here, we will use BASE = %2d\n", base );
  printf ( "\n" );
  printf ( "  Test effect on PI:\n" );
  printf ( "  X = %24.16f\n", x );
  printf ( "\n" );
  printf ( "  NPLACE  XROUND\n" );
  printf ( "\n" );

  for ( i = 0; i <= 20; i++ )
  {
    nplace = i;
    xround = r8_roundb ( base, nplace, x );
    printf ( "  %8d  %24.16f\n", i, xround );
  }

  printf ( "\n" );
  printf ( "  Try with a negative base:\n" );
  x = 121.0;
  base = -3;
  nplace = 3;
  printf ( "\n" );
  printf ( "  Input quantity is X = %24.16f\n", x );
  printf ( "  to be rounded in base %d\n", base );

  for ( nplace = 1; nplace <= 5; nplace++ )
  {
    xround = r8_roundb ( base, nplace, x );

    printf ( "\n" );
    printf ( "  Output value to %d places is %f\n", nplace, xround );
  }

  return;
}
/******************************************************************************/

void r8_roundx_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_ROUNDX_TEST tests R8_ROUNDX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2011

  Author:

    John Burkardt
*/
{
  int i;
  int nplace;
  int seed;
  double x;
  double xround;

  seed = 123456789;
  x = r8_pi ( );

  printf ( "\n" );
  printf ( "R8_ROUNDX_TEST\n" );
  printf ( "  R8_ROUNDX rounds a number to a \n" );
  printf ( "  specified number of decimal digits.\n" );
  printf ( "\n" );
  printf ( "  Test effect on PI:\n" );
  printf ( "  X = %f\n", x );
  printf ( "\n" );
  printf ( "  NPLACE  XROUND\n" );
  printf ( "\n" );

  for ( i = 0; i <= 10; i++ )
  {
    nplace = i;
    xround = r8_roundx ( nplace, x );
    printf ( "  %6d  %20.12f\n", i, xround );
  }

  printf ( "\n" );
  printf ( "  Test effect on random values:\n" );
  printf ( "\n" );
  printf ( "  NPLACE  X     XROUND\n" );
  printf ( "\n" );

  for ( i = 1; i <= 5; i++ )
  {
    x = r8_uniform_01 ( &seed );

    printf ( "\n" );

    for ( nplace = 0; nplace <= 10; nplace = nplace + 2 )
    {
      xround = r8_roundx ( nplace, x );

      printf ( "  %6d  %16f  %20f\n", nplace, x, xround );
    }
  }

  return;
}
/******************************************************************************/

void r8_secd_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_SECD_TEST tests R8_SECD.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2015

  Author:

    John Burkardt
*/
{
  double angle;
  int i;

  printf ( "\n" );
  printf ( "R8_SECD_TEST\n" );
  printf ( "  R8_SECD computes the secant of an angle\n" );
  printf ( "  given in degrees.\n" );
  printf ( "\n" );
  printf ( "  ANGLE    R8_SECD(ANGLE)\n" );
  printf ( "\n" );
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( ( i + 90 ) % 180 == 0 )
    {
      printf ( "  %8.2f    Undefined\n", angle );
    }
    else
    {
      printf ( "  %8.2f  %14.6g\n", angle, r8_secd ( angle ) );
    }
  }
 
  return;
}
/******************************************************************************/

void r8_sign_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN_TEST tests R8_SIGN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  double r8;
  double r8_test[5] = { -1.25, -0.25, 0.0, +0.5, +9.0 };
  double s;
  int test;
  const int test_num = 5;

  printf ( "\n" );
  printf ( "R8_SIGN_TEST\n" );
  printf ( "  R8_SIGN returns the sign of a number.\n" );
  printf ( "\n" );
  printf ( "    X        R8_SIGN(X)\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    r8 = r8_test[test];
    s = r8_sign ( r8 );
    printf ( "  %8.4f  %8.0f\n", r8, s );
  }

  return;
}
/******************************************************************************/

void r8_sign3_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN3_TEST tests R8_SIGN3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 September 2014

  Author:

    John Burkardt
*/
{
  double r8;
  double r8_test[5] = { -1.25, -0.25, 0.0, +0.5, +9.0 };
  double s;
  int test;
  const int test_num = 5;

  printf ( "\n" );
  printf ( "R8_SIGN3_TEST\n" );
  printf ( "  R8_SIGN3 returns the three way sign of a number.\n" );
  printf ( "\n" );
  printf ( "    X        R8_SIGN3(X)\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    r8 = r8_test[test];
    s = r8_sign3 ( r8 );
    printf ( "  %8.4f  %8.0f\n", r8, s );
  }

  return;
}
/******************************************************************************/

void r8_sind_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_SIND_TEST tests R8_SIND.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2015

  Author:

    John Burkardt
*/
{
  double angle;
  int i;

  printf ( "\n" );
  printf ( "R8_SIND_TEST\n" );
  printf ( "  R8_SIND computes the sine of an angle\n" );
  printf ( "  given in degrees.\n" );
  printf ( "\n" );
  printf ( "  ANGLE    R8_SIND(ANGLE)\n" );
  printf ( "\n" );
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    printf ( "  %8.2f  %14.6g\n", angle, r8_sind ( angle ) );
  }
 
  return;
}
/******************************************************************************/

void r8_swap_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_SWAP_TEST tests R8_SWAP.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 June 2011

  Author:

    John Burkardt
*/
{
  double x;
  double y;

  printf ( "\n" );
  printf ( "R8_SWAP_TEST\n" );
  printf ( "  R8_SWAP swaps two reals.\n" );

  x = 1.0;
  y = 3.14159;

  printf ( "\n" );
  printf ( "  Before swapping: \n" );
  printf ( "\n" );
  printf ( "    X = %f\n", x );
  printf ( "    Y = %f\n", y );

  r8_swap ( &x, &y );

  printf ( "\n" );
  printf ( "  After swapping: \n" );
  printf ( "\n" );
  printf ( "    X = %f\n", x );
  printf ( "    Y = %f\n", y );

  return;
}
/******************************************************************************/

void r8_swap3_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_SWAP3_TEST tests R8_SWAP3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2013

  Author:

    John Burkardt
*/
{
  int i;
  double x;
  double y;
  double z;

  printf ( "\n" );
  printf ( "R8_SWAP3_TEST\n" );
  printf ( "  R8_SWAP3 swaps three reals.\n" );

  x = 1.0;
  y = 3.14159;
  z = 1952.0;

  printf ( "\n" );
  printf ( "              X       Y       Z\n" );
  printf ( "\n" );
  printf ( "  Start: %g  %g  %g\n", x, y, z );

  for ( i = 1; i <= 3; i++ )
  {
    r8_swap3 ( &x, &y, &z );
    printf ( "  Swap %d  %g  %g  %g\n", i, x, y, z );
  }

  return;
}
/******************************************************************************/

void r8_tand_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_TAND_TEST tests R8_TAND.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 January 2015

  Author:

    John Burkardt
*/
{
  double angle;
  int i;

  printf ( "\n" );
  printf ( "R8_TAND_TEST\n" );
  printf ( "  R8_TAND computes the tangent of an angle\n" );
  printf ( "  given in degrees.\n" );
  printf ( "\n" );
  printf ( "  ANGLE    R8_TAND(ANGLE)\n" );
  printf ( "\n" );
 
  for ( i = 0; i <= 360; i = i + 15 )
  {
    angle = ( double ) ( i );
    if ( ( i + 90 ) % 180 == 0 )
    {
      printf ( "  %8.2f    Undefined\n", angle );
    }
    else
    {
      printf ( "  %8.2f  %14.6g\n", angle, r8_tand ( angle ) );
    }
  }
 
  return;
}
/******************************************************************************/

void r8_to_r8_discrete_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_TO_R8_DISCRETE_TEST tests R8_TO_R8_DISCRETE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 June 2011

  Author:

    John Burkardt
*/
{
  int ndx = 19;
  double r;
  double rd;
  double rhi = 10.0;
  double rhi2;
  double rlo = 1.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 15;

  printf ( "\n" );
  printf ( "R8_TO_R8_DISCRETE_TEST\n" );
  printf ( "  R8_TO_R8_DISCRETE maps numbers to a discrete set\n" );
  printf ( "  of equally spaced numbers in an interval.\n" );
  printf ( "\n" );
  printf ( "  Number of discrete values = %d\n", ndx );
  printf ( "  Real interval: [%f, %f]\n", rlo, rhi );
  printf ( "\n" );
  printf ( "      R         RD\n" );
  printf ( "\n" );

  seed = 123456789;

  rlo2 = rlo - 2.0;
  rhi2 = rhi + 2.0;

  for ( test = 0; test < test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, &seed );
    rd = r8_to_r8_discrete ( r, rlo, rhi, ndx );
    printf ( "  %14f  %14f\n", r, rd );
  }

  return;
}
/******************************************************************************/

void r8_to_i4_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_TO_I4_TEST tests R8_TO_I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2014

  Author:

    John Burkardt
*/
{
  int ix;
  int ixmax;
  int ixmin;
  double x;
  double xmax;
  double xmin;

  printf ( "\n" );
  printf ( "R8_TO_I4_TEST\n" );
  printf ( "  R8_TO_I4 finds an integer IX in [IXMIN,IXMAX]\n" );
  printf ( "  corresponding to X in [XMIN,XMAX].\n" );

  xmin = 2.5;
  x = 3.5;
  xmax = 5.5;

  ixmin = 10;
  ixmax = 40;

  ix = r8_to_i4 ( xmin, xmax, x, ixmin, ixmax );

  printf ( "\n" );
  printf ( "   XMIN %f,   X = %f,  XMAX = %f\n", xmin, x, xmax );
  printf ( "  IXMIN %d,  IX = %d, IXMAX = %d\n", ixmin, ix, ixmax );

  return;
}
/******************************************************************************/

void r8_uniform_01_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
# define N 1000

  int i;
  double max;
  double mean;
  double min;
  int n;
  int seed = 123456789;
  double x[N];
  double variance;

  printf ( "\n" );
  printf ( "R8_UNIFORM_01_TEST\n" );
  printf ( "  R8_UNIFORM_01 samples a uniform random distribution in [0,1].\n" );
  printf ( "  distributed random numbers.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  for ( i = 0; i < N; i++ )
  {
    x[i] = r8_uniform_01 ( &seed );
  }

  printf ( "\n" );
  printf ( "  First few values:\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    printf ( "  %6d  %14f\n", i, x[i] );
  }
  min = r8vec_min ( N, x );
  max = r8vec_max ( N, x );
  mean = r8vec_mean ( N, x );
  variance = r8vec_variance ( N, x );

  printf ( "\n" );
  printf ( "  Number of samples was %d\n", N );
  printf ( "  Minimum value was %f\n", min );
  printf ( "  Maximum value was %f\n", max );
  printf ( "  Average value was %f\n", mean );
  printf ( "  Variance was      %f\n", variance );

  return;
# undef N
}
/******************************************************************************/

void r8_uniform_ab_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_AB_TEST tests R8_UNIFORM_AB.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  int i;
  int seed;

  b = 10.0;
  c = 25.0;
  seed = 17;

  printf ( "\n" );
  printf ( "R8_UNIFORM_AB_TEST\n" );
  printf ( "  R8_UNIFORM_AB produces a random real in a given range.\n" );
  printf ( "\n" );
  printf ( "  Using range %f <= A <= %f.\n", b, c );
  printf ( "\n" );

  printf ( "\n" );
  printf ( "     I      A\n" );
  printf ( "\n" );
  for ( i = 0; i < 10; i++ )
  {
    a = r8_uniform_ab ( b, c, &seed );
    printf ( "  %4d  %10f\n", i, a );
  }

  return;
}
/******************************************************************************/

void r8_walsh_1d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_WALSH_1D_TEST tests R8_WALSH_1D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
  int i;
  double w0;
  double wm1;
  double wm2;
  double wm3;
  double wp1;
  double wp2;
  double x;

  printf ( "\n" );
  printf ( "R8_WALSH_1D_TEST\n" );
  printf ( "  R8_WALSH_1D evaluates 1D Walsh functions:\n" );
  printf ( "\n" );
  printf ( "  X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)\n" );
  printf ( "\n" );

  for ( i = 0; i <= 32; i++ )
  {
    x = ( double ) ( i ) / 4.0;

    wp2 = r8_walsh_1d ( x,  2 );
    wp1 = r8_walsh_1d ( x,  1 );
    w0  = r8_walsh_1d ( x,  0 );
    wm1 = r8_walsh_1d ( x, -1 );
    wm2 = r8_walsh_1d ( x, -2 );
    wm3 = r8_walsh_1d ( x, -3 );

    printf ( "  %10g  %2g  %2g  %2g  %2g  %2g  %2g\n",
      x, wp2, wp1, w0, wm1, wm2, wm3 );
  }

  return;
}
/******************************************************************************/

void r8_wrap_test ( )

/******************************************************************************/
/*
  Purpose:

    R8_WRAP_TEST tests R8_WRAP;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2011

  Author:

    John Burkardt
*/
{
  double a = - 2.0;
  double b = 12.0;
  double r;
  double r2;
  double rhi = 6.5;
  double rlo = 3.0;
  int seed;
  int test;
  int test_num = 20;

  printf ( "\n" );
  printf ( "R8_WRAP_TEST\n" );
  printf ( "  R8_WRAP \"wraps\" an R8 to lie within an interval:\n" );
  printf ( "\n" );
  printf ( "  Wrapping interval is %g, %g\n", rlo, rhi );
  printf ( "\n" );
  printf ( "      R      R8_WRAP ( R )\n" );
  printf ( "\n" );
  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( a, b, &seed );
    r2 = r8_wrap ( r, rlo, rhi );
    printf ( "  %14g  %14g\n", r, r2 );
  }

  return;
}
/******************************************************************************/

void r82col_print_part_test ( )

/******************************************************************************/
/*
  Purpose:

    R82COL_PRINT_PART_TEST tests R82COL_PRINT_PART.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 April 2015

  Author:

    John Burkardt
*/
{
  int max_print;
  int n = 10;
  double v[10*2] = {
    11.0,  21.0, 31.0, 41.0, 51.0, 61.0, 71.0, 81.0, 91.0, 101.0, 
    12.0,  22.0, 32.0, 42.0, 52.0, 62.0, 72.0, 82.0, 92.0, 102.0 };

  printf ( "\n" );
  printf ( "R82COL_PRINT_PART_TEST\n" );
  printf ( "  R82COL_PRINT_PART prints part of an R82COL.\n" );

  max_print = 2;
  r82col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r82col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r82col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
/******************************************************************************/

void r82poly2_type_test ( )

/******************************************************************************/
/*
  Purpose:

    R82POLY2_TYPE_TEST tests R82POLY2_TYPE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 12

  double a;
  double a_test[TEST_NUM] = {
    9.0, 4.0, 9.0,  1.0, 0.0,
    1.0, 0.0, 0.0,  0.0, 0.0,
    0.0, 0.0 };
  double b;
  double b_test[TEST_NUM] = {
    -4.0, 1.0,  16.0,  1.0,  0.0,
     2.0, 1.0,   1.0,  1.0,  0.0,
     0.0, 0.0 };
  double c;
  double c_test[TEST_NUM] = {
     0.0, -4.0,   0.0,   0.0, 1.0,
     0.0,  0.0,   0.0,  0.0,  0.0,
     0.0,  0.0 };
  double d;
  double r8_test[TEST_NUM] = {
    -36.0,  3.0,  36.0,  -6.0, 3.0,
    -2.0,   0.0,   0.0,  0.0,  2.0,
     0.0, 0.0 };
  double e;
  double e_test[TEST_NUM] = {
    -24.0, -4.0, -32.0, -10.0, -1.0,
     16.0, -6.0, -6.0, -2.0, -1.0,
     0.0, 0.0 };
  double f;
  double f_test[TEST_NUM] = {
    -36.0,  1.0, -92.0, 115.0, -3.0,
     33.0, +8.0, 10.0,  +1.0,  1.0,
      0.0, 1.0 };
  int test;
  int type;

  printf ( "\n" );
  printf ( "R82POLY2_TYPE_TEST\n" );
  printf ( "  R82POLY2_TYPE determines the type of a second order\n" );
  printf ( "  equation in two variables.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    a = a_test[test];
    b = b_test[test];
    c = c_test[test];
    d = r8_test[test];
    e = e_test[test];
    f = f_test[test];

    printf ( "\n" );

    r82poly2_print ( a, b, c, d, e, f );

    type = r82poly2_type ( a, b, c, d, e, f );

    printf ( "  Type = %d\n", type );

    r82poly2_type_print ( type );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r82row_order_type_test ( )

/******************************************************************************/
/*
  Purpose:

    R82ROW_ORDER_TYPE_TEST tests R82ROW_ORDER_TYPE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
# define N 4
# define TEST_NUM 10

  int i;
  int j;
  int order;
  int seed = 123456789;
  int test;
  double *x;

  printf ( "\n" );
  printf ( "R82ROW_ORDER_TYPE_TEST\n" );
  printf ( "  R82ROW_ORDER_TYPE classifies an R8VEC as\n" );
  printf ( "  -1: no order\n" );
  printf ( "   0: all equal;\n" );
  printf ( "   1: ascending;\n" );
  printf ( "   2: strictly ascending;\n" );
  printf ( "   3: descending;\n" );
  printf ( "   4: strictly descending.\n" );
  printf ( "\n" );

  for ( test = 1; test <= TEST_NUM; test++ )
  {
    x = r8mat_uniform_01_new ( 2, N, &seed );

    for ( j = 0; j < N; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        x[i+j*2] = ( double ) ( r8_nint ( 3.0 * x[i+j*2] ) );
      }
    }
    order = r82row_order_type ( N, x );

    printf ( "  Order type = %d\n", order );

    r82row_print ( N, x, " " );

    free ( x );
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void r82row_part_quick_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R82ROW_PART_QUICK_A_TEST tests R82ROW_PART_QUICK_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
# define N 12

  double *a;
  double b = 0.0E+00;
  double c = 2.0E+00;
  int i;
  int l;
  int r;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R82ROW_PART_QUICK_A_TEST\n" );
  printf ( "  R82ROW_PART_QUICK_A reorders an R82ROW\n" );
  printf ( "  as part of a quick sort.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  a = r8mat_uniform_ab_new ( 2, N, b, c, &seed );

  r82row_print ( N, a, "  Before rearrangment:" );

  r82row_part_quick_a ( N, a, &l, &r );

  printf ( "\n" );
  printf ( "  Rearranged array\n" );
  printf ( "  Left index =  %d\n", l );
  printf ( "  Key index =   %d\n", l + 1 );
  printf ( "  Right index = %d\n", r );

  r82row_print ( l,     a,         "  Left half:" );
  r82row_print ( 1,     a+2*l,     "  Key:" );
  r82row_print ( N-l-1, a+2*(l+1), "  Right half:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r82row_print_part_test ( )

/******************************************************************************/
/*
  Purpose:

    R82ROW_PRINT_PART_TEST tests R82ROW_PRINT_PART.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2015

  Author:

    John Burkardt
*/
{
  int max_print;
  int n = 10;
  double v[2*10] = {
     11.0,  21.0, 
     12.0,  22.0, 
     13.0,  23.0, 
     14.0,  24.0, 
     15.0,  25.0, 
     16.0,  26.0, 
     17.0,  27.0, 
     18.0,  28.0, 
     19.0,  29.0, 
     20.0,  30.0 };

  printf ( "\n" );
  printf ( "R82ROW_PRINT_PART_TEST\n" );
  printf ( "  R82ROW_PRINT_PART prints part of an R82ROW.\n" );

  max_print = 2;
  r82row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r82row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r82row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
/******************************************************************************/

void r82row_sort_heap_index_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R82ROW_SORT_HEAP_INDEX_A_TEST tests R82ROW_SORT_HEAP_INDEX_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2012

  Author:

    John Burkardt
*/
{
# define N 12

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int *indx;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R82ROW_SORT_HEAP_INDEX_A_TEST\n" );
  printf ( "  R82ROW_SORT_HEAP_INDEX_A index sorts an R82ROW\n" );
  printf ( "  using heapsort.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  a = r8mat_uniform_ab_new ( 2, N, b, c, &seed );
/*
  Give a few elements the same first component.
*/
  a[0+2*2] = a[0+4*2];
  a[0+3*2] = a[0+11*2];
/*
  Give a few elements the same second component.
*/
  a[1+5*2] = a[1+0*2];
  a[1+1*2] = a[1+8*2];
/*
  Make two entries equal.
*/
  a[0+6*2] = a[0+10*2];
  a[1+6*2] = a[1+10*2];

  r82row_print ( N, a, "  Before rearrangement:" );

  indx = r82row_sort_heap_index_a ( N, a );

  printf ( "\n" );
  printf ( "         I     Index   A(Index)\n" );
  printf ( "\n" );

  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %8d  %12f  %12f\n", 
      i, indx[i], a[0+indx[i]*2], a[1+indx[i]*2] );
  }

  r82row_permute ( N, indx, a );

  r82row_print ( N, a, "  After rearrangement by R82ROW_PERMUTE:" );

  free ( a );
  free ( indx );

  return;
# undef N
}
/******************************************************************************/

void r82row_sort_quick_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R82ROW_SORT_QUICK_A_TEST tests R82ROW_SORT_QUICK_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2005

  Author:

    John Burkardt
*/
{
# define N 12

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R82ROW_SORT_QUICK_A_TEST\n" );
  printf ( "  R82ROW_SORT_QUICK_A sorts an R82ROW\n" );
  printf ( "  as part of a quick sort.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  a = r8mat_uniform_ab_new ( 2, N, b, c, &seed );
/*
  For better testing, give a few elements the same first component.
*/
  a[2*(3-1)+0] = a[2*(5-1)+0];
  a[2*(4-1)+0] = a[2*(12-1)+0];
/*
  Make two entries equal.
*/
  a[2*(7-1)+0] = a[2*(11-1)+0];
  a[2*(7-1)+1] = a[2*(11-1)+1];

  r82row_print ( N, a, "  Before sorting:" );

  r82row_sort_quick_a ( N, a );

  r82row_print ( N, a, "  Sorted array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r83col_print_part_test ( )

/******************************************************************************/
/*
  Purpose:

    R83COL_PRINT_PART_TEST tests R83COL_PRINT_PART.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 April 2015

  Author:

    John Burkardt
*/
{
  int max_print;
  int n = 10;
  double v[10*3] = {
    11.0,  21.0, 31.0, 41.0, 51.0, 61.0, 71.0, 81.0, 91.0, 101.0, 
    12.0,  22.0, 32.0, 42.0, 52.0, 62.0, 72.0, 82.0, 92.0, 102.0,
    13.0,  23.0, 33.0, 43.0, 53.0, 63.0, 73.0, 83.0, 93.0, 103.0 };

  printf ( "\n" );
  printf ( "R83COL_PRINT_PART_TEST\n" );
  printf ( "  R83COL_PRINT_PART prints part of an R83COL.\n" );

  max_print = 2;
  r83col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r83col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r83col_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
/******************************************************************************/

void r83row_print_part_test ( )

/******************************************************************************/
/*
  Purpose:

    R83ROW_PRINT_PART_TEST tests R83ROW_PRINT_PART.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 April 2015

  Author:

    John Burkardt
*/
{
  int max_print;
  int n = 10;
  double v[3*10] = {
     11.0,  21.0,  31.0, 
     12.0,  22.0,  32.0,
     13.0,  23.0,  33.0,
     14.0,  24.0,  34.0,
     15.0,  25.0,  35.0,
     16.0,  26.0,  36.0,
     17.0,  27.0,  37.0,
     18.0,  28.0,  38.0,
     19.0,  29.0,  39.0,
     20.0,  30.0,  40.0, };

  printf ( "\n" );
  printf ( "R83ROW_PRINT_PART_TEST\n" );
  printf ( "  R83ROW_PRINT_PART prints part of an R83ROW.\n" );

  max_print = 2;
  r83row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 2" );

  max_print = 5;
  r83row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 5" );

  max_print = 25;
  r83row_print_part ( n, v, max_print, "  Output with MAX_PRINT = 25" );

  return;
}
/******************************************************************************/

void r8block_expand_linear_test ( )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_EXPAND_LINEAR_TEST tests R8BLOCK_EXPAND_LINEAR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2012

  Author:

    John Burkardt
*/
{
# define L 4
# define M 3
# define N 2

  int l2;
  int lfat = 1;
  int m2;
  int mfat = 2;
  int n2;
  int nfat = 1;
  double x[L*M*N] = {
        1.0,  2.0,  3.0,   4.0,  1.0,
        4.0,  9.0, 16.0,   1.0,  8.0,
       27.0, 64.0,  2.0,   4.0,  6.0,
        8.0,  2.0,  8.0,  18.0, 32.0,
        2.0, 16.0, 54.0, 128.0 };
  double *xfat;

  l2 = ( L - 1 ) * ( lfat + 1 ) + 1;
  m2 = ( M - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( N - 1 ) * ( nfat + 1 ) + 1;

  printf ( "\n" );
  printf ( "R8BLOCK_EXPAND_LINEAR_TEST\n" );
  printf ( "  R8BLOCK_EXPAND_LINEAR linearly interpolates new data\n" );
  printf ( "  between old values in a 3D block.\n" );

  r8block_print ( L, M, N, x, "  Original block:" );

  printf ( "\n" );
  printf ( "  LFAT = %d\n", lfat );
  printf ( "  MFAT = %d\n", mfat );
  printf ( "  NFAT = %d\n", nfat );

  xfat = r8block_expand_linear ( L, M, N, x, lfat, mfat, nfat );

  r8block_print ( l2, m2, n2, xfat, "  Fattened block:" );

  free ( xfat );

  return;
# undef L
# undef M
# undef N
}
/******************************************************************************/

void r8block_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_NEW_TEST tests R8BLOCK_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 March 2012

  Author:

    John Burkardt
*/
{
  double ***a;
  double ***b;
  int i;
  int j;
  int k;
  int l;
  int m;
  int n;

  printf ( "\n" );
  printf ( "R8BLOCK_NEW_TEST:\n" );
  printf ( "  R8BLOCK_NEW dynamically creates a 3D array.\n" );
  printf ( "  Array entries can be addressed using the\n" );
  printf ( "  notation \"a[i][j][k]\".\n" );
/*
  These dimensions could be entered by the user; they could depend on
  some other calculation; or they could be changed repeatedly during this
  computation, as long as old memory is deleted by R8BLOCK_DELETE and new memory
  requested by R8BLOCK_NEW.
*/
  l = 2;
  m = 3;
  n = 2;
/*
  Allocate memory.
*/
  printf ( "\n" );
  printf ( "  Allocating memory for array A of size %d by %d by %d.\n", l, m, n );

  a = r8block_new ( l, m, n );

  printf ( "\n" );
  printf ( "  Assigning values to A.\n" );
/*
  Store values in A.
*/
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        a[i][j][k] = ( double ) ( 100 * i + 10 * j + k );
      }
    }
  }
/*
  Print A.
*/
  printf ( "\n" );
  printf ( "  Dynamically allocated matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        printf ( "  %8g", a[i][j][k] );
      }
      printf ( "\n" );
    }
    printf ( "\n" );
  }
/*
  Free memory.
*/
  r8block_delete ( l, m, n, a );

  return;
}
/******************************************************************************/

void r8block_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_PRINT_TEST tests R8BLOCK_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 June 2012

  Author:

    John Burkardt
*/
{
# define L 4
# define M 3
# define N 2

  double x[L*M*N] = {
        1.0,  2.0,  3.0,   4.0,  1.0,
        4.0,  9.0, 16.0,   1.0,  8.0,
       27.0, 64.0,  2.0,   4.0,  6.0,
        8.0,  2.0,  8.0,  18.0, 32.0,
        2.0, 16.0, 54.0, 128.0 };

  printf ( "\n" );
  printf ( "R8BLOCK_PRINT_TEST\n" );
  printf ( "  R8BLOCK_PRINT prints an R8BLOCK.\n" );

  r8block_print ( L, M, N, x, "  The 3D array:" );

  return;
# undef L
# undef M
# undef N
}
/******************************************************************************/

void r8cmat_to_r8mat_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8CMAT_TO_R8MAT_NEW_TEST tests R8CMAT_TO_R8MAT_NEW;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  printf ( "\n" );
  printf ( "R8CMAT_TO_R8MAT_NEW_TEST\n" );
  printf ( "  R8CMAT_TO_R8MAT_NEW converts an R8CMAT to an R8MAT.\n" );
  printf ( "\n" );
  printf ( "  Data is of order (%d,%d).\n", m, n );
/*
  Set the R8MAT.
*/
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
/*
  Convert.
*/
  b = r8mat_to_r8cmat_new ( m, n, a );
  r8cmat_print ( m, n, b, "  The R8CMAT B:" );
/*
  Recover the matrix.
*/
  c = r8cmat_to_r8mat_new ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
/*
  Free memory.
*/
  free ( a );
  r8cmat_delete ( m, n, b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8col_find_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_FIND_TEST tests R8COL_FIND.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  int col;
  double dtab[M*N];
  double r8vec[M];
  int i;
  int j;
  int k;

  k = 1;

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      dtab[i+j*M] = ( double ) k;
      if ( j == 2 )
      {
        r8vec[i] = ( double ) k;
      }
      k = k + 1;
    }
  }

  col = r8col_find ( M, N, dtab, r8vec );

  printf ( "\n" );
  printf ( "R8COL_FIND_TEST\n" );
  printf ( "  R8COL_FIND finds a column in a table matching\n" );
  printf ( "  a given set of data.\n" );
  printf ( "\n" );
  printf ( "  R8COL_FIND returns COL = %d\n", col );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_insert_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_INSERT_TEST tests R8COL_INSERT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N_MAX 10

  double a[M*N_MAX] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0 };
  int col;
  double r8vec1[M] = { 3.0, 7.0, 11.0 };
  double r8vec2[M] = { 3.0, 4.0, 18.0 };
  int n;

  printf ( "\n" );
  printf ( "R8COL_INSERT_TEST\n" );
  printf ( "  R8COL_INSERT inserts new columns into a sorted R8COL.\n" );

  n = 4;

  r8mat_print ( M, n, a, "  The unsorted matrix:" );

  r8col_sort_heap_a ( M, n, a );

  r8mat_print ( M, n, a, "  The sorted matrix:" );

  r8vec_print ( M, r8vec1, "  New column:" );

  col = r8col_insert ( N_MAX, M, n, a, r8vec1 );

  if ( col < 0 )
  {
    printf ( "\n" );
    printf ( "  The data was already in column %d\n", abs ( col ) );
  }
  else
  {
    r8mat_print ( M, n, a, "  The updated matrix:" );
  }

  r8vec_print ( M, r8vec2, "  New column:" );

  col = r8col_insert ( N_MAX, M, n, a, r8vec2 );

  if ( col < 0 )
  {
    printf ( "\n" );
    printf ( "  The data was already in column %d\n", abs ( col ) );
  }
  else
  {
    r8mat_print ( M, n, a, "  The updated matrix:" );
  }

  return;
# undef M
# undef N_MAX
}
/******************************************************************************/

void r8col_sort_heap_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_HEAP_A_TEST tests R8COL_SORT_HEAP_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N_MAX 10

  double a[M*N_MAX] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,
    0.0,  0.0,  0.0 };
  int col;
  double r8vec1[M] = { 3.0, 7.0, 11.0 };
  double r8vec2[M] = { 3.0, 4.0, 18.0 };
  int n;

  printf ( "\n" );
  printf ( "R8COL_SORT_HEAP_A_TEST\n" );
  printf ( "  R8COL_SORT_HEAP_A ascending heap sorts a table of columns.\n" );

  n = 4;

  r8mat_print ( M, n, a, "  The unsorted matrix:" );

  r8col_sort_heap_a ( M, n, a );

  r8mat_print ( M, n, a, "  The sorted matrix:" );

  return;
# undef M
# undef N_MAX
}
/******************************************************************************/

void r8col_sort_heap_index_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_HEAP_INDEX_A_TEST tests R8COL_SORT_HEAP_INDEX_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 15

  double a[M*N] = {
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    3.0,  4.0, 18.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int i;
  int *indx;
  int j;
  int j2;
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8COL_SORT_HEAP_INDEX_A_TEST\n" );
  printf ( "  R8COL_SORT_HEAP_INDEX_A computes an index vector which\n" );
  printf ( "  ascending sorts an R8COL.\n" );

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  indx = r8col_sort_heap_index_a ( m, n, a );

  printf ( "\n" );
  printf ( "  The implicitly sorted R8COL (transposed)\n" );
  printf ( "\n" );

  for ( j = 0; j < n; j++ )
  {
    j2 = indx[j];
    printf ( "  %4d:", j2 );
    for ( i = 0; i < m; i++ )
    {
      printf ( "  %10f", a[i+j2*m] );
    }
    printf ( "\n" );
  }

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_sort_quick_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_QUICK_A_TEST tests R8COL_SORT_QUICK_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 10

  double *a;
  double b = 0.0;
  double c = 10.0;
  int seed;

  printf ( "\n" );
  printf ( "R8COL_SORT_QUICK_A_TEST\n" );
  printf ( "  R8COL_SORT_QUICK_A sorts a table of columns.\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, &seed );

  r8mat_print ( M, N, a, "  The unsorted matrix:" );

  r8col_sort_quick_a ( M, N, a );

  r8mat_print ( M, N, a, "  The sorted matrix:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_sorted_tol_unique_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_TOL_UNIQUE_TEST tests R8COL_SORTED_TOL_UNIQUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  printf ( "\n" );
  printf ( "R8COL_SORTED_TOL_UNIQUE_TEST\n" );
  printf ( "  R8COL_SORTED_TOL_UNIQUE finds tolerably unique columns \n" );
  printf ( "  in a sorted R8COL.\n" );

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  printf ( "\n" );
  printf ( "  Using tolerance = %g\n", tol );

  unique_num = r8col_sorted_tol_unique ( m, n, a, tol );

  printf ( "\n" );
  printf ( "  Number of tolerably unique columns is %d\n", unique_num );

  r8mat_transpose_print ( m, unique_num, a,
    "  The sorted tolerably unique R8COL (transposed):" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_sorted_unique_count_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_UNIQUE_COUNT_TEST tests R8COL_SORTED_UNIQUE_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  printf ( "\n" );
  printf ( "R8COL_SORTED_UNIQUE_COUNT_TEST\n" );
  printf ( "  R8COL_SORTED_UNIQUE_COUNT counts tolerably unique columns \n" );
  printf ( "  in a sorted R8COL.\n" );

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  printf ( "\n" );
  printf ( "  Using tolerance = %g\n", tol );

  unique_num = r8col_sorted_tol_unique_count ( m, n, a, tol );

  printf ( "\n" );
  printf ( "  Number of tolerably unique columns is %d\n", unique_num );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_sorted_tol_undex_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_TOL_UNDEX_TEST tests R8COL_SORTED_TOL_UNDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2010

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  printf ( "\n" );
  printf ( "R8COL_SORTED_TOL_UNDEX_TEST\n" );
  printf ( "  R8COL_SORTED_TOL_UNDEX produces index vectors which create a sorted\n" );
  printf ( "  list of the tolerably unique columns of a sorted R8COL,\n" );
  printf ( "  and a map from the original R8COL to the (implicit)\n" );
  printf ( "  R8COL of sorted tolerably unique elements.\n" );

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  r8col_sort_heap_a ( m, n, a );

  r8mat_transpose_print ( m, n, a, "  The sorted R8COL (transposed):" );

  tol = 0.25;

  printf ( "\n" );
  printf ( "  Using tolerance = %g\n", tol );

  n_unique = r8col_sorted_tol_unique_count ( m, n, a, tol );

  printf ( "\n" );
  printf ( "  Number of tolerably unique columns is %d\n", n_unique );

  au = ( double * ) malloc ( m * n_unique * sizeof ( double ) );
  undx = ( int * ) malloc ( n_unique * sizeof ( int ) );
  xdnu = ( int * ) malloc ( n * sizeof ( int ) );

  r8col_sorted_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  printf ( "\n" );
  printf ( "  XDNU points to the representative for each item.\n" );
  printf ( "  UNDX selects the representatives.\n" );
  printf ( "\n" );
  printf ( "     I  XDNU  UNDX\n" );
  printf ( "\n" );
  for ( i = 0; i < n_unique; i++ )
  {
    printf ( "  %4d  %4d  %4d\n", i, xdnu[i], undx[i] );
  }
  for ( i = n_unique; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, xdnu[i] );
  }
  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  free ( au );
  free ( undx );
  free ( xdnu );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_max_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_MAX_TEST tests R8COL_MAX;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8COL_MAX_TEST\n" );
  printf ( "  R8COL_MAX computes maximums of an R8COL;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  amax = r8col_max ( M, N, a );

  r8vec_print ( N, amax, "  Column maximums:" );

  free ( amax );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_MEAN_TEST tests R8COL_MEAN and R8COL_SUM;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *mean;

  printf ( "\n" );
  printf ( "R8COL_MEAN_TEST\n" );
  printf ( "  R8COL_MEAN computes means of an R8COL;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  mean = r8col_mean ( M, N, a );

  r8vec_print ( N, mean, "  The column means:" );

  free ( mean );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_min_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_MIN_TEST tests R8COL_MIN;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  double *amin;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8COL_MIN_TEST\n" );
  printf ( "  R8COL_MIN computes minimums of an R8COL;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  amin = r8col_min ( M, N, a );

  r8vec_print ( N, amin, "  Column minimums:" );

  free ( amin );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_permute_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_PERMUTE_TEST tests R8COL_PERMUTE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 5

  double a[M*N] = {
    11.0, 21.0, 31.0,
    12.0, 22.0, 32.0,
    13.0, 23.0, 33.0,
    14.0, 24.0, 34.0,
    15.0, 25.0, 35.0 };
  int perm[N] = { 1, 3, 4, 0, 2 };


  printf ( "\n" );
  printf ( "R8COL_PERMUTE_TEST\n" );
  printf ( "  R8COL_PERMUTE permutes an R8COL in place.\n" );

  r8mat_print ( M, N, a, "  A (unpermuted):" );

  i4vec_print ( N, perm, "  The (column) permutation vector:" );

  r8col_permute ( M, N, perm, a );

  r8mat_print ( M, N, a, "  A (permuted):" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_sortr_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTR_A_TEST tests R8COL_SORTR_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
# define M 10
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int key;
  int seed;

  printf ( "\n" );
  printf ( "R8COL_SORTR_A_TEST\n" );
  printf ( "  R8COL_SORTR_A is given an array, and reorders\n" );
  printf ( "  it so that a particular column is sorted.\n" );

  key = 2;
  printf ( "\n" );
  printf ( "  Here, the special column is %d\n", key );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, &seed );

  r8mat_print ( M, N, a, "  Unsorted array:" );

  r8col_sortr_a ( M, N, a, key );

  r8mat_print ( M, N, a, "  Sorted array:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_sum_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SUM_TEST tests R8COL_SUM;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  double *colsum;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8COL_SUM_TEST\n" );
  printf ( "  R8COL_SUM computes sums of an R8COL;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  colsum = r8col_sum ( M, N, a );

  r8vec_print ( N, colsum, "  The column sums:" );

  free ( colsum );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_swap_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_SWAP_TEST tests R8COL_SWAP;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int icol1;
  int icol2;
  int m = 3;
  int n = 4;

  printf ( "\n" );
  printf ( "R8COL_SWAP_TEST\n" );
  printf ( "  R8COL_SWAP swaps two columns of an R8COL;\n" );

  a = r8mat_indicator_new ( m, n );

  r8mat_print ( m, n, a, "  The array:" );

  icol1 = 1;
  icol2 = 3;

  printf ( "\n" );
  printf ( "  Swap columns %d and %d:\n", icol1, icol2 );

  r8col_swap ( m, n, a, icol1, icol2 );

  r8mat_print ( m, n, a, "  The updated matrix:" );

  free ( a );

  return;
}
/******************************************************************************/

void r8col_to_r8vec_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_TO_R8VEC_TEST tests R8COL_TO_R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  double *x;

  printf ( "\n" );
  printf ( "R8COL_TO_R8VEC_TEST\n" );
  printf ( "  R8COL_TO_R8VEC converts an array of columns to a vector.\n" );
  printf ( "\n" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*M] = ( double ) ( 10 * i + j );
    }
  }

  r8mat_print ( M, N, a, "  The array of columns:" );

  x = r8col_to_r8vec ( M, N, a );

  r8vec_print ( M*N, x, "  The resulting vector of columns:" );

  free ( x );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_tol_undex_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_TOL_UNDEX_TEST tests R8COL_TOL_UNDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  double tol;
  int *undx;
  int unique_num;
  int *xdnu;

  printf ( "\n" );
  printf ( "R8COL_TOL_UNDEX_TEST\n" );
  printf ( "  R8COL_TOL_UNDEX produces index vectors which create a sorted\n" );
  printf ( "  list of the tolerably unique columns of an R8COL,\n" );
  printf ( "  and a map from the original R8COL to the (implicit)\n" );
  printf ( "  R8COL of sorted tolerably unique elements.\n" );

  r8mat_transpose_print ( m, n, a, "  The unsorted R8COL (transposed):" );

  tol = 0.25;

  printf ( "\n" );
  printf ( "  Using tolerance = %f\n", tol );

  n_unique = r8col_tol_unique_count ( m, n, a, tol );

  printf ( "\n" );
  printf ( "  Number of tolerably unique columns is %d\n", n_unique );

  au = ( double * ) malloc ( m * n_unique * sizeof ( double ) );
  undx = ( int * ) malloc ( n_unique * sizeof ( int ) );
  xdnu = ( int * ) malloc ( n * sizeof ( n ) );

  r8col_tol_undex ( m, n, a, n_unique, tol, undx, xdnu );

  printf ( "\n" );
  printf ( "  XDNU points to the representative for each item.\n" );
  printf ( "  UNDX selects the representatives.\n" );
  printf ( "\n" );
  printf ( "     I  XDNU  UNDX\n" );
  printf ( "\n" );
  for ( i = 0; i < n_unique; i++ )
  {
    printf ( "  %4d  %4d  %4d\n", i, xdnu[i], undx[i] );
  }
  for ( i = n_unique; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, xdnu[i] );
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au,
    "  The tolerably unique R8COL (transposed):" );

  free ( au );
  free ( undx );
  free ( xdnu );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_undex_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNDEX_TEST tests R8COL_UNDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  double *au;
  int i;
  int j;
  int j2;
  int m = M;
  int n = N;
  int n_unique;
  int *undx;
  int unique_num;
  int *xdnu;

  printf ( "\n" );
  printf ( "R8COL_UNDEX_TEST\n" );
  printf ( "  R8COL_UNDEX produces index vectors which create a sorted\n" );
  printf ( "  list of the unique columns of an (unsorted) R8COL,\n" );
  printf ( "  and a map from the original R8COL to the (implicit)\n" );
  printf ( "  R8COL of sorted unique elements.\n" );

  r8mat_transpose_print ( m, n, a, "  The R8COL (transposed):" );

  n_unique = r8col_unique_count ( m, n, a );

  printf ( "\n" );
  printf ( "  Number of unique columns is %d\n", n_unique );

  au = ( double * ) malloc ( m * n_unique * sizeof ( double ) );
  undx = ( int * ) malloc ( n_unique * sizeof ( int ) );
  xdnu = ( int * ) malloc ( n * sizeof ( int ) );

  r8col_undex ( m, n, a, n_unique, undx, xdnu );

  printf ( "\n" );
  printf ( "  XDNU points to the representative for each item.\n" );
  printf ( "  UNDX selects the representatives.\n" );
  printf ( "\n" );
  printf ( "     I  XDNU  UNDX\n" );
  printf ( "\n" );
  for ( i = 0; i < n_unique; i++ )
  {
    printf ( "  %4d  %4d  %4d\n", i, xdnu[i], undx[i] );
  }
  for ( i = n_unique; i < n; i++ )
  {
    printf ( "  %4d  %4d\n", i, xdnu[i] );
  }

  for ( j = 0; j < n_unique; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      au[i+j*m] = a[i+undx[j]*m];
    }
  }

  r8mat_transpose_print ( m, n_unique, au, "  The Unique R8COL (transposed):" );

  free ( au );
  free ( undx );
  free ( xdnu );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_unique_count_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIQUE_COUNT_TEST tests R8COL_UNIQUE_COUNT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 22

  double a[M*N] = {
    1.9,  0.0, 10.0,
    2.0,  6.0, 10.0,
    4.0,  8.0, 12.0,
    1.0,  5.0,  9.0,
    3.0,  7.0, 11.0,
    2.0,  6.0,  0.0,
    2.0,  0.0, 10.1,
    2.0,  0.1, 10.0,
    3.0,  4.0, 18.0,
    1.9,  8.0, 10.0,
    0.0,  0.0,  0.0,
    0.0,  6.0, 10.0,
    2.1,  0.0, 10.0,
    2.0,  6.0, 10.0,
    3.0,  7.0, 11.0,
    2.0,  0.0, 10.0,
    2.0,  0.0, 10.0,
    2.0,  6.0, 10.0,
    1.0,  5.0,  9.0,
    2.0,  0.0, 10.1,
    1.0,  5.0,  9.1,
    1.0,  5.1,  9.0 };
  int m = M;
  int n = N;
  double tol;
  int unique_num;

  printf ( "\n" );
  printf ( "R8COL_UNIQUE_COUNT_TEST\n" );
  printf ( "  R8COL_UNIQUE_COUNT counts unique columns.\n" );

  r8mat_transpose_print ( m, n, a, "  The R8COL (transposed):" );

  unique_num = r8col_unique_count ( m, n, a );

  printf ( "\n" );
  printf ( "  Number of unique columns is %d\n", unique_num );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8col_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    R8COL_VARIANCE_TEST tests R8COL_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  printf ( "\n" );
  printf ( "R8COL_VARIANCE_TEST\n" );
  printf ( "  R8COL_VARIANCE computes variances of an R8COL;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( k );
    }
  }

  r8mat_print ( M, N, a, "  The array:" );

  variance = r8col_variance ( M, N, a );

  printf ( "\n" );
  printf ( "  Column  variance:\n" );
  printf ( "\n" );

  for ( j = 0; j < N; j++ )
  {
    printf ( "  %6d  %10g\n", j, variance[j] );
  }

  free ( variance );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8r8vec_index_insert_unique_test ( )

/******************************************************************************/
/*
  Purpose:

    R8R8VEC_INDEX_INSERT_UNIQUE_TEST tests R8R8VEC_INDEX_INSERT_UNIQUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 30

  int i;
  int ierror;
  int indx[N_MAX];
  int ival;
  int n;
  int seed;
  double x[N_MAX];
  double x_max = 4.0;
  double x_min = 1.0;
  double xval;
  double y[N_MAX];
  double y_max = 3.0;
  double y_min = 1.0;
  double yval;

  n = 0;

  printf ( "\n" );
  printf ( "R8R8VEC_INDEX_INSERT_UNIQUE_TEST\n" );
  printf ( "  R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n" );
  printf ( "  index sorted array.\n" );
  printf ( "\n" );
  printf ( "  Generate %d random values:\n", N_MAX );
  printf ( "\n" );
  printf ( "    XVAL    YVAL   Index\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( x_min, x_max, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    yval = r8_uniform_ab ( y_min, y_max, &seed );
    yval = ( double ) ( r8_nint ( yval ) );

    r8r8vec_index_insert_unique ( N_MAX, &n, x, y, indx, xval, yval,
      &ival, &ierror );

    printf ( "  %6d  %12f  %12f\n", ival, xval, yval );
  }

  printf ( "\n" );
  printf ( "  Vector of unique X Y values:\n" );
  printf ( "\n" );
  printf ( "  I  X(I)   Y(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %12f  %12f\n", i+1, x[i], y[i] );
  }

  printf ( "\n" );
  printf (  "  X, Y sorted by index\n" );
  printf (  "\n" );
  printf (  "  I  INDX(I)  X(INDX(I))  Y(INDX(I))\n" );
  printf (  "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %6d  %12f  %12f\n", 
      i+1, indx[i], x[indx[i]-1], y[indx[i]-1] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8r8r8vec_index_insert_unique_test ( )

/******************************************************************************/
/*
  Purpose:

    R8R8R8VEC_INDEX_INSERT_UNIQUE_TEST tests R8R8R8VEC_INDEX_INSERT_UNIQUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 30

  int i;
  int ierror;
  int indx[N_MAX];
  int ival;
  int n;
  int seed;
  double x[N_MAX];
  double xval;
  double y[N_MAX];
  double yval;
  double z[N_MAX];
  double zval;

  n = 0;

  printf ( "\n" );
  printf ( "R8R8R8VEC_INDEX_INSERT_UNIQUE_TEST\n" );
  printf ( "  R8R8R8VEC_INDEX_INSERT_UNIQUE inserts unique values into\n" );
  printf ( "  an index sorted array.\n" );
  printf ( "\n" );
  printf ( "  Number of random values to generate = %d\n", N_MAX );
  printf ( "\n" );
  printf ( "    XVAL    YVAL  ZVAL  Index\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( 1.0, 4.0, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    yval = r8_uniform_ab ( 1.0, 3.0, &seed );
    yval = ( double ) ( r8_nint ( yval ) );
    zval = r8_uniform_ab ( 1.0, 4.0, &seed );
    zval = ( double ) ( r8_nint ( zval ) );

    r8r8r8vec_index_insert_unique ( N_MAX, &n, x, y, z, indx,
      xval, yval, zval, &ival, &ierror );

    printf ( "  %6g  %6g  %6g  %6d\n", xval, yval, zval, ival );
  }

  printf ( "\n" );
  printf ( "  Vector of unique X Y Z values:\n" );
  printf ( "\n" );
  printf ( "  I  X(I)   Y(I)    Z(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %6g  %6g  %6g\n", i+1, x[i], y[i], z[i] );
  }

  printf ( "\n" );
  printf ( "  X Y Z sorted by index:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %6d  %6g  %6g  %6g\n", 
      i+1, indx[i], x[indx[i]-1], y[indx[i]-1], z[indx[i]-1] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8int_to_i4int_test ( )

/******************************************************************************/
/*
  Purpose:

    R8INT_TO_I4INT_TEST tests R8INT_TO_I4INT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2012

  Author:

    John Burkardt
*/
{
  int i;
  int ihi = 11;
  int ilo = 1;
  int ir;
  double r;
  double r2;
  double rhi = 200.0;
  double rhi2;
  double rlo = 100.0;
  double rlo2;
  int seed;
  int test;
  int test_num = 10;

  printf ( "\n" );
  printf ( "R8INT_TO_I4INT_TEST\n" );
  printf ( "  For data in an interval,\n" );
  printf ( "  R8INT_TO_I4INT converts a real to an integer.\n" );
  printf ( "\n" );
  printf ( "  Integer interval: [%d, %d]\n", ilo, ihi );
  printf ( "  Real interval:    [%g, %g]\n", rlo, rhi );
  printf ( "\n" );
  printf ( "  R   I(R)  R(I(R))\n" );
  printf ( "\n" );

  seed = 123456789;

  rlo2 = rlo - 15.0;
  rhi2 = rhi + 15.0;

  for ( test = 1; test <= test_num; test++ )
  {
    r = r8_uniform_ab ( rlo2, rhi2, &seed );
    ir = r8int_to_i4int ( rlo, rhi, r, ilo, ihi );
    r2 = i4int_to_r8int ( ilo, ihi, ir, rlo, rhi );
    printf ( "  %12g  %6d  %12g\n", r, ir, r2 );
  }

  return;
}
/******************************************************************************/

void r8mat_cholesky_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHOLESKY_INVERSE_TEST tests R8MAT_CHOLESKY_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *a2;
  double *a3;
  int i;
  int j;
  int n = 5;
 
  printf ( "\n" );
  printf ( "R8MAT_CHOLESKY_INVERSE_TEST\n" );
  printf ( "  For a positive definite symmetric matrix,\n" );
  printf ( "  R8MAT_CHOLESKY_INVERSE computes the inverse.\n" );

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*n] = -1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Matrix to be inverted:" );

  a2 = r8mat_copy_new ( n, n, a );

  r8mat_cholesky_inverse ( n, a2 );

  r8mat_print ( n, n, a2, "  Inverse matrix:" );

  a3 = r8mat_mm_new ( n, n, n, a2, a );
  
  r8mat_print ( n, n, a3, "  Product inv(A) * A:" );

  free ( a );
  free ( a2 );
  free ( a3 );

  return;
}
/******************************************************************************/

void r8mat_cholesky_solve_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHOLESKY_SOLVE_TEST tests R8MAT_CHOLESKY_SOLVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *d;
  int flag;
  int i;
  int j;
  double *l;
  int n = 5;
  double *r;
  int test;
  double *x;

  printf ( "\n" );
  printf ( "R8MAT_CHOLESKY_SOLVE_TEST\n" );
  printf ( "  For a positive definite symmetric matrix,\n" );
  printf ( "  R8MAT_CHOLESKY_SOLVE solves a linear system\n" );
  printf ( "  using the lower Cholesky factorization.\n" );

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*n] = -1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }
  r8mat_print ( n, n, a, "  Matrix to be factored:" );
/*
  Compute the Cholesky factor.
*/
  l = r8mat_cholesky_factor ( n, a, &flag );

  if ( flag != 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  R8MAT_CHOLESKY_FACTOR failed.\n" );
    return;
  }
  r8mat_print ( n, n, l, "  Cholesky factor L:" );
  d = r8mat_mmt_new ( n, n, n, l, l );
  r8mat_print ( n, n, d, "  Product L * L':" );
/*
  Solve a system.
*/
  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = 0.0;
  }
  b[n-1] = ( double ) ( n + 1 );

  r8vec_print ( n, b, "  Right hand side:" );

  x = r8mat_cholesky_solve ( n, l, b );

  r8vec_print ( n, x, "  Computed solution:" );

  free ( a );
  free ( b );
  free ( d );
  free ( l );
  free ( x );

  return;
}
/******************************************************************************/

void r8mat_cholesky_solve_upper_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHOLESKY_SOLVE_UPPER_TEST tests R8MAT_CHOLESKY_SOLVE_UPPER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 October 2013

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *d;
  int flag;
  int i;
  int j;
  double *l;
  int n = 5;
  double *r;
  int test;
  double *x;

  printf ( "\n" );
  printf ( "R8MAT_CHOLESKY_SOLVE_UPPER_TEST\n" );
  printf ( "  For a positive definite symmetric matrix,\n" );
  printf ( "  R8MAT_CHOLESKY_SOLVE_UPPER solves a linear system\n" );
  printf ( "  using the upper Cholesky factorization.\n" );

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a[i+j*n] = -1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Matrix to be factored:" );
/*
  Compute the Cholesky factor.
*/

  r = r8mat_cholesky_factor_upper ( n, a, &flag );
  if ( flag != 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "  R8MAT_CHOLESKY_FACTOR_UPPER failed.\n" );
    return;
  }
  r8mat_print ( n, n, r, "  Cholesky factor R:" );
  d = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, d, "  Product R' * R:" );
/*
  Solve a system.
*/
  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = 0.0;
  }
  b[n-1] = ( double ) ( n + 1 );

  r8vec_print ( n, b, "  Right hand side:" );

  x = r8mat_cholesky_solve_upper ( n, r, b );

  r8vec_print ( n, x, "  Computed solution:" );

  free ( a );
  free ( b );
  free ( d );
  free ( r );
  free ( x );

  return;
}
/******************************************************************************/

void r8mat_det_2d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_2D_TEST tests R8MAT_DET_2D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 June 2012

  Author:

    John Burkardt
*/
{
# define N 2

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0 };

  printf ( "\n" );
  printf ( "R8MAT_DET_2D_TEST\n" );
  printf ( "  R8MAT_DET_2D: determinant of a 2 by 2 matrix;\n" );

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_2d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  printf ( "\n" );
  printf ( "  R8MAT_DET_2D computes determinant: %g\n", det );
/*
  Special formula for the determinant of a Vandermonde matrix:
*/
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  printf ( "  Exact determinant is %g\n", det );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8mat_det_3d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_3D_TEST tests R8MAT_DET_3D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2012

  Author:

    John Burkardt
*/
{
# define N 3

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0 };

  printf ( "\n" );
  printf ( "R8MAT_DET_3D_TEST\n" );
  printf ( "  R8MAT_DET_3D: determinant of a 3 by 3 matrix;\n" );

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_3d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  printf ( "\n" );
  printf ( "  R8MAT_DET_3D computes determinant: %g\n", det );
/*
  Special formula for the determinant of a Vandermonde matrix:
*/
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  printf ( "  Exact determinant is %g\n", det );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8mat_det_4d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_4D_TEST tests R8MAT_DET_4D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2012

  Author:

    John Burkardt
*/
{
# define N 4

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0, 2.0 };

  printf ( "\n" );
  printf ( "R8MAT_DET_4D\n" );
  printf ( "  R8MAT_DET_4D determinant of a 4 by 4 matrix;\n" );

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_4d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  printf ( "\n" );
  printf ( "  R8MAT_DET_4D computes determinant:%g\n", det );
/*
  Special formula for the determinant of a Vandermonde matrix:
*/
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  printf ( "  Exact determinant is %g\n", det );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8mat_det_5d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_5D_TEST tests R8MAT_DET_5D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double det;
  int i;
  int j;
  double x[N] = { 1.0, 10.0, 4.0, 2.0, 3.0 };

  printf ( "\n" );
  printf ( "R8MAT_DET_5D_TEST\n" );
  printf ( "  R8MAT_DET_5D determinant of a 5 by 5 matrix;\n" );

  a = r8mat_vand2 ( N, x );
  det = r8mat_det_5d ( a );

  r8mat_print ( N, N, a, "  Matrix:" );

  printf ( "\n" );
  printf ( "  R8MAT_DET_5D computes determinant:%g\n", det );
/*
  Special formula for the determinant of a Vandermonde matrix:
*/
  det = 1.0;
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      det = det * ( x[i] - x[j] );
    }
  }
  printf ( "  Exact determinant is %g\n", det );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8mat_expand_linear_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_EXPAND_LINEAR_TEST tests R8MAT_EXPAND_LINEAR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2012

  Author:

    John Burkardt
*/
{
# define M 4
# define N 3

  int m2;
  int mfat = 2;
  int n2;
  int nfat = 1;
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double x[M*N] = {
    1.0, 2.0, 3.0, 4.0, 1.0,
    4.0, 9.0, 16.0, 1.0, 8.0,
    27.0, 64.0 };
  double *xfat;

  m2 = ( M - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( N - 1 ) * ( nfat + 1 ) + 1;

  printf ( "\n" );
  printf ( "R8MAT_EXPAND_LINEAR_TEST\n" );
  printf ( "  R8MAT_EXPAND_LINEAR linearly interpolates new data\n" );
  printf ( "  between old values in a matrix.\n" );

  r8mat_print ( M, N, x, "  Original matrix:" );

  printf ( "\n" );
  printf ( "  MFAT = %d\n", mfat );
  printf ( "  NFAT = %d\n", nfat );

  xfat = r8mat_expand_linear ( M, N, x, mfat, nfat );

  r8mat_print ( m2, n2, xfat, "  Fattened matrix:" );

  free ( xfat );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_expand_linear2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_EXPAND_LINEAR2_TEST tests R8MAT_EXPAND_LINEAR2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 2

  double a[M*N];
  double *a2;
  int i;
  int j;
  int m2 = 10;
  int n2 = 5;

  printf ( "\n" );
  printf ( "R8MAT_EXPAND_LINEAR2_TEST\n" );
  printf ( "  R8MAT_EXPAND_LINEAR2 fills in a large array by\n" );
  printf ( "  interpolating data from a small array.\n" );
  printf ( "\n" );
  printf ( "  Original matrix has dimensions:\n" );
  printf ( "\n" );
  printf ( "  M = %d, N = %d\n", M, N );
  printf ( "\n" );
  printf ( "  Expanded matrix has dimensions:\n" );
  printf ( "\n" );
  printf ( "  M2 = %d, N2 = %d\n", m2, n2 );

  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      a[i-1+(j-1)*M] = 10.0 * ( double ) ( i ) + ( double ) ( j );
    }
  }

  r8mat_print ( M, N, a, "  The little matrix A:" );

  a2 = r8mat_expand_linear2 ( M, N, a, m2, n2 );

  r8mat_print ( m2, n2, a2, "  Expanded array A2:" );

  free ( a2 );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_fs_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FS_NEW_TEST tests R8MAT_FS_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 January 2013

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "R8MAT_FS_NEW_TEST\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8MAT_FS_NEW factors and solves a linear system.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
/*
  Set the matrix.
*/
  a = r8mat_uniform_01_new ( n, n, &seed );
/*
  Set the desired solutions.
*/
  b = ( double * ) malloc ( n * sizeof ( double ) );

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*n] * x[j];
    }
  }
/*
  Factor and solve the system.
*/
  free ( x );

  x = r8mat_fs_new ( n, a, b );
  
  r8vec_print ( n, x, "  Solution:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
# undef NB
}
/******************************************************************************/

void r8mat_fss_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FSS_NEW_TEST tests R8MAT_FSS_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt
*/
{
# define N 10
# define NB 3

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int nb = NB;
  int seed = 123456789;
  double *x;

  printf ( "\n" );
  printf ( "R8MAT_FSS_NEW_TEST\n" );
  printf ( "  For a matrix in general storage,\n" );
  printf ( "  R8MAT_FSS_NEW factors and solves multiple linear systems.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", n );
  printf ( "  Number of systems NB = %d\n", nb );
/*
  Set the matrix.
*/
  a = r8mat_uniform_01_new ( n, n, &seed );
/*
  Set the desired solutions.
*/
  b = ( double * ) malloc ( n * nb * sizeof ( double ) );

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
  }
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  k = 1;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( i % 3 ) + 1;
  }
  k = 2;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
/*
  Factor and solve the system.
*/
  free ( x );

  x = r8mat_fss_new ( n, a, nb, b );
  
  r8mat_print ( n, nb, x, "  Solutions:" );

  free ( a );
  free ( b );
  free ( x );

  return;
# undef N
# undef NB
}
/******************************************************************************/

void r8mat_givens_post_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_GIVENS_POST_TEST tests R8MAT_GIVENS_POST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2012

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N*N];
  double *ag;
  int col;
  double *g;
  int i;
  int j;
  int row;

  printf ( "\n" );
  printf ( "R8MAT_GIVENS_POST_TEST\n" );
  printf ( "  R8MAT_GIVENS_POST computes a Givens postmultiplier rotation matrix.\n" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_power ( i + 1, j );
    }
  }

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 3;
  col = 2;

  printf ( "\n" );
  printf ( "  I = %d  J = %d\n", row, col );

  g = r8mat_givens_post ( N, a, row, col );

  r8mat_print ( N, N, g, "  G" );

  ag = r8mat_mm_new ( N, N, N, a, g );

  r8mat_print ( N, N, ag, "  A*G" );

  free ( ag );
  free ( g );

  return;
# undef N
}
/******************************************************************************/

void r8mat_givens_pre_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_GIVENS_PRE_TEST tests R8MAT_GIVENS_PRE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N*N];
  int col;
  double *g;
  double *ga;
  int i;
  int j;
  int row;

  printf ( "\n" );
  printf ( "R8MAT_GIVENS_PRE_TEST\n" );
  printf ( "  R8MAT_GIVENS_PRE computes a Givens premultiplier rotation matrix.\n" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_power ( i + 1, j );
    }
  }

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 3;
  col = 2;

  printf ( "\n" );
  printf ( "  I = %d  J = %d\n", row, col );

  g = r8mat_givens_pre ( N, a, row, col );

  r8mat_print ( N, N, g, "  G" );

  ga = r8mat_mm_new ( N, N, N, g, a );

  r8mat_print ( N, N, ga, "  G*A" );

  free ( g );
  free ( ga );

  return;
# undef N
}
/******************************************************************************/

void r8mat_hess_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HESS_TEST tests R8MAT_HESS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define N 3

  double *h;
  double x[N] = { 1.0, 2.0, 3.0 };

  printf ( "\n" );
  printf ( "R8MAT_HESS_TEST\n" );
  printf ( "  R8MAT_HESS estimates the Hessian matrix\n" );
  printf ( "  of a scalar function.\n" );

  h = r8mat_hess ( r8mat_hess_f, N, x );

  r8mat_print ( N, N, h, "  Estimated Hessian:" );

  free ( h );

  h = r8mat_hess_exact ( N, x );

  r8mat_print ( N, N, h, "  Exact Hessian:" );

  free ( h );

  return;
# undef N
}
/******************************************************************************/

double r8mat_hess_f ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HESS_F is a sample nonlinear function for treatment by R8MAT_HESS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of parameters.

    Input, double X[N], the parameter values.

    Output, double R8MAT_HESS_F, the function value.
*/
{
  double f;

  f = x[0] * x[0] + x[0] * x[1] + x[1] * cos ( 10.0 * x[2] );

  return f;
}
/******************************************************************************/

double *r8mat_hess_exact ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HESS_EXACT is the exact Hessian of R8MAT_HESS_F.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of parameters.

    Input, double X[N], the parameter values.

    Output, double R8MAT_HESS_EXACT[N*N], the Hessian values.
*/
{
  double *h;

  h = ( double * ) malloc ( n * n * sizeof ( double ) );

  h[0+0*3] = 2.0;
  h[0+1*3] = 1.0;
  h[0+2*3] = 0.0;

  h[1+0*3] = 1.0;
  h[1+1*3] = 0.0;
  h[1+2*3] = -10.0 * sin ( 10.0 * x[2] );

  h[2+0*3] = 0.0;
  h[2+1*3] = -10.0 * sin ( 10.0 * x[2] );
  h[2+2*3] = -100.0 * x[1] * cos ( 10.0 * x[2] );

  return h;
}
/******************************************************************************/

void r8mat_house_axh_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_AXH_TEST tests R8MAT_HOUSE_AXH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 March 2015

  Author:

    John Burkardt
*/
{
  double *a;
  double *a_col;
  double *ah;
  double *h;
  double *ha;
  int k;
  int n = 5;
  double r8_hi;
  double r8_lo;
  int seed;
  double *v;

  printf ( "\n" );
  printf ( "R8MAT_HOUSE_AXH_TEST\n" );
  printf ( "  R8MAT_HOUSE_AXH multiplies a matrix A times a\n" );
  printf ( "  compact Householder matrix.\n" );

  r8_lo = -5.0;
  r8_hi = +5.0;
  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, r8_lo, r8_hi, &seed );

  r8mat_print ( n, n, a, "  Matrix A:" );
/*
  Request V, the compact form of the Householder matrix H
  such that H*A packs column 3 of A.

  Look at me using a pointer!
*/
  k = 3;
  a_col = ( a + ( k - 1 ) * n );
  v = r8vec_house_column ( n, a_col, k );

  r8vec_print ( n, v, "  Compact vector V so H*A packs column 3:" );

  h = r8mat_house_form ( n, v );

  r8mat_print ( n, n, h, "  Householder matrix H:" );
/*
  Compute A*H.
*/
  ah = r8mat_house_axh_new ( n, a, v );

  r8mat_print ( n, n, ah, "  Indirect product A*H:" );

  free ( ah );
/*
  Compare with a direct calculation.
*/
  ah = r8mat_mm_new ( n, n, n, a, h );

  r8mat_print ( n, n, ah, "  Direct product A*H:" );
/*
  Compute H*A to verify packing.
*/
  ha = r8mat_mm_new ( n, n, n, h, a );

  r8mat_print ( n, n, ha, "  Product H*A has packed column 3:" );
/*
  Free memory.
*/
  free ( a );
  free ( ah );
  free ( h );
  free ( ha );
  free ( v );

  return;
# undef N
}
/******************************************************************************/

void r8mat_house_form_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_FORM_TEST tests R8MAT_HOUSE_FORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *h;
  double v[N] = { 0.0, 0.0, 1.0, 2.0, 3.0 };

  printf ( "\n" );
  printf ( "R8MAT_HOUSE_FORM_TEST\n" );
  printf ( "  R8MAT_HOUSE_FORM forms a Householder\n" );
  printf ( "  matrix from its compact form.\n" );

  r8vec_print ( N, v, "  Compact vector form V:" ) ;

  h = r8mat_house_form ( N, v );

  r8mat_print ( N, N, h, "  Householder matrix H:" );

  free ( h );

  return;
# undef N
}
/******************************************************************************/

void r8mat_house_post_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_POST_TEST tests R8MAT_HOUSE_POST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *ah;
  double b = 0.0;
  double c = 5.0;
  int col;
  double *h;
  int n = 5;
  int row;
  int seed;

  printf ( "\n" );
  printf ( "R8MAT_HOUSE_POST_TEST\n" );
  printf ( "  R8MAT_HOUSE_POST computes a Householder postmultiplier;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, b, c, &seed );

  r8mat_print ( n, n, a, "  Matrix A:" );

  row = 2;
  col = 3;

  printf ( "\n" );
  printf ( "  I = %d  J = %d\n", row, col );

  h = r8mat_house_post ( n, a, row, col );

  r8mat_print ( n, n, h, "  Householder matrix H:" );

  ah = r8mat_mm_new ( n, n, n, a, h );

  r8mat_print ( n, n, ah, "  Product A*H:" );

  free ( a );
  free ( ah );
  free ( h );

  return;
}
/******************************************************************************/

void r8mat_house_pre_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_PRE_TEST tests R8MAT_HOUSE_PRE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double b = 0.0;
  double c = 5.0;
  int col;
  double *h;
  double *ha;
  int row;
  int seed;

  printf ( "\n" );
  printf ( "R8MAT_HOUSE_PRE_TEST\n" );
  printf ( "  R8MAT_HOUSE_PRE computes a Householder premultiplier;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( N, N, b, c, &seed );

  r8mat_print ( N, N, a, "  Matrix A:" );

  row = 2;
  col = 3;

  printf ( "\n" );
  printf ( "  I = %d  J = %d\n", row, col );

  h = r8mat_house_pre ( N, a, row, col );

  r8mat_print ( N, N, h, "  Householder matrix H:" );

  ha = r8mat_mm_new ( N, N, N, h, a );

  r8mat_print ( N, N, ha, "  Product H*A:" );

  free ( a );
  free ( h );
  free ( ha );

  return;
# undef N
}
/******************************************************************************/

void r8mat_indicator_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INDICATOR_NEW_TEST tests R8MAT_INDICATOR_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int m = 5;
  int n = 4;

  printf ( "\n" );
  printf ( "R8MAT_INDICATOR_NEW_TEST\n" );
  printf ( "  R8MAT_INDICATOR_NEW returns an indicator matrix;\n" );

  a = r8mat_indicator_new ( m, n );

  r8mat_print ( m, n, a, "  Indicator matrix:" );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_inverse_2d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INVERSE_2D_TEST tests R8MAT_INVERSE_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define N 2

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8MAT_INVERSE_2D_TEST\n" );
  printf ( "  R8MAT_INVERSE_2D inverts a 2 by 2 matrix.\n" );

  a[0+0*N] = 1.0;
  a[0+1*N] = 2.0;

  a[1+0*N] = 3.0;
  a[1+1*N] = 4.0;

  r8mat_print ( 2, 2, a, "  Matrix A:" );

  b = r8mat_inverse_2d ( a );

  r8mat_print ( 2, 2, b, "  Inverse matrix A:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 2, 2, c, "  Product C = A * B:" );

  free ( b );

  return;

# undef N
}
/******************************************************************************/

void r8mat_inverse_3d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INVERSE_3D_TEST tests R8MAT_INVERSE_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8MAT_INVERSE_3D_TEST\n" );
  printf ( "  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.\n" );

  a[0+0*N] = 3.0;
  a[0+1*N] = 2.0;
  a[0+2*N] = 1.0;

  a[1+0*N] = 2.0;
  a[1+1*N] = 2.0;
  a[1+2*N] = 1.0;

  a[2+0*N] = 0.0;
  a[2+1*N] = 1.0;
  a[2+2*N] = 1.0;

  r8mat_print ( 3, 3, a, "  Matrix A:" );

  b = r8mat_inverse_3d ( a );

  r8mat_print ( 3, 3, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 3, 3, c, "  C = A * B:" );

  free ( b );

  return;

# undef N
}
/******************************************************************************/

void r8mat_inverse_4d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INVERSE_4D_TEST tests R8MAT_INVERSE_4D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N];
  double *b;
  double c[N*N];
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8MAT_INVERSE_4D_TEST\n" );
  printf ( "  R8MAT_INVERSE_4D inverts a 4 x 4 matrix.\n" );


  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i <= j )
      {
        a[i+j*N] = ( double ) ( N - j );
      }
      else if ( j == i - 1 )
      {
        a[i+j*N] = ( double ) ( N - j - 1 );
      }
      else {
        a[i+j*N] = 0.0;
      }
    }
  }

  r8mat_print ( 4, 4, a, "  Matrix A:" );

  b = r8mat_inverse_4d ( a );

  r8mat_print ( 4, 4, b, "  Inverse matrix B:" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        c[i+j*N] = c[i+j*N] + a[i+k*N] * b[k+j*N];
      }
    }
  }

  r8mat_print ( 4, 4, c, "  C = A * B:" );

  free ( b );

  return;

# undef N
}
/******************************************************************************/

void r8mat_jac_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_JAC_TEST tests R8MAT_JAC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N 4

  double eps = 0.00001;
  double *fprime;
  int m = 3;
  double x[N] = { 1.0, 2.0, 3.0, 4.0 };

  printf ( "\n" );
  printf ( "R8MAT_JAC_TEST\n" );
  printf ( "  R8MAT_JAC estimates the M by N jacobian matrix\n" );
  printf ( "  of a nonlinear function.\n" );

  fprime = r8mat_jac ( m, N, eps, r8mat_jac_f, x );

  r8mat_print ( m, N, fprime, "  Estimated jacobian:" );

  free ( fprime );

  fprime = r8mat_jac_exact ( m, N, x );

  r8mat_print (  m, N, fprime, "  Exact jacobian:" );

  free ( fprime );

  return;
# undef N
}
/******************************************************************************/

double *r8mat_jac_f ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_JAC_F is a sample nonlinear function for treatment by R8MAT_JAC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of functions.

    Input, int N, the number of parameters.

    Input, double X[N], the parameter values.

    Output, double R8MAT_JAC_F[M], the function values.
*/
{
  double *f;

  f = ( double * ) malloc ( m * sizeof ( double ) );

  f[0] = sin ( x[0] * x[1] );
  f[1] = sqrt ( 1.0 + x[0] * x[0] ) + x[2];
  f[2] = x[0] + 2.0 * x[1] + 3.0 * x[2] + 4.0 * x[3];

  return f;
}
/******************************************************************************/

double *r8mat_jac_exact ( int m, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_JAC_EXACT is the exact jacobian of R8MAT_F.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of functions.

    Input, int N, the number of parameters.

    Input, double X[N], the parameter values.

    Output, double R8MAT_JAC_EXACT[M*N], the jacobian values.
*/
{
  double *fprime;

  fprime = ( double * ) malloc ( m * n * sizeof ( double ) );

  fprime[0+0*3] = cos ( x[0] * x[1] ) * x[1];
  fprime[0+1*3] = cos ( x[0] * x[1] ) * x[0];
  fprime[0+2*3] = 0.0;
  fprime[0+3*3] = 0.0;

  fprime[1+0*3] = x[0] / sqrt ( 1.0 + x[0] * x[0] );
  fprime[1+1*3] = 0.0;
  fprime[1+2*3] = 1.0;
  fprime[1+3*3] = 0.0;

  fprime[2+0*3] = 1.0;
  fprime[2+1*3] = 2.0;
  fprime[2+2*3] = 3.0;
  fprime[2+3*3] = 4.0;

  return fprime;
}
/******************************************************************************/

void r8mat_kronecker_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_KRONECKER_TEST tests R8MAT_KRONECKER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2013

  Author:

    John Burkardt
*/
{
  int m1 = 2;
  int m2 = 3;
  int m = m1 * m2;
  int n1 = 3;
  int n2 = 2;
  int n = n1 * n2;

  double a[2*3] = {
    1.0, 4.0, 
    2.0, 5.0, 
    3.0, 6.0 };
  double b[3*2] = {
    7.0,  9.0, 11.0, 
    8.0, 10.0, 12.0 };
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_KRONECKER_TEST\n" );
  printf ( "  R8MAT_KRONECKER computes the Kronecker product\n" );
  printf ( "  of two matrices.\n" );

  r8mat_print ( m1, n1, a, "  Factor matrix A:" );
  r8mat_print ( m2, n2, b, "  Factor matrix B:" );

  c = r8mat_kronecker ( m1, n1, a, m2, n2, b );

  r8mat_print ( m, n, c, "  Kronecker product C = kron(A,B)" );

  free ( c );

  return;
}
/******************************************************************************/

void r8mat_l_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L_INVERSE_TEST tests R8MAT_L_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N 4
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N*N] = {
    1.0, 2.0, 4.0,  7.0,
    0.0, 3.0, 5.0,  8.0,
    0.0, 0.0, 6.0,  9.0,
    0.0, 0.0, 0.0, 10.0 };
  double *b;
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_L_INVERSE_TEST\n" );
  printf ( "  R8MAT_L_INVERSE inverts a lower triangular matrix.\n" );

  r8mat_print ( N, N, a, "  Matrix A to be inverted:" );

  b = r8mat_l_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8mat_l_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L_PRINT_TEST tests R8MAT_L_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double a1[28] = {
    11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          22.0, 32.0, 42.0, 52.0, 62.0, 72.0,
                33.0, 43.0, 53.0, 63.0, 73.0,
                      44.0, 54.0, 64.0, 74.0,
                            55.0, 65.0, 75.0,
                                  66.0, 76.0,
                                        77.0 };
  double a2[18] = {
    11.0, 21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          22.0, 32.0, 42.0, 52.0, 62.0, 72.0,
                33.0, 43.0, 53.0, 63.0, 73.0 };
  double a3[10] = {
    11.0, 21.0, 31.0, 41.0,
          22.0, 32.0, 42.0,
                33.0, 43.0,
                      44.0 };
  int m1 = 7;
  int m2 = 7;
  int m3 = 4;
  int n1 = 7;
  int n2 = 3;
  int n3 = 7;

  printf ( "\n" );
  printf ( "R8MAT_L_PRINT_TEST\n" );
  printf ( "  R8MAT_L_PRINT prints a lower triangular matrix\n" );
  printf ( "  stored compactly.  Only the (possibly) nonzero\n" );
  printf ( "  elements are printed.\n" );

  r8mat_l_print ( m1, n1, a1, "  A 7 by 7 matrix." );

  r8mat_l_print ( m2, n2, a2, "  A 7 by 3 matrix." );

  r8mat_l_print ( m3, n3, a3, "  A 4 by 7 matrix." );

  return;
}
/******************************************************************************/

void r8mat_l1_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L1_INVERSE_TEST tests R8MAT_L1_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N 6
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N*N] = {
     1.0, 2.0, 0.0, 5.0, 0.0, 75.0,
     0.0, 1.0, 0.0, 0.0, 0.0,  0.0,
     0.0, 0.0, 1.0, 3.0, 0.0,  0.0,
     0.0, 0.0, 0.0, 1.0, 0.0,  6.0,
     0.0, 0.0, 0.0, 0.0, 1.0,  4.0,
     0.0, 0.0, 0.0, 0.0, 0.0,  1.0 };
  double *b;
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_L1_INVERSE_TEST\n" );
  printf ( "  R8MAT_L1_INVERSE inverts a unit lower triangular matrix.\n" );

  r8mat_print ( N, N, a, "  Matrix A to be inverted:" );

  b = r8mat_l1_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8mat_lu_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_LU_TEST tests R8MAT_LU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define M 5
# define N 5

  double *a;
  double l[M*M];
  double *lu;
  double p[M*M];
  double *plu;
  double u[M*N];
  double x[N] = { 1.0, 10.0, 4.0, 2.0, 3.0 };

  printf ( "\n" );
  printf ( "R8MAT_LU_TEST\n" );
  printf ( "  R8MAT_LU computes the LU factors of a matrix.\n" );

  a = r8mat_vand2 ( N, x );

  r8mat_print ( M, N, a, "  Matrix to be factored:" );

  r8mat_lu ( M, N, a, l, p, u );

  r8mat_print ( M, M, p, "  P factor:" );

  r8mat_print ( M, M, l, "  L factor:" );

  r8mat_print ( M, N, u, "  U factor:" );

  lu = r8mat_mm_new ( M, M, N, l, u );

  plu = r8mat_mm_new ( M, M, N, p, lu );

  r8mat_print ( M, N, plu, "  P*L*U:" );

  free ( a );
  free ( lu );
  free ( plu );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_max_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAX_TEST tests R8MAT_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp;

  printf ( "\n" );
  printf ( "R8MAT_MAX_TEST\n" );
  printf ( "  For a real matrix,\n" );
  printf ( "  R8MAT_MAX computes the maximum value;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, &seed );

  r8mat_print ( m, n, a, "  Random array:" );

  temp = r8mat_max ( m, n, a );

  printf ( "\n" );
  printf ( "  Maximum value = %g\n", temp );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_max_index_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAX_INDEX_TEST tests R8MAT_MAX_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define M 5
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int j;
  int seed;

  printf ( "\n" );
  printf ( "R8MAT_MAX_INDEX\n" );
  printf ( "  R8MAT_MAX_INDEX locates the maximum entry of an R8MAT;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, &seed );

  r8mat_print ( M, N, a, "  Random array:" );

  r8mat_max_index ( M, N, a, &i, &j );

  printf ( "\n" );
  printf ( "  Maximum I,J indices            %d  %d\n", i, j );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_maxcol_minrow_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAXCOL_MINROW_TEST tests R8MAT_MAXCOL_MINROW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  printf ( "\n" );
  printf ( "R8MAT_MAXCOL_MINROW_TEST\n" );
  printf ( "  R8MAT_MAXCOL_MINROW computes the maximum over\n" );
  printf ( "  columns of the mininum over rows;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, &seed );

  r8mat_print ( m, n, a, "  Random array:" );

  printf ( "  MAXCOL_MINROW = %g\n", r8mat_maxcol_minrow ( m, n, a ) );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_maxrow_mincol_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAXROW_MINCOL_TEST tests R8MAT_MAXROW_MINCOL_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  printf ( "\n" );
  printf ( "R8MAT_MAXROW_MINCOL_TEST\n" );
  printf ( "  R8MAT_MAXROW_MINCOL computes the maximum over\n" );
  printf ( "  rows of the mininum over columns;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, &seed );

  r8mat_print ( m, n, a, "  Random array:" );

  printf ( "  MAXROW_MINCOL = %g\n", r8mat_maxrow_mincol ( m, n, a ) );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_min_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MIN_TEST tests R8MAT_MAX and R8MAT_MIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp;

  printf ( "\n" );
  printf ( "R8MAT_MIN_TEST\n" );
  printf ( "  For a real matrix,\n" );
  printf ( "  R8MAT_MIN computes the minimum value;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, &seed );

  r8mat_print ( m, n, a, "  Random array:" );

  temp = r8mat_min ( m, n, a );

  printf ( "\n" );
  printf ( "  Minimum value = %g\n", temp );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_min_index_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MIN_INDEX_TEST tests R8MAT_MIN_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
# define M 5
# define N 3

  double *a;
  double b = 0.0;
  double c = 10.0;
  int i;
  int j;
  int seed;

  printf ( "\n" );
  printf ( "R8MAT_MIN_INDEX_TEST\n" );
  printf ( "  R8MAT_MIN_INDEX locates the minimum entry of an R8MAT;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( M, N, b, c, &seed );

  r8mat_print ( M, N, a, "  Random array:" );

  printf ( "\n" );
  r8mat_min_index ( M, N, a, &i, &j );
  printf ( "  Minimum I,J indices            %d  %d\n", i, j );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_mincol_maxrow_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINCOL_MAXROW_TEST tests R8MAT_MINCOL_MAXROW_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  printf ( "\n" );
  printf ( "R8MAT_MINCOL_MAXROW_TEST\n" );
  printf ( "  R8MAT_MINCOL_MAXROW computes the minimum over\n" );
  printf ( "  columns of the maxinum over rows;\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, &seed );

  r8mat_print ( m, n, a, "  Random array:" );

  printf ( "  MINCOL_MAXROW = %g\n", r8mat_mincol_maxrow ( m, n, a ) );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_minrow_maxcol_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINROW_MAXCOL_TEST tests R8MAT_MINROW_MAXCOL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 10.0;
  int m = 5;
  int n = 3;
  int seed;
  double temp1;
  double temp2;

  printf ( "\n" );
  printf ( "R8MAT_MINROW_MAXCOL_TEST\n" );
  printf ( "  R8MAT_MINROW_MAXCOL computes the minimum over\n" );
  printf ( "  rows of the maxinum over columns;\n" );
  printf ( "\n" );

  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, b, c, &seed );

  r8mat_print ( m, n, a, "  Random array:" );

  printf ( "  MINROW_MAXCOL = %g\n", r8mat_minrow_maxcol ( m, n, a ) );
 
  free ( a );

  return;
}
/******************************************************************************/

void r8mat_mm_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM_NEW_TEST tests R8MAT_MM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N1 2
# define N2 3
# define N3 4
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2*N3] = {
     1.0,  2.0,  3.0,
     4.0,  5.0,  6.0,
     0.0,  0.0,  1.0,
    -1.0,  2.0, -1.0 };
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_MM_NEW_TEST\n" );
  printf ( "  R8MAT_MM_NEW multiplies two (rectangular) matrices\n" );
  printf ( "  and returns the result as the function value.\n" );

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8mat_print ( N2, N3, b, "  Matrix B:" );

  c = r8mat_mm_new ( N1, N2, N3, a, b );

  r8mat_print ( N1, N3, c, "  Product C = A * B:" );

  free ( c );

  return;
# undef N1
# undef N2
# undef N3
}
/******************************************************************************/

void r8mat_mm_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM_TEST tests R8MAT_MM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N1 2
# define N2 3
# define N3 4
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2*N3] = {
     1.0,  2.0,  3.0,
     4.0,  5.0,  6.0,
     0.0,  0.0,  1.0,
    -1.0,  2.0, -1.0 };
  double c[N1*N3];

  printf ( "\n" );
  printf ( "R8MAT_MM_TEST\n" );
  printf ( "  R8MAT_MM multiplies two (rectangular) matrices\n" );
  printf ( "  and returns the result as an argument.\n" );

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8mat_print ( N2, N3, b, "  Matrix B:" );

  r8mat_mm ( N1, N2, N3, a, b, c );

  r8mat_print ( N1, N3, c, "  Product C = A * B:" );

  return;
# undef N1
# undef N2
# undef N3
}
/******************************************************************************/

void r8mat_mv_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MV_NEW_TEST tests R8MAT_MV_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N1 2
# define N2 3
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2] = {
     1.0,  2.0,  3.0 };
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_MV_NEW_TEST\n" );
  printf ( "  R8MAT_MV_NEW multiplies a (rectangular) matrix times a vector,\n" );
  printf ( "  and returns the result as the function value.\n" );

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N2, b, "  Vector B:" );

  c = r8mat_mv_new ( N1, N2, a, b );

  r8vec_print ( N1, c, "  Product C = A * B:" );

  free ( c );

  return;
# undef N1
# undef N2
}
/******************************************************************************/

void r8mat_mv_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MV_TEST tests R8MAT_MV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N1 2
# define N2 3
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N2] = {
     1.0,  2.0,  3.0 };
  double c[N1];

  printf ( "\n" );
  printf ( "R8MAT_MV_TEST\n" );
  printf ( "  R8MAT_MV multiplies a (rectangular) matrix times a vector,\n" );
  printf ( "  and returns the result as an argument.\n" );

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N2, b, "  Vector B:" );

  r8mat_mv ( N1, N2, a, b, c );

  r8vec_print ( N1, c, "  Product C = A * B:" );

  return;
# undef N1
# undef N2
}
/******************************************************************************/

void r8mat_mtv_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MTV_NEW_TEST tests R8MAT_MTV_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N1 2
# define N2 3
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N1] = {
     1.0,  2.0 };
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_MTV_NEW_TEST\n" );
  printf ( "  R8MAT_MTV_NEW multiplies a transposed matrix times a vector,\n" );
  printf ( "  and returns the result as the function value.\n" );

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N1, b, "  Vector B:" );

  c = r8mat_mtv_new ( N1, N2, a, b );

  r8vec_print ( N2, c, "  Product C = A' * B:" );

  free ( c );

  return;
# undef N1
# undef N2
}
/******************************************************************************/

void r8mat_mtv_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MTV_TEST tests R8MAT_MTV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2012

  Author:

    John Burkardt
*/
{
# define N1 2
# define N2 3
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N1*N2] = {
     1.0, 4.0,
     2.0, 5.0,
     3.0, 6.0 };
  double b[N1] = {
     1.0,  2.0 };
  double c[N2];

  printf ( "\n" );
  printf ( "R8MAT_MTV_TEST\n" );
  printf ( "  R8MAT_MTV multiplies a transposed matrix times a vector,\n" );
  printf ( "  and returns the result as an argument.\n" );

  r8mat_print ( N1, N2, a, "  Matrix A:" );

  r8vec_print ( N1, b, "  Vector B:" );

  r8mat_mtv ( N1, N2, a, b, c );

  r8vec_print ( N2, c, "  Product C = A' * B:" );

  return;
# undef N1
# undef N2
}
/******************************************************************************/

void r8mat_nint_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NINT_TEST tests R8MAT_NINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int m;
  int n;
  int seed;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "R8MAT_NINT_TEST\n" );
  printf ( "  R8MAT_NINT rounds an R8MAT.\n" );

  m = 5;
  n = 4;
  x1 = -5.0;
  x2 = +5.0;
  seed = 123456789;
  a = r8mat_uniform_ab_new ( m, n, x1, x2, &seed );
  r8mat_print ( m, n, a, "  Matrix A:" );
  r8mat_nint ( m, n, a );
  r8mat_print ( m, n, a, "  Rounded matrix A:" );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_nonzeros_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NONZEROS_TEST tests R8MAT_NONZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int j;
  int m = 5;
  int n = 4;
  int c1;
  int c2;

  printf ( "\n" );
  printf ( "R8MAT_NONZEROS_TEST\n" );
  printf ( "  R8MAT_NONZEROS counts nonzeros in an R8MAT.\n" );

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  c1 = 0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( ( i % 2 ) == 0 && ( j % 2 ) == 0 )
      {
        a[i+j*m] = 1;
        c1 = c1 + 1;
      }
      else
      {
        a[i+j*m] = 0;
      }
    }
  }

  r8mat_print ( m, n, a, "  Matrix A:" );

  c2 = r8mat_nonzeros ( m, n, a );

  printf ( "\n" );
  printf ( "  Expected nonzeros = %d\n", c1 );
  printf ( "  Computed nonzeros = %d\n", c2 );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_norm_fro_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_FRO_TEST tests R8MAT_NORM_FRO.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int j;
  int k;
  int m = 5;
  int n = 4;
  double t1;
  double t2;

  printf ( "\n" );
  printf ( "R8MAT_NORM_FRO_TEST\n" );
  printf ( "  R8MAT_NORM_FRO computes the Frobenius norm of a matrix.\n" );

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  k = 0;
  t1 = 0.0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      k = k + 1;
      a[i+j*m] = ( double ) ( k );
      t1 = t1 + k * k;
    }
  }
  t1 = sqrt ( t1 );

  r8mat_print ( m, n, a, "  Matrix A:" );

  t2 = r8mat_norm_fro ( m, n, a );

  printf ( "\n" );
  printf ( "  Expected Frobenius norm = %g\n", t1 );
  printf ( "  Computed Frobenius norm = %g\n", t2 );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_norm_l1_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_L1_TEST tests R8MAT_NORM_L1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int m;
  int n;
  int seed;
  double t;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "R8MAT_NORM_L1_TEST\n" );
  printf ( "  R8MAT_NORM_L1 computes the L1 norm of a matrix.\n" );

  m = 5;
  n = 4;
  x1 = -5.0;
  x2 = +5.0;
  seed = 123456789;

  a = r8mat_uniform_ab_new ( m, n, x1, x2, &seed );
  r8mat_nint ( m, n, a );

  r8mat_print ( m, n, a, "  Matrix A:" );

  t = r8mat_norm_l1 ( m, n, a );

  printf ( "\n" );
  printf ( "  L1 norm = %g\n", t );

  free ( a );

  return;
}
/******************************************************************************/

void r8mat_nullspace_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NULLSPACE_TEST tests R8MAT_NULLSPACE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  double *ax;
  int m = M;
  int n = N;
  double *nullspace;
  int nullspace_size;

  printf ( "\n" );
  printf ( "R8MAT_NULLSPACE_TEST\n" );
  printf ( "  R8MAT_NULLSPACE computes the nullspace of a matrix.\n" );

  r8mat_print ( m, n, a, "  Input A:" );

  nullspace_size = r8mat_nullspace_size ( m, n, a );

  printf ( "\n" );
  printf ( "  Nullspace size is %d\n", nullspace_size );

  nullspace = r8mat_nullspace ( m, n, a, nullspace_size );

  r8mat_print ( n, nullspace_size, nullspace, "  Nullspace vectors:" );

  ax = r8mat_mm_new ( m, n, nullspace_size, a, nullspace );

  r8mat_print ( m, nullspace_size, ax, "  Product A * Nullspace vectors:" );

  free ( ax );
  free ( nullspace );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_nullspace_size_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NULLSPACE_SIZE_TEST tests R8MAT_NULLSPACE_SIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;
  int nullspace_size;

  printf ( "\n" );
  printf ( "R8MAT_NULLSPACE_SIZE_TEST\n" );
  printf ( "  R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.\n" );

  r8mat_print ( m, n, a, "  Input A:" );

  nullspace_size = r8mat_nullspace_size ( m, n, a );

  printf ( "\n" );
  printf ( "  Nullspace size is %d\n", nullspace_size );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_orth_uniform_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ORTH_UNIFORM_NEW_TEST tests R8MAT_ORTH_UNIFORM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *at;
  double *ata;
  int n = 5;
  int seed;

  printf ( "\n" );
  printf ( "R8MAT_ORTH_UNIFORM_NEW_TEST\n" );
  printf ( "  R8MAT_ORTH_UNIFORM_NEW computes a random orthogonal matrix.\n" );

  seed = 123456789;

  a = r8mat_orth_uniform_new ( n, &seed );

  r8mat_print ( n, n, a, "  Random orthogonal matrix A" );

  at = r8mat_transpose_new ( n, n, a );

  ata = r8mat_mm_new ( n, n, n, at, a );

  r8mat_print ( n, n, ata, "  AT*A" );

  free ( a );
  free ( at );
  free ( ata );

  return;
}
/******************************************************************************/

void r8mat_plot_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PLOT_TEST tests R8MAT_PLOT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define M 10
# define N 100

  double a[M*N];
  int i;
  int im1;
  int ip1;

  r8mat_zero ( M, N, a );

  for ( i = 0; i < M; i++ )
  {
    a[i+i*M] = -2.0;

    if ( i+1 < N )
    {
      ip1 = i+1;
    }
    else
    {
      ip1 = 0;
    }

    a[i+ip1*M] = 1.0;

    if ( 0 <= i-1 )
    {
      im1 = i-1;
    }
    else
    {
      im1 = N-1;
    }
    a[i+im1*M] = 1.0;
  }

  printf ( "\n" );
  printf ( "R8MAT_PLOT_TEST\n" );
  printf ( "  R8MAT_PLOT prints a symbolic picture of a matrix.\n" );
  printf ( "  Typically,\n" );
  printf ( "\n" );
  printf ( "    - for negative, \n" );
  printf ( "    0 for zero, and\n" );
  printf ( "    + for positive entries\n" );
  printf ( "\n" );
  printf ( "  or\n" );
  printf ( "\n" );
  printf ( "    X for nonzero and\n" );
  printf ( "    0 for zero.\n" );
  printf ( "\n" );

  r8mat_plot ( M, N, a, "  A plot of the matrix:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_power_method_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_POWER_METHOD_TEST tests R8MAT_POWER_METHOD.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double a[N*N];
  double *av;
  int i;
  int j;
  double r;
  double v[N];

  printf ( "\n" );
  printf ( "R8MAT_POWER_METHOD_TEST\n" );
  printf ( "  R8MAT_POWER_METHOD applies the power method\n" );
  printf ( "  to a matrix.\n" );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( j == i - 1 || j == i + 1 )
      {
        a[i+j*N] = -1.0;
      }
      else if ( j == i )
      {
        a[i+j*N] = 2.0;
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }
  r8vec_zero ( N, v );

  r8mat_power_method ( N, a, &r, v );

  printf ( "\n" );
  printf ( "  Estimated eigenvalue = %g\n", r );

  r8vec_print ( N, v, "  Estimated eigenvector V:" );

  av = r8mat_mv_new ( N, N, a, v );

  r8vec_print ( N, av, "  Value of A*V:" );

  free ( av );

  return;
# undef N
}
/******************************************************************************/

void r8mat_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_TEST tests R8MAT_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8MAT_PRINT_TEST\n" );
  printf ( "  R8MAT_PRINT prints an R8MAT.\n" );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print ( m, n, a, "  The matrix:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_print_some_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
# define M 6
# define N 4

  double a[M*N];
  int i;
  int j;
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8MAT_PRINT_SOME_TEST\n" );
  printf ( "  R8MAT_PRINT_SOME prints some of an R8MAT.\n" );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( ( i + 1 ) * 10 + ( j + 1 ) );
    }
  }
  r8mat_print_some ( m, n, a, 2, 1, 4, 2, "  Rows 2:4, Cols 1:2:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_ref_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_REF_TEST tests R8MAT_REF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8MAT_REF_TEST\n" );
  printf ( "  R8MAT_REF computes the row echelon form of a matrix.\n" );

  r8mat_print ( m, n, a, "  Input A:" );

  r8mat_ref ( m, n, a );

  r8mat_print ( m, n, a, "  REF form:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_rref_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_RREF_TEST tests R8MAT_RREF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define M 4
# define N 7

  double a[M*N] = {
    1.0, -2.0, 3.0, -1.0,
    3.0, -6.0, 9.0, -3.0,
    0.0,  0.0, 0.0,  0.0,
    2.0, -2.0, 0.0,  1.0,
    6.0, -8.0, 6.0,  0.0,
    3.0,  3.0, 6.0,  9.0,
    1.0,  1.0, 2.0,  3.0 };
  int m = M;
  int n = N;

  printf ( "\n" );
  printf ( "R8MAT_RREF_TEST\n" );
  printf ( "  R8MAT_RREF computes the reduced row echelon form of a matrix.\n" );

  r8mat_print ( m, n, a, "  Input A:" );

  r8mat_rref ( m, n, a );

  r8mat_print ( m, n, a, "  REF form:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_solve_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE_TEST tests R8MAT_SOLVE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 3
# define RHS_NUM 2

  double a[N*(N+RHS_NUM)] = {
     1.0,  4.0,  7.0,
     2.0,  5.0,  8.0,
     3.0,  6.0,  0.0,
    14.0, 32.0, 23.0,
     7.0, 16.0,  7.0 };
  int i;
  int info;
  int j;

  printf ( "\n" );
  printf ( "R8MAT_SOLVE_TEST\n" );
  printf ( "  R8MAT_SOLVE solves linear systems.\n" );
/*
  Print out the matrix to be inverted.
*/
  r8mat_print ( N, N+RHS_NUM, a, "  The linear system:" );
/*
  Solve the systems.
*/
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "  The input matrix was singular.\n" );
    printf ( "  The solutions could not be computed.\n" );
    return;
  }

  printf ( "\n" );
  printf ( "  The computed solutions:\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    for ( j = N; j < N+RHS_NUM; j++ )
    {
      printf ( "%10g  \n", a[i+j*N] );
    }
    printf ( "\n" );
  }

  return;
# undef N
# undef RHS_NUM
}
/******************************************************************************/

void r8mat_solve_2d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE_2D_TEST tests R8MAT_SOLVE_2D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 2;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  printf ( "\n" );
  printf ( "R8MAT_SOLVE_2D_TEST\n" );
  printf ( "  R8MAT_SOLVE_2D solves 2D linear systems.\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, &seed );
    x = r8vec_uniform_01_new ( n, &seed );
    b = r8mat_mv_new ( n, n, a, x );

    x2 = r8mat_solve_2d ( a, b, &det );

    printf ( "\n" );
    printf ( "  Solution / Computed:\n" );
    printf ( "\n" );

    for ( i = 0; i < n; i++ )
    {
      printf ( "  %14g  %14g\n", x[i], x2[i] );
    }

    free ( a );
    free ( b );
    free ( x );
    free ( x2 );
  }

  return;
}
/******************************************************************************/

void r8mat_solve_3d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE_3D_TEST tests R8MAT_SOLVE_3D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double det;
  int i;
  int n = 3;
  int seed;
  int test;
  int test_num = 5;
  double *x;
  double *x2;

  printf ( "\n" );
  printf ( "R8MAT_SOLVE_3D_TEST\n" );
  printf ( "  R8MAT_SOLVE_3D solves 3D linear systems.\n" );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8mat_uniform_01_new ( n, n, &seed );
    x = r8vec_uniform_01_new ( n, &seed );
    b = r8mat_mv_new ( n, n, a, x );

    x2 = r8mat_solve_3d ( a, b, &det );

    printf ( "\n" );
    printf ( "  Solution / Computed:\n" );
    printf ( "\n" );

    for ( i = 0; i < n; i++ )
    {
      printf ( "  %14g  %14g\n", x[i], x2[i] );
    }

    free ( a );
    free ( b );
    free ( x );
    free ( x2 );
  }

  return;
}
/******************************************************************************/

void r8mat_solve2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE2_TEST tests R8MAT_SOLVE2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 4

  double *a;
  double a1[2*2] = {
    1.0, 3.0,
    2.0, 4.0 };
  double a2[3*3] = {
    2.0, 1.0, 1.0,
    1.0, 1.0, 0.0,
    1.0, 0.0, 1.0 };
  double a3[4*4] = {
    1.0, 2.0, 1.0, 3.0,
    0.0, 1.0, 2.0, 1.0,
    0.0, 0.0, 3.0, 2.0,
    1.0, 3.0, 0.0, 1.0 };
  double a4[3*3] = {
    2.0, 1.0, 3.0,
    4.0, 2.0, 6.0,
    1.0, 4.0, 5.0 };
  double *b;
  double b1[2] = { 5.0, 11.0 };
  double b2[3] = { 4.0, 2.0, 2.0 };
  double b3[4] = { 5.0, 11.0, 16.0, 15.0 };
  double b4[3] = { 13.0, 17.0, 20.0 };
  int ierror;
  int n;
  int n_test[TEST_NUM] = { 2, 3, 4, 3 };
  int test;
  double *x;

  printf ( "\n" );
  printf ( "R8MAT_SOLVE2_TEST\n" );
  printf ( "  R8MAT_SOLVE2 is a linear solver.\n" );
  printf ( "\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = n_test[test];

    if ( test == 0 )
    {
      a = a1;
      b = b1;
    }
    else if ( test == 1 )
    {
      a = a2;
      b = b2;
    }
    else if ( test == 2 )
    {
      a = a3;
      b = b3;
    }
    else if ( test == 3 )
    {
      a = a4;
      b = b4;
    }

    r8vec_print ( n, b, "  Right hand side:" );

    x = r8mat_solve2 ( n, a, b, &ierror );

    printf ( "\n" );
    if ( ierror == 0 )
    {
      printf ( "  The system is nonsingular.\n" );
    }
    else if ( ierror == 1 )
    {
      printf ( "  The system is singular, but consistent.\n" );
    }
    else if ( ierror == 2 )
    {
      printf ( "  The system is singular and inconsistent.\n" );
    }

    r8vec_print ( n, x, "  Computed solution:" );

    free ( x );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8mat_sub_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SUB_NEW_TEST tests R8MAT_SUB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2015

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *c;
  int m = 4;
  int n = 4;

  printf ( "\n" );
  printf ( "R8MAT_SUB_NEW_TEST\n" );
  printf ( "  R8MAT_SUB_NEW computes C = A - B for R8MAT's\n" );

  a = r8mat_indicator_new ( m, n );

  b = r8mat_transpose_new ( m, n, a );

  c = r8mat_sub_new ( m, n, a, b );

  r8mat_print ( m, n, a, "  A:" );
  r8mat_print ( m, n, b, "  B:" );
  r8mat_print ( m, n, c, "  C = A-B:" );

  free ( a );
  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void r8mat_symm_jacobi_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SYMM_JACOBI_TEST tests R8MAT_SYMM_JACOBI;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int n = 5;
  double *q;
  int seed;
  double *x;

  printf ( "\n" );
  printf ( "R8MAT_SYMM_JACOBI_TEST\n" );
  printf ( "  For a symmetric R8MAT:\n" );
  printf ( "  R8MAT_SYMM_JACOBI diagonalizes;\n" );
/*
  Choose eigenvalues.
*/
  x = r8vec_indicator1_new ( n );
/*
  Choose eigenvectors.
*/
  seed = 123456789;

  q = r8mat_orth_uniform_new ( n, &seed );
/*
  Now get A = Q*X*Q.
*/
  a = r8mat_symm_eigen ( n, x, q );

  r8mat_print ( n, n, a, "  Matrix to diagonalize:" );

  r8mat_symm_jacobi ( n, a );

  printf ( "\n" );
  printf ( "  Computed Eigenvalues:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %12g\n", a[i+i*n] );
  }

  free ( a );
  free ( q );
  free ( x );

  return;
}
/******************************************************************************/

void r8mat_to_r8cmat_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TO_R8CMAT_NEW_TEST tests R8MAT_TO_R8CMAT_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  printf ( "\n" );
  printf ( "R8MAT_TO_R8CMAT_NEW_TEST\n" );
  printf ( "  R8MAT_TO_R8CMAT_NEW converts an R8MAT to an R8CMAT.\n" );
  printf ( "\n" );
  printf ( "  Data is of order (%d,%d).\n", m, n );
/*
  Set the R8MAT.
*/
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
/*
  Convert.
*/
  b = r8mat_to_r8cmat_new ( m, n, a );
  r8cmat_print ( m, n, b, "  The R8CMAT B:" );
/*
  Recover the matrix.
*/
  c = r8cmat_to_r8mat_new ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
/*
  Free memory.
*/
  free ( a );
  r8cmat_delete ( m, n, b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8mat_to_r8plu_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TO_R8PLU_TEST tests R8MAT_TO_R8PLU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double a2[N*N];
  double b = 0.0;
  double c = 1.0;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8MAT_TO_R8PLU_TEST\n" );
  printf ( "  R8MAT_TO_R8PLU determines the compressed PLU factors\n" );
  printf ( "  of a real general matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8mat_uniform_ab_new ( N, N, b, c, &seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
/*
  Factor the matrix.
*/
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "Warning!\n" );
    printf ( "  R8MAT_TO_R8PLU declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
  }
/*
  Display the gory details.
*/
  i4vec_print ( N, pivot, "  The pivot vector P:" );

  r8mat_print ( N, N, lu, "  The compressed LU factors:" );
/*
  Recover the matrix from the PLU factors.
*/
  r8plu_to_r8mat ( N, pivot, lu, a2 );

  r8mat_print ( N, N, a2, "  The recovered matrix A2:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8mat_to_r8rmat_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TO_R8RMAT_TEST tests R8MAT_TO_R8RMAT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  printf ( "\n" );
  printf ( "R8MAT_TO_R8RMAT_TEST\n" );
  printf ( "  R8MAT_TO_R8RMAT converts an R8MAT to an R8RMAT.\n" );
  printf ( "\n" );
  printf ( "  Data is of order (%d,%d).\n", m, n );
/*
  Set the R8MAT.
*/
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
/*
  Convert.
*/
  b = r8mat_to_r8rmat ( m, n, a );
  r8rmat_print ( m, n, b, "  The R8RMAT B:" );
/*
  Recover the matrix.
*/
  c = r8rmat_to_r8mat ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
/*
  Free memory.
*/
  free ( a );
  r8rmat_delete ( m, n, b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8mat_trace_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRACE_TEST tests R8MAT_TRACE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 4

  double a[N*N];
  int i;
  int j;
  double trace;

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i <= j )
      {
        a[i+j*N] = ( double ) ( N - j );
      }
      else if ( j == i - 1 )
      {
        a[i+j*N] = ( double ) ( N - j - 1 );
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  printf ( "\n" );
  printf ( "R8MAT_TRACE_TEST\n" );
  printf ( "  R8MAT_TRACE computes the trace of a matrix\n" );

  r8mat_print ( N, N, a, "  Matrix:" );

  trace = r8mat_trace ( N, a );

  printf ( "\n" );
  printf ( "  Trace is %g\n", trace );

  return;
# undef N
}
/******************************************************************************/

void r8mat_transpose_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_NEW_TEST tests R8MAT_TRANSPOSE_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2015

  Author:

    John Burkardt
*/
{
  double *a;
  double *at;
  int m = 5;
  int n = 4;

  printf ( "\n" );
  printf ( "R8MAT_TRANSPOSE_NEW_TEST\n" );
  printf ( "  R8MAT_TRANSPOSE_NEW transposes an R8MAT.\n" );

  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  Matrix A:" );

  at = r8mat_transpose_new ( m, n, a );
  r8mat_print ( n, m, at, "  Transposed matrix At:" );

  free ( a );
  free ( at );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_TEST tests R8MAT_TRANSPOSE_PRINT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define M 7
# define N 12

  double a[M*N];
  int i;
  int j;

  printf ( "\n" );
  printf ( "R8MAT_TRANSPOSE_PRINT_TEST\n" );
  printf ( "  R8MAT_TRANSPOSE_PRINT prints an R8MAT,\n" );
  printf ( "  transposed.\n" );
  printf ( "\n" );
  printf ( "  Matrix row order M =    %d\n", M );
  printf ( "  Matrix column order N = %d\n", N );
/*
  Set the matrix.
*/
  for ( i = 1; i <= M; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      a[i-1+(j-1)*M] = ( double ) ( i * 100 + j );
    }
  }

  r8mat_transpose_print ( M, N, a, "  The transposed matrix A:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8mat_u_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_U_INVERSE_TEST tests R8MAT_U_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 4
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N*N] = {
    1.0, 0.0, 0.0,  0.0,
    2.0, 3.0, 0.0,  0.0,
    4.0, 5.0, 6.0,  0.0,
    7.0, 8.0, 9.0, 10.0 };
  double *b;
  double *c;
  int i;

  printf ( "\n" );
  printf ( "R8MAT_U_INVERSE_TEST\n" );
  printf ( "  R8MAT_U_INVERSE inverts an upper triangular matrix.\n" );

  r8mat_print ( N, N, a, "  Input matrix A" );

  b = r8mat_u_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8mat_u1_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_U1_INVERSE_TEST tests R8MAT_U1_INVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 6
/*
  Each row of this definition is a COLUMN of the matrix.
*/
  double a[N*N] = {
    1.0, 0.0, 0.0, 0.0, 0.0,  0.0,
    2.0, 1.0, 0.0, 0.0, 0.0,  0.0,
    0.0, 0.0, 1.0, 0.0, 0.0,  0.0,
    5.0, 0.0, 3.0, 1.0, 0.0,  0.0,
    0.0, 0.0, 0.0, 0.0, 1.0,  0.0,
   75.0, 0.0, 0.0, 6.0, 4.0,  1.0 };
  double *b;
  double *c;

  printf ( "\n" );
  printf ( "R8MAT_U1_INVERSE_TEST\n" );
  printf ( "  R8MAT_U1_INVERSE inverts a unit upper triangular matrix.\n" );

  r8mat_print ( N, N, a, "  Input matrix A" );

  b = r8mat_u1_inverse ( N, a );

  r8mat_print ( N, N, b, "  Inverse matrix B:" );

  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  free ( b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8mat_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB_NEW_TEST tests R8MAT_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2005

  Author:

    John Burkardt
*/
{
# define M 5
# define N 4

  double *a;
  double b = 2.0E+00;
  double c = 10.0E+00;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8MAT_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R8MAT_UNIFORM_AB_NEW sets an R8MAT to random values in [A,B].\n" );
  printf ( "\n" );

  a = r8mat_uniform_ab_new ( M, N, b, c, &seed );

  r8mat_print ( M, N, a, "  The random R8MAT:" );

  free ( a );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8plu_det_test ( )

/******************************************************************************/
/*
  Purpose:

    R8PLU_DET_TEST tests R8PLU_DET;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double b = 0.0;
  double c = 1.0;
  double det;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8PLU_DET_TEST\n" );
  printf ( "  R8PLU_DET determines the determinant of a matrix from its\n" );
  printf ( "  compressed PLU factors.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8mat_uniform_ab_new ( N, N, b, c, &seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
/*
  Factor the matrix.
*/
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "Warning!\n" );
    printf ( "  R8MAT_TO_R8PLU declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
  }
/*
  Compute the determinant.
*/
  det = r8plu_det ( N, pivot, lu );

  printf ( "\n" );
  printf ( "  The determinant = %g\n", det );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8plu_inverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8PLU_INVERSE_TEST tests R8PLU_INVERSE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double b[N*N];
  double *c;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8PLU_INVERSE_TEST\n" );
  printf ( "  R8PLU_INVERSE determines the inverse of a matrix from its\n" );
  printf ( "  compressed PLU factors.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8mat_uniform_01_new ( N, N, &seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
/*
  Factor the matrix.
*/
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "Warning!\n" );
    printf ( "  R8MAT_TO_R8PLU declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
  }
/*
  Compute the inverse.
*/
  r8plu_inverse ( N, pivot, lu, b );

  r8mat_print ( N, N, b, "  The inverse B:" );
/*
  Compute the product C = A * B.
*/
  c = r8mat_mm_new ( N, N, N, a, b );

  r8mat_print ( N, N, c, "  Product C = A * B:" );

  free ( a );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8plu_mul_test ( )

/******************************************************************************/
/*
  Purpose:

    R8PLU_MUL_TEST tests R8PLU_MUL;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;
  double x[N];

  printf ( "\n" );
  printf ( "R8PLU_MUL_TEST\n" );
  printf ( "  R8PLU_MUL computes the product A*x\n" );
  printf ( "  using the compressed PLU factors of A.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8mat_uniform_01_new ( N, N, &seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
/*
  Set the right hand side B1.
*/
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) (i+1);
  }

  b = r8mat_mv_new ( N, N, a, x );

  r8vec_print ( N, b, "  The right hand side B (computed from A):" );
/*
  Factor the matrix.
*/
  info = r8mat_to_r8plu ( N, a, pivot, lu );
/*
  Solve the system.
*/
  r8plu_mul ( N, pivot, lu, x, b );

  r8vec_print ( N, b, "  The right hand side B (computed from PLU):" );

  free ( a );
  free ( b );

  return;
# undef N
}
/******************************************************************************/

void r8plu_sol_test ( )

/******************************************************************************/
/*
  Purpose:

    R8PLU_SOL_TEST tests R8PLU_SOL;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;
  double x[N];

  printf ( "\n" );
  printf ( "R8PLU_SOL_TEST\n" );
  printf ( "  R8PLU_SOL solves the linear system A*x=b\n" );
  printf ( "  using the compressed PLU factors of A.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8mat_uniform_01_new ( N, N, &seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
/*
  Set the desired solution.
*/
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) (i+1);
  }
/*
  Set the right hand side.
*/
  b = r8mat_mv_new ( N, N, a, x );

  r8vec_print ( N, b, "  The right hand side B (computed from A):" );
/*
  Destroy the desired solution (no cheating now!)
*/
  r8vec_zero ( N, x );
/*
  Factor the matrix.
*/
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "Fatal error!\n" );
    printf ( "  R8MAT_TO_R8PLU declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
    return;
  }
/*
  Solve the system.
*/
  r8plu_sol ( N, pivot, lu, b, x );

  r8vec_print ( N, x, "  The computed solution X:" );

  free ( a );
  free ( b );

  return;
# undef N
}
/******************************************************************************/

void r8plu_to_r8mat_test ( )

/******************************************************************************/
/*
  Purpose:

    R8PLU_TO_R8MAT_TEST tests R8MAT_TO_R8PLU and R8PLU_TO_R8MAT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;
  double a2[N*N];
  double b = 0.0;
  double c = 1.0;
  int info;
  double lu[N*N];
  int pivot[N];
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8PLU_TO_R8MAT_TEST\n" );
  printf ( "  R8PLU_TO_R8MAT determines the original matrix from\n" );
  printf ( "  the compressed PLU factors.\n" );
  printf ( "\n" );
  printf ( "  Matrix order N = %d\n", N );
/*
  Set the matrix.
*/
  a = r8mat_uniform_ab_new ( N, N, b, c, &seed );

  r8mat_print ( N, N, a, "  The matrix A:" );
/*
  Factor the matrix.
*/
  info = r8mat_to_r8plu ( N, a, pivot, lu );

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "Warning!\n" );
    printf ( "  R8MAT_TO_R8PLU declares the matrix is singular!\n" );
    printf ( "  The value of INFO is %d\n", info );
  }
/*
  Display the gory details.
*/
  i4vec_print ( N, pivot, "  The pivot vector P:" );

  r8mat_print ( N, N, lu, "  The compressed LU factors:" );
/*
  Recover the matrix from the PLU factors.
*/
  r8plu_to_r8mat ( N, pivot, lu, a2 );

  r8mat_print ( N, N, a2, "  The recovered matrix A2:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8poly_degree_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_DEGREE_TEST tests R8POLY_DEGREE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2015

  Author:

    John Burkardt
*/
{
  double c1[4] = { 1.0, 2.0, 3.0, 4.0 }; 
  double c2[4] = { 1.0, 2.0, 3.0, 0.0 };
  double c3[4] = { 1.0, 2.0, 0.0, 4.0 };
  double c4[4] = { 1.0, 0.0, 0.0, 0.0 };
  double c5[4] = { 0.0, 0.0, 0.0, 0.0 };
  int d;
  int m;
 
  printf ( "\n" );
  printf ( "R8POLY_DEGREE_TEST\n" );
  printf ( "  R8POLY_DEGREE determines the degree of an R8POLY.\n" );

  m = 3;

  r8poly_print ( m, c1, "  The R8POLY:" );
  d = r8poly_degree ( m, c1 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c2, "  The R8POLY:" );
  d = r8poly_degree ( m, c2 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c3, "  The R8POLY:" );
  d = r8poly_degree ( m, c3 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c4, "  The R8POLY:" );
  d = r8poly_degree ( m, c4 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  r8poly_print ( m, c5, "  The R8POLY:" );
  d = r8poly_degree ( m, c5 );
  printf ( "  Dimensioned degree = %d,  Actual degree = %d\n", m, d );

  return;
}
/******************************************************************************/

void r8poly_deriv_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_DERIV_TEST tests R8POLY_DERIV.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define N 4

  double *c;
  double *cp;
  int d;
  double *x;

  printf ( "\n" );
  printf ( "R8POLY_DERIV_TEST\n" );
  printf ( "  R8POLY_DERIV computes the coefficients of\n" );
  printf ( "  the derivative of a polynomial.\n" );

  x = r8vec_indicator1_new ( N );

  c = roots_to_r8poly ( N, x );

  r8poly_print ( N, c, "  The initial polynomial" );

  for ( d = 0; d <= N; d++ )
  {
    cp = r8poly_deriv ( N, c, d );
    printf ( "\n" );
    printf ( "  The derivative of order %d\n", d );
    printf ( "\n" );
    r8poly_print ( N-d, cp, " " );
    free ( cp );
  }

  free ( c );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void r8poly_lagrange_coef_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_COEF_TEST tests R8POLY_LAGRANGE_COEF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define NPOL 3

  int i;
  int ipol;
  double *pcof;
  double *xpol;

  printf ( "\n" );
  printf ( "R8POLY_LAGRANGE_COEF_TEST\n" );
  printf ( "  R8POLY_LAGRANGE_COEF returns the coefficients for\n" );
  printf ( "  a Lagrange basis polynomial.\n" );

  xpol = r8vec_indicator1_new ( NPOL );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );

  for ( ipol = 1; ipol <= NPOL; ipol++ )
  {
    pcof = r8poly_lagrange_coef ( NPOL, ipol, xpol );

    printf ( "\n" );
    printf ( "  Lagrange basis polynomial %4d:\n", ipol );
    printf ( "\n" );

    for ( i = 0; i < NPOL; i++ )
    {
      printf ( "%10g  %4d\n", pcof[i], i );
    }
    free ( pcof );

  }

  free ( xpol );

  return;
# undef NPOL
}
/******************************************************************************/

void r8poly_lagrange_0_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_0_TEST tests R8POLY_LAGRANGE_0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define NPOL 5

  int ival;
  int nx;
  double wval;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  printf ( "\n" );
  printf ( "R8POLY_LAGRANGE_0_TEST\n" );
  printf ( "  R8POLY_LAGRANGE_0 evaluates the Lagrange\n" );
  printf ( "  factor W(X) at a point.\n" );
  printf ( "\n" );
  printf ( "  The number of data points is %d\n", NPOL );
/*
  Set the abscissas of the polynomials.
*/
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
/*
  Evaluate W(X).
*/
  printf ( "\n" );
  printf ( "      X          W(X)\n" );
  printf ( "\n" );

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    wval = r8poly_lagrange_0 ( NPOL, xpol, xval );

    printf ( "%12g  %12g\n", xval, wval );
  }

  free ( xpol );

  return;
# undef NPOL
}
/******************************************************************************/

void r8poly_lagrange_1_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_1_TEST tests R8POLY_LAGRANGE_1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define NPOL 5

  double dwdx;
  int ival;
  int nx;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  printf ( "\n" );
  printf ( "R8POLY_LAGRANGE_1_TEST\n" );
  printf ( "  R8POLY_LAGRANGE_1 evaluates the Lagrange\n" );
  printf ( "  factor W'(X) at a point.\n" );
  printf ( "\n" );
  printf ( "  The number of data points is %d\n", NPOL );
/*
  Set the abscissas of the polynomials.
*/
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
/*
  Evaluate W'(X).
*/
  printf ( "\n" );
  printf ( "      X          W'(X)\n" );
  printf ( "\n" );

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    dwdx = r8poly_lagrange_1 ( NPOL, xpol, xval );

    printf ( "%12g  %12g\n", xval, dwdx );
  }

  free ( xpol );

  return;
# undef NPOL
}
/******************************************************************************/

void r8poly_lagrange_2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_2_TEST tests R8POLY_LAGRANGE_2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define NPOL 5

  double dw2dx2;
  int ival;
  int nx;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  printf ( "\n" );
  printf ( "R8POLY_LAGRANGE_2_TEST\n" );
  printf ( "  R8POLY_LAGRANGE_2 evaluates the Lagrange\n" );
  printf ( "  factor W''(X) at a point.\n" );
  printf ( "\n" );
  printf ( "  The number of data points is %d\n", NPOL );
/*
  Set the abscissas of the polynomials.
*/
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );

  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
/*
  Evaluate W(X), W'(X), W''.
*/
  printf ( "\n" );
  printf ( "      X          W''(X)\n" );
  printf ( "\n" );

  nx = 4 * NPOL - 1;

  for ( ival = 1; ival <= nx; ival++ )
  {
    xval = r8vec_even_select ( nx, xlo, xhi, ival );

    dw2dx2 = r8poly_lagrange_2 ( NPOL, xpol, xval );

    printf ( "%12g  %12g\n", xval, dw2dx2 );
  }

  free ( xpol );

  return;
# undef NPOL
}
/******************************************************************************/

void r8poly_lagrange_factor_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_FACTOR_TEST tests R8POLY_LAGRANGE_FACTOR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define NPOL 5

  double dwdx;
  int i;
  double wval;
  double xhi;
  double xlo;
  double xpol[NPOL];
  double xval;

  printf ( "\n" );
  printf ( "R8POLY_LAGRANGE_FACTOR_TEST\n" );
  printf ( "  R8POLY_LAGRANGE_FACTOR evaluates the Lagrange\n" );
  printf ( "  factor W(X) at a point.\n" );
  printf ( "\n" );
  printf ( "  For this test, we use %d functions.\n", NPOL );
/*
  Set the abscissas of the polynomials.
*/
  xlo = 0.0;
  xhi = ( double ) ( NPOL - 1 );
  for ( i = 0; i < NPOL; i++ )
  {
    xpol[i] = ( ( double ) ( NPOL - i ) * xlo + ( double ) i * xhi )
      / ( double ) ( NPOL );
  }

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
/*
  Evaluate W(X) and W'(X).
*/
  printf ( "\n" );
  printf ( "      X          W(X)          W'(X)\n" );
  printf ( "\n" );

  for ( i = 0; i < 2 * NPOL - 2; i++ )
  {
    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, i );

    r8poly_lagrange_factor ( NPOL, xpol, xval, &wval, &dwdx );

    printf ( "%10g  %10g  %10g\n", xval, wval, dwdx );
  }

  return;
# undef NPOL
}
/******************************************************************************/

void r8poly_lagrange_val_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_VAL_TEST tests R8POLY_LAGRANGE_VAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
# define NPOL 5

  int i;
  int ipol;
  int ival;
  double p1;
  double p2;
  double xhi;
  double xlo;
  double *xpol;
  double xval;

  printf ( "\n" );
  printf ( "R8POLY_LAGRANGE_VAL_TEST\n" );
  printf ( "  R8POLY_LAGRANGE_VAL evaluates a Lagrange\n" );
  printf ( "  interpolating polynomial at a point.\n" );
  printf ( "\n" );
  printf ( "  For this test, we use %d functions.\n", NPOL );
/*
  Set the abscissas of the polynomials.
*/
  xlo = 0.0E+00;
  xhi = ( double ) ( NPOL - 1 );
  xpol = r8vec_even_new ( NPOL, xlo, xhi );

  r8vec_print ( NPOL, xpol, "  Abscissas:" );
/*
  Evaluate the polynomials.
*/
  printf ( "\n" );
  printf ( "  Here are the values of the functions at\n" );
  printf ( "  several points:\n" );
  printf ( "\n" );
  printf ( "      X          L1          L2          L3      L4          L5\n" );
  printf ( "\n" );

  for ( ival = 0; ival < 2 * NPOL - 1; ival++ )
  {

    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, ival );
    printf ( "%10g  ", xval );

    for ( ipol = 0; ipol < NPOL; ipol++ )
    {
      r8poly_lagrange_val ( NPOL, ipol, xpol, xval, &p1, &p2 );
      printf ( "%10g  ", p1 );
    }
    printf ( "\n" );
  }
  printf ( "\n" );
  printf ( "  And the derivatives:\n" );
  printf ( "\n" );
  printf ( "      X          L'1         L'2         L'3     L'4         L'5\n" );
  printf ( "\n" );

  for ( ival = 0; ival < 2 * NPOL - 1; ival++ )
  {
    xval = r8vec_even_select ( 2 * NPOL - 1, xhi, xlo, ival );
    printf ( "%10g  ", xval );

    for ( ipol = 0; ipol < NPOL; ipol++ )
    {
      r8poly_lagrange_val ( NPOL, ipol, xpol, xval, &p1, &p2 );
      printf ( "%10g  ", p2 );
    }
    printf ( "\n" );
  }

  free ( xpol );

  return;
# undef NPOL
}
/******************************************************************************/

void r8poly_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_PRINT_TEST tests R8POLY_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 January 2015

  Author:

    John Burkardt
*/
{
  double c[6] = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
  int m = 5;

  printf ( "\n" );
  printf ( "R8POLY_PRINT_TEST\n" );
  printf ( "  R8POLY_PRINT prints an R8POLY.\n" );

  r8poly_print ( m, c, "  The R8POLY:" );

  return;
}
/******************************************************************************/

void r8poly_value_horner_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2015

  Author:

    John Burkardt
*/
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double p;
  double *x;
  double x_hi;
  double x_lo;

  printf ( "\n" );
  printf ( "R8POLY_VALUE_HORNER_TEST\n" );
  printf ( "  R8POLY_VALUE_HORNER evaluates a polynomial at\n" );
  printf ( "  one point, using Horner's method.\n" );

  r8poly_print ( m, c, "  The polynomial coefficients:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  printf ( "\n" );
  printf ( "   I    X    P(X)\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    p = r8poly_value_horner ( m, c, x[i] );
    printf ( "  %2d  %8.4f  %14.6g\n", i, x[i], p );
  }

  free ( x );

  return;
}
/******************************************************************************/

void r8poly_values_horner_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUES_HORNER_TEST tests R8POLY_VALUE_HORNER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 December 2013

  Author:

    John Burkardt
*/
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double *p;
  double *x;
  double x_hi;
  double x_lo;

  printf ( "\n" );
  printf ( "R8POLY_VALUES_HORNER_TEST\n" );
  printf ( "  R8POLY_VALUES_HORNER evaluates a polynomial at a\n" );
  printf ( "  point, using Horner's method.\n" );

  r8poly_print ( m, c, "  The polynomial:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  p = r8poly_values_horner ( m, c, n, x );

  r8vec2_print ( n, x, p, "  X, P(X):" );

  free ( p );
  free ( x );

  return;
# undef N
}
/******************************************************************************/

void r8poly2_ex_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_EX_TEST tests R8POLY2_EX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  int ierror;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  printf ( "\n" );
  printf ( "R8POLY2_EX_TEST\n" );
  printf ( "  R8POLY2_EX finds the extreme value\n" );
  printf ( "  of a parabola determined by three points.\n" );

  a =  2.0;
  b = -4.0;
  c = 10.0;

  x1 = 1.0;
  y1 = a * x1 * x1 + b * x1 + c;
  x2 = 2.0;
  y2 = a * x2 * x2 + b * x2 + c;
  x3 = 3.0;
  y3 = a * x3 * x3 + b * x3 + c;

  printf ( "\n" );
  printf ( "  Parabolic coefficients A = %g, B = %g, c = %g\n", a, b, c );
  printf ( "\n" );
  printf ( "  X, Y data:\n" );
  printf ( "\n" );
  printf ( "  %10g.4  %10g.4\n", x1, y1 );
  printf ( "  %10g.4  %10g.4\n", x2, y2 );
  printf ( "  %10g.4  %10g.4\n", x3, y3 );

  a = 0.0;
  b = 0.0;
  c = 0.0;

  ierror = r8poly2_ex ( x1, y1, x2, y2, x3, y3, &xmin, &ymin );

  if ( ierror == 0 )
  {
    printf ( "\n" );
    printf ( "  R8POLY2_EX returns XMIN = %g, YMIN = %g\n", xmin, ymin );
  }
  else
  {
    printf ( "\n" );
    printf ( "  R8POLY2_EX returns error code %d.\n", ierror );
  }

  return;
}
/******************************************************************************/

void r8poly2_ex2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_EX2_TEST tests R8POLY2_EX2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 June 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  int ierror;
  double x1;
  double x2;
  double x3;
  double xmin;
  double y1;
  double y2;
  double y3;
  double ymin;

  printf ( "\n" );
  printf ( "R8POLY2_EX2_TEST\n" );;
  printf ( "  R8POLY2_EX2 finds the extreme value\n" );
  printf ( "  of a parabola determined by three points.\n" );

  a =  2.0;
  b = -4.0;
  c = 10.0;

  x1 = 1.0;
  y1 = a * x1 * x1 + b * x1 + c;
  x2 = 2.0;
  y2 = a * x2 * x2 + b * x2 + c;
  x3 = 3.0;
  y3 = a * x3 * x3 + b * x3 + c;

  printf ( "\n" );
  printf ( "  Parabolic coefficients A = %g, B = %g, c = %g\n", a, b, c );
  printf ( "\n" );
  printf ( "  X, Y data:\n" );
  printf ( "\n" );
  printf ( "  %10g.4  %10g.4\n", x1, y1 );
  printf ( "  %10g.4  %10g.4\n", x2, y2 );
  printf ( "  %10g.4  %10g.4\n", x3, y3 );

  a = 0.0;
  b = 0.0;
  c = 0.0;

  ierror = r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, &xmin, &ymin, &a, &b, &c );

  if ( ierror == 0 )
  {
    printf ( "\n" );
    printf ( "  R8POLY2_EX2 returns XMIN = %d, YMIN = %g\n", xmin, ymin );
    printf ( "  and A = %g, B = %g, c = %g\n", a, b, c );
  }
  else
  {
    printf ( "\n" );
    printf ( "  R8POLY2_EX2 returns error code %d.\n", ierror );
  }

  return;
}
/******************************************************************************/

void r8poly2_val_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_VAL_TEST tests R8POLY2_VAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  int i;
  double x;
  double x1;
  double x2;
  double x3;
  double y;
  double y1;
  double y2;
  double y3;
  double yp;
  double ypp;

  printf ( "\n" );
  printf ( "R8POLY2_VAL_TEST\n" );
  printf ( "  R8POLY2_VAL evaluates a parabola given\n" );
  printf ( "  3 data points.\n" );
  printf ( "\n" );
  printf ( "  Our parabola will be 2*x^2 + 3 * x + 1.\n" );
  printf ( "\n" );
  printf ( "  Case 1: 3 distinct data points:\n" );
  printf ( "\n" );

  x1 = -1.0;
  x2 = 1.0;
  x3 = 3.0;

  r8poly2_f ( x1, &y1, &yp, &ypp );
  r8poly2_f ( x2, &y2, &yp, &ypp );
  r8poly2_f ( x3, &y3, &yp, &ypp );

  printf ( "  %g  %g\n", x1, y1 );
  printf ( "  %g  %g\n", x2, y2 );
  printf ( "  %g  %g\n", x3, y3 );

  printf ( "\n" );
  printf ( "  Sampled data:\n" );
  printf ( "\n" );
  printf ( "  X, Y, Y', Y''\n" );
  printf ( "\n" );
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    printf ( "  %g  %g  %g  %g\n", x, y, yp, ypp );
  }

  printf ( "\n" );
  printf ( "  Case 2: X1=X2, X3 distinct:\n" );
  printf ( "\n" );

  x1 = -1.0;
  x2 = -1.0;
  x3 = 3.0;

  r8poly2_f ( x1, &y1, &y2, &ypp );
  r8poly2_f ( x3, &y3, &yp, &ypp );

  printf ( "  %g  %g\n", x1, y1 );
  printf ( "  %g  %g\n", x2, y2 );
  printf ( "  %g  %g\n", x3, y3 );

  printf ( "\n" );
  printf ( "  Sampled data:\n" );
  printf ( "\n" );
  printf ( "   X, Y, Y', Y''\n" );
  printf ( "\n" );
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    printf ( "  %g  %g  %g  %g\n", x, y, yp, ypp );
  }

  printf ( "\n" );
  printf ( "  Case 3: X1=X2=X3:\n" );
  printf ( "\n" );

  x1 = -1.0;
  x2 = -1.0;
  x3 = -1.0;

  r8poly2_f ( x1, &y1, &y2, &y3 );

  printf ( "  %g  %g\n", x1, y1 );
  printf ( "  %g  %g\n", x2, y2 );
  printf ( "  %g  %g\n", x3, y3 );

  printf ( "\n" );
  printf ( "  Sampled data:\n" );
  printf ( "\n" );
  printf ( "  X, Y, Y', Y''\n" );
  printf ( "\n" );
  for ( i = 0; i < 4; i++ )
  {
    x = ( double ) i;
    r8poly2_val ( x1, y1, x2, y2, x3, y3, x, &y, &yp, &ypp );
    printf ( "  %g  %g  %g  %g\n", x, y, yp, ypp );
  }

  return;
}
/******************************************************************************/

void r8poly2_f ( double x, double *y, double *yp, double *ypp )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_F evaluates a parabola for us.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
  *y = 2.0 * x * x + 3.0 * x + 1.0;
  *yp = 4.0 * x + 3.0;
  *ypp = 4.0;

  return;
}
/******************************************************************************/

void r8poly2_val2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_VAL2_TEST tests R8POLY2_VAL2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define NDATA 5

  int i;
  int left;
  double xdata[NDATA];
  double xval;
  double ydata[NDATA];
  double yval;
  double zdata[NDATA];
  double zval;

  printf ( "\n" );
  printf ( "R8POLY2_VAL2_TEST\n" );
  printf ( "  R8POLY2_VAL2 evaluates parabolas through\n" );
  printf ( "  3 points in a table\n" );
  printf ( "\n" );
  printf ( "  Our data tables will actually be parabolas:\n" );
  printf ( "    A: 2*x^2 + 3 * x + 1.\n" );
  printf ( "    B: 4*x^2 - 2 * x + 5.\n" );
  printf ( "\n" );

  for ( i = 0; i < NDATA; i++ )
  {
    xval = 2.0 * ( double ) ( i + 1 );
    xdata[i] = xval;
    ydata[i] = 2.0 * xval * xval + 3.0 * xval + 1.0;
    zdata[i] = 4.0 * xval * xval - 2.0 * xval + 5.0;
    printf ( "%6d  %10g  %10g  %10g\n", i, xdata[i], ydata[i], zdata[i] );
  }

  printf ( "\n" );
  printf ( "  Interpolated data:\n" );
  printf ( "\n" );
  printf ( "  LEFT, X, Y1, Y2\n" );
  printf ( "\n" );

  for ( i = 0; i <= 4; i++ )
  {
    xval = ( double ) ( 2 * i + 1 );
    left = i;
    if ( NDATA - 3 < left )
    {
      left = NDATA - 3;
    }
    if ( left < 0 )
    {
      left = 0;
    }
    r8poly2_val2 ( NDATA, xdata, ydata, left, xval, &yval );
    r8poly2_val2 ( NDATA, xdata, zdata, left, xval, &zval );

    printf ( "%6d  %10g  %10g  %10g\n", left, xval, yval, zval );
  }

  return;
# undef NDATA
}
/******************************************************************************/

void r8rmat_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_NEW_TEST tests R8RMAT_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 September 2013

  Author:

    John Burkardt
*/
{
  double **a;
  double **b;
  int i;
  int j;
  int k;
  int m;
  int n;

  printf ( "\n" );
  printf ( "R8RMAT_NEW_TEST:\n" );
  printf ( "  R8RMAT_NEW dynamically creates a 2D row major array.\n" );
  printf ( "  Array entries can be addressed using the\n" );
  printf ( "  notation \"a[i][j]\".\n" );
/*
  These dimensions could be entered by the user; they could depend on
  some other calculation; or they could be changed repeatedly during this
  computation, as long as old memory is deleted by R8RMAT_DELETE and new memory
  requested by R8RMAT_NEW.
*/
  m = 4;
  n = 5;
/*
  Allocate memory.
*/
  printf ( "\n" );
  printf ( "  Allocating memory for array A of size %d by %d.\n", m, n );

  a = r8rmat_new ( m, n );

  printf ( "\n" );
  printf ( "  Assigning values to A.\n" );
/*
  Store values in A.
*/
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = ( double ) ( 10 * i + j );
    }
  }
/*
  Print A.
*/
  printf ( "\n" );
  printf ( "  Dynamically allocated matrix A:\n" );
  printf ( "\n" );
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %8g", a[i][j] );
    }
    printf ( "\n" );
  }
/*
  Create a new matrix B to store A' * A.
*/
  b = r8rmat_new ( n, n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        b[i][j] = b[i][j] + a[k][i] * a[k][j];
      }
    }
  }
/*
  Print the matrix.
*/
  printf ( "\n" );
  printf ( "  Dynamically allocated matrix B = A' * A:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      printf ( "  %8g", b[i][j] );
    }
    printf ( "\n" );
  }
/*
  Free memory.
*/
  r8rmat_delete ( m, n, a );
  r8rmat_delete ( n, n, b );

  return;
}
/******************************************************************************/

void r8rmat_to_r8mat_test ( )

/******************************************************************************/
/*
  Purpose:

    R8RMAT_TO_R8MAT_TEST tests R8RMAT_TO_R8MAT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 January 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double **b;
  double *c;
  int m = 5;
  int n = 4;

  printf ( "\n" );
  printf ( "R8RMAT_TO_R8MAT_TEST\n" );
  printf ( "  R8RMAT_TO_R8MAT converts an R8RMAT to an R8MAT.\n" );
  printf ( "\n" );
  printf ( "  Data is of order (%d,%d).\n", m, n );
/*
  Set the R8MAT.
*/
  a = r8mat_indicator_new ( m, n );
  r8mat_print ( m, n, a, "  The R8MAT A:" );
/*
  Convert.
*/
  b = r8mat_to_r8rmat ( m, n, a );
  r8rmat_print ( m, n, b, "  The R8RMAT B:" );
/*
  Recover the matrix.
*/
  c = r8rmat_to_r8mat ( m, n, b );
  r8mat_print ( m, n, c, "  The R8MAT C:" );
/*
  Free memory.
*/
  free ( a );
  r8rmat_delete ( m, n, b );
  free ( c );

  return;
# undef N
}
/******************************************************************************/

void r8row_max_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_MAX_TEST tests R8ROW_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8ROW_MAX_TEST\n" );
  printf ( "  For an R8ROW (a matrix regarded as rows):\n" );
  printf ( "  R8ROW_MAX computes maximums;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  amax = r8row_max ( M, N, a );

  r8vec_print ( M, amax, "  Row maximums:" );

  free ( amax );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8row_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_MEAN_TEST tests R8ROW_MEAN ;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *mean;

  printf ( "\n" );
  printf ( "R8ROW_MEAN_TEST\n" );
  printf ( "  For an R8ROW (a matrix regarded as rows):\n" );
  printf ( "  R8ROW_MEAN computes means;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  mean = r8row_mean ( M, N, a );

  r8vec_print ( M, mean, "  The row means:" );

  free ( mean );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8row_min_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_MIN_TEST tests R8ROW_MIN;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  double *amax;
  double *amin;
  int i;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8ROW_MIN_TEST\n" );
  printf ( "  For an R8ROW (a matrix regarded as rows):\n" );
  printf ( "  R8ROW_MIN computes minimums;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  amin = r8row_min ( M, N, a );

  r8vec_print ( M, amin, "  Row minimums:" );

  free ( amin );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8row_sum_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_SUM_TEST tests R8ROW_SUM;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *rowsum;

  printf ( "\n" );
  printf ( "R8ROW_SUM_TEST\n" );
  printf ( "  For an R8ROW (a matrix regarded as rows):\n" );
  printf ( "  R8ROW_SUM computes sums;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  rowsum = r8row_sum ( M, N, a );

  r8vec_print ( M, rowsum, "  The row sums:" );

  free ( rowsum );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8row_swap_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_SWAP_TEST tests R8ROW_SWAP;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int row1;
  int row2;
  int j;
  int k;

  printf ( "\n" );
  printf ( "R8ROW_SWAP_TEST\n" );
  printf ( "  For an R8ROW (a matrix regarded as rows):\n" );
  printf ( "  R8ROW_SWAP swaps two rows;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  row1 = 1;
  row2 = 3;

  printf ( "\n" );
  printf ( "  Swap rows %d and %d\n", row1, row2 );

  r8row_swap ( M, N, a, row1, row2 );

  r8mat_print ( M, N, a, "  The modified matrix:" );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8row_to_r8vec_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_TO_R8VEC_TEST tests R8ROW_TO_R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *x;

  printf ( "\n" );
  printf ( "R8ROW_TO_R8VEC_TEST\n" );
  printf ( "  R8ROW_TO_R8VEC converts an array of rows into a vector.\n" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) ( 10 * i + j );
    }
  }

  r8mat_print ( M, N, a, "  The array of rows:" );

  x = r8row_to_r8vec ( M, N, a );

  r8vec_print ( M*N, x, "  The resulting vector of rows:" );

  free ( x );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8row_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    R8ROW_VARIANCE_TEST tests R8ROW_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define M 3
# define N 4

  double a[M*N];
  int i;
  int j;
  int k;
  double *variance;

  printf ( "\n" );
  printf ( "R8ROW_VARIANCE_TEST\n" );
  printf ( "  For an R8ROW (a matrix regarded as rows):\n" );
  printf ( "  R8ROW_VARIANCE computes variances;\n" );

  k = 0;
  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      k = k + 1;
      a[i+j*M] = ( double ) k;
    }
  }

  r8mat_print ( M, N, a, "  The original matrix:" );

  variance = r8row_variance ( M, N, a );

  printf ( "\n" );
  printf ( "  Row variances:\n" );
  printf ( "\n" );

  for ( i = 0; i < M; i++ )
  {
    printf ( "  %3d  %10g\n", i+1, variance[i] );
  }

  free ( variance );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8slmat_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8SLMAT_PRINT_TEST tests R8SLMAT_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define TEST_NUM 3

  double *a;
  double a1[21] = {
    21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          32.0, 42.0, 52.0, 62.0, 72.0,
                43.0, 53.0, 63.0, 73.0,
                      54.0, 64.0, 74.0,
                            65.0, 75.0,
                                  76.0 };
  double a2[15] = {
    21.0, 31.0, 41.0, 51.0, 61.0, 71.0,
          32.0, 42.0, 52.0, 62.0, 72.0,
                43.0, 53.0, 63.0, 73.0 };
  double a3[6] = {
    21.0, 31.0, 41.0,
          32.0, 42.0,
                43.0 };
  int m;
  int m_test[TEST_NUM] = { 7, 7, 4 };
  int n;
  int n_test[TEST_NUM] = { 7, 3, 7 };
  int size;
  int size_test[TEST_NUM] = { 21, 15, 6 };
  int test;

  printf ( "\n" );
  printf ( "R8SLMAT_PRINT_TEST\n" );
  printf ( "  R8SLMAT_PRINT prints a strictly lower triangular matrix\n" );
  printf ( "  stored compactly.  Only the (possibly) nonzero \n" );
  printf ( "  elements are printed.\n" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    m = m_test[test];
    n = n_test[test];
    size = size_test[test];
    a = ( double * ) malloc ( size * sizeof ( double ) );

    if ( test == 0 )
    {
      r8vec_copy ( size, a1, a );
    }
    else if ( test == 1 )
    {
      r8vec_copy ( size, a2, a );
    }
    else if ( test == 2 )
    {
      r8vec_copy ( size, a3, a );
    }

    r8slmat_print ( m, n, a, "  R8SLMAT:" );

    free ( a );
  }

  return;
# undef TEST_NUM
}
/******************************************************************************/

void r8vec_amax_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMAX_TEST tests R8VEC_AMAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_AMAX_TEST\n" );
  printf ( "  For an R8VEC:\n" );
  printf ( "  R8VEC_AMAX:      maximum magnitude entry;\n" );

  b = - ( double ) N;
  c =  ( double ) N;

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );

  aval = r8vec_amax ( N, a );
  printf ( "  Maximum absolute:         %g\n", aval );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_amin_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMIN_TEST tests R8VEC_AMIN;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_AMIN_TEST\n" );
  printf ( "  For an R8VEC:\n" );
  printf ( "  R8VEC_AMIN:      minimum magnitude entry.\n" );

  b = - ( double ) N;
  c =  ( double ) N;

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );

  aval = r8vec_amin ( N, a );
  printf ( "  Minimum absolute:         %g\n", aval );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_bracket_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET_TEST tests R8VEC_BRACKET.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define N 10
# define TEST_NUM 6

  int i;
  int left;
  int right;
  int test;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  printf ( "\n" );
  printf ( "R8VEC_BRACKET_TEST\n" );
  printf ( "  R8VEC_BRACKET finds a pair of entries in a\n" );
  printf ( "  sorted real array which bracket a value.\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!)" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    xval = xtest[test];

    printf ( "\n" );
    printf ( "  Search for XVAL = %g\n", xval );

    r8vec_bracket ( N, x, xval, &left, &right );

    printf ( "  X[%d-1] = %g\n", left, x[left-1] );
    printf ( "  X[%d-1] = %g\n", right, x[right-1] );
  }

  return;

# undef N
# undef TEST_NUM
}
/******************************************************************************/

void r8vec_bracket2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET2_TEST tests R8VEC_BRACKET2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define N 10
# define TEST_NUM 6

  int i;
  int left;
  int right;
  int start;
  int test;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  printf ( "\n" );
  printf ( "R8VEC_BRACKET2_TEST\n" );
  printf ( "  R8VEC_BRACKET2 finds a pair of entries in a\n" );
  printf ( "  sorted R8VEC which bracket a value.\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!)" );

  for ( test = 0; test < TEST_NUM; test++ )
  {

    xval = xtest[test];

    printf ( "\n" );
    printf ( "  Search for XVAL = %g\n", xval );

    if ( 0 < left )
    {
      start = left;
    }
    else
    {
      start = ( N + 1 ) / 2;
    }

    printf ( "  Start = %d\n", start );

    r8vec_bracket2 ( N, x, xval, start, &left, &right );

    printf ( "  Left =  %d\n", left );
    printf ( "  Right = %d\n", right );

    if ( 1 <= left )
    {
      printf ( "  X[%d-1] = %g\n", left, x[left-1] );
    }

    if ( 1 <= right )
    {
      printf ( "  X[%d-1] = %g\n", right, x[right-1] );
    }
  }
  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void r8vec_bracket3_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET3_TEST tests R8VEC_BRACKET3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 June 2012

  Author:

    John Burkardt
*/
{
# define N 10
# define TEST_NUM 6

  int i;
  int itest;
  int left;
  double x[N];
  double xtest[TEST_NUM] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  printf ( "\n" );
  printf ( "R8VEC_BRACKET3_TEST\n" );
  printf ( "  R8VEC_BRACKET3 finds a pair of entries in a\n" );
  printf ( "  sorted real array which bracket a value.\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  The array (must be in ascending order!):" );

  left = ( N - 1 ) / 2;

  for ( itest = 0; itest < TEST_NUM; itest++ )
  {
    xval = xtest[itest];

    printf ( "\n" );
    printf ( "  Search for XVAL = %g\n", xval );

    printf ( "  Starting guess for interval is = %d\n", left );

    r8vec_bracket3 ( N, x, xval, &left );

    printf ( "  Nearest interval:\n" );
    printf ( "   X[%d]= %g\n", left, x[left] );
    printf ( "   X[%d]= %g\n", left+1, x[left+1]  );
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void r8vec_bracket5_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET5_TEST tests R8VEC_BRACKET5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 October 2012

  Author:

    John Burkardt
*/
{
  int left;
  int n = 10;
  int right;
  int test;
  int test_num = 6;
  double *x;
  double xtest[6] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  printf ( "\n" );
  printf ( "R8VEC_BRACKET5_TEST\n" );
  printf ( "  R8VEC_BRACKET5 finds a pair of entries in a\n" );
  printf ( "  sorted R8VEC which bracket a value.\n" );

  x = r8vec_indicator1_new ( n );
  x[5] = x[4];

  r8vec_print ( n, x, "  Sorted array:" );

  printf ( "\n" );
  printf ( "        LEFT                   RIGHT\n" );
  printf ( "      X(LEFT)       XVAL     X(RIGHT)\n" );
  printf ( "\n" );

  for ( test = 0; test < test_num; test++ )
  {
    xval = xtest[test];

    left = r8vec_bracket5 ( n, x, xval );

    if ( left == -1 )
    {
      printf ( "  %10d\n", left );
      printf ( "              %10.4f  (Not bracketed!)\n", xval );
    }
    else
    {
      right = left + 1;
      printf ( "  %10d              %10d\n", left, right );
      printf ( "  %10.4f  %10.4f  %10.4f\n", x[left], xval, x[right] );
    }
  }

  free ( x );

  return;
}
/******************************************************************************/

void r8vec_chebyspace_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CHEBYSPACE_NEW_TEST tests R8VEC_CHEBYSPACE_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2011

  Author:

    John Burkardt
*/
{
  int n;
  double *r;
  double r1;
  double r2;

  printf ( "\n" );
  printf ( "R8VEC_CHEBYSPACE_NEW_TEST\n" );
  printf ( "  R8VEC_CHEBYSPACE_NEW computes N Chebyshev points in [R1,R2].\n" );

  r1 = -1.0;
  r2 = +1.0;
  n = 5;

  r = r8vec_chebyspace_new ( n, r1, r2 );

  printf ( "\n" );
  printf ( "  N = %d, R1 = %f, R2 = %f\n", n, r1, r2 );

  r8vec_print ( n, r, "  Chebyshev points:" );

  free ( r );

  r1 =   0.0;
  r2 = +10.0;
  n = 7;

  r = r8vec_chebyspace_new ( n, r1, r2 );

  printf ( "\n" );
  printf ( "  N = %d, R1 = %f, R2 = %f\n", n, r1, r2 );

  r8vec_print ( n, r, "  Chebyshev points:" );

  free ( r );

  return;
}
/******************************************************************************/

void r8vec_concatenate_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CONCATENATE_NEW_TEST tests R8VEC_CONCATENATE_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 October 2014

  Author:

    John Burkardt
*/
{
  int n1 = 5;
  int n2 = 3;
  int n3 = n1 + n2;

  double a1[5] = { 91.1, 31.2, 71.3, 51.4, 31.5 };
  double a2[3] = { 42.6, 22.7, 12.8 };
  double *a3;

  printf ( "\n" );
  printf ( "R8VEC_CONCATENATE_NEW_TEST\n" );
  printf ( "  R8VEC_CONCATENATE_NEW concatenates two R8VECs\n" );

  r8vec_print ( n1, a1, "  Array 1:" );
  r8vec_print ( n2, a2, "  Array 2:" );
  a3 = r8vec_concatenate_new ( n1, a1, n2, a2 );
  r8vec_print ( n3, a3, "  Array 3 = Array 1 + Array 2:" );

  free ( a3 );

  return;
}
/******************************************************************************/

void r8vec_convolution_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CONVOLUTION_TEST tests R8VEC_CONVOLUTION

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt
*/
{
# define M 4
# define N 3

  int m = M;
  int n = N;

  double x[M] = { 1.0, 2.0, 3.0, 4.0 };
  double y[N] = { -1.0, 5.0, 3.0 };
  double *z;
  double z_correct[M+N-1] = { -1.0, 3.0, 10.0, 17.0, 29.0, 12.0 };

  printf ( "\n" );
  printf ( "R8VEC_CONVOLUTION_TEST\n" );
  printf ( "  R8VEC_CONVOLUTION computes the convolution\n" );
  printf ( "  of two vectors.\n" );

  r8vec_print ( m, x, "  The factor X:" );
  r8vec_print ( n, y, "  The factor Y:" );

  z = r8vec_convolution ( m, x, n, y );

  r8vec_print ( m + n - 1, z, "  The convolution z = x star y:" );

  r8vec_print ( m + n - 1, z_correct, "  Correct answer:" );

  free ( z );

  return;
# undef M
# undef N
}
/******************************************************************************/

void r8vec_convolution_circ_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CONVOLUTION_CIRC_TEST tests R8VEC_CONVOLUTION_CIRC

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 August 2011

  Author:

    John Burkardt
*/
{
# define N 4

  double x[N] = { 1.0, 2.0, 3.0, 4.0 };
  double y[N] = { 1.0, 2.0, 4.0, 8.0 };
  double *z;
  double z_correct[N] = { 37.0, 44.0, 43.0, 26.0 };

  printf ( "\n" );
  printf ( "R8VEC_CONVOLUTION_CIRC_TEST\n" );
  printf ( "  R8VEC_CONVOLUTION_CIRC computes the circular convolution\n" );
  printf ( "  of two vectors.\n" );

  r8vec_print ( N, x, "  The factor X:" );
  r8vec_print ( N, y, "  The factor Y:" );

  z = r8vec_convolution_circ ( N, x, y );

  r8vec_print ( N, z, "  The circular convolution z = xCCy:" );

  r8vec_print ( N, z_correct, "  Correct answer:" );

  free ( z );

  return;
# undef N
}
/******************************************************************************/

void r8vec_dif_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIF_TEST tests R8VEC_DIF.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  double *cof;
  double fdif;
  double h = 0.01;
  int i;
  int n = 4;
  double x = 1.0;
  double xi;

  printf ( "\n" );
  printf ( "R8VEC_DIF_TEST\n" );
  printf ( "  R8VEC_DIF estimates derivatives.\n" );
  printf ( "\n" );
  printf ( "  Estimate the derivative of order N = %d\n", n );
  printf ( "  Using H = %g\n", h );
  printf ( "  at argument X = %g\n", x );
/*
  Get the coefficients.
*/
  cof = r8vec_dif ( n, h );

  r8vec_print ( n+1, cof, "  The difference coefficients:" );

  fdif = 0.0;
  for ( i = 0; i <= n; i++ )
  {
    xi = x + ( double ) ( 2 * i - n ) * h;
    fdif = fdif + cof[i] * r8vec_dif_f ( xi );
  }

  printf ( "\n" );
  printf ( "  Estimate is FDIF = %g\n", fdif );

  free ( cof );

  return;
}
/******************************************************************************/

double r8vec_dif_f ( double x )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIF_F evaluates the function used in R8VEC_DIF_TEST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  double value;

  value = exp ( x );

  return value;
}
/******************************************************************************/

void r8vec_direct_product_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIRECT_PRODUCT_TEST tests R8VEC_DIRECT_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  int factor_num = 3;
  int point_num = 24;

  int factor_index;
  int factor_order;
  double *factor_value;
  int i;
  int j;
  double x[factor_num*point_num];

  printf ( "\n" );
  printf ( "R8VEC_DIRECT_PRODUCT_TEST\n" );
  printf ( "  R8VEC_DIRECT_PRODUCT forms the entries of a\n" );
  printf ( "  direct product of a given number of R8VEC factors.\n" );

  for ( j = 0; j < point_num; j++ )
  {
    for ( i = 0; i < factor_num; i++ )
    {
      x[i+j*factor_num] = 0.0;
    }
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = ( double * ) malloc ( factor_order * sizeof ( double ) );
      factor_value[0] = 1.0;
      factor_value[1] = 2.0;
      factor_value[2] = 3.0;
      factor_value[3] = 4.0;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = ( double * ) malloc ( factor_order * sizeof ( double ) );
      factor_value[0] = 50.0;
      factor_value[1] = 60.0;
      factor_value[2] = 70.0;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = ( double * ) malloc ( factor_order * sizeof ( double ) );
      factor_value[0] = 800.0;
      factor_value[1] = 900.0;
    }

    r8vec_direct_product ( factor_index, factor_order, factor_value,
      factor_num, point_num, x );

    free ( factor_value );
  }

  printf ( "\n" );
  printf ( "     J         X(1)      X(2)      X(3)\n" );
  printf ( "\n" );

  for ( j = 0; j < point_num; j++ )
  {
    printf ( "  %4d  ", j );
    for ( i = 0; i < factor_num; i++ )
    {
      printf ( "  %8g", x[i+j*factor_num] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void r8vec_direct_product2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIRECT_PRODUCT2_TEST tests R8VEC_DIRECT_PRODUCT2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
  int factor_num = 3;
  int point_num = 24;

  int factor_index;
  int factor_order;
  double *factor_value;
  int i;
  int j;
  double w[point_num];

  printf ( "\n" );
  printf ( "R8VEC_DIRECT_PRODUCT2_TEST\n" );
  printf ( "  R8VEC_DIRECT_PRODUCT2 forms the entries of a\n" );
  printf ( "  direct product of a given number of R8VEC factors.\n" );

  for ( j = 0; j  < point_num; j++ )
  {
    w[j] = 1.0;
  }

  for ( factor_index = 0; factor_index < factor_num; factor_index++ )
  {
    if ( factor_index == 0 )
    {
      factor_order = 4;
      factor_value = ( double * ) malloc ( factor_order * sizeof ( double ) );
      factor_value[0] = 2.0;
      factor_value[1] = 3.0;
      factor_value[2] = 5.0;
      factor_value[3] = 7.0;
    }
    else if ( factor_index == 1 )
    {
      factor_order = 3;
      factor_value = ( double * ) malloc ( factor_order * sizeof ( double ) );
      factor_value[0] = 11.0;
      factor_value[1] = 13.0;
      factor_value[2] = 17.0;
    }
    else if ( factor_index == 2 )
    {
      factor_order = 2;
      factor_value = ( double * ) malloc ( factor_order * sizeof ( double ) );
      factor_value[0] = 19.0;
      factor_value[1] = 21.0;
    }

    r8vec_direct_product2 ( factor_index, factor_order, factor_value,
      factor_num, point_num, w );

    free ( factor_value );
  }

  printf ( "\n" );
  printf ( "     J         W(J)\n" );
  printf ( "\n" );

  for ( j = 0; j < point_num; j++ )
  {
    printf ( "  %4d    %8g\n", j, w[j] );
  }

  return;
}
/******************************************************************************/

void r8vec_even_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN_TEST tests R8VEC_EVEN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *x;
  double xhi = 99.0;
  double xlo = 0.0;

  printf ( "\n" );
  printf ( "R8VEC_EVEN_TEST\n" );
  printf ( "  R8VEC_EVEN computes N evenly spaced values\n" );
  printf ( "  between XLO and XHI.\n" );
  printf ( "\n" );
  printf ( "  XLO = %g\n", xlo );
  printf ( "  XHI = %g\n", xhi );
  printf ( "  while N = %d\n", N );

  x = r8vec_even_new ( N, xlo, xhi );

  r8vec_print ( N, x, "  Resulting array:" );

  free ( x );

  return;
# undef N
}
/******************************************************************************/

void r8vec_even2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN2_TEST tests R8VEC_EVEN2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define NOLD 5
# define MAXVAL 20

  int i;
  int istar;
  int jstar;
  int nfill[NOLD-1] = { 4, 3, 5, 0 };
  int nval;
  double xold[NOLD] = { 0.0, 1.0, 5.0, 2.0, 0.0 };
  double xval[MAXVAL];

  printf ( "\n" );
  printf ( "R8VEC_EVEN2_TEST:\n" );
  printf ( "  R8VEC_EVEN2 interpolates a specified number of\n" );
  printf ( "  points pairs of values in a vector.\n" );
  printf ( "\n" );
  printf ( "  Input data:\n" );
  printf ( "\n" );
  for ( i = 1; i <= NOLD; i++ )
  {
    printf ( "  %12g\n", xold[i-1] );
    if ( i < NOLD )
    {
      printf ( "  (%d)\n", nfill[i-1] );
    }
  }

  r8vec_even2 ( MAXVAL, nfill, NOLD, xold, &nval, xval );

  printf ( "\n" );
  printf ( "  Resulting vector:\n" );
  printf ( "\n" );

  istar = 1;
  jstar = 1;
  for ( i = 1; i <= nval; i++ )
  {
    if ( i == istar )
    {
      printf ( "  *  %12g\n", xval[i-1] );

      if ( jstar < NOLD )
      {
        istar = istar + nfill[jstar-1] + 1;
        jstar = jstar + 1;
      }
    }
    else
    {
      printf ( "     %12g\n", xval[i-1] );
    }
  }

  return;
# undef MAXVAL
# undef NOLD
}
/******************************************************************************/

void r8vec_expand_linear_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EXPAND_LINEAR_TEST tests R8VEC_EXPAND_LINEAR.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N 6

  int fat = 3;
  int nfat;
  double x[N] = { 16.0, 4.0, 0.0, 4.0, 16.0, 36.0 };
  double *xfat;

  printf ( "\n" );
  printf ( "R8VEC_EXPAND_LINEAR_TEST\n" );
  printf ( "  R8VEC_EXPAND_LINEAR linearly interpolates new data\n" );
  printf ( "  between old values.\n" );
  printf ( "\n" );

  r8vec_print ( N, x, "  Original vector:" );

  printf ( "\n" );
  printf ( "  Expansion factor is %d\n", fat );

  xfat = r8vec_expand_linear ( N, x, fat );

  nfat = ( N - 1 ) * ( fat + 1 ) + 1;

  r8vec_print ( nfat, xfat, "  Fattened vector:" );

  free ( xfat );

  return;
# undef N
}
/******************************************************************************/

void r8vec_frac_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_FRAC_TEST tests R8VEC_FRAC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double afrac;
  int k;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_FRAC_TEST\n" );
  printf ( "  R8VEC_FRAC: K-th smallest R8VEC entry;\n" );

  seed = 123456789;

  a = r8vec_uniform_01_new ( N, &seed );

  r8vec_print ( N, a, "  Array to search:" );

  printf ( "\n" );
  printf ( "  Fractile  Value\n" );
  printf ( "\n" );

  for ( k = 1; k < N; k = k + N/2 )
  {
    afrac = r8vec_frac ( N, a, k );
    printf ( "  %6d  %14g\n", k, afrac );
  }

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_histogram_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_HISTOGRAM_TEST tests R8VEC_HISTOGRAM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define HISTO_NUM 20
# define N 1000

  double *a;
  double a_hi;
  double a_lo;
  double bin_hi;
  double bin_lo;
  int *histo_gram;
  int i;
  int seed = 123456789;
  int test;
  int test_num = 2;

  printf ( "\n" );
  printf ( "R8VEC_HISTOGRAM_TEST\n" );
  printf ( "  R8VEC_HISTOGRAM histograms a real vector.\n" );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      printf ( "\n" );
      printf ( "  Uniform data:\n" );

      a_lo =  0.0;
      a_hi = +1.0;
      a = r8vec_uniform_01_new ( N, &seed );
    }
    else if ( test == 2 )
    {
      printf ( "\n" );
      printf ( "  Normal data:\n" );
      a_lo = -3.0;
      a_hi = +3.0;
      a = r8vec_normal_01_new ( N, &seed );
    }

    histo_gram = r8vec_histogram ( N, a, a_lo, a_hi, HISTO_NUM );

    printf ( "\n" );
    printf ( "  Histogram of data:\n" );
    printf ( "\n" );

    for ( i = 0; i <= HISTO_NUM+1; i++ )
    {
      if ( i == 0 )
      {
        printf ( "              %10g  %6d\n", a_lo, histo_gram[i] );
      }
      else if ( i <= HISTO_NUM )
      {
        bin_lo = ( ( double ) ( HISTO_NUM - i + 1 ) * a_lo
                 + ( double ) (             i - 1 ) * a_hi )
                 / ( double ) ( HISTO_NUM         );

        bin_hi = ( ( double ) ( HISTO_NUM - i     ) * a_lo
                 + ( double ) (             i     ) * a_hi )
                 / ( double ) ( HISTO_NUM         );

        printf ( "  %10g  %10g  %6d\n", bin_lo, bin_hi, histo_gram[i] );
      }
      else if ( i == HISTO_NUM+1 )
      {
        printf ( "  %10g               %6d\n", a_hi, histo_gram[i] );
      }
    }
    free ( a );
    free ( histo_gram );
  }

  return;
# undef HISTO_NUM
# undef N
}
/******************************************************************************/

void r8vec_house_column_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_HOUSE_COLUMN_TEST tests R8VEC_HOUSE_COLUMN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b = 0.0;
  double c = 5.0;
  double *h;
  double *ha;
  int i;
  int j;
  int k;
  int n = 4;
  int seed;
  double *v;

  printf ( "\n" );
  printf ( "R8VEC_HOUSE_COLUMN_TEST\n" );
  printf ( "  R8VEC_HOUSE_COLUMN returns the compact form of\n" );
  printf ( "  a Householder matrix that packs a column\n" );
  printf ( "  of a matrix.\n" );
/*
  Get a random matrix.
*/
  seed = 123456789;

  a = r8mat_uniform_ab_new ( n, n, b, c, &seed );

  r8mat_print ( n, n, a, "  Matrix A:" );

  for ( k = 1; k <= n-1; k++ )
  {
    printf ( "\n" );
    printf ( "  Working on column K = %d\n", k );

    v = r8vec_house_column ( n, a+(k-1)*n, k );

    h = r8mat_house_form ( n, v );

    r8mat_print ( n, n, h, "  Householder matrix H:" );

    ha = r8mat_mm_new ( n, n, n, h, a );

    r8mat_print ( n, n, ha, "  Product H*A:" );
/*
  If we set A := HA, then we can successively convert A to upper
  triangular form.
*/
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        a[i+j*n] = ha[i+j*n];
      }
    }
    free ( h );
    free ( ha );
    free ( v );
  }
  free ( a );

  return;
}
/******************************************************************************/

void r8vec_index_delete_all_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_DELETE_ALL_TEST tests R8VEC_INDEX_DELETE_ALL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_DELETE_ALL_TEST\n" );
  printf ( "  R8VEC_INDEX_DELETE_ALL deletes all copies of a\n" );
  printf ( "  particular value.\n" );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "  %g\n", xval );
    r8vec_index_insert ( &n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  printf ( "\n" );
  printf ( "  Call R8VEC_INDEX_DELETE_ALL to delete all values of 7:\n" );

  xval = 7.0;
  r8vec_index_delete_all ( n, x, indx, xval, &n, x, indx );

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_delete_dupes_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_DELETE_DUPES_TEST tests R8VEC_INDEX_DELETE_DUPES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_DELETE_DUPES_TEST\n" );
  printf ( "  R8VEC_INDEX_DELETE_DUPES deletes duplicates.\n" );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "  %g\n", xval );
    r8vec_index_insert ( &n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  printf ( "\n" );
  printf ( "  Call R8VEC_INDEX_DELETE_DUPES to delete duplicates:\n" );

  r8vec_index_delete_dupes ( n, x, indx, &n, x, indx );

  printf ( "\n" );
  printf ( "  Indexed list of unique entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g\n", i+1, indx[i], x[i] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_delete_one_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_DELETE_ONE_TEST tests R8VEC__INDEX_DELETE_ONE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_DELETE_ONE_TEST\n" );
  printf ( "  R8VEC_INDEX_DELETE_ONE deletes one copies of a\n" );
  printf ( "  particular value.\n" );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "  %g\n", xval );
    r8vec_index_insert ( &n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  printf ( "\n" );
  printf ( "  Call R8VEC_INDEX_DELETE_ONE to delete one value of 8:\n" );

  xval = 8.0;
  r8vec_index_delete_one ( n, x, indx, xval, &n, x, indx );

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_insert_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_INSERT_TEST tests R8VEC_INDEX_INSERT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 25

  int i;
  int indx[N_MAX];
  int n;
  int n2;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_INSERT_TEST\n" );
  printf ( "  R8VEC_INDEX_INSERT inserts values into an\n" );
  printf ( "  index sorted array.\n" );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  seed = 123456789;

  for ( i = 1; i <= 20; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "  %g\n", xval );
    r8vec_index_insert ( &n, x, indx, xval );
  }

  xval = 7.0;
  r8vec_index_insert ( &n, x, indx, xval );

  xval = 8.0;
  r8vec_index_insert ( &n, x, indx, xval );

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_insert_unique_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_INSERT_UNIQUE_TEST tests R8VEC_INDEX_INSERT_UNIQUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double b;
  double c;
  int equal;
  int i;
  int indx[N_MAX];
  int less;
  int more;
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_INSERT_UNIQUE_TEST\n" );
  printf ( "  R8VEC_INDEX_INSERT_UNIQUE inserts unique values into an\n" );
  printf ( "  index sorted array.\n" );

  b = 0.0;
  c = 20.0;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( b, c, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "    %6g\n", xval );
    r8vec_index_insert_unique ( &n, x, indx, xval );
  }

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 1; i <= n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_order_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_ORDER_TEST tests R8VEC_INDEX_ORDER.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  int i;
  int indx[N_MAX];
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_ORDER_TEST\n" );
  printf ( "  R8VEC_INDEX_ORDER sorts an index sorted array.\n" );
  printf ( "\n" );
  printf ( "  Generate some random values:\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( 0.0, 20.0, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "  %6g\n", xval );
    r8vec_index_insert_unique ( &n, x, indx, xval );
  }

  printf ( "\n" );
  printf ( "  Indexed list of unique entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  printf ( "\n" );
  printf ( "  Now call R8VEC_INDEX_ORDER to carry out the sorting:\n" );

  r8vec_index_order ( n, x, indx );

  r8vec_print ( n, x, "  X:" );

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_search_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_SEARCH_TEST tests R8VEC_INDEX_SEARCH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 20

  double b;
  double c;
  int equal;
  int i;
  int indx[N_MAX];
  int less;
  int more;
  int n;
  int seed;
  double x[N_MAX];
  double xval;

  n = 0;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_SEARCH_TEST\n" );
  printf ( "  R8VEC_INDEX_SEARCH searches for an entry \n" );
  printf ( "  with a given value.\n" );

  b = 0.0;
  c = 20.0;
  seed = 123456789;

  for ( i = 1; i <= N_MAX; i++ )
  {
    xval = r8_uniform_ab ( b, c, &seed );
    xval = ( double ) ( r8_nint ( xval ) );
    printf ( "    %6g\n", xval );
    r8vec_index_insert_unique ( &n, x, indx, xval );
  }

  printf ( "\n" );
  printf ( "  Indexed list of entries:\n" );
  printf ( "\n" );
  printf ( "  I  INDX(I)  X(I)  X(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 1; i <= n; i++ )
  {
    printf ( "  %3d  %3d  %6g  %6g\n", i+1, indx[i], x[i], x[indx[i]-1] );
  }

  printf ( "\n" );
  printf ( "  Results of search for given XVAL:\n" );
  printf ( "\n" );
  printf ( "  XVAL  Less Equal More\n" );
  printf ( "\n" );

  for ( i = 0; i <= N_MAX; i++ )
  {
    xval = ( double ) ( i );
    r8vec_index_search ( n, x, indx, xval, &less, &equal, &more );
    printf ( "  %6g  %3d  %3d  %3d\n", xval, less, equal, more );
  }

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_index_sorted_range_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_SORTED_RANGE_TEST tests R8VEC_INDEX_SORTED_RANGE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2010

  Author:

    John Burkardt
*/
{
  int i;
  int i_hi;
  int i_lo;
  int *indx;
  int n = 20;
  double r[20];
  double r_lo;
  double r_hi;
  int seed;
  double t;
  int test;

  printf ( "\n" );
  printf ( "R8VEC_INDEX_SORTED_RANGE_TEST\n" );
  printf ( "  R8VEC_INDEX_SORTED_RANGE seeks the range I_LO:I_HI\n" );
  printf ( "  of entries of sorted indexed R so that\n" );
  printf ( "  R_LO <= R(INDX(I)) <= R_HI for I_LO <= I <= I_HI.\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r8vec_uniform_01 ( n, &seed, r );

    r8vec_print ( n, r, "  Array" );

    indx = r8vec_sort_heap_index_a_new ( n, r );

    printf ( "\n" );
    printf ( "     I  INDX    R(INDX(I))\n" );
    printf ( "\n" );
    for ( i = 0; i < n; i++ )
    {
      printf ( "  %4d  %4d  %14f\n", i, indx[i], r[indx[i]] );
    }

    r_lo = r8_uniform_01 ( &seed );
    r_hi = r8_uniform_01 ( &seed );

    if ( r_hi < r_lo )
    {
      t = r_lo;
      r_lo = r_hi;
      r_hi = t;
    }

    r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, &i_lo, &i_hi );

    printf ( "\n" );
    if ( i_hi < i_lo )
    {
      printf ( "  R_LO        %14f\n", r_lo );
      printf ( "  R_HI        %14f\n", r_hi );
      printf ( "  Empty range in R.\n" );
    }
    else
    {
      printf ( "  R_LO        %14f\n", r_lo );
      for ( i = i_lo; i <= i_hi; i++ )
      {
        printf ( "  %4d  %4d  %14f\n", i, indx[i], r[indx[i]] );
      }
      printf ( "  R_HI        %14f\n", r_hi );
    }
    free ( indx );
  }

  return;
}
/******************************************************************************/

void r8vec_indexed_heap_d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_TEST tests R8VEC_INDEXED_HEAP_D;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  double a[20] = {
    101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0,
    111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0 };
  int i;
  int indx[10] = {
    0, 10, 16, 4, 6, 12, 14, 2, 18, 8 };
  int m = 20;
  int n = 10;

  printf ( "\n" );
  printf ( "R8VEC_INDEXED_HEAP_D_TEST\n" );
  printf ( "  R8VEC_INDEXED_HEAP_D creates a descending heap\n" );
  printf ( "  from an indexed vector.\n" );
/*
  Print before.
*/
  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  printf ( "\n" );
  printf ( "  A(INDX):\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Heap the data.
*/
  r8vec_indexed_heap_d ( n, a, indx );
/*
  Print afterwards.  Only INDX should change.
*/
  r8vec_print ( m, a, "  The data vector (should NOT change):" );
  i4vec_print ( n, indx, "  The index vector (may change):" );
  printf ( "\n" );
  printf ( "  A(INDX) is now a descending heap:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }

  return;
}
/******************************************************************************/

void r8vec_indexed_heap_d_extract_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_EXTRACT_TEST tests R8VEC_INDEXED_HEAP_D_EXTRACT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  printf ( "\n" );
  printf ( "R8VEC_INDEXED_HEAP_D_EXTRACT_TEST\n" );
  printf ( "  For an indexed R8VEC,\n" );
  printf ( "  R8VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;\n" );
/*
  Set the data array.  To keep things easy, we will use the indicator vector.
*/
  a = r8vec_indicator1_new ( m );
/*
  The index array will initially be a random subset of the numbers 1 to M,
  in random order.
*/
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  printf ( "\n" );
  printf ( "  A(INDX):\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Create a descending heap from the indexed array.
*/
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  printf ( "\n" );
  printf ( "  A(INDX) after heaping:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Insert five entries, and monitor the maximum.
*/
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    printf ( "\n" );
    printf ( "  Inserting value %f\n", a[indx_insert] );

    r8vec_indexed_heap_d_insert ( &n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    printf ( "  Current maximum is %f\n", a[indx_max] );
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  printf ( "\n" );
  printf ( "  A(INDX) after insertions:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Extract the first 5 largest elements.
*/
  printf ( "\n" );
  printf ( "  Now extract the maximum several times.\n" );
  printf ( "\n" );

  for ( i = 0; i < 5; i++ )
  {
    indx_extract = r8vec_indexed_heap_d_extract ( &n, a, indx );
    printf ( "  Extracting maximum element A[%d] = %f\n",
      indx_extract, a[indx_extract] );
  }
  r8vec_print ( m, a, "  The data vector after extractions:" );
  i4vec_print ( n, indx, "  The index vector after extractions:" );
  printf ( "\n" );
  printf ( "  A(INDX) after extractions:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_indexed_heap_d_insert_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_INSERT_TEST tests R8VEC_INDEXED_HEAP_D_INSERT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  printf ( "\n" );
  printf ( "R8VEC_INDEXED_HEAP_D_INSERT_TEST\n" );
  printf ( "  For an indexed R8VEC,\n" );
  printf ( "  R8VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.\n" );
/*
  Set the data array.  To keep things easy, we will use the indicator vector.
*/
  a = r8vec_indicator1_new ( m );
/*
  The index array will initially be a random subset of the numbers 1 to M,
  in random order.
*/
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  printf ( "\n" );
  printf ( "  A(INDX):\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Create a descending heap from the indexed array.
*/
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  printf ( "\n" );
  printf ( "  A(INDX) after heaping:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Insert five entries, and monitor the maximum.
*/
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    printf ( "\n" );
    printf ( "  Inserting value %f\n", a[indx_insert] );

    r8vec_indexed_heap_d_insert ( &n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    printf ( "  Current maximum is %f\n", a[indx_max] );
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  printf ( "\n" );
  printf ( "  A(INDX) after insertions:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_indexed_heap_d_max_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_MAX_TEST tests R8VEC_INDEXED_HEAP_D_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int indx[20];
  int indx_extract;
  int indx_insert;
  int indx_max;
  int m = 20;
  int n;
  int n_max = 20;

  printf ( "\n" );
  printf ( "R8VEC_INDEXED_HEAP_D_MAX_TEST\n" );
  printf ( "  For an indexed R8VEC,\n" );;
  printf ( "  R8VEC_INDEXED_HEAP_D_MAX reports the maximum value.\n" );
/*
  Set the data array.  To keep things easy, we will use the indicator vector.
*/
  a = r8vec_indicator1_new ( m );
/*
  The index array will initially be a random subset of the numbers 1 to M,
  in random order.
*/
  n = 5;
  indx[0]  =  8;
  indx[1]  =  1;
  indx[2]  =  7;
  indx[3]  = 13;
  indx[4]  =  4;
  indx[5]  =  6;
  indx[6]  = 14;
  indx[7]  =  0;
  indx[8]  = 18;
  indx[9]  = 19;
  indx[10] =  2;

  r8vec_print ( m, a, "  The data vector:" );
  i4vec_print ( n, indx, "  The index vector:" );
  printf ( "\n" );
  printf ( "  A(INDX):\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Create a descending heap from the indexed array.
*/
  r8vec_indexed_heap_d ( n, a, indx );

  i4vec_print ( n, indx, "  The index vector after heaping:" );
  printf ( "\n" );
  printf ( "  A(INDX) after heaping:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }
/*
  Insert five entries, and monitor the maximum.
*/
  for ( i = 0; i < 5; i++ )
  {
    indx_insert = indx[n];

    printf ( "\n" );
    printf ( "  Inserting value %f\n", a[indx_insert] );

    r8vec_indexed_heap_d_insert ( &n, a, indx, indx_insert );

    indx_max = r8vec_indexed_heap_d_max ( n, a, indx );

    printf ( "  Current maximum is %f\n", a[indx_max] );
  }
  r8vec_print ( m, a, "  The data vector after insertions:" );
  i4vec_print ( n, indx, "  The index vector after insertions:" );
  printf ( "\n" );
  printf ( "  A(INDX) after insertions:\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10f\n", i, a[indx[i]] );
  }

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_indicator0_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDICATOR0_NEW_TEST tests R8VEC_INDICATOR0_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 September 2014

  Author:

    John Burkardt
*/
{
  int n;
  double *v;

  printf ( "\n" );
  printf ( "R8VEC_INDICATOR0_NEW_TEST\n" );
  printf ( "  R8VEC_INDICATOR0_NEW returns an indicator vector.\n" );

  n = 10;
  v = r8vec_indicator0_new ( n );
  r8vec_print ( n, v, "  Indicator0 vector:" );
  free ( v );

  return;
# undef N
}
/******************************************************************************/

void r8vec_legendre_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LEGENDRE_TEST tests R8VEC_LEGENDRE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2011

  Author:

    John Burkardt
*/
{
  int n;
  double *r;
  double r1;
  double r2;

  printf ( "\n" );
  printf ( "R8VEC_LEGENDRE_TEST\n" );
  printf ( "  R8VEC_LEGENDRE_NEW computes N Legendre points in [R1,R2].\n" );

  r1 = -1.0;
  r2 = +1.0;
  n = 5;

  r = r8vec_legendre_new ( n, r1, r2 );

  printf ( "\n" );
  printf ( "  N = %d,   R1 = %g,  R2 = %g\n", n, r1, r2 );

  r8vec_print ( n, r, "  Legendre points:" );

  free ( r );

  r1 =   0.0;
  r2 = +10.0;
  n = 7;

  r = r8vec_legendre_new ( n, r1, r2 );

  printf ( "\n" );
  printf ( "  N = %d,   R1 = %g,  R2 = %g\n", n, r1, r2 );

  r8vec_print ( n, r, "  Legendre points:" );

  free ( r );

  return;
}
/******************************************************************************/

void r8vec_linspace_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW_TEST tests R8VEC_LINSPACE_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int n = 5;
  double *x;

  printf ( "\n" );
  printf ( "R8VEC_LINSPACE_NEW_TEST\n" );
  printf ( "  For a R8VEC:\n" );
  printf ( "  R8VEC_LINSPACE_NEW: evenly spaced points between A and B;\n" );

  a = 10.0;
  b = 20.0;

  x = r8vec_linspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_linspace ( 5, 10, 20 )" );
  free ( x );

  return;
}
/******************************************************************************/

void r8vec_max_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX_TEST tests R8VEC_MAX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int  n;
  double rmax;
  double rmin;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_MAX_TEST\n" );
  printf ( "  R8VEC_MAX produces the maximum entry in a real array.\n" );

  n = 10;
  seed = 123456789;

  a = r8vec_uniform_01_new ( n, &seed );

  r8vec_print ( n, a, "  The array:" );

  rmax = r8vec_max ( n, a );

  printf ( "\n" );
  printf ( "  R8VEC_MAX reports the maximum value is %g.\n", rmax );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_max_index_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX_INDEX_TEST tests R8VEC_MAX_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int ival;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_MAX_INDEX_TEST\n" );
  printf ( "  For an R8VEC:\n" );
  printf ( "  R8VEC_MAX_INDEX: index of maximum entry;\n" );

  b = - ( double ) ( N );
  c =   ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );

  ival = r8vec_max_index ( N, a );
  printf ( "  Maximum index:           %d\n", ival );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_mean_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MEAN_TEST tests R8VEC_MEAN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2015

  Author:

    John Burkardt
*/
{
  double *a;
  double mean;
  int n = 10;
  double r8_hi;
  double r8_lo;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_MEAN_TEST\n" );
  printf ( "  R8VEC_MEAN computes the mean of an R8VEC.\n" );

  r8_lo = - 5.0;
  r8_hi = + 5.0;
  seed = 123456789;
  a = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );

  r8vec_print ( n, a, "  Input vector:" );

  mean = r8vec_mean ( n, a );

  printf ( "\n" );
  printf ( "  Mean:    %g\n", mean );
 
  free ( a );

  return;
}
/******************************************************************************/

void r8vec_median_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MEDIAN_TEST tests R8VEC_MEDIAN;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double b;
  double c;
  double median;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_MEDIAN_TEST\n" );
  printf ( "  For an R8VEC:\n" );
  printf ( "  R8VEC_MEAN:      mean value;\n" );

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );

  median = r8vec_median ( N, a );
  printf ( "  Median:  %g\n", median );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_midspace_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIDSPACE_NEW_TEST tests R8VEC_MIDSPACE_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int n = 5;
  double *x;

  printf ( "\n" );
  printf ( "R8VEC_MIDSPACE_NEW_TEST\n" );
  printf ( "  For a R8VEC:\n" );
  printf ( "  R8VEC_MIDSPACE_NEW: evenly spaced midpoints between A and B\n" );

  a = 10.0;
  b = 20.0;

  x = r8vec_midspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_midspace ( 5, 10, 20 )" );
  free ( x );

  return;
}
/******************************************************************************/

void r8vec_min_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN_TEST tests R8VEC_MIN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int  n;
  double rmax;
  double rmin;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_MIN_TEST\n" );
  printf ( "  R8VEC_MIN produces the minimum entry.\n" );

  n = 10;
  seed = 123456789;

  a = r8vec_uniform_01_new ( n, &seed );

  r8vec_print ( n, a, "  The array:" );

  rmin = r8vec_min ( n, a );

  printf ( "\n" );
  printf ( "  R8VEC_MIN reports the minimum value is %g.\n", rmin );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_min_index_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN_INDEX_TEST tests R8VEC_MIN_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double aval;
  double b;
  double c;
  int ival;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_MIN_INDEX_TEST\n" );
  printf ( "  For an R8VEC:\n" );
  printf ( "  R8VEC_MIN_INDEX: index of minimum entry;\n" );

  b = - ( double ) ( N );
  c =   ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );

  ival = r8vec_min_index ( N, a );
  printf ( "  Minimum index:           %d\n", ival );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_nint_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NINT_TEST tests R8VEC_NINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  int n;
  int seed;
  double x1;
  double x2;

  printf ( "\n" );
  printf ( "R8VEC_NINT_TEST\n" );
  printf ( "  R8VEC_NINT rounds an R8VEC.\n" );

  n = 5;
  x1 = -5.0;
  x2 = +5.0;
  seed = 123456789;
  a = r8vec_uniform_ab_new ( n, x1, x2, &seed );
  r8vec_print ( n, a, "  Vector A:" );
  r8vec_nint ( n, a );
  r8vec_print ( n, a, "  Rounded vector A:" );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_norm_l0_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L0_TEST tests R8VEC_NORM_L0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 January 2015

  Author:

    John Burkardt
*/
{
  double *a;
  double a_hi;
  double a_lo;
  int n;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_NORM_L0_TEST\n" );
  printf ( "  R8VEC_NORM_L0 computes the L0 'norm' of an R8VEC.\n" );

  n = 10;
  a_lo = - 2.0;
  a_hi = + 2.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, a_lo, a_hi, &seed );
  r8vec_nint ( n, a );

  r8vec_print ( n, a, "  Input vector:" );

  printf ( "\n" );
  printf ( "  L0 norm:           %g\n", r8vec_norm_l0 ( n, a ) );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_norm_l1_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L1_TEST tests R8VEC_NORM_L1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_NORM_L1_TEST\n" );
  printf ( "  R8VEC_NORM_L1 computes the L1 norm of an R8VEC.\n" );

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );
  printf ( "  L1 norm:           %g\n", r8vec_norm_l1 ( N, a ) );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_norm_l2_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L2_TEST tests R8VEC_NORM_L2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_NORM_L2_TEST\n" );
  printf ( "  R8VEC_NORM_L2 computes the L2 norm of an R8VEC.\n" );

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );
  printf ( "  L2 norm:           %g\n", r8vec_norm_l2 ( N, a ) );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_norm_li_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_LI_TEST tests R8VEC_NORM_LI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_NORM_LI_TEST\n" );
  printf ( "  R8VEC_NORM_LI computes the Loo norm of an R8VEC.\n" );

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  printf ( "\n" );
  printf ( "  L-Infinity norm:   %g\n", r8vec_norm_li ( N, a ) );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_normal_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMAL_01_NEW_TEST tests R8VEC_NORMAL_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N_MAX 1000

  int i;
  int n;
  int seed = 123456789;
  double *x;
  double x_max;
  double x_mean;
  double x_min;
  double x_var;

  printf ( "\n" );
  printf ( "R8VEC_NORMAL_01_NEW_TEST\n" );
  printf ( "  R8VEC_NORMAL_01 computes a vector of normally\n" );
  printf ( "  distributed random numbers.\n" );
  printf ( "  Using initial random number seed = %d\n", seed );
/*
  Test 1:
  Simply call 5 times for 1 value, and print.
*/
  printf ( "\n" );
  printf ( "  Test #1: Call 5 times, 1 value each time.\n" );
  printf ( "\n" );

  n = 1;
  for ( i = 0; i < 5; i++ )
  {
    x = r8vec_normal_01_new ( n, &seed );
    printf ( "  %6d  %14g\n", i, x[0] );
    free ( x );
  }
/*
  Test 2:
  Restore the random number seed, and repeat.
*/
  printf ( "\n" );
  printf ( "  Test #2: Restore the random number seed.\n" );
  printf ( "  Call 5 times, 1 value each time.\n" );
  printf ( "  The results should be identical.\n" );
  printf ( "\n" );

  seed = 123456789;

  n = 1;
  for ( i = 0; i < 5; i++ )
  {
    x = r8vec_normal_01_new ( n, &seed );
    printf ( "  %6d  %14g\n", i, x[0] );
    free ( x );
  }
/*
  Test 3:
  Restore the random number seed, compute all 5 values at once.
*/
  printf ( "\n" );
  printf ( "  Test #3: Restore the random number seed.\n" );
  printf ( "  Call 1 time for 5 values.\n" );
  printf ( "  The results should be identical.\n" );
  printf ( "\n" );

  seed = 123456789;

  n = 5;
  x = r8vec_normal_01_new ( n, &seed );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14g\n", i, x[i] );
  }
  free ( x );
/*
  Test 4:
  Restore the random number seed, compute all 5 values at once.
*/
  printf ( "\n" );
  printf ( "  Test #4: Restore the random number seed.\n" );
  printf ( "  Call for 2, 1, and 2 values.\n" );
  printf ( "  The results should be identical.\n" );
  printf ( "\n" );

  seed = 123456789;

  n = 2;
  x = r8vec_normal_01_new ( n, &seed );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14g\n", i, x[i] );
  }
  free ( x );

  n = 1;
  x = r8vec_normal_01_new ( n, &seed );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14g\n", i, x[i] );
  }
  free ( x );

  n = 2;
  x = r8vec_normal_01_new ( n, &seed );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %6d  %14g\n", i, x[i] );
  }
  free ( x );
/*
  Test 5:
  Determine the minimum, maximum, mean and variance.
*/
  n = N_MAX;
  x = r8vec_normal_01_new ( n, &seed );
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  x_mean = r8vec_mean ( n, x );
  x_var = r8vec_variance ( n, x );
  free ( x );

  printf ( "\n" );
  printf ( "  Test #5:\n" );
  printf ( "  Number of samples was %d\n", n );
  printf ( "  Minimum value was %g\n", x_min );
  printf ( "  Maximum value was %g\n", x_max );
  printf ( "  Average value was %g\n", x_mean );
  printf ( "  Variance was      %g\n", x_var );
  printf ( "  Expected average  = 0.0\n" );
  printf ( "  Expected variance = 1.0\n" );

  return;
# undef N_MAX
}
/******************************************************************************/

void r8vec_normalize_l1_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMALIZE_L1_TEST tests R8VEC_NORMALIZE_L1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_NORMALIZE_L1_TEST\n" );
  printf ( "  For an R8VEC:\n" );
  printf ( "  R8VEC_NORMALIZE_L1:  make unit sum;\n" );

  b = - ( double ) ( N );
  c = ( double ) ( N );

  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Input vector:" );

  r8vec_normalize_l1 ( N, a );

  r8vec_print ( N, a, "  After calling R8VEC_NORMALIZE_L1:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_order_type_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ORDER_TYPE_TEST tests R8VEC_ORDER_TYPE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 4
# define TEST_NUM 6

  int itest;
  int j;
  int order;
  double x[N];

  printf ( "\n" );
  printf ( "R8VEC_ORDER_TYPE_TEST\n" );
  printf ( "  R8VEC_ORDER_TYPE classifies a real vector as\n" );
  printf ( "  -1: no order\n" );
  printf ( "   0: all equal;\n" );
  printf ( "   1: ascending;\n" );
  printf ( "   2: strictly ascending;\n" );
  printf ( "   3: descending;\n" );
  printf ( "   4: strictly descending.\n" );
  printf ( "\n" );

  for ( itest = 1; itest <= TEST_NUM; itest++ )
  {
    if ( itest == 1 )
    {
      x[0] = 1.0;
      x[1] = 3.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 2 )
    {
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 2.0;
    }
    else if ( itest == 3 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 4 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
    }
    else if ( itest == 5 )
    {
      x[0] = 4.0;
      x[1] = 4.0;
      x[2] = 3.0;
      x[3] = 1.0;
    }
    else if ( itest == 6 )
    {
      x[0] = 9.0;
      x[1] = 7.0;
      x[2] = 3.0;
      x[3] = 0.0;
    }

    order = r8vec_order_type ( N, x );

    printf ( "\n" );
    printf ( "  The following vector has order type %d.\n", order );

    r8vec_print ( N, x, "" );
  }

  return;
# undef N
# undef TEST_NUM
}
/******************************************************************************/

void r8vec_permute_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE_TEST tests R8VEC_PERMUTE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 October 2014

  Author:

    John Burkardt
*/
{
  int i;
  int n = 5;
  int p[5] = { 1, 3, 4, 0, 2 };
  double x[5] = { 1.1, 2.2, 3.3, 4.4, 5.5 };

  printf ( "\n" );
  printf ( "R8VEC_PERMUTE_TEST\n" );
  printf ( "  R8VEC_PERMUTE permutes an R8VEC.\n" );

  r8vec_print ( n, x, "  Original array X[]:" );

  i4vec_print ( n, p, "  Permutation vector P[]:" );

  r8vec_permute ( n, p, x );

  r8vec_print ( n, x, "  Permuted array X[P[]]:" );

  return;
}
/******************************************************************************/

void r8vec_permute_uniform_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE_UNIFORM_TEST tests R8VEC_PERMUTE_UNIFORM.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 May 2015

  Author:

    John Burkardt
*/
{
  double *a;
  int i;
  int n = 10;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_PERMUTE_UNIFORM_TEST\n" );
  printf ( "  R8VEC_PERMUTE_UNIFORM randomly reorders an R8VEC.\n" );

  a = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( double ) ( 101 + i );
  }
  seed = 123456789;

  r8vec_print ( n, a, "  A, before rearrangement:" );

  r8vec_permute_uniform ( n, a, &seed );

  r8vec_print ( n, a, "  A, after random permutation:" );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_polarize_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_POLARIZE_TEST tests R8VEC_POLARIZE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 3

  double a[N] = { 1.0, 2.0,  3.0 };
  double a2[N];
  double a_normal[N];
  double a_parallel[N];
  double ap_norm;
  int i;
  double p[N] = { 3.0, 1.0, -2.0 };
  double p_norm;
  double pan;
  double pap;

  printf ( "\n" );
  printf ( "R8VEC_POLARIZE_TEST\n" );
  printf ( "  R8VEC_POLARIZE decomposes an R8VEC into\n" );
  printf ( "  components parallel and normal to a direction.\n" );

  r8vec_print ( N, a, "  Original vector:" );

  r8vec_print ( N, p, "  Direction vector:" );

  r8vec_polarize ( N, a, p, a_normal, a_parallel );

  r8vec_print ( N, a_normal, "  Normal component:" );

  r8vec_print ( N, a_parallel, "  Parallel component:" );

  pan = r8vec_dot_product ( N, p, a_normal );
  p_norm = r8vec_norm ( N, p );
  ap_norm = r8vec_norm ( N, a_parallel );

  pap = r8vec_dot_product ( N, p, a_parallel ) / ( p_norm * ap_norm );

  printf ( "\n" );
  printf ( "  Dot product of P and A_normal (should be 0) %g\n", pan );
  printf ( "  Cosine of angle between P and A_parallel (should be 1 or -1) = %g\n", pap );

  for ( i = 0; i < N; i++ )
  {
    a2[i] = a_normal[i] + a_parallel[i];
  }

  r8vec_print ( N, a2, "  Sum of components (should equal A):" );

  return;
# undef N
}
/******************************************************************************/

void r8vec_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_TEST tests R8VEC_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 August 2014

  Author:

    John Burkardt
*/
{
  double a[4] = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
  int n = 4;

  printf ( "\n" );
  printf ( "R8VEC_PRINT_TEST\n" );
  printf ( "  R8VEC_PRINT prints an R8VEC.\n" );

  r8vec_print ( n, a, "  The R8VEC:" );

  return;
}
/******************************************************************************/

void r8vec_rotate_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ROTATE_TEST tests R8VEC_ROTATE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double a[N] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  int m = 2;

  printf ( "\n" );
  printf ( "R8VEC_ROTATE_TEST\n" );
  printf ( "  R8VEC_ROTATE rotates an R8VEC in place.\n" );
  printf ( "\n" );
  printf ( "  Rotate entries %d places to the right.\n", m );

  r8vec_print ( N, a, "  Original array:" );

  r8vec_rotate ( N, a, m );

  r8vec_print ( N, a, "  Rotated array:" );

  return;
# undef N
}
/******************************************************************************/

void r8vec_reverse_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_REVERSE_TEST tests R8VEC_REVERSE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 5

  double *a;

  printf ( "\n" );
  printf ( "R8VEC_REVERSE_TEST\n" );
  printf ( "  R8VEC_REVERSE reverses an R8VEC.\n" );

  a = r8vec_indicator1_new ( N );

  r8vec_print ( N, a, "  Original array:" );

  r8vec_reverse ( N, a );

  r8vec_print ( N, a, "  Reversed array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_search_binary_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SEARCH_BINARY_A_TEST tests R8VEC_SEARCH_BINARY_A;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a;
  int index;
  int nc;
  double search_val;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8VEC_SEARCH_BINARY_A_TEST\n" );
  printf ( "  For ascending order:\n" );
  printf ( "  R8VEC_SEARCH_BINARY_A searches a sorted R8VEC;\n" );
  printf ( "\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  a = r8vec_uniform_01_new ( N, &seed );

  search_val = a[0];

  r8vec_sort_heap_a ( N, a );

  r8vec_print ( N, a, "  Sorted vector A:" );
/*
  Now search the sorted array for a given value.
*/
  printf ( "\n" );
  printf ( "  Search the array for the value %g\n", search_val );

  index = r8vec_search_binary_a ( N, a, search_val );

  printf ( "\n" );
  printf ( "  SEARCH RESULT:\n" );
  printf ( "\n" );

  if ( 0 < index )
  {
    printf ( "    The value occurs in index %d\n", index );
  }
  else
  {
    printf ( "    The value does not occur in the array.\n" );
  }

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_bubble_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_BUBBLE_A_TEST tests R8VEC_SORT_BUBBLE_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  int i;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_BUBBLE_A_TEST\n" );
  printf ( "  R8VEC_SORT_BUBBLE_A sorts an R8VEC.\n" );

  seed = 123456789;

  printf ( "\n" );
  printf ( "  Using random seed %d.\n", seed );

  a = r8vec_uniform_01_new ( N, &seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_bubble_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_heap_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_A_TEST tests R8VEC_SORT_HEAP_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_HEAP_A_TEST\n" );
  printf ( "  R8VEC_SORT_HEAP_A ascending sorts an R8VEC.\n" );
  printf ( "\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print_some ( N, a, 1, 10, "  Original array:" );

  r8vec_sort_heap_a ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Ascending sorted array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_heap_d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_D_TEST tests R8VEC_SORT_HEAP_D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_HEAP_D_TEST\n" );
  printf ( "  R8VEC_SORT_HEAP_D descending sorts an R8VEC.\n" );
  printf ( "\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print_some ( N, a, 1, 10, "  Original array:" );

  r8vec_sort_heap_d ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Descending sorted array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_heap_index_a_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_A_NEW_TEST tests R8VEC_SORT_HEAP_INDEX_A_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_HEAP_INDEX_A_NEW_TEST\n" );
  printf ( "  R8VEC_SORT_HEAP_INDEX_A_NEW creates an ascending\n" );
  printf ( "  sort index for an R8VEC.\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  indx = r8vec_sort_heap_index_a_new ( N, a );

  printf ( "\n" );
  printf ( "  After indexed ascending sort:\n" );
  printf ( "\n" );
  printf ( "         I INDX(I)       A(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %8d  %12g\n", i, indx[i], a[i] );
  }

  printf ( "\n" );
  printf ( "  Now use the index array to carry out the\n" );
  printf ( "  permutation implicitly.\n" );
  printf ( "\n" );
  printf ( "   INDX(I)  A(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %12g\n", indx[i], a[indx[i]] );
  }
  printf ( "\n" );
  printf ( "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n" );
  printf ( "\n" );

  r8vec_permute ( N, indx, a );

  r8vec_print ( N, a, "  I, A(I)" );

  free ( a );
  free ( indx );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_heap_index_d_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_D_NEW_TEST tests R8VEC_SORT_HEAP_INDEX_D_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_HEAP_INDEX_D_NEW_TEST\n" );
  printf ( "  R8VEC_SORT_HEAP_INDEX_D_NEW creates a descending\n" );
  printf ( "  sort index for an R8VEC.\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  indx = r8vec_sort_heap_index_d_new ( N, a );

  printf ( "\n" );
  printf ( "  After indexed descending sort:\n" );
  printf ( "\n" );
  printf ( "         I  INDX(I)  A(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %8d  %12g\n", i, indx[i], a[i] );
  }

  printf ( "\n" );
  printf ( "  Now use the index array to carry out the\n" );
  printf ( "  permutation implicitly.\n" );
  printf ( "\n" );
  printf ( "   INDX(I)  ARRAY(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %12g\n", indx[i], a[indx[i]] );
  }

  free ( a );
  free ( indx );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_heap_mask_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_MASK_A_TEST tests R8VEC_SORT_HEAP_MASK_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define MASK_NUM 10
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int mask[MASK_NUM] = { 1, 3, 6, 7, 8, 11, 12, 15, 17, 18 };
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_HEAP_MASK_A_TEST\n" );
  printf ( "  R8VEC_SORT_HEAP_MASK_A creates an ascending\n" );
  printf ( "  sort index for a masked R8VEC.\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  i4vec_print ( MASK_NUM, mask, "  The mask array:" );

  r8vec_mask_print ( N, a, MASK_NUM, mask, "  The masked unsorted array:" );

  indx = r8vec_sort_heap_mask_a ( N, a, MASK_NUM, mask );

  printf ( "\n" );
  printf ( "  After masked indexed ascending sort:\n" );
  printf ( "\n" );
  printf ( "  I, INDX(I), MASK(INDX(I)), A(MASK(INDX(I)))\n" );
  printf ( "\n" );
  for ( i = 0; i < MASK_NUM; i++ )
  {
    printf ( "  %6d  %6d  %6d  %14g\n", 
      i, indx[i], mask[indx[i]], a[mask[indx[i]]] );
  }

  printf ( "\n" );
  printf ( "  Call I4VEC_PERMUTE to carry out the index permutation\n" );
  printf ( "  explicitly on the MASK vector.\n" );
  printf ( "\n" );

  i4vec_permute ( MASK_NUM, indx, mask );
/*
  Essentially, INDX becomes the identity vector now.
*/
  free ( indx );

  indx = i4vec_indicator1_new ( MASK_NUM );

  i4vec_print ( MASK_NUM, mask, "  The reordered mask array:" );

  r8vec_mask_print ( N, a, MASK_NUM, mask,
    "  The reordered masked sorted array:" );

  free ( a );
  free ( indx );

  return;
# undef MASK_NUM
# undef N
}
/******************************************************************************/

void r8vec_sort_insert_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_INSERT_A_TEST tests R8VEC_SORT_INSERT_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_INSERT_A_TEST\n" );
  printf ( "  R8VEC_SORT_INSERT_A ascending sorts an R8VEC.\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print_some ( N, a, 1, 10, "  Unsorted array:" );

  r8vec_sort_insert_a ( N, a );

  r8vec_print_some ( N, a, 1, 10, "  Sorted array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_insert_index_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_INSERT_INDEX_A_TEST tests R8VEC_SORT_INSERT_INDEX_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int i;
  int *indx;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_INSERT_INDEX_A_TEST\n" );
  printf ( "  R8VEC_SORT_INSERT_INDEX_A creates an ascending\n" );
  printf ( "  sort index for an R8VEC.\n" );

  b = 0.0;
  c = 3.0 * ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print_some ( N, a, 1, 10, "  Unsorted array:" );

  indx = r8vec_sort_insert_index_a ( N, a );

  printf ( "\n" );
  printf ( "  After indexed ascending sort:\n" );
  printf ( "\n" );
  printf ( "  I, INDX(I), A(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6d  %12g\n", i+1, indx[i], a[i] );
  }

  printf ( "\n" );
  printf ( "  Now use the index array to carry out the\n" );
  printf ( "  permutation implicitly.\n" );
  printf ( "\n" );
  printf ( "  I, INDX(I), A(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %6d  %6d  %12g\n", i+1, indx[i], a[indx[i]] );
  }

  printf ( "\n" );
  printf ( "  Call R8VEC_PERMUTE to carry out the permutation explicitly.\n" );
  printf ( "\n" );

  r8vec_permute ( N, indx, a );

  r8vec_print_some ( N, a, 1, 10, "  Permuted data" );

  free ( a );
  free ( indx );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sort_quick_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_QUICK_A_TEST tests R8VEC_SORT_QUICK_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  double *a;
  double b;
  double c;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_SORT_QUICK_A_TEST\n" );
  printf ( "  R8VEC_SORT_QUICK_A sorts an R8VEC\n" );
  printf ( "  using quick sort.\n" );

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_quick_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sorted_merge_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_MERGE_A_TEST tests R8VEC_SORTED_MERGE_A;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double *b;
  double *c;
  int na = 10;
  int nb = 10;
  int nc;
  int seed = 123456789;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_MERGE_A_TEST\n" );
  printf ( "  For ascending order:\n" );
  printf ( "  R8VEC_SORTED_MERGE_A merges two sorted R8VEC's;\n" );
  printf ( "\n" );
  printf ( "  Using initial random number seed = %d\n", seed );

  a = r8vec_uniform_01_new ( na, &seed );
  b = r8vec_uniform_01_new ( nb, &seed );

  r8vec_sort_heap_a ( na, a );

  r8vec_sort_heap_a ( nb, b );

  r8vec_print ( na, a, "  Sorted vector A:" );

  r8vec_print ( nb, b, "  Sorted vector B:" );

  c = r8vec_sorted_merge_a ( na, a, nb, b, &nc );

  r8vec_print ( nc, c, "  Merged vector C:" );

  free ( a );
  free ( b );
  free ( c );

  return;
}
/******************************************************************************/

void r8vec_sorted_nearest_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_NEAREST_TEST tests R8VEC_SORTED_NEAREST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double b;
  double c;
  int i;
  int j;
  int seed = 123456789;
  double *x;
  double xval;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_NEAREST_TEST\n" );
  printf ( "  R8VEC_SORTED_NEAREST finds the nearest entry\n" );
  printf ( "  in a sorted R8VEC.\n" );

  b = 0.0;
  c = 10.0;

  x = r8vec_uniform_ab_new ( N, b, c, &seed );

  r8vec_sort_heap_a ( N, x );

  r8vec_print ( N, x, "  Sorted array:" );

  printf ( "\n" );
  printf ( "     Test        Nearest\n" );
  printf ( "     Value    Index   Value\n" );
  printf ( "\n" );
  for ( i = 1; i <= 10; i++ )
  {
    xval = r8_uniform_01 ( &seed );

    j = r8vec_sorted_nearest ( N, x, xval );

    printf ( "  %8g  %6d  %8g\n", xval, j, x[j-1] );
  }

  free ( x );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sorted_range_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_RANGE_TEST tests R8VEC_SORTED_RANGE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2010

  Author:

    John Burkardt
*/
{
  int i;
  int i_hi;
  int i_lo;
  int n = 10;
  double r[10];
  double r_lo;
  double r_hi;
  int seed;
  double t;
  int test;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_RANGE_TEST\n" );
  printf ( "  R8VEC_SORTED_RANGE seeks the range of indices\n" );
  printf ( "  in a sorted R8VEC so that\n" );
  printf ( "  R_LO <= R(I_LO:I_HI) <= R_HI.\n" );

  seed = 123456789;

  for ( test = 1; test <= 5; test++ )
  {
    r8vec_uniform_01 ( n, &seed, r );

    r8vec_sort_heap_a ( n, r );

    r8vec_print ( n, r, "  Sorted array R:" );

    r_lo = r8_uniform_01 ( &seed );
    r_hi = r8_uniform_01 ( &seed );

    if ( r_hi < r_lo )
    {
      t = r_lo;
      r_lo = r_hi;
      r_hi = t;
    }

    r8vec_sorted_range ( n, r, r_lo, r_hi, &i_lo, &i_hi );

    printf ( "\n" );
    if ( i_hi < i_lo )
    {
      printf ( "  R_LO  %14f\n", r_lo );
      printf ( "  R_HI  %14f\n", r_hi );
      printf ( "  Empty range in R.\n" );
    }
    else
    {
      printf ( "  R_LO  %14f\n", r_lo );
      for ( i = i_lo; i <= i_hi; i++)
      {
        printf ( "  %4d  %14f\n", i, r[i] );
      }
      printf ( "  R_HI  %14f\n", r_hi );
    }
  }
  return;
}
/******************************************************************************/

void r8vec_sorted_split_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_SPLIT_TEST tests R8VEC_SORTED_SPLIT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b;
  double c;
  int i;
  int i_gt;
  int i_lt;
  int isplit;
  int n = 25;
  int seed;
  double split;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_SPLIT_TEST\n" );
  printf ( "  R8VEC_SORTED_SPLIT splits a sorted R8VEC into\n" );
  printf ( "  entries less than and greater than a\n" );
  printf ( "  splitting value.\n" );
  printf ( "\n" );

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, &seed );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.5 * ( double ) ( r8_nint ( a[i] ) );
  }

  r8vec_sort_heap_a ( n, a );

  split = 0.5 * ( a[0] + a[n-1] );

  r8vec_print ( n, a, "  The sorted array:" );

  printf ( "\n" );
  printf ( "  Splitting value is %g\n", split );
  printf ( "\n" );

  r8vec_sorted_split ( n, a, split, &i_lt, &i_gt );

  printf ( "  Lower index I_LT = %d\n", i_lt );
  printf ( "  Upper index I_GT = %d\n", i_gt );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_sorted_undex_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNDEX_TEST tests R8VEC_SORTED_UNDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define X_NUM 9

  int i;
  double tol;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  double x_val[X_NUM] = { 11.0, 11.0, 11.0, 22.0, 22.0, 33.0, 33.0, 55.0, 55.0 };
  int *xdnu;
  double *xu_val;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_UNDEX_TEST\n" );
  printf ( "  R8VEC_SORTED_UNDEX produces index vectors which create a sorted\n" );
  printf ( "  list of the unique elements of a sorted R8VEC,\n" );
  printf ( "  and a map from the original vector to the (implicit)\n" );
  printf ( "  vector of sorted unique elements.\n" );

  r8vec_print ( x_num, x_val, "  The vector X:" );

  tol = r8_epsilon ( );
  x_unique_num = r8vec_sorted_unique_count ( x_num, x_val, tol );

  undx = ( int * ) malloc ( x_unique_num * sizeof ( int ) );
  xu_val = ( double * ) malloc ( x_unique_num * sizeof ( double ) );

  xdnu = ( int * ) malloc ( x_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "  Number of unique entries in X is %d\n", x_unique_num );

  r8vec_sorted_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu );

  printf ( "\n" );
  printf ( "  UNDX can be used to list the unique elements of X\n" );
  printf ( "  in sorted order.\n" );
  printf ( "\n" );
  printf ( "     I  UNDX   X(UNDX)\n" );
  printf ( "\n" );

  for ( i = 0; i < x_unique_num; i++ )
  {
    printf ( "  %4d  %4d  %8g\n", i, undx[i], x_val[undx[i]] );
  }

  for ( i = 0; i < x_unique_num; i++ )
  {
    xu_val[i] = x_val[undx[i]];
  }
  printf ( "\n" );
  printf ( "  UNDX can be used to created XU, a copy of X\n" );
  printf ( "  containing only the unique elements, in sorted order.\n" );
  printf ( "\n" );
  printf ( "     I  UNDX     XU(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < x_unique_num; i++ )
  {
    printf ( "  %4d  %4d  %8g\n", i, undx[i], xu_val[i] );
  }

  printf ( "\n" );
  printf ( "  XDNU can be used to match each element of X with one of the\n" );
  printf ( "  unique elements\n" );
  printf ( "\n" );
  printf ( "     I  XDNU  X(I)   XU(XDNU(I))\n" );
  printf ( "\n" );

  for ( i = 0; i < x_num; i++ )
  {
    printf ( "  %4d  %4d  %6g  %12g", i, xdnu[i], x_val[i], xu_val[xdnu[i]] );
  }

  free ( undx );
  free ( xdnu );
  free ( xu_val );

  return;
# undef X_NUM
}
/******************************************************************************/

void r8vec_sorted_unique_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNIQUE_TEST tests R8VEC_SORTED_UNIQUE;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 December 2014

  Author:

    John Burkardt
*/
{
  double *a;
  double *a_unique;
  double b;
  double c;
  int i;
  int n = 20;
  int seed;
  double tol = 0.25;
  int unique_num;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_UNIQUE_TEST\n" );
  printf ( "  R8VEC_SORTED_UNIQUE finds unique entries in a sorted R8VEC;\n" );

  b = 0.0;
  c = ( double ) ( n );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, &seed );

  r8vec_nint ( n, a );

  r8vec_print_some ( n, a, 1, 10, "  Unsorted array:" );

  r8vec_sort_heap_a ( n, a );

  a_unique = r8vec_sorted_unique ( n, a, tol, &unique_num );

  r8vec_print ( unique_num, a_unique, "  Unique entries" );

  free ( a );
  free ( a_unique );

  return;
}
/******************************************************************************/

void r8vec_sorted_unique_count_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNIQUE_COUNT_TEST tests R8VEC_SORTED_UNIQUE_COUNT;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 30

  double *a;
  double b;
  double c;
  int i;
  int unique_num;
  int seed;
  double tol = 0.25;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_UNIQUE_COUNT_TEST\n" );
  printf ( "  R8VEC_SORTED_UNIQUE_COUNT counts unique entries in a sorted R8VEC;\n" );

  b = 0.0;
  c = ( double ) ( N );
  seed = 123456789;

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  for ( i = 0; i < N; i++ )
  {
    a[i] = ( double ) r8_nint ( a[i] );
  }

  unique_num = r8vec_sorted_unique_count ( N, a, tol );

  printf ( "\n" );
  printf ( "  Using a tolerance of %g\n", tol );
  printf ( "  R8VEC_SORTED_UNIQUE_COUNT counts %d unique entries in A.\n", unique_num );

  free ( a );

  return;
# undef N
}
/******************************************************************************/

void r8vec_sorted_unique_hist_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNIQUE_HIST_TEST tests R8VEC_SORTED_UNIQUE_HIST.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define MAXUNIQ 30
# define N 30

  double *a;
  int acount[MAXUNIQ];
  double auniq[MAXUNIQ];
  double b;
  double c;
  int i;
  int unique_num;
  int seed;
  double tol = 0.25;

  printf ( "\n" );
  printf ( "R8VEC_SORTED_UNIQUE_HIST_TEST\n" );
  printf ( "  R8VEC_SORTED_UNIQUE_HIST stores the unique entries\n" );
  printf ( "  and their multiplicities.\n" );

  b = 0.0;
  c = ( double ) N;
  seed = 123456789;

  printf ( "\n" );
  printf ( "  Using random seed %d.\n", seed );

  a = r8vec_uniform_ab_new ( N, b, c, &seed );

  for ( i = 0; i < N; i++ )
  {
    a[i] = ( ( double ) ( ( int ) a[i] ) ) + 0.5;
  }

  r8vec_print ( N, a, "  Unsorted array:" );

  r8vec_sort_bubble_a ( N, a );

  r8vec_print ( N, a, "  Sorted array:" );

  r8vec_sorted_unique_hist ( N, a, tol, MAXUNIQ, &unique_num, auniq, acount );

  printf ( "\n" );
  printf ( "  R8VEC_SORTED_UNIQUE_HIST counts %d unique entries.\n", unique_num );
  printf ( "\n" );
  printf ( "  Value  Multiplicity\n" );
  printf ( "\n" );
  for ( i = 0; i < unique_num; i++ )
  {
    printf ( "%6d  %12g  %6d\n",  i, auniq[i], acount[i] );
  }

  free ( a );

  return;
# undef MAXUNIQ
# undef N
}
/******************************************************************************/

void r8vec_split_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SPLIT_TEST tests R8VEC_SPLIT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
  double *a;
  double b;
  double c;
  int i;
  int i_gt;
  int i_lt;
  int isplit;
  int n = 25;
  int seed;
  double split;

  printf ( "\n" );
  printf ( "R8VEC_SPLIT_TEST\n" );
  printf ( "  R8VEC_SPLIT splits an R8VEC into\n" );
  printf ( "  entries less than and greater than a\n" );
  printf ( "  splitting value.\n" );
  printf ( "\n" );

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a = r8vec_uniform_ab_new ( n, b, c, &seed );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.5 * ( double ) ( r8_nint ( a[i] ) );
  }

  split = 0.5 * ( a[0] + a[n-1] );

  r8vec_print ( n, a, "  The array:" );

  isplit = r8vec_split ( n, a, split );

  r8vec_print ( n, a, "  The split array:" );

  printf ( "\n" );
  printf ( "  Array entries <= SPLIT up to index %d\n", isplit );

  free ( a );

  return;
}
/******************************************************************************/

void r8vec_transpose_print_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_TRANSPOSE_PRINT_TEST tests R8VEC_TRANSPOSE_PRINT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt
*/
{
  int n = 12;
  int seed;
  double *x;

  seed = 123456789;

  printf ( "\n" );
  printf ( "R8VEC_TRANSPOSE_PRINT_TEST\n" );
  printf ( "  R8VEC_TRANSPOSE_PRINT prints an R8VEC \"tranposed\",\n" );
  printf ( "  that is, placing multiple entries on a line.\n" );

  x = r8vec_uniform_01_new ( n, &seed );

  r8vec_transpose_print ( n, x, "  The vector X:" );

  free ( x );

  return;
}
/******************************************************************************/

void r8vec_undex_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNDEX_TEST tests R8VEC_UNDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define X_NUM 9

  int i;
  double tol;
  int *undx;
  int x_num = X_NUM;
  int x_unique_num;
  double x_val[X_NUM] = { 33.0, 55.0, 11.0, 11.0, 55.0, 33.0, 22.0, 22.0, 11.0 };
  int *xdnu;
  double *xu_val;

  printf ( "\n" );
  printf ( "R8VEC_UNDEX_TEST\n" );
  printf ( "  R8VEC_UNDEX produces index vectors which create a sorted\n" );
  printf ( "  list of the unique elements of an (unsorted) R8VEC,\n" );
  printf ( "  and a map from the original vector to the (implicit)\n" );
  printf ( "  vector of sorted unique elements.\n" );

  r8vec_print ( x_num, x_val, "  The vector X:" );

  tol = r8_epsilon ( );
  x_unique_num = r8vec_unique_count ( x_num, x_val, tol );

  undx = ( int * ) malloc ( x_unique_num * sizeof ( int ) );
  xu_val = ( double * ) malloc ( x_unique_num * sizeof ( double ) );

  xdnu = ( int * ) malloc ( x_num * sizeof ( int ) );

  printf ( "\n" );
  printf ( "  Number of unique entries in X is %d\n", x_unique_num );

  r8vec_undex ( x_num, x_val, x_unique_num, tol, undx, xdnu );

  printf ( "\n" );
  printf ( "  UNDX can be used to list the unique elements of X\n" );
  printf ( "  in sorted order.\n" );
  printf ( "\n" );
  printf ( "     I  UNDX   X(UNDX)\n" );
  printf ( "\n" );

  for ( i = 0; i < x_unique_num; i++ )
  {
    printf ( "  %4d  %4d  %8g\n", i, undx[i], x_val[undx[i]] );
  }

  for ( i = 0; i < x_unique_num; i++ )
  {
    xu_val[i] = x_val[undx[i]];
  }
  printf ( "\n" );
  printf ( "  UNDX can be used to created XU, a copy of X\n" );
  printf ( "  containing only the unique elements, in sorted order.\n" );
  printf ( "\n" );
  printf ( "     I  UNDX     XU(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < x_unique_num; i++ )
  {
    printf ( "  %4d  %4d  %8g\n", i, undx[i], xu_val[i] );
  }

  printf ( "\n" );
  printf ( "  XDNU can be used to match each element of X with one of the\n" );
  printf ( "  unique elements\n" );
  printf ( "\n" );
  printf ( "     I  XDNU  X(I)   XU(XDNU(I))\n" );
  printf ( "\n" );

  for ( i = 0; i < x_num; i++ )
  {
    printf ( "  %4d  %4d  %4g  %12g\n", i, xdnu[i], x_val[i], xu_val[xdnu[i]] );
  }

  free ( undx );
  free ( xdnu );
  free ( xu_val );

  return;
# undef X_NUM
}
/******************************************************************************/

void r8vec_uniform_01_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW_TEST tests R8VEC_UNIFORM_01_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  int j;
  double *r;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_UNIFORM_01_NEW_TEST\n" );
  printf ( "  R8VEC_UNIFORM_01_NEW returns a random R8VEC\n" );
  printf ( "  with entries in a given range [ 0.0, 1.0 ]\n" );
  printf ( "\n" );

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    printf ( "\n" );
    printf ( "  Input SEED = %d\n", seed );
    printf ( "\n" );

    r = r8vec_uniform_01_new ( N, &seed );

    r8vec_print ( N, r, "  Random R8VEC:" );

    free ( r );
  }

  return;
# undef N
}
/******************************************************************************/

void r8vec_uniform_ab_new_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB_NEW_TEST tests R8VEC_UNIFORM_AB_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double a = 10.0;
  double b = 20.0;
  int j;
  double *r;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC_UNIFORM_AB_NEW_TEST\n" );
  printf ( "  R8VEC_UNIFORM_AB_NEW returns a random R8VEC\n" );
  printf ( "  with entries in a given range [ A, B ]\n" );
  printf ( "\n" );
  printf ( "  For this problem:\n" );
  printf ( "  A = %g\n", a );
  printf ( "  B = %g\n", b );
  printf ( "\n" );

  seed = 123456789;

  for ( j = 1; j <= 3; j++ )
  {
    printf ( "\n" );
    printf ( "  Input SEED = %d\n", seed );
    printf ( "\n" );

    r = r8vec_uniform_ab_new ( N, a, b, &seed );

    r8vec_print ( N, r, "  Random R8VEC:" );

    free ( r );
  }

  return;
# undef N
}
/******************************************************************************/

void r8vec_variance_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC_VARIANCE_TEST tests R8VEC_VARIANCE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2015

  Author:

    John Burkardt
*/
{
  double *a;
  int n = 10;
  double r8_hi;
  double r8_lo;
  int seed;
  double variance;

  printf ( "\n" );
  printf ( "R8VEC_VARIANCE_TEST\n" );
  printf ( "  R8VEC_VARIANCE computes the variance of an R8VEC.\n" );

  r8_lo = - 5.0;
  r8_hi = + 5.0;
  seed = 123456789;
  a = r8vec_uniform_ab_new ( n, r8_lo, r8_hi, &seed );

  r8vec_print ( n, a, "  Input vector:" );

  variance = r8vec_variance ( n, a );

  printf ( "\n" );
  printf ( "  Variance:    %g\n", variance );
 
  free ( a );

  return;
}
/******************************************************************************/

void r8vec2_sort_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORT_A_TEST tests R8VEC2_SORT_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC2_SORT_A_TEST\n" );
  printf ( "  For a pair of R8VEC's:\n" );
  printf ( "  R8VEC2_SORT_A ascending sorts;\n" );

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, &seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, &seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_a ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after ascending sort:" );

  free ( a1 );
  free ( a2 );

  return;
}
/******************************************************************************/

void r8vec2_sort_d_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORT_D_TEST tests R8VEC2_SORT_D.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC2_SORT_D_TEST\n" );
  printf ( "  For a pair of R8VEC's:\n" );
  printf ( "  R8VEC2_SORT_D descending sorts;\n" );

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, &seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, &seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_d ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after descending sort:" );

  free ( a1 );
  free ( a2 );

  return;
}
/******************************************************************************/

void r8vec2_sort_heap_index_a_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORT_HEAP_INDEX_A_TEST tests R8VEC2_SORT_HEAP_INDEX_A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 20

  int i;
  int *indx;
  int seed = 123456789;
  double x[N];
  double y[N];

  printf ( "\n" );
  printf ( "R8VEC2_SORT_HEAP_INDEX_A_TEST\n" );
  printf ( "  R8VEC2_SORT_HEAP_INDEX_A creates a sort index\n" );
  printf ( "  for an (X,Y) array.\n" );

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) i4_uniform_ab ( 0, N, &seed ) / ( double ) N;
    y[i] = ( double ) i4_uniform_ab ( 0, N, &seed ) / ( double ) N;
  }

  printf ( "\n" );
  printf ( "  The unsorted array:\n" );
  printf ( "\n" );
  printf ( "         I  X(I), Y(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %12g  %12g", i, x[i], y[i] );
  }

  indx = r8vec2_sort_heap_index_a ( N, x, y );

  printf ( "\n" );
  printf ( "  After sorting:\n" );
  printf ( "\n" );
  printf ( "         I  INDX(I), X(I), Y(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %8d  %12g  %12g\n", i, indx[i], x[i], y[i] );
  }

  printf ( "\n" );
  printf ( "  Now use the index array to carry out the\n" );
  printf ( "  permutation implicitly.\n" );
  printf ( "\n" );
  printf ( "         I  INDX(I), X(INDX(I)), Y(INDX(I))\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %8d  %12g  %12g\n", i, indx[i], x[indx[i]], y[indx[i]] );
  }

  printf ( "\n" );
  printf ( "  R8VEC_PERMUTE carries out the permutation.\n" );

  r8vec_permute ( N, indx, x );
  r8vec_permute ( N, indx, y );

  printf ( "\n" );
  printf ( "         I X(I), Y(I)\n" );
  printf ( "\n" );
  for ( i = 0; i < N; i++ )
  {
    printf ( "  %8d  %12g  %12g\n", i, x[i], y[i] );
  }

  return;
# undef N
}
/******************************************************************************/

void r8vec2_sorted_unique_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORTED_UNIQUE_TEST tests R8VEC2_SORTED_UNIQUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
  double *a1;
  double *a2;
  double b;
  double c;
  int n = 10;
  int seed;
  int unique_num;

  printf ( "\n" );
  printf ( "R8VEC2_SORTED_UNIQUE_TEST\n" );
  printf ( "  For a pair of R8VEC's:\n" );
  printf ( "  R8VEC2_SORTED_UNIQUE counts unique entries.\n" );

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, &seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, &seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( n, a1, a2, "  The pair of arrays:" );

  r8vec2_sort_a ( n, a1, a2 );

  r8vec2_print ( n, a1, a2, "  Arrays after ascending sort:" );

  r8vec2_sorted_unique ( n, a1, a2, &unique_num );

  r8vec2_print ( unique_num, a1, a2, "  UNIQed array:" );

  free ( a1 );
  free ( a2 );

  return;
}
/******************************************************************************/

void r8vec2_sorted_unique_index_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORTED_UNIQUE_INDEX_TEST tests R8VEC2_SORTED_UNIQUE_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2012

  Author:

    John Burkardt
*/
{
# define N 10

  double *a1;
  double *a2;
  double b;
  double c;
  int indx[N];
  int seed;
  int unique_num;

  printf ( "\n" );
  printf ( "R8VEC2_SORTED_UNIQUE_INDEX_TEST\n" );
  printf ( "  For a pair of R8VEC's:\n" );
  printf ( "  R8VEC2_SORTED_UNIQUE_INDEX indexes unique entries.\n" );

  b = 1.0;
  c = 3.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( N, b, c, &seed );

  b = 5.0;
  c = 10.0;

  a2 = r8vec_uniform_ab_new ( N, b, c, &seed );

  a1[2] = a1[0];
  a2[2] = a2[0];

  a1[5] = a1[1];
  a2[5] = a2[1];

  a1[8] = a1[0];
  a2[8] = a2[0];

  r8vec2_print ( N, a1, a2, "  The pair of arrays:" );

  r8vec2_sorted_unique_index ( N, a1, a2, &unique_num, indx );

  printf ( "\n" );
  printf ( "  The number of unique elements is %d\n", unique_num );

  i4vec_print ( unique_num, indx, "  Index of Unique Elements:" );

  r8vec_index_order ( unique_num, a1, indx );
  r8vec_index_order ( unique_num, a2, indx );

  r8vec2_print ( unique_num, a1, a2, "  After Indexed Nonunique Deletion." );

  free ( a1 );
  free ( a2 );

  return;
# undef N
}
/******************************************************************************/

void r8vec2_sum_max_index_test ( )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SUM_MAX_INDEX_TEST tests R8VEC2_SUM_MAX_INDEX.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt
*/
{
  double *a1;
  double *a2;
  double b;
  double c;
  int ival;
  int n = 10;
  int seed;

  printf ( "\n" );
  printf ( "R8VEC2_SUM_MAX_INDEX_TEST\n" );
  printf ( "  For a pair of R8VEC's:\n" );
  printf ( "  R8VEC2_SUM_MAX_INDEX: index of the sum vector\n" );
  printf ( "  with maximum value.\n" );

  b = 0.0;
  c = 10.0;
  seed = 123456789;

  a1 = r8vec_uniform_ab_new ( n, b, c, &seed );

  b = 0.0;
  c = 5.0;

  a2 = r8vec_uniform_ab_new ( n, b, c, &seed );

  r8vec2_print ( n, a1, a2, "  The pair of vectors:" );

  ival = r8vec2_sum_max_index ( n, a1, a2 );

  printf ( "\n" );
  printf ( "  Index of maximum in A+B: %d\n", ival );

  free ( a1 );
  free ( a2 );

  return;
}
/******************************************************************************/

void roots_to_r8poly_test ( )

/******************************************************************************/
/*
  Purpose:

    ROOTS_TO_R8POLY_TEST tests ROOTS_TO_R8POLY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 March 2015

  Author:

    John Burkardt
*/
{
  double *c;
  int n = 5;
  double x[5] = { 1.0, -4.0, 3.0, 0.0, 3.0 };

  printf ( "\n" );
  printf ( "ROOTS_TO_R8POLY_TEST:\n" );
  printf ( "  ROOTS_TO_R8POLY is given N real roots,\n" );
  printf ( "  and constructs the coefficient vector\n" );
  printf ( "  of the corresponding polynomial.\n" );

  r8vec_print ( n, x, "  N real roots:" );

  c = roots_to_r8poly ( n, x );

  r8poly_print ( n, c, "  Corresponding polynomial:" );

  free ( c );

  return;
}
