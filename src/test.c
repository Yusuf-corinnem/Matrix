#include <check.h>

#include "s21_matrix.h"

START_TEST(s21_create_matrix_test) {
  matrix_t A = {0};
  int reza = s21_create_matrix(2, 2, &A);
  ck_assert_int_eq(0, reza);
  s21_remove_matrix(&A);

  matrix_t B = {0};
  int reza2 = s21_create_matrix(-2, -2, &B);
  ck_assert_int_eq(1, reza2);
  s21_remove_matrix(&B);

  matrix_t C = {0};
  int reza3 = s21_create_matrix(0, 2, &C);
  ck_assert_int_eq(1, reza3);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(s21_remove_matrix_test) {
  matrix_t A = {0};
  s21_create_matrix(3, 3, &A);
  s21_remove_matrix(&A);
  ck_assert_int_eq(A.rows, 0);
  ck_assert_int_eq(A.columns, 0);
  ck_assert_ptr_eq(NULL, A.matrix);
}
END_TEST

START_TEST(s21_eq_matrix_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_test_not) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(2, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);

  int reza2 = 0;
  matrix_t C = {0};
  matrix_t D = {0};
  s21_create_matrix(3, 3, &C);
  s21_create_matrix(3, 3, &D);
  C.matrix[0][0] = 2;
  C.matrix[0][1] = 4;
  C.matrix[0][2] = 6;
  C.matrix[1][0] = 8;
  C.matrix[1][1] = 10;
  C.matrix[1][2] = 12;
  C.matrix[2][0] = 14;
  C.matrix[2][1] = 16;
  C.matrix[2][2] = 18;
  s21_initialization(&D);
  reza2 = s21_eq_matrix(&C, &D);
  ck_assert_int_eq(0, reza2);
  s21_remove_matrix(&C);
  s21_remove_matrix(&D);
}
END_TEST

START_TEST(s21_sum_matrix_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &R2);
  R2.matrix[0][0] = 2;
  R2.matrix[0][1] = 4;
  R2.matrix[0][2] = 6;
  R2.matrix[1][0] = 8;
  R2.matrix[1][1] = 10;
  R2.matrix[1][2] = 12;
  R2.matrix[2][0] = 14;
  R2.matrix[2][1] = 16;
  R2.matrix[2][2] = 18;
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_sum_matrix(&A, &B, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_sum_matrix_test_incorrect) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  s21_create_matrix(-2, -2, &A);
  s21_create_matrix(3, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_sum_matrix(&A, &B, &R);
  ck_assert_int_eq(1, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
}
END_TEST

START_TEST(s21_sum_matrix_test_calc) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_sum_matrix(&A, &B, &R);
  ck_assert_int_eq(2, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
}
END_TEST

START_TEST(s21_sub_matrix_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &R2);
  R2.matrix[0][0] = 0;
  R2.matrix[0][1] = 0;
  R2.matrix[0][2] = 0;
  R2.matrix[1][0] = 0;
  R2.matrix[1][1] = 0;
  R2.matrix[1][2] = 0;
  R2.matrix[2][0] = 0;
  R2.matrix[2][1] = 0;
  R2.matrix[2][2] = 0;
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_sub_matrix(&A, &B, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_sub_matrix_test_incorrect) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  s21_create_matrix(-2, -2, &A);
  s21_create_matrix(3, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_sub_matrix(&A, &B, &R);
  ck_assert_int_eq(1, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
}
END_TEST

START_TEST(s21_sub_matrix_test_calc) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_sub_matrix(&A, &B, &R);
  ck_assert_int_eq(2, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
}
END_TEST

START_TEST(s21_mult_number_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &R2);
  s21_initialization(&A);
  s21_initialization(&R2);
  reza = s21_mult_number(&A, 1, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_mult_matrix_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &R2);
  R2.matrix[0][0] = 30;
  R2.matrix[0][1] = 36;
  R2.matrix[0][2] = 42;
  R2.matrix[1][0] = 66;
  R2.matrix[1][1] = 81;
  R2.matrix[1][2] = 96;
  R2.matrix[2][0] = 102;
  R2.matrix[2][1] = 126;
  R2.matrix[2][2] = 150;
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_mult_matrix(&A, &B, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_mult_matrix_test_calc) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t B = {0};
  matrix_t R = {0};
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(3, 3, &B);
  s21_initialization(&A);
  s21_initialization(&B);
  reza = s21_mult_matrix(&A, &B, &R);
  ck_assert_int_eq(2, reza);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&R);
}
END_TEST

START_TEST(s21_transpose_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &R2);
  R2.matrix[0][0] = 1;
  R2.matrix[0][1] = 2;
  R2.matrix[0][2] = 3;
  R2.matrix[1][0] = 4;
  R2.matrix[1][1] = 5;
  R2.matrix[1][2] = 6;
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 4;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 5;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 6;
  reza = s21_transpose(&A, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_calc_complements_test_1x1) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &R2);
  s21_initialization(&R2);
  s21_initialization(&A);
  reza = s21_calc_complements(&A, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_calc_complements_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &R2);
  R2.matrix[0][0] = 0;
  R2.matrix[0][1] = 10;
  R2.matrix[0][2] = -20;
  R2.matrix[1][0] = 4;
  R2.matrix[1][1] = -14;
  R2.matrix[1][2] = 8;
  R2.matrix[2][0] = -8;
  R2.matrix[2][1] = -2;
  R2.matrix[2][2] = 4;
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;
  reza = s21_calc_complements(&A, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

START_TEST(s21_determinant_test) {
  int reza = 0;
  double r = 0;
  matrix_t A = {0};
  s21_create_matrix(3, 3, &A);
  s21_initialization(&A);
  reza = s21_determinant(&A, &r);
  ck_assert_int_eq(0, reza);
  ck_assert_double_eq(0, r);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_test_1x1) {
  int reza = 0;
  double r = 0;
  matrix_t A = {0};
  s21_create_matrix(1, 1, &A);
  s21_initialization(&A);
  reza = s21_determinant(&A, &r);
  ck_assert_int_eq(0, reza);
  ck_assert_double_eq(1, r);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_test_calc) {
  int reza = 0;
  double r = 0;
  matrix_t A = {0};
  s21_create_matrix(2, 1, &A);
  s21_initialization(&A);
  reza = s21_determinant(&A, &r);
  ck_assert_int_eq(2, reza);
  ck_assert_double_eq(0, r);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_test) {
  int reza = 0;
  matrix_t A = {0};
  matrix_t R = {0};
  matrix_t R2 = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &R2);
  R2.matrix[0][0] = 1;
  R2.matrix[0][1] = -1;
  R2.matrix[0][2] = 1;
  R2.matrix[1][0] = -38;
  R2.matrix[1][1] = 41;
  R2.matrix[1][2] = -34;
  R2.matrix[2][0] = 27;
  R2.matrix[2][1] = -29;
  R2.matrix[2][2] = 24;
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;
  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;
  reza = s21_inverse_matrix(&A, &R);
  ck_assert_int_eq(0, reza);
  ck_assert_int_eq(1, s21_eq_matrix(&R2, &R));
  s21_remove_matrix(&A);
  s21_remove_matrix(&R);
  s21_remove_matrix(&R2);
}
END_TEST

int main() {
  Suite *s1 = suite_create("s21_matrix");
  TCase *tc1_1 = tcase_create("s21_matrix");
  SRunner *sr = srunner_create(s1);
  int nf = 0;
  suite_add_tcase(s1, tc1_1);

  tcase_add_test(tc1_1, s21_create_matrix_test);
  tcase_add_test(tc1_1, s21_remove_matrix_test);
  tcase_add_test(tc1_1, s21_eq_matrix_test);
  tcase_add_test(tc1_1, s21_eq_matrix_test_not);
  tcase_add_test(tc1_1, s21_sum_matrix_test);
  tcase_add_test(tc1_1, s21_sum_matrix_test_incorrect);
  tcase_add_test(tc1_1, s21_sum_matrix_test_calc);
  tcase_add_test(tc1_1, s21_sub_matrix_test);
  tcase_add_test(tc1_1, s21_sub_matrix_test_incorrect);
  tcase_add_test(tc1_1, s21_sub_matrix_test_calc);
  tcase_add_test(tc1_1, s21_mult_number_test);
  tcase_add_test(tc1_1, s21_mult_matrix_test);
  tcase_add_test(tc1_1, s21_mult_matrix_test_calc);
  tcase_add_test(tc1_1, s21_transpose_test);
  tcase_add_test(tc1_1, s21_calc_complements_test_1x1);
  tcase_add_test(tc1_1, s21_calc_complements_test);
  tcase_add_test(tc1_1, s21_determinant_test);
  tcase_add_test(tc1_1, s21_determinant_test_1x1);
  tcase_add_test(tc1_1, s21_determinant_test_calc);
  tcase_add_test(tc1_1, s21_inverse_matrix_test);

  srunner_run_all(sr, CK_ENV);
  nf = srunner_ntests_failed(sr);
  srunner_free(sr);

  return nf == 0 ? 0 : 1;
}
