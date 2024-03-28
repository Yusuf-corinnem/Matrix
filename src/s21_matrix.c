#include "s21_matrix.h"

int s21_check_matrix(matrix_t *check) {
  int res = 0;
  if (check == NULL || check->matrix == NULL || check->columns < 1 ||
      check->rows < 1) {
    res = 1;
  }
  return res;
}

void s21_initialization(matrix_t *A) {
  int count = 1;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      A->matrix[i][j] = count++;
    }
  }
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (result == NULL || rows < 1 || columns < 1) return 1;
  int res = 0;
  result->rows = rows;
  result->columns = columns;
  result->matrix =
      malloc(rows * columns * sizeof(double) + rows * sizeof(double));
  if (result->matrix != NULL) {
    double *pointer = (double *)(result->matrix + rows);
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = pointer + columns * i;
    }
  } else
    res = 1;
  return res;
}

void s21_remove_matrix(matrix_t *A) {
  free(A->matrix);
  A->rows = 0, A->columns = 0, A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = 1;
  if (A->rows == B->rows && A->columns == B->columns &&
      s21_check_matrix(A) == 0 && s21_check_matrix(B) == 0) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
          res = 0;
        }
      }
    }
  } else
    res = 0;
  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) || s21_check_matrix(B)) return 1;
  if (A->rows == B->rows && A->columns == B->columns) {
    if (s21_create_matrix(A->rows, A->columns, result) == 0) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else
      res = 1;
  } else
    res = 2;
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) || s21_check_matrix(B)) return 1;
  if (A->rows == B->rows && A->columns == B->columns) {
    if (s21_create_matrix(A->rows, A->columns, result) == 0) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else
      res = 1;
  } else
    res = 2;
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (s21_check_matrix(A) != 0) return 1;
  int res = 0;
  if (s21_create_matrix(A->rows, A->columns, result) == 0) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  } else
    res = 1;
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0, m = A->rows, q = B->columns, n = A->columns;
  if (s21_check_matrix(A) || s21_check_matrix(B))
    res = 1;
  else if (s21_create_matrix(A->rows, B->columns, result) != 0)
    res = 1;
  else if (A->columns == B->rows) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < q; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < n; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  } else
    res = 2;

  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A))
    res = 1;
  else if (s21_create_matrix(A->columns, A->rows, result) != 0)
    res = 1;
  else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (s21_check_matrix(A) || result == NULL) return 1;
  int res = 0;
  if (A->columns != A->rows) {
    res = 1;
  } else if (A->rows == 1) {
    if (s21_create_matrix(A->rows, A->columns, result) == 0) {
      result->matrix[0][0] = A->matrix[0][0];
    }
  } else if (s21_create_matrix(A->rows, A->columns, result) == 0) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        matrix_t temp = {0};
        double det = 0;
        s21_minor(A, &temp, i, j);
        if (s21_determinant(&temp, &det) == 0) {
          result->matrix[i][j] = pow(-1, (i + j)) * det;
          s21_remove_matrix(&temp);
        } else
          res = 2;
      }
    }
  }
  return res;
}

int s21_minor(matrix_t *A, matrix_t *matrix_minor, int index_rows,
              int index_columns) {
  int res = 0;
  if (s21_create_matrix(A->rows - 1, A->columns - 1, matrix_minor) == 0) {
    int minor_row = 0;
    for (int i = 0; i < A->rows; i++) {
      if (i == index_rows) continue;
      int minor_columns = 0;
      for (int j = 0; j < A->columns; j++) {
        if (j == index_columns) continue;
        matrix_minor->matrix[minor_row][minor_columns] = A->matrix[i][j];
        minor_columns++;
      }
      minor_row++;
    }
  } else
    res = 1;
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  if (s21_check_matrix(A) || result == NULL) return 1;
  int res = 0;
  *result = 0;
  if (A->columns != A->rows) {
    res = 2;
  } else if (A->rows == 1) {
    *result = A->matrix[0][0];
  } else if (A->rows == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    int sign = 1;
    for (int i = 0; i < A->columns; i++) {
      matrix_t temp = {0};
      int check_minor = s21_minor(A, &temp, 0, i);
      double temp_d = 0;
      if (s21_determinant(&temp, &temp_d) == 0 && check_minor == 0) {
        *result += A->matrix[0][i] * temp_d * sign;
        sign = -sign;
        s21_remove_matrix(&temp);
      } else
        res = 2;
    }
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) != 0)
    res = 1;
  else if (A->rows != A->columns)
    res = 2;
  matrix_t calc = {0}, transpons = {0};
  double determinant = 0;
  s21_determinant(A, &determinant);
  if (fabs(determinant) < 1e-7)
    res = 2;
  else
    s21_calc_complements(A, &calc);
  s21_transpose(&calc, &transpons);
  s21_mult_number(&transpons, 1 / determinant, result);
  s21_remove_matrix(&calc);
  s21_remove_matrix(&transpons);
  return res;
}