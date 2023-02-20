#ifndef MEM_EKF_COMPONENTS_HELPERS_HELPER_FUNCTIONS_HPP_
#define MEM_EKF_COMPONENTS_HELPERS_HELPER_FUNCTIONS_HPP_

#include <Eigen/Dense>

#include "ellipse.hpp"

namespace eot {
  /* B = 0.5 * (A + A') */
  template <size_t matrix_size>
  Eigen::Matrix<double, matrix_size, matrix_size> MakeMatrixSymetric(Eigen::Matrix<double, matrix_size, matrix_size> & matrix) {
    Eigen::Matrix<double, matrix_size, matrix_size> output;
    
    output = 0.5 * (matrix + matrix.transpose());

    return output;
  }

  /* Kronecker tensor product */
  template <size_t a_rows, size_t a_cols, size_t b_rows, size_t b_cols>
  const Eigen::Matrix<double, a_rows * b_rows, a_cols * b_cols> & KroneckerProduct(const Eigen::Matrix<double, a_rows, a_cols> & a, const Eigen::Matrix<double, b_rows, b_cols> & b) {
    static Eigen::Matrix<double, a_rows * b_rows, a_cols * b_cols> output;

    for (size_t row = 0u; row < a_rows; row++) {
      for (size_t col = 0u; col < a_cols; col++) {
        const auto output_row_index = b_rows * row;
        const auto output_col_index = b_cols * col;
        const auto new_value = a(row, col) * b;
        for (auto updated_row = 0; updated_row < b_rows; updated_row++) {
          for (auto updated_col = 0; updated_col < b_cols; updated_col++) {
            output(updated_row + output_row_index, updated_col + output_col_index) = new_value(updated_row, updated_col);
          }
        }
      }
    }

    return output;
  }

  /* Convert ellipse to vector */
  const Eigen::Vector<double, 3u> & ConvertEllipseToVector(const Ellipse & ellipse) {
    static Eigen::Vector<double, 3u> vector;

    vector(0u) = ellipse.alpha;
    vector(1u) = ellipse.l1;
    vector(2u) = ellipse.l2;

    return vector;
  }

  /* Convert vector to ellipse */
  const Ellipse & ConvertVectorToEllipse(const Eigen::Vector<double, 3u> & vector) {
    static Ellipse ellipse;

    ellipse.alpha = vector(0u);
    ellipse.l1 = vector(1u);
    ellipse.l2 = vector(2u);

    return ellipse;
  }

  /* Diagonal to matrix */
  template <size_t size>
  Eigen::Matrix<double, size, size> ConvertDiagonalToMatrix(const std::array<double, size> & diagonal) {
    static Eigen::Matrix<double, size, size> output = Eigen::Matrix<double, size, size>::Zero();

    for (size_t index = 0u; index < size; index++)
      output(index, index) = diagonal.at(index);

    return output;
  }
} //  namespace eot

#endif  //  MEM_EKF_COMPONENTS_HELPERS_HELPER_FUNCTIONS_HPP_
