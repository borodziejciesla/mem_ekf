#ifndef MEM_EKF_COMPONENTS_HELPERS_HELPER_FUNCTIONS_HPP_
#define MEM_EKF_COMPONENTS_HELPERS_HELPER_FUNCTIONS_HPP_

#include <Eigen/Dense>

#include "ellipse.hpp"

namespace eot {
  /* B = 0.5 * (A + A') */
  template <size_t matrix_size>
  void MakeMatrixSymetric(Eigen::Matrix<double, matrix_size, matrix_size> & matrix) {
    matrix +=  matrix.transpose();
    matrix *= 0.5;
  }

  /* Kronecker tensor product */
  template <size_t a_rows, size_t a_cols, size_t b_rows, size_t b_cols, typename matrix_a, typename matrix_b>
  const Eigen::Matrix<double, a_rows * b_rows, a_cols * b_cols> & KroneckerProduct(const matrix_a & a, const matrix_b & b) {
    static Eigen::Matrix<double, a_rows * b_rows, a_cols * b_cols> output;

    for (size_t row = 0u; row < a_rows; row++) {
      for (size_t col = 0u; col < a_cols; col++) {
        const auto output_row_index = b_rows * (row - 1u);
        const auto output_col_index = b_cols * (col - 1u);
        output.block<b_rows, b_cols>(output_row_index, output_col_index) = a(row, col) * b;
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

    ellipse.alpha =  vector(0u);
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
