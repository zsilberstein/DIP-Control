/**
 * @file LQR.h
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Class to represent a Linear Quadratic Regulator (LQR)
 * @date 2024-06-05
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "LQR.h"
#include "Eigen/Dense"

/**
 * @brief Solves the continuous algebraic riccati equation (CARE) with a
 * matrix sign algorithm.
 * @cite J. D. GARDINER AND A. J. LAUB, A generalization of the
 * matrix-sign-function solution for algebraic Riccati equations, Internat.
 * J. Control, 44 (1986), pp. 823–832.
 *
 * @return Eigen::MatrixXd Returns the solution X (n by n) to the matrix
 * equation A.T*X + X*A − X*B*R^(−1)*B.T*X + Q = 0
 */
Eigen::MatrixXd LQR::_solveCARE() {
    Eigen::MatrixXd H_curr(2 * _n, 2 * _n);
    H_curr << _A, -_B * _R.inverse() * _B.transpose(), -_Q, -_A.transpose();

    Eigen::MatrixXd H_prev =
        Eigen::MatrixXd::Identity(2 * _n, 2 * _n) * INFINITY;

    while (abs((H_curr - H_prev).maxCoeff() / H_curr.maxCoeff()) > 0.0001) {
        H_prev = H_curr;
        H_curr = 1.0 / 2.0 * (H_prev + H_prev.inverse());
    }

    Eigen::MatrixXd W11 =
        H_curr.block(0, 0, _n, _n) + Eigen::MatrixXd::Identity(_n, _n);
    Eigen::MatrixXd W12 = H_curr.block(0, _n, _n, _n);

    return -W12.inverse() * W11;
}

/**
 * @brief Construct a new LQR object
 *
 * @param A State matrix (n by n)
 * @param B Input matrix (n by m)
 * @param Q State cost matrix (n by n)
 * @param R Input cost matrix (m by m)
 */
LQR::LQR(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
         const Eigen::MatrixXd &Q, const Eigen::MatrixXd &R)
    : _A(A), _B(B), _Q(Q), _R(R), _n(_B.rows()), _m(_B.cols()) {
    _K = R.inverse() * B.transpose() * _solveCARE();
}

/**
 * @brief Get the LQR state-feedback control u (m by 1).
 *
 * @param state Current state of the system (n by 1).
 * @return Eigen::MatrixXd
 */
Eigen::MatrixXd LQR::getFeedbackControl(Eigen::VectorXd state) {
    return -_K * state;
}
