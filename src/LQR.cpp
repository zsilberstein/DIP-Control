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

    // Compute matrix sign
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
 * @brief Solves the discrete algebraic riccati equation (DARE) with a
 * structure-preserving doubling algorithm.
 * @cite E. K.-W. Chu, H.-Y. Fan, W.-W. Lin & C.-S. Wang (2004)
 * Structure-Preserving Algorithms for Periodic Discrete-Time Algebraic
 * Riccati Equations, International Journal of Control, 77:8, 767-788,
 * DOI: 10.1080/00207170410001714988
 *
 * @return Eigen::MatrixXd Returns the solution X (n by n) to the matrix
 * equation A.T*X*A − X − (A.T*X*B)(B.T*X*B+R)^(-1)*(B.T*X*A) + Q = 0
 */
Eigen::MatrixXd LQR::_solveDARE() {
    Eigen::MatrixXd A_curr = _A;
    Eigen::MatrixXd G_curr = _B * _R.inverse() * _B.transpose();
    Eigen::MatrixXd H_curr = _Q;

    Eigen::MatrixXd A_prev(_n, _n);
    Eigen::MatrixXd G_prev(_n, _n);
    Eigen::MatrixXd inv(_n, _n);
    Eigen::MatrixXd H_prev = Eigen::MatrixXd::Identity(_n, _n) * INFINITY;

    while (abs((H_curr - H_prev).maxCoeff() / H_curr.maxCoeff()) > 0.0001) {
        A_prev = A_curr;
        G_prev = G_curr;
        H_prev = H_curr;

        inv = (Eigen::MatrixXd::Identity(_n, _n) + G_prev * H_prev).inverse();
        A_curr = A_prev * inv * A_prev;
        G_curr = G_prev + A_prev * inv * G_prev * A_prev.transpose();
        H_curr = H_prev + A_prev.transpose() * H_prev * inv * A_prev;
    }
    return H_curr;
}

/**
 * @brief Construct a new LQR object
 *
 * @param A State matrix (n by n)
 * @param B Input matrix (n by m)
 * @param Q State cost matrix (n by n)
 * @param R Input cost matrix (m by m)
 * @param continuous true for a continuous system, false for a discrete
 * system
 */
LQR::LQR(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd Q,
         Eigen::MatrixXd R, bool continuous)
    : _A(A), _B(B), _Q(Q), _R(R), _n(_B.rows()), _m(_B.cols()) {
    if (continuous) {
        _K = R.inverse() * B.transpose() * _solveCARE();
    }
    else {
        Eigen::MatrixXd S = _solveDARE();
        _K = (_B.transpose() * S * _B + R).inverse() * (_B.transpose() * S * _A);
    }
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
