/**
 * @file LQR.h
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Class to represent a Linear Quadratic Regulator (LQR)
 * @date 2024-06-05
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "Eigen/Dense"

#ifndef LQR_H
#define LQR_H

/**
 * @brief Class to represent a Linear Quadratic Regulator (LQR)
 *
 */
class LQR {
   private:
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
    Eigen::MatrixXd _solveCARE();

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
    Eigen::MatrixXd _solveDARE();

    // State matrix (n by n)
    const Eigen::MatrixXd &_A;
    // Input matrix (n by m)
    const Eigen::MatrixXd &_B;
    // State cost matrix (n by n)
    const Eigen::MatrixXd &_Q;
    // Input cost matrix (m by m)
    const Eigen::MatrixXd &_R;
    // Optimal gain (m by n)
    Eigen::MatrixXd _K;
    // Number of states
    const int _n;
    // Number of inputs
    const int _m;

   public:
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
    LQR(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd Q,
        Eigen::MatrixXd R, bool continuous);

    /**
     * @brief Get the LQR state-feedback control u (m by 1).
     *
     * @param state Current state of the system (n by 1).
     * @return Eigen::MatrixXd
     */
    Eigen::MatrixXd getFeedbackControl(Eigen::VectorXd state);
};

#endif