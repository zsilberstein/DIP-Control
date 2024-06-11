/**
 * @file LQR_test.cpp
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Tests for the LQR class
 * @date 2024-06-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <gtest/gtest.h>
#include "../src/LQR.h"
#include "../src/DIP.h"

namespace {

TEST(LQRTest, ContinuousLQR) {
    // Test with 1 state and 1 input
    Eigen::Matrix<double, 1, 1> A1{0.87};
    Eigen::Matrix<double, 1, 1> B1{-1.43};
    Eigen::Matrix<double, 1, 1> Q1{1.0};
    Eigen::Matrix<double, 1, 1> R1{1.0};
    Eigen::Vector<double, 1> state1{0.48};
    LQR lqr1(A1, B1, Q1, R1, true);

    // Ideal values are computed with control toolbox
    Eigen::Matrix<double, 1, 1> ideal1{0.85388};

    EXPECT_TRUE(ideal1.isApprox(lqr1.getFeedbackControl(state1), 0.00001));

    // Test with 2 states and 2 inputs
    Eigen::Matrix2d A2{{0, 1}, {1, 0}};
    Eigen::Matrix2d B2{{1, 0}, {0, 1}};
    Eigen::Matrix2d Q2{{1, 0}, {0, 1}};
    Eigen::Matrix2d R2{{1, 0}, {0, 1}};
    Eigen::Vector2d state2{1, -1};
    LQR lqr2(A2, B2, Q2, R2, true);

    // Ideal values are computed with control toolbox
    Eigen::Matrix<double, 2, 1> ideal2{-0.41421, 0.41421};

    EXPECT_TRUE(ideal2.isApprox(lqr2.getFeedbackControl(state2), 0.00001));

    // Test with 6 states and 1 input with the DIP system
    Eigen::Vector<double, 6> state3{3.3, M_PI / 8, M_PI / 12, 2.2, 0.6, 0.5};
    DIP dipRandom{3.4, 2.8, 0.7, 1.3, 2.6, 1.4, 0.5, 0.3, 0.02, state3};
    Eigen::Matrix<double, 6, 6> Q3 = Eigen::MatrixXd::Identity(6, 6);
    Eigen::Matrix<double, 1, 1> R3{1.0};

    LQR lqr3(dipRandom.getAc(), dipRandom.getBc(), Q3, R3, true);

    // Ideal values are computed with control toolbox
    Eigen::Matrix<double, 1, 1> ideal3{-371.38286};

    EXPECT_TRUE(ideal3.isApprox(lqr3.getFeedbackControl(state3), 0.00001));
}

TEST(LQRTest, DiscreteLQR) {
    // Test with 1 state and 1 input
    Eigen::Matrix<double, 1, 1> A1{0.87};
    Eigen::Matrix<double, 1, 1> B1{-1.43};
    Eigen::Matrix<double, 1, 1> Q1{1.0};
    Eigen::Matrix<double, 1, 1> R1{1.0};
    Eigen::Vector<double, 1> state1{0.48};
    LQR lqr1(A1, B1, Q1, R1, false);

    // Ideal values are computed with control toolbox
    Eigen::Matrix<double, 1, 1> ideal1{0.210706};

    EXPECT_TRUE(ideal1.isApprox(lqr1.getFeedbackControl(state1), 0.00001));

    // Test with 2 states and 2 inputs
    Eigen::Matrix2d A2{{0, 1}, {1, 0}};
    Eigen::Matrix2d B2{{1, 0}, {0, 1}};
    Eigen::Matrix2d Q2{{1, 0}, {0, 1}};
    Eigen::Matrix2d R2{{1, 0}, {0, 1}};
    Eigen::Vector2d state2{1, -1};
    LQR lqr2(A2, B2, Q2, R2, false);

    // Ideal values are computed with control toolbox
    Eigen::Matrix<double, 2, 1> ideal2{0.61803, -0.61803};

    EXPECT_TRUE(ideal2.isApprox(lqr2.getFeedbackControl(state2), 0.00001));

    // Test with 6 states and 1 input with the DIP system
    Eigen::Vector<double, 6> state3{3.3, M_PI / 8, M_PI / 12, 2.2, 0.6, 0.5};
    DIP dipRandom{3.4, 2.8, 0.7, 1.3, 2.6, 1.4, 0.5, 0.3, 0.02, state3};
    Eigen::Matrix<double, 6, 6> Q3 = Eigen::MatrixXd::Identity(6, 6);
    Eigen::Matrix<double, 1, 1> R3{1.0};

    LQR lqr3(dipRandom.getAc(), dipRandom.getBc(), Q3, R3, false);

    // Ideal values are computed with control toolbox
    Eigen::Matrix<double, 1, 1> ideal3{21.92422};

    EXPECT_TRUE(ideal3.isApprox(lqr3.getFeedbackControl(state3), 0.00001));
}
}