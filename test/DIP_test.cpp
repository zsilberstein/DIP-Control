/**
 * @file DIP_test.cpp
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Tests for the DIP class
 * @date 2024-06-10
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <gtest/gtest.h>
#include "../src/DIP.h"
namespace {
  /**
   * @brief Test fixture for the DIP class
   * 
   */
class DIPTest : public testing::Test {
   protected:
    // Initialize dipRandom to have some arbitrary values
    DIPTest()
        : dipRandom(3.4, 2.8, 0.7, 1.3, 2.6, 1.4, 0.5, 0.3, 0.02,
                    Eigen::Vector<double, 6>{3.3, M_PI / 8, M_PI / 12, 2.2, 0.6, 0.5}) {}
    // Default constructor will initialize dipDefault
    DIP dipDefault;
    DIP dipRandom;
};

// Test that basic getters work as expected and that the DIP objects were
// initialized with expected values
TEST_F(DIPTest, GetCartMass) {
    EXPECT_EQ(dipDefault.getCartMass(), 1.0);
    EXPECT_EQ(dipRandom.getCartMass(), 3.4);
}

TEST_F(DIPTest, GetMassOne) {
    EXPECT_EQ(dipDefault.getMassOne(), 1.0);
    EXPECT_EQ(dipRandom.getMassOne(), 2.8);
}

TEST_F(DIPTest, GetMassTwo) {
    EXPECT_EQ(dipDefault.getMassTwo(), 1.0);
    EXPECT_EQ(dipRandom.getMassTwo(), 0.7);
}

TEST_F(DIPTest, GetLinkOneLength) {
    EXPECT_EQ(dipDefault.getLinkOneLength(), 1.0);
    EXPECT_EQ(dipRandom.getLinkOneLength(), 1.3);
}

TEST_F(DIPTest, GetLinkTwoLength) {
    EXPECT_EQ(dipDefault.getLinkTwoLength(), 1.0);
    EXPECT_EQ(dipRandom.getLinkTwoLength(), 2.6);
}

TEST_F(DIPTest, GetState) {
    EXPECT_EQ(dipDefault.getState(), Eigen::VectorXd::Zero(6));
    Eigen::Vector<double, 6> randomState{3.3, M_PI / 8, M_PI / 12, 2.2, 0.6, 0.5};
    EXPECT_EQ(dipRandom.getState(), randomState);
}

TEST_F(DIPTest, GetMassOnePos) {
    // Should be straight up since angles are zero
    Eigen::Vector2d defaultMassOnePos{0, 1};
    EXPECT_EQ(dipDefault.getMassOnePos(), defaultMassOnePos);

    // X position is shifted over by cart location (3.3m)
    Eigen::Vector2d randomMassOnePos{3.3 + 1.3 * sin(M_PI / 8),
                                     1.3 * cos(M_PI / 8)};
    EXPECT_EQ(dipRandom.getMassOnePos(), randomMassOnePos);
}

TEST_F(DIPTest, GetMassTwoPos) {
    // Should be straight up since angles are zero
    Eigen::Vector2d defaultMassTwoPos{0, 2};
    EXPECT_EQ(dipDefault.getMassTwoPos(), defaultMassTwoPos);

    // X position is shifted over by cart location (3.3m)
    Eigen::Vector2d randomMassTwoPos{3.3 + 1.3 * sin(M_PI / 8) + 2.6 * sin(M_PI/8 + M_PI/12),
                                     1.3 * cos(M_PI / 8) + 2.6 * cos(M_PI/8 + M_PI/12)};
    EXPECT_EQ(dipRandom.getMassTwoPos(), randomMassTwoPos);
}

TEST_F(DIPTest, GetPotentialEnergy) {
    double m1Pot = 1 * 9.81 * 1; // Mass of 1kg and 1m up
    double m2Pot = 1 * 9.81 * 2; // Mass of 1kg and 2m up
    EXPECT_EQ(dipDefault.getPotentialEnergy(), m1Pot + m2Pot);

    m1Pot = 2.8 * 9.81 * 1.3 * cos(M_PI / 8); // M * g * h
    m2Pot = 0.7 * 9.81 * (1.3 * cos(M_PI / 8) + 2.6 * cos(M_PI/8 + M_PI/12)); // M * g * h
    EXPECT_EQ(dipRandom.getPotentialEnergy(), m1Pot + m2Pot);
}

TEST_F(DIPTest, GetKineticEnergy) {
    EXPECT_EQ(dipDefault.getKineticEnergy(), 0); // Zero velocity for all members

    double cartKin = 0.5 * 3.4 * pow(2.2, 2);  // 1/2*m*v^2
    // See EOMDerivation.py
    double m1Kin = 0.5 * 2.8 *
                   (pow(1.3 * 0.6 * cos(M_PI / 8) + 2.2, 2) +
                    pow(1.3 * 0.6 * sin(M_PI / 8), 2));
    double m2Kin = 0.5 * 0.7 *
                   (pow(1.3 * 0.6 * cos(M_PI / 8) + 2.2 +
                            2.6 * (0.5 + 0.6) * cos(M_PI / 8 + M_PI / 12),
                        2) +
                    pow(1.3 * 0.6 * sin(M_PI / 8) +
                            2.6 * (0.5 + 0.6) * sin(M_PI / 8 + M_PI / 12),
                        2));
    EXPECT_EQ(dipRandom.getKineticEnergy(), cartKin + m1Kin + m2Kin);
}

// See EOMDerivation.py
const Eigen::Matrix<double, 6, 6> defaultAc{{0, 0, 0, 1, 0, 0},
                                            {0, 0, 0, 0, 1, 0},
                                            {0, 0, 0, 0, 0, 1},
                                            {0, -2.0 * 9.81, 0, 0, 0, 0},
                                            {0, 3 * 9.81, -9.81, 0, 0, 0},
                                            {0, -3 * 9.81, 3 * 9.81, 0, 0, 0}};

const Eigen::Matrix<double, 6, 6> randomAc{
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 1},
    {0, -9.81 * (2.8 + 0.7) / 3.4, 0, -1.4 / 3.4, 0.5 / (1.3 * 3.4),
     -0.3 / (1.3 * 3.4)},
    {0, 9.81 * (3.4 + 2.8 + 0.7) / (3.4 * 1.3), -9.81 * 0.7 / (2.8 * 1.3),
     1.4 / (3.4 * 1.3), -0.5 * (3.4 + 2.8) / (1.3 * 1.3 * 3.4 * 2.8),
     0.3 * (1.3 * 3.4 + 2.6 * 3.4 + 2.6 * 2.8) / (1.3 * 1.3 * 2.6 * 3.4 * 2.8)},
    {0, -9.81 * (3.4 + 2.8 + 0.7) / (3.4 * 1.3),
     9.81 * (1.3 * 2.8 + 1.3 * 0.7 + 2.6 * 0.7) / (1.3 * 2.6 * 2.8),
     -1.4 / (1.3 * 3.4),
     0.5 * (1.3 * 3.4 + 2.6 * 3.4 + 2.6 * 2.8) / (1.3 * 1.3 * 2.6 * 3.4 * 2.8),
     -0.3 * (1 / (2.6 * 2.6 * 0.7) + 1 / (2.6 * 2.6 * 2.8) +
             2 / (1.3 * 2.6 * 2.8) + 1 / (1.3 * 1.3 * 3.4) +
             1 / (1.3 * 1.3 * 2.8))}};

TEST_F(DIPTest, GetAc) {
    EXPECT_EQ(dipDefault.getAc(), defaultAc);
    // Use isApprox() due to double rounding
    EXPECT_TRUE(randomAc.isApprox(dipRandom.getAc()));
}

const Eigen::Vector<double, 6> defaultBc{0, 0, 0, 1, -1, 1};

const Eigen::Vector<double, 6> randomBc{
    0, 0, 0, 1 / 3.4, -1 / (3.4 * 1.3), 1 / (3.4 * 1.3)};

TEST_F(DIPTest, GetBc) {
    EXPECT_EQ(defaultBc, dipDefault.getBc());
    EXPECT_EQ(randomBc, dipRandom.getBc());
}

TEST_F(DIPTest, GetAd) {
    Eigen::Matrix<double, 6, 6> defaultAd =
        Eigen::MatrixXd::Identity(6, 6) + 0.01 * defaultAc +
        pow(0.01, 2) / 2 * defaultAc * defaultAc +
        pow(0.01, 3) / 6 * defaultAc * defaultAc * defaultAc +
        pow(0.01, 4) / 24 * defaultAc * defaultAc * defaultAc * defaultAc;
    Eigen::Matrix<double, 6, 6> randomAd =
        Eigen::MatrixXd::Identity(6, 6) + 0.02 * randomAc +
        pow(0.02, 2) / 2 * randomAc * randomAc +
        pow(0.02, 3) / 6 * randomAc * randomAc * randomAc +
        pow(0.02, 4) / 24 * randomAc * randomAc * randomAc * randomAc;

    EXPECT_EQ(defaultAd, dipDefault.getAd());
    // Use isApprox() due to double rounding
    EXPECT_TRUE(randomAd.isApprox(dipRandom.getAd()));
}

TEST_F(DIPTest, GetBd) {
    Eigen::Vector<double, 6> defaultBd =
        (Eigen::MatrixXd::Identity(6, 6) * 0.01 + pow(0.01, 2) / 2 * defaultAc +
         pow(0.01, 3) / 6 * defaultAc * defaultAc +
         pow(0.01, 4) / 24 * defaultAc * defaultAc * defaultAc) *
        defaultBc;
    Eigen::Vector<double, 6> randomBd =
        (Eigen::MatrixXd::Identity(6, 6) * 0.02 + pow(0.02, 2) / 2 * randomAc +
         pow(0.02, 3) / 6 * randomAc * randomAc +
         pow(0.02, 4) / 24 * randomAc * randomAc * randomAc) *
        randomBc;

    EXPECT_EQ(defaultBd, dipDefault.getBd());
    EXPECT_EQ(randomBd, dipRandom.getBd());
}

TEST_F(DIPTest, UpdateState) {
    // Confirm updateState() works with zero input force
    dipDefault.updateState(0);
    dipRandom.updateState(0);

    // DIP will remain in equilibrium if no external force acts on it
    Eigen::Vector<double, 6> defaultIdealNextState{0, 0, 0, 0, 0, 0};
    EXPECT_EQ(defaultIdealNextState, dipDefault.getState());

    // Confirmed by performing Rk4
    Eigen::Vector<double, 6> randomIdealNextState{3.34331, 0.40573, 0.27093,
                                                  2.13151, 0.70298, 0.41297};
    EXPECT_TRUE(randomIdealNextState.isApprox(dipRandom.getState(), 0.00001));

    // Confirm updateState() works with non-zero input force
    dipDefault.updateState(76.34);
    dipRandom.updateState(-43.6);

    defaultIdealNextState = {0.00382, -0.00382, 0.00382,
                             0.76363, -0.76388, 0.76411};
    EXPECT_TRUE(defaultIdealNextState.isApprox(dipDefault.getState(), 0.00001));

    randomIdealNextState = {3.38308, 0.42241, 0.27661,
                            1.84617, 0.96444, 0.15610};
    EXPECT_TRUE(randomIdealNextState.isApprox(dipRandom.getState(), 0.00001));
}
}