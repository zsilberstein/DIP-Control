/**
 * @file DIP.h
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Class to represent a double inverted pendulum on a cart
 * @date 2024-06-03
 *
 * @copyright Copyright (c) 2024
 *
 */
#include "Eigen/Dense"

#ifndef DIP_H
#define DIP_H

/**
 * @brief Class to represent a double inverted pendulum on a cart
 *
 */
class DIP {
   private:
    /**
     * @brief Calculates the time derivative of the state.
     * See the EOMDerivation.py file for a derivation of the equations of
     * motion.
     *
     * @param currentState The current state of the DIP to calculate the
     * derivative at
     * @param u The input force on the cart
     * @return Eigen::Vector<double, 6> The time derivative of the state
     */
    Eigen::Vector<double, 6> _getDerivative(
        Eigen::Vector<double, 6> &currentState, double u);

    /**
     * @brief Performs 4th order Runge–Kutta integration with zero-order hold on
     * the input force
     *
     * @param prevState The current state of the DIP to calculate the next state
     * from
     * @param u The input force on the cart
     * @return Eigen::Vector<double, 6> Estimated next state at a time dt in the
     * future
     */
    Eigen::Vector<double, 6> _Rk4(Eigen::Vector<double, 6> &prevState,
                                  double u);

    const double _mc;        // Cart mass (kg)
    const double _m1;        // Mass one (kg)
    const double _m2;        // Mass two (kg)
    const double _L1;        // Link one length (m)
    const double _L2;        // Link two length (m)
    const double _Dc;        // Cart damping (N*s/m)
    const double _D1;        // Theta one damping (N*m*s/rad)
    const double _D2;        // Theta two damping (N*m*s/rad)
    const double _dt;        // Time between states (s)
    const double _g = 9.81;  // Gravity (m/s^2)
    // Current state of the DIP system as [cart position (m), theta one (rads),
    // theta two (rads), cart velocity (m/s), Omega one (rads/s), Omega two
    // (rads/s)]
    Eigen::Vector<double, 6> _state;
    Eigen::Matrix<double, 6, 6> _Ac;  // Continuous state-space matrix A
    Eigen::Vector<double, 6> _Bc;     // Continuous state-space matrix B
    Eigen::Matrix<double, 6, 6> _Ad;  // Discrete state-space matrix A
    Eigen::Vector<double, 6> _Bd;     // Discrete state-space matrix B

   public:
    /**
     * @brief Construct a new DIP object with default system parameters
     *
     */
    DIP();

    /**
     * @brief Construct a new DIP object
     *
     * @param cartMass Cart mass (kg)
     * @param massOne Mass one (kg)
     * @param massTwo Mass two (kg)
     * @param linkOneLength Link one length (m)
     * @param linkTwoLength Link two length (m)
     * @param cartDamping Cart damping (N*s/m)
     * @param thetaOneDamping Theta one damping (N*m*s/rad)
     * @param thetaTwoDamping Theta two damping (N*m*s/rad)
     * @param dt Time between states (s)
     * @param initalState Initial state of the DIP system as [cart position (m),
     * theta one (rads), theta two (rads), cart velocity (m/s), Omega one
     * (rads/s), Omega two (rads/s)]
     */
    DIP(double cartMass, double massOne, double massTwo, double linkOneLength,
        double linkTwoLength, double cartDamping, double thetaOneDamping,
        double thetaTwoDamping, double dt,
        Eigen::Vector<double, 6> initialState);

    double getCartMass();
    double getMassOne();
    double getMassTwo();
    double getLinkOneLength();
    double getLinkTwoLength();

    /**
     * @brief Gets the (x,y) position of mass one in meters
     *
     * @return Eigen::Vector2d
     */
    Eigen::Vector2d getMassOnePos();

    /**
     * @brief Gets the (x,y) position of mass two in meters
     *
     * @return Eigen::Vector2d
     */
    Eigen::Vector2d getMassTwoPos();

    /**
     * @brief Gets current state of the DIP system as [cart position (m), theta
     * one (rads), theta two (rads), cart velocity (m/s), Omega one (rads/s),
     * Omega two (rads/s)]
     *
     * @return Eigen::Vector<double, 6>
     */
    Eigen::Vector<double, 6> getState();

    /**
     * @brief Updates the state of the DIP system to a time dt away by
     * integrating with RK4
     *
     * @param u The input force on the cart
     * @return Eigen::Vector<double, 6> The new state of the DIP system
     */
    Eigen::Vector<double, 6> updateState(double u);

    /**
     * @brief Get the potential energy of the DIP at the current state
     *
     * @return double
     */
    double getPotentialEnergy();

    /**
     * @brief Get the kinetic energy of the DIP at the current state
     *
     * @return double
     */
    double getKineticEnergy();

    /**
     * @brief Get the continuous state-space matrix A
     *
     * @return Eigen::Matrix<double, 6, 6>
     */
    Eigen::Matrix<double, 6, 6> getAc();

    /**
     * @brief Get the continuous state-space matrix B
     *
     * @return Eigen::Vector<double, 6>
     */
    Eigen::Vector<double, 6> getBc();

    /**
     * @brief Get the discrete state-space matrix A
     *
     * @return Eigen::Matrix<double, 6, 6>
     */
    Eigen::Matrix<double, 6, 6> getAd();

    /**
     * @brief Get the discrete state-space matrix B
     *
     * @return Eigen::Vector<double, 6>
     */
    Eigen::Vector<double, 6> getBd();
};

#endif