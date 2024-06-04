/**
 * @file DIP.h
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Class to represent a double inverted pendulum on a cart
 * @date 2024-06-03
 *
 * @copyright Copyright (c) 2024
 *
 */
#include "DIP.h"

#include "Eigen/Dense"

DIP::DIP(double cartMass, double massOne, double massTwo, double linkOneLength,
         double linkTwoLength, double cartDamping, double thetaOneDamping,
         double thetaTwoDamping, double dt,
         Eigen::Vector<double, 6> initialState)
    : _mc(cartMass),
      _m1(massOne),
      _m2(massTwo),
      _L1(linkOneLength),
      _L2(linkTwoLength),
      _Dc(cartDamping),
      _D1(thetaOneDamping),
      _D2(thetaTwoDamping),
      _dt(dt),
      _state(initialState) {

    // See EOMDerivation.py for a derivation of _Ac and _Bc
    _Ac = Eigen::Matrix<double, 6, 6>{
        {0, 0, 0, 1.0, 0, 0},
        {0, 0, 0, 0, 1.0, 0},
        {0, 0, 0, 0, 0, 1.0},
        {0, -_g * (_m1 + _m2) / _mc, 0,
         -_Dc / _mc, _D1 / (_L1 * _mc), -_D2 / (_L1 * _mc)},
        {0, _g * (_mc + _m1 + _m2) / (_L1 * _mc),
         -_m2 * _g / (_L1 * _m1),
         _Dc / (_L1 * _mc), 
         -_D1 * (_mc + _m1) / (_L1 * _L1 * _mc * _m1),
         _D2 * (_L1 * _mc + _L2 * _m1 + _L2 * _mc) / (_L1 * _L1 * _L2 * _mc * _m1)},
        {0, -_g * (_mc + _m1 + _m2) / (_L1 * _mc),
         _g * (1 / _L2 + _m2 / _m1 * (1 / _L1 + 1 / _L2)),
         -_Dc / (_L1 * _mc),
         _D1 * (_L1 * _mc + _L2 * _m1 + _L2 * _mc) / (_L1 * _L1 * _L2 * _mc * _m1),
         -_D2 * (1 / (_L2 * _L2 * _m2) + 1 / (_L2 * _L2 * _m1) +
                 2 / (_L1 * _L2 * _m1) + 1 / (_L1 * _L1 * _mc) +
                 1 / (_L1 * _L1 * _m1))}};

    _Bc = Eigen::Vector<double, 6>{
        0, 0, 0, 1 / _mc, -1 / (_L1 * _mc), 1 / (_L1 * _mc)};

    // _Ad and _Bd are derived from using Rk4 on Ẋ = _Ac*X + _Bc*U
    _Ad = Eigen::Matrix<double, 6, 6>::Identity() + _dt * _Ac +
          pow(_dt, 2) / 2 * _Ac * _Ac + pow(_dt, 3) / 6 * _Ac * _Ac * _Ac +
          pow(_dt, 4) / 24 * _Ac * _Ac * _Ac * _Ac;

    _Bd =
        (_dt * Eigen::Matrix<double, 6, 6>::Identity() + pow(_dt, 2) / 2 * _Ac +
         pow(_dt, 3) / 6 * _Ac * _Ac + pow(_dt, 4) / 24 * _Ac * _Ac * _Ac) *
        _Bc;
}

/**
 * @brief Construct a new DIP object with default system parameters
 *
 */
DIP::DIP()
    : DIP(1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.01,
          Eigen::Vector<double, 6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)) {}

double DIP::getCartMass() { return _mc; }
double DIP::getMassOne() { return _m1; }
double DIP::getMassTwo() { return _m2; }
double DIP::getLinkOneLength() { return _L1; }
double DIP::getLinkTwoLength() { return _L2; }

/**
 * @brief Gets the (x,y) position of mass one in meters
 *
 * @return Eigen::Vector2d
 */
Eigen::Vector2d DIP::getMassOnePos() {
    return Eigen::Vector2d(_state[0] + _L1 * sin(_state[1]),
                           _L1 * cos(_state[1]));
}

/**
 * @brief Gets the (x,y) position of mass two in meters
 *
 * @return Eigen::Vector2d
 */
Eigen::Vector2d DIP::getMassTwoPos() {
    return getMassOnePos() + _L2 * Eigen::Vector2d(sin(_state[1] + _state[2]),
                                                   cos(_state[1] + _state[2]));
}

/**
 * @brief Gets current state of the DIP system as [cart position (m), theta
 * one (rads), theta two (rads), cart velocity (m/s), Omega one (rads/s),
 * Omega two (rads/s)]
 *
 * @return Eigen::Vector<double, 6>
 */
Eigen::Vector<double, 6> DIP::getState() { return _state; }

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
Eigen::Vector<double, 6> DIP::_getDerivative(
    const Eigen::Vector<double, 6> &currentState, double u) {
    // Unpack state
    double th1 = currentState[1];
    double th2 = currentState[2];
    double cartVel = currentState[3];
    double om1 = currentState[4];
    double om2 = currentState[5];
    double c12 = cos(th1 + th2);
    double s12 = sin(th1 + th2);
    double s1 = sin(th1);
    double s2 = sin(th2);
    double c1 = cos(th1);
    double c2 = cos(th2);
    // MẊ = Gg + C - D + U
    // See EOMDerivation.py for a derivation of M, G, C, D, and U
    Eigen::Matrix3d M{
        {_mc + _m1 + _m2, _L1 * c1 * (_m1 + _m2) + _m2 * _L2 * c12,
         _m2 * _L2 * c12},
        {_L1 * c1 * (_m1 + _m2) + _m2 * _L2 * c12,
         _L1 * _L1 * (_m1 + _m2) + 2.0 * _m2 * _L1 * _L2 * c2 + _L2 * _L2 * _m2,
         _m2 * _L1 * _L2 * c2 + _m2 * _L2 * _L2},
        {_L2 * _m2 * c12, _L2 * _m2 * (_L1 * c2 + _L2), _m2 * _L2 * _L2}};
    Eigen::Vector3d G{0.0, _L1 * s1 * (_m1 + _m2) + _L2 * s12 * _m2,
                      _L2 * _m2 * s12};
    Eigen::Vector3d C{_m1 * _L1 * s1 * om1 * om1 + _m2 *
                        (_L1 * s1 * om1 * om1 + _L2 * (om1 + om2) * (om1 + om2) * s12),
                     _m2 * (2.0 * _L1 * _L2 * s2 * om1 * om2 + _L1 * _L2 * s2 * om2 * om2),
                     -_m2 * _L1 * _L2 * s2 * om1 * om1};
    Eigen::Vector3d D{_Dc * cartVel, _D1 * om1, _D2 * om2};
    Eigen::Vector3d U{u, 0, 0};

    Eigen::Vector<double, 6> derivative;
    // First half of state's derivative is the second half
    derivative << cartVel, om1, om2, M.inverse() * (G * _g + C - D + U);

    return derivative;
}

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
Eigen::Vector<double, 6> DIP::_Rk4(Eigen::Vector<double, 6> &prevState,
                                   double u) {
    Eigen::Vector<double, 6> k1 = _getDerivative(prevState, u);
    Eigen::Vector<double, 6> k2 = _getDerivative(prevState + k1 / 2.0 * _dt, u);
    Eigen::Vector<double, 6> k3 = _getDerivative(prevState + k2 / 2.0 * _dt, u);
    Eigen::Vector<double, 6> k4 = _getDerivative(prevState + k3 * _dt, u);

    return prevState + _dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

/**
 * @brief Updates the state of the DIP system to a time dt away by
 * integrating with RK4
 *
 * @param u The input force on the cart
 * @return Eigen::Vector<double, 6> The new state of the DIP system
 */
Eigen::Vector<double, 6> DIP::updateState(double u) {
    _state = _Rk4(_state, u);
    // Convert angles to be within [-Pi, Pi]
    _state[1] = std::fmod(_state[1] + M_PI, 2 * M_PI) - M_PI;
    _state[2] = std::fmod(_state[2] + M_PI, 2 * M_PI) - M_PI;
    return _state;
}

/**
 * @brief Get the potential energy of the DIP at the current state
 *
 * @return double
 */
double DIP::getPotentialEnergy() {
    // Unpack state
    double th1 = _state[1];
    double th2 = _state[2];
    double c12 = cos(th1 + th2);
    double c1 = cos(th1);

    // Potential energy is measured from cart track. This means that when masses
    // are below cart, potential energy is negative. Also means cart has no
    // potential energy.
    double m1Pot = _L1 * _m1 * _g * c1;
    double m2Pot = _m2 * _g * (_L1 * c1 + _L2 * c12);

    return m1Pot + m2Pot;
}

/**
 * @brief Get the kinetic energy of the DIP at the current state
 *
 * @return double
 */
double DIP::getKineticEnergy() {
    double th1 = _state[1];
    double th2 = _state[2];
    double cartVel = _state[3];
    double om1 = _state[4];
    double om2 = _state[5];
    double c12 = cos(th1 + th2);
    double s12 = sin(th1 + th2);
    double s1 = sin(th1);
    double c1 = cos(th1);

    // See EOMDerivation.py for the derivation of the kinetic energies
    double cartKin = 0.5 * _mc * cartVel * cartVel;
    double m1Kin = 0.5 * _m1 *
                   (_L1 * _L1 * om1 * om1 + 2.0 * _L1 * c1 * om1 * cartVel +
                    cartVel * cartVel);
    double m2kin = 0.5 * _m2 *
                   ((_L1 * s1 * om1 + _L2 * (om1 + om2) * s12) *
                        (_L1 * s1 * om1 + _L2 * (om1 + om2) * s12) +
                    (_L1 * c1 * om1 + _L2 * (om1 + om2) * c12 + cartVel) *
                        (_L1 * c1 * om1 + _L2 * (om1 + om2) * c12 + cartVel));

    return cartKin + m1Kin + m2kin;
}

/**
 * @brief Get the continuous state-space matrix A
 *
 * @return Eigen::Matrix<double, 6, 6>
 */
Eigen::Matrix<double, 6, 6> DIP::getAc() { return _Ac; }

/**
 * @brief Get the continuous state-space matrix B
 *
 * @return Eigen::Vector<double, 6>
 */
Eigen::Vector<double, 6> DIP::getBc() { return _Bc; }

/**
 * @brief Get the discrete state-space matrix A
 *
 * @return Eigen::Matrix<double, 6, 6>
 */
Eigen::Matrix<double, 6, 6> DIP::getAd() { return _Ad; }

/**
 * @brief Get the discrete state-space matrix B
 *
 * @return Eigen::Vector<double, 6>
 */
Eigen::Vector<double, 6> DIP::getBd() { return _Bd; }
