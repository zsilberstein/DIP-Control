import sympy as sp

sp.init_printing()
# System constants and time
mc, m1, m2, L1, L2, g, Dc, D1, D2, dt, t = sp.symbols(
    'mc m1 m2 L1 L2 g Dc D1 D2 dt t', positive=True)
u = sp.symbols('u')

# x as a function of time and it's derivatives
xt = sp.Function('x')(t)
x_d = sp.diff(xt, t)
x_d_d = sp.diff(xt, t, 2)

# Theta 1 as a function of time and it's derivatives
th1t = sp.Function('th1')(t)
th1_d = sp.diff(th1t, t)
th1_d_d = sp.diff(th1t, t, 2)

# Theta 2 as a function of time and it's derivatives
th2t = sp.Function('th2')(t)
th2_d = sp.diff(th2t, t)
th2_d_d = sp.diff(th2t, t, 2)

# Cart kinetic and potential energy
Xc = xt
Xc_dot = x_d
Tc = 1/2*mc*Xc_dot**2
Vc = 0

# Mass one kinetic and potential energy
X1 = xt + L1*sp.sin(th1t)
X1_dot = sp.diff(X1, t)
Y1 = L1*sp.cos(th1t)
Y1_dot = sp.diff(Y1, t)
T1 = 1/2 * m1 * (X1_dot**2+Y1_dot**2)
V1 = m1 * g * Y1

# Mass two kinetic and potential energy
X2 = xt + L1*sp.sin(th1t) + L2*sp.sin(th1t+th2t)
X2_dot = sp.diff(X2, t)
Y2 = L1*sp.cos(th1t) + L2*sp.cos(th1t+th2t)
Y2_dot = sp.diff(Y2, t)
T2 = 1/2 * m2 * (X2_dot**2 + Y2_dot**2)
V2 = m2 * g * Y2

# Lagrangian
L = Tc + T1 + T2 - Vc - V1 - V2

# Set up equations of motion
cart_EOM = sp.simplify(sp.diff(sp.diff(L, x_d), t) -
                       sp.diff(L, xt) - u + x_d*Dc)
massOne_EOM = sp.simplify(sp.diff(sp.diff(L, th1_d), t) -
                          sp.diff(L, th1t) + th1_d*D1)
massTwo_EOM = sp.simplify(sp.diff(sp.diff(L, th2_d), t) -
                          sp.diff(L, th2t) + th2_d*D2)

# State variables
x, th1, th2, v, omega1, omega2, v_dot, omega1_dot, omega2_dot = sp.symbols(
    'x theta1 theta2 v omega1 omega2 vdot omegadot1 omegadot2')

# Reduce EOMs to be first order with state variables
cart_EOM = cart_EOM.subs([(x_d_d, v_dot), (th1_d_d, omega1_dot), (th2_d_d, omega2_dot),
                          (x_d, v), (th1_d, omega1), (th2_d, omega2), (x_d, v),
                          (th1t, th1), (th2t, th2)])
massOne_EOM = massOne_EOM.subs([(x_d_d, v_dot), (th1_d_d, omega1_dot), (th2_d_d, omega2_dot),
                                (x_d, v), (th1_d, omega1), (th2_d, omega2), (x_d, v),
                                (th1t, th1), (th2t, th2)])
massTwo_EOM = massTwo_EOM.subs([(x_d_d, v_dot), (th1_d_d, omega1_dot), (th2_d_d, omega2_dot),
                                (x_d, v), (th1_d, omega1), (th2_d, omega2), (x_d, v),
                                (th1t, th1), (th2t, th2)])

print('Cart EOM')
print(sp.pretty(cart_EOM))
print()
print('Mass One EOM')
print(sp.pretty(massOne_EOM))
print()
print('Mass Two EOM')
print(sp.pretty(massTwo_EOM))
print()
