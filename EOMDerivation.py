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
# Print EOMs
print('Cart EOM:')
print(sp.pretty(cart_EOM))
print()
print('Mass One EOM:')
print(sp.pretty(massOne_EOM))
print()
print('Mass Two EOM:')
print(sp.pretty(massTwo_EOM))

# Collect like terms to simplify building matrices
cartEOM_dict = sp.collect(cart_EOM.expand(),
                          [v_dot, omega1_dot, omega2_dot, g, Dc, u],
                          evaluate=False)
Mass1EOM_dict = sp.collect(massOne_EOM.expand(),
                           [v_dot, omega1_dot, omega2_dot, g, D1, u],
                           evaluate=False)
Mass2EOM_dict = sp.collect(massTwo_EOM.expand(),
                           [v_dot, omega1_dot, omega2_dot, g, D2, u],
                           evaluate=False)

# Build EOMs into matrix form as MẊ - C - G + D - U = 0
M = sp.Matrix([[cartEOM_dict[v_dot],
                sp.collect(cartEOM_dict[omega1_dot], L1*sp.cos(th1)),
                cartEOM_dict[omega2_dot]],
               [sp.collect(Mass1EOM_dict[v_dot], L1*sp.cos(th1)),
                sp.collect(Mass1EOM_dict[omega1_dot], L1*L1),
                sp.collect(Mass1EOM_dict[omega2_dot], m2*L2)],
               [Mass2EOM_dict[v_dot],
                sp.collect(Mass2EOM_dict[omega1_dot], L2*m2),
                Mass2EOM_dict[omega2_dot]]])

C = sp.Matrix([sp.collect(-cartEOM_dict[1], [m1, m2]),
               sp.collect(-Mass1EOM_dict[1], m2),
               sp.collect(-Mass2EOM_dict[1], m2)])

G = sp.Matrix([0,
               sp.collect(-Mass1EOM_dict[g], L1*sp.sin(th1)),
               -Mass2EOM_dict[g]])

D = sp.Matrix([Dc*cartEOM_dict[Dc],
               D1*Mass1EOM_dict[D1],
               D2*Mass2EOM_dict[D2]])

U = sp.Matrix([-u*cartEOM_dict[u], 0, 0])

X_dot = sp.Matrix([v_dot, omega1_dot, omega2_dot])

system = M*X_dot - C - G*g + D - U

# Confirm that matrix representation is equivalent to EOMs found earlier
if not sp.simplify(cart_EOM - system[0]) == 0:
    print('Warning: Difference found in matrix representation')
if not sp.simplify(massOne_EOM - system[1]) == 0:
    print('Warning: Difference found in matrix representation')
if not sp.simplify(massTwo_EOM - system[2]) == 0:
    print('Warning: Difference found in matrix representation')

# Print matrices
print()
print('-' * 150)
print()
print('Matrix representation of equations of motion: MẊ - C - G*g + D - U = 0')
print()
print('M matrix:')
print(sp.pretty(M))
print()
print('C matrix:')
print(sp.pretty(C))
print()
print('G matrix:')
print(sp.pretty(G))
print()
print('D matrix:')
print(sp.pretty(D))
print()
print('U matrix:')
print(sp.pretty(U))

# Calculate linearized state space matrices A and B such that Ẋ = AX + BU at the equilibrium point.
M_inv = M.adjugate() / M.det()
X_Dot = sp.Matrix([v, omega1, omega2, M_inv * (G*g + C - D + U)])

A = sp.simplify(X_Dot.jacobian([x, th1, th2, v, omega1, omega2]).subs(
    [(th1, 0), (th2, 0), (omega1, 0), (omega2, 0), (v, 0)]))

B = sp.simplify(X_Dot.jacobian([u]).subs(
    [(th1, 0), (th2, 0), (omega1, 0), (omega2, 0), (v, 0)]))
print()
print('-' * 150)
print()
print('Linearized state-space representation: Ẋ = AX + BU')
print('A matrix:')
print(sp.pretty(A))
print()
print(f'B matrix:')
print(sp.pretty(B))
