# DIP-Control
Dynamical simulation and control of the unstable double inverted pendulum (DIP) on a cart.
<p align="center">
  <img src="https://github.com/zsilberstein/DIP-Control/blob/main/images/DIP.gif?raw=true" width="640" height="360"/>
</p>

## Features
- Adjustable double inverted pendulum system on cart
- Adjustable initial conditions
- Open-loop simulation of system
- Continuous and Discrete LQR control of system

## Dependencies
All dependencies in this project work on Linux, macOS, and Windows. In addition to C++ 20 and a c++ complier, this project depends on the following:
- [CMake](https://cmake.org/) - The build system used for this project. See [here](https://cmake.org/download/) for getting started with CMake.
- [raylib](https://www.raylib.com/index.html) - A library built on GLFW for video game programming. If not already installed, the CMakeLists.txt file will install it for you.
- [raygui](https://github.com/raysan5/raygui) - A header only file to help design simple GUIs. This file is already included in the /src directory. 
- [Eigen](http://eigen.tuxfamily.org/) - A C++ template library for linear algebra. Eigen must be installed prior to running this project, see [here](https://eigen.tuxfamily.org/dox/GettingStarted.html) for simple  instructions on how to get started with Eigen.
- [GoogleTest](https://google.github.io/googletest/) - A C++ testing framework by Google. This is required to run tests and, if not already installed, the CMakeLists.txt file will install it for you.


## Installation
The following instructions are designed for Linux. Instructions for other operating systems may slightly vary.
Clone this project with git:
```
git clone https://github.com/zsilberstein/DIP-Control.git
```
Navigate into the project:
```
cd DIP-Control
```
Create a build directory:
```
mkdir build
```
Build the project with CMake:
```
cmake -S . -B ./build/
```
Move to the build directory and make the project:
```
cd build && make
```
Run the executable:
```
./DIP-Control
```
To run the tests:
```
ctest
```

## The DIP System
The double inverted pendulum (DIP) on a cart system consists of cart and two masses connected by rigid, massless rods as illustrated below. In addition, there is an horizontal external force u on the cart and damping Dc, D1, D2 on the cart velocity, angular velocity of mass one, and angular velocity of mass two respectively.
<p align="center">
  <img src="https://github.com/zsilberstein/DIP-Control/blob/main/images/DIP_FBD.png?raw=true" width="512" height="288"/>
</p>

The EOMDerivation.py file derives the equations of motion of the system with lagrangian dynamics. For simplicity, the results are shown below.

Cart EOM:  

<img src="https://latex.codecogs.com/svg.image?D_c%20v+m_c\dot{v}+m_1\left(\dot{v}+L_1(\dot{\omega_1}\cos{\theta_1}-\omega_1^2\sin{\theta_1})\right)+m_2\left(\dot{v}+L_1(\dot{\omega_1}\cos{\theta_1}-\omega_1^2\sin{\theta_1})+L_2\left((\dot{\omega_1}+\dot{\omega_2})\cos{(\theta_1+\theta_2)}-(\omega_1+\omega_2)^2\sin{(\theta_1+\theta_2)}\right)\right)-u=0" />


<!-- $$
D_c v + m_c \dot{v} + m_1 \left(\dot{v} + L_1 (\dot{\omega_1} \cos{\theta_1} - \omega_1^2  \sin{\theta_1})\right) + 
m_2 \left(\dot{v} + L_1 (\dot{\omega_1} \cos{\theta_1} - \omega_1^2  \sin{\theta_1}) 
+ L_2 \left((\dot{\omega_1} + \dot{\omega_2})\cos{(\theta_1+\theta_2)} - (\omega_1 + \omega_2)^2\sin{(\theta_1+\theta_2)}\right)\right) - u = 0
$$ -->
Mass One EOM:

<img src="https://latex.codecogs.com/svg.image?D_1\omega_1+m_1%20L_1\left(L_1\dot{\omega_1}-g\sin{\theta_1}+\dot{v}\cos{\theta_1}\right)+\\m_2\left(L_1\left(L_1\dot{\omega_1}-g\sin{\theta_1}+\dot{v}\cos{\theta_1}\right)+L_2\left(-2%20L_1\omega_1\omega_2\sin{\theta_2}-L_1\omega_2^2\sin{\theta_2}+2%20L_1\dot{\omega_1}\cos{\theta_2}+L_1\dot{\omega_2}\cos{\theta_2}+L_2\dot{\omega_1}+L_2\dot{\omega_2}-g\sin{\left(\theta_1+\theta_2\right)}+\dot{v}\cos{\left(\theta_1+\theta_2\right)}\right)\right)=0" />

<!-- $$
D_1 \omega_1 + m_1 L_1 \left(L_1 \dot{\omega_1} - g \sin{\theta_1} + \dot{v} \cos{\theta_1} \right) + \\
m_2 \left(L_1 \left(L_1 \dot{\omega_1} - g \sin{\theta_1} + \dot{v} \cos{\theta_1} \right) + 
L_2 \left(- 2 L_1 \omega_1 \omega_2 \sin{\theta_2} - L_1 \omega_2^2 \sin{\theta_2} + 2 L_1 \dot{\omega_1} \cos{\theta_2} + L_1 \dot{\omega_2} \cos{\theta_2} + L_2 \dot{\omega_1} + L_2 \dot{\omega_2} - g \sin{\left(\theta_1 + \theta_2 \right)} + \dot{v} \cos{\left(\theta_1 + \theta_2 \right)}\right)\right)=0
$$ -->

Mass Two EOM:

<img src="https://latex.codecogs.com/svg.image?D_2\omega_2+m_2%20L_2\left(L_1\omega_1^2\sin{\theta_2}+L_1\dot{\omega_1}\cos{\theta_2}+L_2\dot{\omega_1}+L_2\dot{\omega_2}-g\sin{\left(\theta_1+\theta_2\right)}+\dot{v}\cos{\left(\theta_1+\theta_2\right)}\right)=0" />
<!-- $$
D_2 \omega_2 + m_2 L_2 \left(L_1 \omega_1^2 \sin{\theta_2} + L_1 \dot{\omega_1} \cos{\theta_2} + L_2 \dot{\omega_1} + L_2 \dot{\omega_2} - g \sin{\left(\theta_1 + \theta_2 \right)} + \dot{v} \cos{\left(\theta_1 + \theta_2 \right)}\right)=0
$$ -->

## LQR Control
A linear–quadratic regulator (LQR) was created to stabilize the DIP system about its unstable equilibrium with both rods straight up. Since LQR is a linear control method, the initial conditions when using the LQR have been restricted to enforce the small angle approximation which is necessary for the nonlinear DIP system to be approximated as linear. That said, there are still a few cases where the system will leave the small angle range and the controller will fail. The LQR controller developed is capable of handling both continuous and discrete time systems. 

## Future Improvements
- Add additional controllers
- Add plotting of state over time
- Improve the GUI

### Have your own improvement in mind? 
Feel free to fork the project and then put in a pull request with the desired improvement. Please be sure to include relevant tests in any pull request.

## License
This project is licensed under the [MIT License](https://github.com/zsilberstein/robot-gait-vis/blob/master/LICENSE). See the LICENSE page for more information. Be sure to also read the license of any project dependency. 

## Author
**Zach Silberstein**  
Email: **zach.silberstein@gmail.com**   
Github: [@zsilberstein](https://github.com/zsilberstein)