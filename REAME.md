# Dynamic Low-Rank Approximation for Large-Scale Models

A modern model-order reduction approach to solving large-scale ODE systems. 
Two main examples are provided in this repository to demonstrate the methods:
* **Lyapunov Equation** (Script: `Lyapunov.m`)
* **Allen-Cahn Equation** (Script: `Allen-Cahn.m`)

## 🧠 Core Functions & Algorithms

### 1. Basic DLR Approximation
* `dlr.m` : Implementation of the basic approach to dynamical low-rank approximation.
  > **Reference:** Koch, O., & Lubich, C. (2007). *Dynamical low-rank approximation*. SIAM Journal on Matrix Analysis and Applications, 29(2), 434-454.

### 2. BUG Integrators
* `[rk_method]BUG.m` : The BUG integrator implemented with respective Runge-Kutta methods. In this work, the following methods are proposed:
  * Heun
  * Euler
  * Midpoint
  
  > **Reference:** Nobile, F., & Riffaud, S. (2025). *Robust high-order low-rank BUG integrators based on explicit Runge-Kutta methods*. arXiv preprint arXiv:2502.07040.

#### Optimized Implementations
Specific, optimized versions of the HeunBUG integrator for the Allen-Cahn and Lyapunov equations are given in the following functions:
* `HeunBUGLyapunov.m`
* `HeunBUGAllenCahn.m`

## 🛠️ Requirements 
* **MATLAB R2024a** or newer.

## 🚀 How to Run
1. Clone the repository and set it as your Current Folder in MATLAB.
2. Open and run either `Lyapunov.m` or `Allen-Cahn.m` to start the simulations and view the results.
