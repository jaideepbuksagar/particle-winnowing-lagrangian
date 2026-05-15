# Particle Winnowing — Lagrangian Dynamics Simulation

Numerical simulation of the winnowing grain-separation process using Lagrangian particle tracking.  
Implements and benchmarks **Euler** and **4th-order Runge–Kutta (RK4)** solvers for a coupled ODE system governing particle motion in a spatially varying planar air jet.

> Submitted as an advanced optional exercise for *Simulations of Mechanical Processes*, OVGU Magdeburg (SoSe 2025).

---

## Physics

Winnowing separates heavy grains from light chaffs by passing a mixture through a horizontal air jet. Particles are released from a point $(x_0, y_0)$ with no initial velocity and fall under gravity while being deflected by the jet.

The motion of each particle is governed by:

$$\frac{dx_p}{dt} = u_p, \qquad \rho_p V_p \frac{du_p}{dt} = F_\text{gravity} + F_\text{buoyancy} + F_\text{drag}$$

**Fluid velocity field** (planar jet):

$$u_f(x,y) = 6.2\, u_0 \sqrt{\frac{h}{x}} \exp\!\left(-50\,\frac{y^2}{x^2}\right), \qquad v_f = 0$$

**Drag** (Schiller–Naumann):

$$C_D = \begin{cases} \dfrac{24}{Re_p}\!\left(1 + 0.15\,Re_p^{0.687}\right) & Re_p < 800 \\ 0.44 & \text{otherwise} \end{cases}$$

**Separation criterion:** heavy grains land left of $(x_c, y_c)$ → Bin 1; light chaffs land right → Bin 2.

---

## Tasks Solved

| Task | Description |
|---|---|
| 1 | Implement Euler and RK4 solvers; plot trajectories of grain and chaff |
| 2 | Verify order-of-accuracy via convergence study (Euler ≈ 1, RK4 ≈ 4) |
| 3 | Find minimum jet velocity $u_{0,\min}$ for separation to 3 s.f. via bisection |
| 4 | Analyse separation robustness under non-zero initial velocities (terminal velocity, extreme angles) |

---

## Key Results

**Task 1 — Trajectories**

Grain (ρ = 750 kg/m³, d = 2.5 mm) falls nearly vertically → Bin 1.  
Chaff (ρ = 50 kg/m³, d = 3.25 mm) is swept by the jet → Bin 2.  
Correct separation confirmed at $u_0 = 20$ cm/s.

**Task 2 — Convergence**

Error measured at a fixed intermediate time $t^*$ to isolate integration error from endpoint interpolation noise.

| Scheme | Measured order | Expected |
|---|---|---|
| Euler | ≈ 1.0 | 1 |
| RK4 | ≈ 4.0 | 4 |

**Task 3 — Minimum jet velocity**

Bisection on $u_0 \in [0, 20]$ cm/s, tolerance $10^{-6}$ m/s:

$$u_{0,\min} = \text{[value to 3 s.f.]} \text{ cm/s}$$

**Task 4 — Non-zero initial velocities**

Terminal velocity computed via Newton–Raphson on the implicit force-balance equation.  
Separation tested at extreme angles (−72°, −108°) for both particle types.

---

## Repository Structure

```
particle-winnowing-lagrangian/
├── src/
│   └── winnowing.m          # Main script — solvers, convergence, bisection, Task 4
├── figures/
│   ├── trajectories_euler.png
│   ├── trajectories_rk4.png
│   ├── convergence_study.png
│   └── bisection_convergence.png
├── report/
│   └── winnowing_report.pdf  # Full write-up with analysis
├── .gitignore
└── README.md
```

---

## Usage

**Requirements:** MATLAB R2018b or later (no additional toolboxes required).

```matlab
% Clone and run
git clone https://github.com/jaideepbuksagar/particle-winnowing-lagrangian
cd particle-winnowing-lagrangian/src
% Open MATLAB, navigate to src/, and run:
winnowing
```

All figures are generated and saved automatically. Console output reports landing positions, bin classifications, convergence orders, and minimum jet velocity.

**Parameters** are defined at the top of `winnowing.m` and can be freely modified:

```matlab
uf0   = 0.20;   % jet velocity (m/s)
dt    = 1e-3;   % time step (s)
```

---

## Physical Parameters

| Parameter | Grain | Chaff |
|---|---|---|
| Density (kg/m³) | 750 | 50 |
| Diameter (mm) | 2.5 | 3.25 |

| Flow parameter | Value |
|---|---|
| Jet velocity $u_0$ | 20 cm/s |
| Slot height $h$ | 10 cm |
| Air density | 1.2 kg/m³ |
| Air viscosity | 1.8 × 10⁻⁵ Pa·s |
| Release point $(x_0, y_0)$ | (50 cm, 50 cm) |
| Bin boundary $(x_c, y_c)$ | (55 cm, −50 cm) |

---

## Numerical Methods

**Euler (1st order):**
$$Y^{n+1} = Y^n + \Delta t \cdot f(Y^n)$$

**RK4 (4th order):**
$$Y^{n+1} = Y^n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

**Bisection** used for Task 3 (minimum $u_0$) — bracket $[0, 20]$ cm/s, tolerance $10^{-6}$ m/s, ~20 iterations.

**Newton–Raphson** used for Task 4 (terminal velocity) — converges in ~5 iterations from Stokes initial guess.

---

## License

MIT — free to use and adapt with attribution.

---

## Author

**Jaideep Buksagarmath**  
MSc Chemical & Energy Engineering, OVGU Magdeburg  
B.Tech. Chemistry, IIT Guwahati | GTRE-DRDO Research Intern  
[GitHub](https://github.com/jaideepbuksagar) · [LinkedIn](https://linkedin.com/in/jaideep-buksagarmath/)
