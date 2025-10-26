# Gas–Gas Ejector — Quasi-1D MATLAB Model

**Author:** Thomas Abraham (22JE1018)  
**Institution:** IIT (ISM) Dhanbad  
**Specialization:** Computational Fluid Dynamics (CFD)  
**Language:** MATLAB  

---

## Project Overview

I developed a quasi-one-dimensional MATLAB solver to model compressible flow in gas–gas ejectors, coding the governing equations entirely from first principles.

Unlike standard Fanno-only formulations, this model explicitly couples **wall friction** and **variable cross-sectional area**, allowing it to capture both **viscous losses** and **area-driven acceleration** effects along the ejector axis.

The solver predicts **Mach number**, **pressure**, and **temperature** variations across each ejector segment (converging nozzle, mixing duct, and diffuser). It provides a compact analytical–numerical framework suitable for **parametric studies**, **geometry optimization**, and **pre-CFD design validation**.

This implementation represents the **quasi-1D analytical stage** and will be validated against **3D ANSYS Fluent simulations** in the next phase.

---

## ⚙️ Numerical Methodology

| Region | Method | Description |
|--------|---------|-------------|
| Converging/Diverging Nozzles | `ode45` | Solves quasi-1D compressible flow equations with friction and area variation |
| Mixing Duct | `fsolve` (Optimization Toolbox) | Simultaneously solves coupled mass and energy conservation equations for the unknown mixture temperature (`T₅`) and velocity (`V₅`) |
| Outlet Section | `fzero` | Determines downstream Mach number under frictional losses |
| Solver Configuration | `optimoptions` | Controls tolerance and iteration limits for nonlinear solvers |

---

## 📄 Files Included
| File | Description |
|------|--------------|
| `GasEjector_Model_Code.m` | MATLAB script implementing the quasi-1D ejector model |
| `GasEjector_Presentation.pptx` | Presentation explaining governing equations, implementation, and results |

---

## 📈 Outputs
- Mach number variation in converging and diverging sections  
- Pressure and temperature distributions along the axis  
- Solved mixture properties at the outlet  
- Generated plots for all major flow variables  

---

## Key Concepts Demonstrated
- Quasi-1D compressible flow modelling  
- Coupling of friction and variable-area effects  
- MATLAB solver integration (`ode45`, `fsolve`, `fzero`)  
- Implementation of conservation laws in ejector physics  
- Pre-CFD analytical validation workflow  

---

## Future Work
- Validation against full 3D ANSYS Fluent simulations  
---

## Author
**Thomas Abraham (22JE1018)**  
B.Tech in Mechanical Engineering, IIT (ISM) Dhanbad  
*(Batch of 2026)*  
Focus: Computational Fluid Dynamics (CFD)

---

## Acknowledgment
This project was carried out under the guidance of **Prof. S. K. Das**, Department of Mechanical Engineering, IIT (ISM) Dhanbad.  
I would also like to thank **Shandharb Singh** for his valuable assistance and insights during the development and implementation of this solver.

---

## 🧾 License
Released under the **MIT License** — feel free to use or adapt with attribution.

