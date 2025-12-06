# Chemical Plant – Decentralized, Distributed, and Centralized Control

This repository contains MATLAB code for analyzing and controlling a **12-state, 6-input chemical plant** using:

- Subsystem decomposition  
- Continuous‐time and discrete‐time modeling  
- Fixed-mode analysis under different information structures  
- LMI-based controller synthesis (CT and DT)  
- Closed-loop simulations for several control architectures  

The project follows the structure of the **Networked and Distributed Control Systems** course.

---

## 📌 1. System Description

The chemical plant model (`MAT08ChemicalPlant.m`) provides:

- Continuous-time dynamics  
  \[
  \dot{x} = A x + B u
  \]
- State matrix `A` (12×12)  
- Input matrix `B` (12×6)  
- Output matrix `C = I₁₂`  

The states correspond to **three chemical reactors**, each with 4 states.

---

## 📌 2. Subsystem Decomposition

The system is decomposed into **3 subsystems**:

- **Subsystem 1:** states 1–4, inputs 1–2  
- **Subsystem 2:** states 5–8, inputs 3–4  
- **Subsystem 3:** states 9–12, inputs 5–6  

This creates decentralized block structures:

- `Bdec{i}` – input matrices for each subsystem  
- `Cdec{i}` – output matrices for each subsystem  

These are later used for fixed-mode analysis and structured LMI design.

---

## 📌 3. CT & DT Models

### Continuous-Time
- Eigenvalues of `A` are computed  
- Spectral abscissa determines stability  

### Discretization
The system is sampled with **h = 0.1 s**:
\[
F = e^{Ah}, \quad G = A^{-1}(F - I)B
\]

We also build:

- `Gdec{i}` – discrete input matrices per subsystem  
- `Hdec{i}` – outputs (same as `Cdec{i}`)

---

## 📌 4. Fixed Modes Analysis

We check fixed modes under several **information structures** using:

```matlab
di_fixed_modes(A, Bdec, Cdec, N, ContStruc, rounding_n)
