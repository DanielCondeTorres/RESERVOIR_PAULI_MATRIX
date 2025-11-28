# RESERVOIR_PAULI_MATRIX

# Experiment: Information Propagation in a Noisy Ising Chain

This script simulates the dynamics of a 6-qubit quantum spin chain governed by the **Transverse Field Ising Hamiltonian**. The goal is to observe how information injected into the first qubit propagates through the system and how the magnetization $\langle X \rangle$ evolves as the input state sweeps across the Bloch sphere.

## 1. Physical System & Parameters

The system is modeled with the following characteristics:

* **System Size:** $N=6$ qubits.
* **Hamiltonian:** Transverse Field Ising Model with random couplings $J_{ij}$ and a strong constant transverse field $h$.
    $$H = \sum_{i<j} J_{ij} Z_i Z_j + h \sum_i X_i$$
* **Regime:** Strong Field ($h = 10J$) and Long Evolution Time ($dt = 10/J$).
* **Noise:** A global dephasing channel with strength $g=0.3$ is applied at each step to simulate decoherence.

### Key Parameters
| Parameter | Value | Description |
| :--- | :--- | :--- |
| `n_qubits` | 6 | Number of spins in the chain. |
| `MAX_PAULI_STRINGS` | 1024 | **Truncation limit.** Only the 1024 strongest correlations are kept. |
| `h_const` | 10.0 | Strong transverse magnetic field. |
| `dt` | 10.0 | Time duration between inputs (allows for relaxation). |
| `g_strength` | 0.3 | Strength of the dephasing noise (Z-axis decay). |

---

## 2. Experimental Protocol

The simulation runs for `k_max = 50` steps. In each step $k$, the following sequence occurs:

### A. Input Encoding (Bloch Sweep)
The script generates a sequence of input states that trace a path on the Bloch sphere from the Z-axis ($|0\rangle$) towards the X-axis ($|+\rangle$).
* **Rx component:** $\sqrt{k / k_{max}}$ (Increases from 0 to 1).
* **Rz component:** $\sqrt{1 - k / k_{max}}$ (Decreases from 1 to 0).

### B. State Injection (`inject_state`)
The state of **Qubit 0** is effectively reset and overwritten with the new Bloch vector:
$$\rho_{new} = \text{Tr}_0(\rho) \otimes \frac{1}{2}(I + r_x X_0 + r_z Z_0)$$
* This uses the unified `inject_state` function.
* **Note:** $r_y$ is kept at 0.

### C. Global Dephasing
A noise channel is applied to the entire system to simulate information loss and enforce the "fading memory" property required for Reservoir Computing.
$$\rho \to (1-g)\rho + g \sum Z \rho Z$$
*(Note: Depending on the exact `apply_global_dephasing` implementation, this typically dampens off-diagonal terms).*

### D. Measurement
The expectation value of the magnetization in the X-basis, $\langle X_i \rangle$, is calculated for every qubit $i=0 \dots 5$.
$$\langle X_i \rangle = \text{Tr}(\rho \cdot X_i)$$

### E. Time Evolution (Stabilized RK4)
The system evolves under the Hamiltonian for a total time `dt = 10.0`.
* **Sub-stepping:** Because the field $h=10.0$ is strong, the dynamics are fast. The time step is split into **2500 sub-steps** to ensure the numerical stability of the Runge-Kutta 4 integrator.
* **Truncation:** Inside this intense loop, `truncate_operator!` is called to discard negligible Pauli strings. This keeps the memory usage constant (at 1024 terms) rather than growing exponentially ($4^6$).

---

## 3. Visual Output

The script concludes by plotting the evolution of $\langle X_i \rangle$ for all qubits over the simulation steps.

**Expected Result:**
* You should observe the "input" qubit (Q0) following the injection sweep closely.
* Subsequent qubits (Q1, Q2...) should show delayed or dampened responses depending on how well the information propagates through the noisy, truncated chain.
* The truncation ($M=1024$) simulates an approximation method (like MPO), allowing us to study how limited quantum resources affect the simulation accuracy compared to an exact full-state simulation.





Here is the explanation for the `run_online_protocol` script. This script implements a classic **Online Quantum Reservoir Computing (QRC)** task.

While the previous experiment (the Bloch sweep) was testing the *spatial* properties of the injection (how $X$ and $Z$ relate), this experiment tests the **temporal** properties: how the system remembers and processes a time-series input.

-----

# Experiment: Online QRC with Square Wave Input

This script runs a full **time-series processing simulation**. It injects a changing input signal into the first qubit of a spin chain and tracks how that signal propagates and influences the rest of the qubits (the "Reservoir") over time.

## 1\. The Input Signal: Square Wave

Unlike the previous smooth sweep, this script uses a **Square Wave** input signal:

  * **Pattern:** Alternates between $+0.9$ and $-0.9$ (High/Low).
  * **Purpose:** This tests the **reactivity** and **memory** of the reservoir.
      * Can the system jump quickly when the input changes?
      * Do the other qubits (Q1-Q5) "echo" this jump with a delay?

<!-- end list -->

```julia
inputs = [isodd(i) ? 0.9 : -0.9 for i in 1:n_inputs]
```

## 2\. The Protocol Loop

The simulation iterates through `n_inputs` (time steps). In every step, three distinct phases occur:

### Phase A: Injection (Erase & Write)

We take the current value of the square wave, $u_k$ (either 0.9 or -0.9), and inject it into **Qubit 0**.

  * **Function:** `rho = inject_state(rho, 0, input_u)`
  * **Physics:**
    $$\rho_{new} = \text{Tr}_0(\rho_{old}) \otimes \frac{1}{2}(I + u_k Z_0)$$
  * **Effect:** This forces Qubit 0 to the state corresponding to the input, erasing its previous history, but keeping the entanglement of the rest of the chain.

### Phase B: Reservoir Mixing (Evolution)

The system evolves for a duration `steps_per_input`. This is where the "computing" happens.

1.  **Hamiltonian Evolution (`step_rk4`):** The information in Qubit 0 spreads to Qubit 1, then Qubit 2, etc., via the interaction terms $J Z_i Z_{i+1}$.
2.  **Fading Memory (`apply_global_dephasing`):** We apply noise to dampen old correlations. This ensures the reservoir doesn't get saturated with infinite history.

### Phase C: Readout (Measurement)

We measure the expectation value of **Z** ($\langle Z_i \rangle$) for every qubit.

  * **Why Z?** Because we injected the information into the Z-axis. We want to see how this Z-magnetization travels down the chain.

## 3\. Key Differences from the Previous Experiment

| Feature | Previous Experiment (Bloch Sweep) | **This Experiment (Online Protocol)** |
| :--- | :--- | :--- |
| **Input** | Gradient vector $(\sqrt{k}, 0, \sqrt{1-k})$ | Scalar Time Series $(+0.9, -0.9, \dots)$ |
| **Goal** | Study state injection physics | Study **Time-Series Memory** & Dynamics |
| **Measurement** | Measured $\langle X \rangle$ (Coherence) | Measures $\langle Z \rangle$ (Population/Bit) |
| **Visual** | A curve sweeping up/down | Oscillating waves (Square pulses) |

## 4\. Expected Outcome

When you run this, the plot will show:

1.  **Qubit 0 (Input):** A perfect (or near-perfect) square wave, flipping between 0.9 and -0.9 immediately.
2.  **Qubit 1:** A similar square wave, but slightly **delayed** and smoother (filtered).
3.  **Qubit 2, 3, 4, 5:** Progressively smaller, smoother, and more delayed waves. This represents the **fading memory** of the reservoir.

This demonstrates that the quantum system acts as a **non-linear fading memory filter**, which is the core requirement for Reservoir Computing.
