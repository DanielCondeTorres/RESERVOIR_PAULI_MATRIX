
***

# Quantum State Injection & Noise Simulation (Schrödinger Picture) injection.jl

This module implements key operations for evolving a quantum system in the **Schrödinger picture** (evolving the density matrix $\rho$). It focuses on two mechanisms essential for Quantum Reservoir Computing:
1.  **Dephasing Noise:** Simulating information loss.
2.  **Input Injection:** The "Reset & Write" mechanism used to feed data into the quantum reservoir.

## 1. Pauli Representation Logic

The code utilizes a symplectic bit-mask representation for Pauli strings. For a specific qubit $k$, the operator $P_k$ is defined by two bits: $x_k$ (`x_mask`) and $z_k$ (`z_mask`).

### `get_pauli_type`
This helper function identifies the Pauli operator acting on a specific qubit index.

**Mapping Logic:**
$$\text{Type} = x_k + 2z_k$$

| Operator | $x_k$ | $z_k$ | Value | Description |
| :--- | :---: | :---: | :---: | :--- |
| $I$ (Identity) | 0 | 0 | **0** | Commutes with everything. |
| $X$ | 1 | 0 | **1** | Bit-flip. |
| $Z$ | 0 | 1 | **2** | Phase-flip. |
| $Y$ ($iXZ$) | 1 | 1 | **3** | Combined flip. |

---

## 2. Global Dephasing Channel

### `apply_global_dephasing`
This function simulates the decoherence of the quantum state. Based on the implementation, it preserves operators that do not have the $z$-bit set (i.e., $I$ and $X$) and decays operators that do ($Z$ and $Y$).

**Mathematical Definition:**
This corresponds to a **dephasing channel along the X-basis**. The channel acts on the density matrix $\rho$ as:

$$\mathcal{E}(\rho) = (1-g)\rho + g X \rho X$$

Under this operation, the Pauli basis transforms as:
* $X \to X$ (Eigenstate, no decay)
* $Z \to (1-2g)Z$
* $Y \to (1-2g)Y$

**Algorithmic Implementation:**
1.  **Decay Factor:** $\lambda = (1 - 2g)$.
2.  **Commutation Check:** It counts how many $Z$ or $Y$ operators exist in a Pauli string `p` using `count_ones(p.z_mask)`.
3.  **Coefficient Update:**
    $$c_p' = c_p \cdot \lambda^{N_{\text{non-commuting}}}$$
    Where $c_p$ is the coefficient of the Pauli string in the density matrix.

---

## 3. General Bloch State Injection

### `inject_general_bloch_state`
This function implements the **"Reset and Write"** operation. It mathematically equates to performing a partial trace over a specific qubit $k$ (discarding its old state) and tensor-producting the remaining system with a new state defined by a Bloch vector $\vec{r} = (r_x, r_y, r_z)$.

**Mathematical Formulation:**
Given an input density matrix $\rho_{in}$, the new state $\rho_{out}$ is:

$$\rho_{out} = \text{Tr}_k(\rho_{in}) \otimes \rho_{new}^{(k)}$$

Where the new single-qubit state is:
$$\rho_{new}^{(k)} = \frac{1}{2} \left( I + r_x X_k + r_y Y_k + r_z Z_k \right)$$

**Algorithmic Steps:**

1.  **Partial Trace (Filtering):**
    Iterate through all Pauli strings $P$ in $\rho_{in}$.
    * If $P$ acts as $X, Y, \text{or } Z$ on qubit $k$, $\text{Tr}(P) = 0$. These terms are discarded.
    * If $P$ acts as $I$ on qubit $k$, it survives.

2.  **Tensor Product (Expansion):**
    For every surviving term $P_{sys} \otimes I_k$ with coefficient $c_p$:
    * **Keep Identity part:** Add $c_p (P_{sys} \otimes I_k)$ to $\rho_{out}$.
    * **Add X part:** Add $c_p \cdot r_x (P_{sys} \otimes X_k)$ to $\rho_{out}$.
    * **Add Y part:** Add $c_p \cdot r_y (P_{sys} \otimes Y_k)$ to $\rho_{out}$.
    * **Add Z part:** Add $c_p \cdot r_z (P_{sys} \otimes Z_k)$ to $\rho_{out}$.

3.  **Fresh Injection:**
    If there was an implicit Identity term in the system (trace of the whole system is usually 1), explicitly add the Bloch vector components for the new qubit state.

---

## 4. Optimized Schrödinger Input Injection

### `inject_input_schrodinger`
This is a specialized version of the general injection optimized for **Quantum Reservoir Computing (QRC)**. In QRC, the input scalar $u$ is typically encoded into the $Z$-polarization of the qubit.

**Configuration:**
* $r_x = 0$
* $r_y = 0$
* $r_z = u$

**Mathematical Operation:**
$$\rho_{out} = \text{Tr}_k(\rho_{in}) \otimes \frac{1}{2}(I + u Z_k)$$

**Key Logic:**
1.  **Memory Retention:** Terms where qubit $k$ is $I$ are kept (representing the state of the rest of the reservoir).
2.  **Non-Linear Interaction:** The existing state of the reservoir (correlations) is multiplied by $u$ and attached to $Z_k$. This creates terms like $u (P_{others} \otimes Z_k)$.
3.  **Fresh Input:** A pure term $u Z_k$ is added.

**Why is this important?**
In Reservoir Computing, this operation represents the injection of the time-series data $u_t$ into the quantum substrate, creating a non-linear mapping between the previous reservoir state and the current input.
