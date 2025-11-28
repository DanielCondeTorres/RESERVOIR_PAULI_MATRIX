

# ⚛️ Quantum State Injection & Noise Simulation (Schrödinger Picture) injection.jl

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

## 2. Global Dephasing Channel (apply_global_dephasing)

Here is the detailed explanation of the `apply_global_dephasing` function in English. This breakdown connects the Julia code logic directly to the physics of the **X-Dephasing Channel** using the Pauli matrix formalism.

-----

## Mathematical & Code Breakdown: Global Dephasing

This function simulates the interaction of the quantum system with a noisy environment. Specifically, it applies a **Phase Flip Channel along the X-basis**.

### 1\. The Physics (The "Why")

We model the noise as a process where, with probability $g$, an $X$ operator is applied to the system. In the density matrix formalism, this channel $\mathcal{E}$ is defined as:

$$
\mathcal{E}(\rho) = (1-g)\rho + g X \rho X
$$

To understand how this affects the state, we look at how the Pauli matrices transform under this operation. Recall the **Commutation Relations**:

  * $X$ commutes with $I$ and $X$ ($[X, X] = 0$).
  * $X$ anti-commutes with $Y$ and $Z$ ($\{X, Z\} = 0$, meaning $XZ = -ZX$).

Therefore, the coefficients $c_i$ of the Pauli string evolve as follows:

| Operator ($P$) | Transformation $X P X$ | Resulting Coefficient Decay |
| :--- | :--- | :--- |
| **I** | $X I X = I$ | $(1-g) + g(1) = \mathbf{1}$ (No Change) |
| **X** | $X X X = X$ | $(1-g) + g(1) = \mathbf{1}$ (No Change) |
| **Y** | $X Y X = -Y$ | $(1-g) + g(-1) = \mathbf{1 - 2g}$ |
| **Z** | $X Z X = -Z$ | $(1-g) + g(-1) = \mathbf{1 - 2g}$ |

**Conclusion:** The operators $Y$ and $Z$ are dampened by a factor $\lambda = (1-2g)$, while $I$ and $X$ remain untouched.

-----

### 2\. Code Walkthrough (The "How")

#### Step 1: Defining the Decay Factor

```julia
decay_factor = (1.0 - 2.0 * g)
```

This line calculates the dampening eigenvalue $\lambda$ derived in the table above.

  * If $g=0$ (no noise), `decay_factor = 1.0` (Identity).
  * If $g=0.5$ (maximum noise), `decay_factor = 0.0` (Complete destruction of quantum information in $Y/Z$).

#### Step 2: Iterating over the Pauli Strings

```julia
for (p, coeff) in rho
```

The density matrix `rho` is stored as a sparse list of Pauli strings (e.g., $X_1 Z_2$) and their complex coefficients. We process them one by one.

#### Step 3: Identifying Non-Commuting Operators

```julia
n_non_commuting = count_ones(p.z_mask)
```

This is the clever bitwise logic. We need to identify if an operator is $Y$ or $Z$.
In the symplectic bit-mask representation:

| Pauli Matrix | Binary Representation | `z_mask` Bit |
| :--- | :--- | :---: |
| $$I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$$ | `x=0, z=0` | **0** |
| $$X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$$ | `x=1, z=0` | **0** |
| $$Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$$ | `x=0, z=1` | **1** |
| $$Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}$$ | `x=1, z=1` | **1** |

  * **Logic:** As you can see, only $Z$ and $Y$ have the `z_mask` bit set to 1.
  * **Result:** `count_ones(p.z_mask)` counts exactly how many operators in the string anti-commute with the noise source ($X$).

#### Step 4: Calculating Total Decay

```julia
total_decay = decay_factor ^ n_non_commuting
```

Since the noise acts on each qubit independently, the decay is multiplicative.

  * If we have a term like $Z_1 \otimes Z_2$:
      * Qubit 1 contributes $(1-2g)$.
      * Qubit 2 contributes $(1-2g)$.
      * Total decay: $(1-2g)^2$.
  * If the term is $X_1 \otimes I_2$:
      * `n_non_commuting` is 0.
      * Total decay: $(1-2g)^0 = 1$. The term is preserved.

#### Step 5: Updating the Matrix

```julia
if abs(coeff * total_decay) > 1e-20
    new_rho[p] = get(new_rho, p, 0.0im) + (coeff * total_decay)
end
```

We multiply the old coefficient by the calculated decay. The `if` statement is a numerical threshold to clean up "dead" terms (floating point underflow) to keep the sparse representation efficient.

-----

### 3\. Visual Example: The Effect on the Bloch Sphere

Let's visualize a single qubit state $\rho$ composed of all Pauli matrices.

**Before Dephasing:**

$$
\rho = \frac{1}{2} \left( I + \color{blue}{0.8 X} + \color{red}{0.5 Y} + \color{red}{0.3 Z} \right)
$$

**Applying the Function** (with $g=0.1 \rightarrow \lambda = 0.8$):

1.  **$I$**: `z_mask` is 0. Factor = $1$. $\rightarrow$ Stays $I$.
2.  **$X$**: `z_mask` is 0. Factor = $1$. $\rightarrow$ Stays $\color{blue}{0.8 X}$.
3.  **$Y$**: `z_mask` is 1. Factor = $0.8$. $\rightarrow$ Becomes $0.5 \times 0.8 = \color{red}{0.4 Y}$.
4.  **$Z$**: `z_mask` is 1. Factor = $0.8$. $\rightarrow$ Becomes $0.3 \times 0.8 = \color{red}{0.24 Z}$.

**Resulting Matrix:**
The components aligned with the noise ($X$) survive. The perpendicular components ($Y, Z$) shrink.

$$
\rho_{new} = \begin{pmatrix} 
\dots & \text{Decayed Coherence} \\
\text{Decayed Coherence} & \dots
\end{pmatrix}
$$










Here is the complete explanation, integrating the matrix visualization with the Hamiltonian dynamics to explain the physical reasoning behind this specific noise channel.

***

# Quantum State, Hamiltonian & Selective Decay

## 1. Constructing the Qubit State

Any qubit state can be constructed as a linear combination of Pauli matrices:

$$
\rho = \frac{1}{2} \left( I + r_x X + r_y Y + r_z Z \right)
$$

If we perform the matrix summation, each coefficient maps to a specific location within the matrix:

$$
\rho = \frac{1}{2} 
\begin{pmatrix} 
1 + r_z & r_x - i r_y \\
r_x + i r_y & 1 - r_z
\end{pmatrix}
$$

* **$r_z$ (The Z coefficient):** Lives on the **DIAGONAL**. It controls the probability of the qubit being 0 or 1 (Population).
* **$r_x$ and $r_y$ (The X and Y coefficients):** Live **OFF-DIAGONAL**. They form the "coherence" (the complex number $r_x \pm i r_y$).

---

## 2. Numerical Example (Noise in X)

Let's recall the values from the previous example after the **dephasing** process (where noise was applied along the X-axis):

* $r_x = 0.8$ (Did **not** decay, because the noise source is aligned with X).
* $r_y = 0.4$ (Decayed, originally 0.5).
* $r_z = 0.24$ (Decayed, originally 0.3).

Now, let's substitute these numbers into the matrix:

$$
\rho_{new} = \frac{1}{2} 
\begin{pmatrix} 
1 + 0.24 & 0.8 - i(0.4) \\
0.8 + i(0.4) & 1 - 0.24
\end{pmatrix}
$$

Which simplifies to:

$$
\rho_{new} = \frac{1}{2} 
\begin{pmatrix} 
1.24 & 0.8 - 0.4i \\
0.8 + 0.4i & 0.76
\end{pmatrix}
$$

---

## 3. The Physical Context: The Hamiltonian

To understand **why** we chose this specific decay pattern, we must look at the system's "Engine": the Hamiltonian. In Quantum Reservoir Computing, we typically use the **Transverse Field Ising Model**:

$$
H = \underbrace{\sum_{i,j} J_{ij} Z_i Z_j}_{\text{Coupling (Memory)}} + \underbrace{\sum_{i} h_i X_i}_{\text{Transverse Field (Quantum Driver)}}
$$



This Hamiltonian has two competing parts:
1.  **The $Z$ terms:** These align the qubits (like magnets) and store the input information.
2.  **The $X$ terms:** This is a magnetic field perpendicular to the spins. It forces them to rotate and creates quantum superposition.

---

## 4. Why do we want $Z$ to decay? (The "Leaky Bucket")

The input $u_t$ is injected directly into the $Z$ component: `inject_input_schrodinger(..., u)`.

**If $Z$ did not decay:**
* The system would have perfect memory.
* The input $u_1$ injected at the start would persist forever in the diagonal elements.
* After 100 steps, the reservoir would be "full" or saturated with old information, making it impossible to distinguish the new input $u_{100}$.

**By forcing $Z$ to decay:**
We implement the **Fading Memory Property**. It acts like a leak in a bucket. The system holds the water (information) long enough to process it, but eventually, it drains away to make room for new water. This ensures the reservoir state depends only on the *recent* history of inputs.

---

## 5. Why do we want $X$ to stay? (The "Quantum Engine")

This is the subtle part. Why not destroy everything?

**If $X$ decayed too:**
* The term $\sum h_i X_i$ in the Hamiltonian would lose its effect.
* The system would effectively become a **Classical** Ising model (just $0$s and $1$s flipping).
* We would lose the quantum advantage (entanglement and superposition).

**By preserving $X$:**
We protect the "Quantum Driver" of the system.
* We allow the Hamiltonian to keep rotating the state around the X-axis.
* This rotation is what mixes the information non-linearly.
* **Result:** We have a system that "forgets" data ($Z$ decays) but remains dynamically active and quantum ($X$ survives).



> **Summary:**
> We kill $Z$ to clean the **Memory**.
> We save $X$ to keep the **Processor** running.
























---

## 3. General Bloch State Injection `inject_general_bloch_state`

Here is the comprehensive, step-by-step mathematical explanation of the `inject_state` function.

This function implements a quantum operation known physically as **"Reset and Write"**. It effectively erases the current state of a specific qubit and replaces it with a new state defined by a Bloch vector $\vec{r} = (r_x, r_y, r_z)$, while preserving the entanglement and correlations of the rest of the system.

-----

# Mathematical Foundation

To understand the code, we must look at the operation in terms of Density Matrices.

### The Objective

We want to transform the system state $\rho_{old}$ into $\rho_{new}$ by resetting qubit $k$:

$$
\rho_{new} = \underbrace{\text{Tr}_k(\rho_{old})}_{\text{Reset (Erase)}} \otimes \underbrace{\rho_{input}^{(k)}}_{\text{Write (Inject)}}
$$

1.  **Partial Trace ($\text{Tr}_k$):** Averages out qubit $k$, effectively deleting its information from the system.
2.  **The New State ($\rho_{input}$):** The new state for qubit $k$ is defined by the input coordinates:
    $$\rho_{input}^{(k)} = \frac{1}{2} (I + r_x X_k + r_y Y_k + r_z Z_k)$$

The function implements the distribution of this tensor product over the existing Pauli strings.

-----

# Code Walkthrough

### 1\. Function Signature & Optimization

```julia
function inject_state(rho::Operator, qubit_idx::Int, rz::Float64; rx::Float64=0.0, ry::Float64=0.0)
    # ... setup ...
    has_x = abs(rx) > tol
    has_y = abs(ry) > tol
    has_z = abs(rz) > tol
```

  * **Math:** The target vector is $\vec{r} = (r_x, r_y, r_z)$.
  * **Logic:** In Reservoir Computing, we typically only use the Z-axis ($r_z$). Calculating Pauli strings for $X$ and $Y$ is computationally expensive. If `rx` or `ry` are 0, the `has_x/y` flags allow the code to skip those calculations completely.

### 2\. The Filter (The "Erase" Step)

```julia
for (p, coeff) in rho
    # get_pauli_type: 0=I, 1=X, 2=Z, 3=Y
    if get_pauli_type(p, qubit_idx) == 0
```

**Mathematical Logic:**
We iterate through every Pauli string $P$ in the density matrix $\rho_{old} = \sum c_i P_i$.
When we calculate the partial trace $\text{Tr}_k(P)$ of a Pauli string:

  * If $P$ acts as $X, Y, \text{or } Z$ on qubit $k$, the trace is **0** (Pauli matrices are traceless). **These terms are discarded.**
  * If $P$ acts as Identity ($I$) on qubit $k$, the trace is **non-zero**. **These terms survive.**

The condition `if get_pauli_type(p, qubit_idx) == 0` implements this filter. It selects only the terms describing the "memory" of the *other* qubits ($P_{env} \otimes I_k$).

### 3\. Tensor Product Expansion (The "Write" Step)

Inside the loop, we have a surviving term with coefficient $c_p$. We now multiply it by the new state:

$$
c_p (P_{env} \otimes I_k) \times (1 \cdot I_k + r_z Z_k + r_x X_k + r_y Y_k)
$$

This expands into four potential new terms:

#### A. The Identity Term (Retention)

```julia
new_rho[p] = get(new_rho, p, 0.0im) + coeff
```

  * **Math:** $c_p (P_{env} \otimes I_k)$
  * **Meaning:** This preserves the existing correlations of the environment.

#### B. The Z Term (Input Injection)

```julia
if has_z
    p_z = PauliString(p.x_mask, p.z_mask | bit_loc)
    new_rho[p_z] = get(new_rho, p_z, 0.0im) + (coeff * rz)
end
```

  * **Math:** $c_p \cdot r_z (P_{env} \otimes Z_k)$
  * **Operation:** The code uses `| bit_loc` (bitwise OR) on the `z_mask` to turn the operator at index $k$ into a $Z$.
  * **Significance:** This couples the existing reservoir state ($P_{env}$) with the new input ($Z_k$). If $P_{env}$ was $Z_1$, the new term is $Z_1 Z_k$. This creates the **non-linear memory** required for computation.

#### C & D. The X and Y Terms (Optional Control)

```julia
if has_x
    p_x = PauliString(p.x_mask | bit_loc, p.z_mask)
    # ... adds (coeff * rx)
end
```

  * **Math:** $c_p \cdot r_x (P_{env} \otimes X_k)$
  * **Operation:** Sets the bit in `x_mask`. If `rx > 0`, this adds quantum coherence or noise control.

### 4\. Fresh Injection (The Vacuum State)

```julia
# 2. Inyectar el término fresco...
if has_z
    p_fresh_z = PauliString(0, bit_loc)
    new_rho[p_fresh_z] = get(new_rho, p_fresh_z, 0.0im) + rz
end
```

**Why is this outside the loop?**
The density matrix always implicitly contains the global identity term $I \otimes I \dots \otimes I$ (usually with coefficient 1.0). This term is not always stored explicitly in the sparse `rho` map to save memory.

We must explicitly calculate the product of the "Vacuum" Identity with the new state:

$$
1 \cdot (I_{env} \otimes I_k) \times r_z Z_k \rightarrow r_z (I_{env} \otimes Z_k)
$$

This adds the "pure" input state to the qubit, independent of previous correlations.

-----

# Summary of the Transformation

This function performs the following linear mapping on the coefficients:

1.  **Delete:** Remove any term where qubit $k$ is not $I$.
2.  **Keep:** Keep terms where qubit $k$ is $I$.
3.  **Clone & Modify:** Take the kept terms, clone them, attach $Z_k$ (or $X_k, Y_k$), and scale them by the input value $u$ (or $r_x, r_y$).

This is the mathematically correct way to **erase old data** and **write new data** into a specific qubit within a highly entangled quantum system.







# ⚛️  Dynamics & Time Integration `dynamics.jl`

This module  handles the **time evolution of a quantum system**. It solves the differential equations governing the system using the **4th-order Runge-Kutta method (RK4)**.

---

## 1. The Differential Equation (derivative)

The core of the simulation is solving the **equation of motion** for an operator $O$ (which represents the density matrix $\rho$ in the **Schrödinger picture** or an observable in the **Heisenberg picture**).

The code implements the **Heisenberg equation of motion** (or the **Liouville–von Neumann equation**):

$$
\frac{dO}{dt} = \mathcal{L}(O) = i[H, O]
$$

Where:

- $O$ is the operator being evolved (e.g., the density matrix $\rho$).
- $H$ is the **Hamiltonian** of the system.
- $[H, O] = HO - OH$ is the **commutator**.
- $i$ is the imaginary unit.

In the code, the function `derivative(O, H)` computes this rate of change. It iterates over the **Pauli strings** of $O$ and $H$, calculates their commutator using symplectic algebra, and returns the gradient:

$$
\frac{dO}{dt}
$$

---

## 2. Numerical Integration: Runge-Kutta 4 (step_rk4)

Since the differential equation cannot always be solved analytically for complex many-body systems, we use **numerical integration** to advance the system from time $t$ to $t + dt$.

The **classic Runge-Kutta method (RK4)** is used. It is a **fourth-order iterative method**, meaning the local error is of order $O(dt^5)$, providing much higher accuracy than the standard Euler method ($O(dt^2)$).

---

### Mathematical Formulation

To solve 

$$
\frac{dO}{dt} = f(t, O)
$$

for a time step $dt$, RK4 computes four "slopes" ($k_1$ through $k_4$) and takes a weighted average to estimate the next state.

1. **Initial slope ($k_1$):**

$$
k_1 = f(t, O_t)
$$

2. **Midpoint slope A ($k_2$):**

$$
k_2 = f\Big(t + \frac{dt}{2}, \, O_t + \frac{dt}{2} k_1 \Big)
$$

3. **Midpoint slope B ($k_3$):**

$$
k_3 = f\Big(t + \frac{dt}{2}, \, O_t + \frac{dt}{2} k_2 \Big)
$$

4. **End slope ($k_4$):**

$$
k_4 = f(t + dt, \, O_t + dt \cdot k_3)
$$

---

### Update Step

The final state $O_{t+dt}$ is calculated as a weighted average of these slopes:

$$
O_{t+dt} = O_t + \frac{dt}{6} \, (k_1 + 2k_2 + 2k_3 + k_4)
$$

In the code, the function `step_rk4` performs these exact operations using the `derivative` function as $f(O)$ and `add_ops` to perform intermediate vector additions.

---

### Why RK4?

- **Stability:** RK4 handles oscillatory quantum dynamics (unitary evolution) much better than lower-order methods, conserving probability (trace) and energy more accurately over long simulations.  
- **Efficiency:** It allows for larger time steps $dt$ while maintaining high precision, reducing the total number of iterations required.





# ⚛️  Pauli Algebra Utilities (`src/utils/pauli_algebra.jl`)

This module provides **fundamental linear algebra operations** on `PauliString` objects and quantum operators. It handles the low-level logic for multiplication, commutation, addition, and memory management (truncation) of operators.

-----

## Constants

```julia
const TOLERANCE = 1e-15
```

  * **Description:** A numerical threshold used to discard very small coefficients during operator addition. This prevents the accumulation of "floating-point noise" (e.g., `1.0e-18`) which should theoretically be zero.

-----

## Functions

### 1\. Multiplying Pauli Strings

```julia
multiply_paulis(p1::PauliString, p2::PauliString) -> (phase, p_new)
```

**Description:**
Computes the product of two Pauli strings $P_1$ and $P_2$.

**Returns:**

  * `phase`: A complex scalar ($1, -1, i, -i$) representing the phase factor resulting from the Pauli multiplication rules.
  * `p_new`: A new `PauliString` object corresponding to the product of the two inputs.

**Implementation Details:**

1.  The `x_mask` and `z_mask` of the two strings are combined using bitwise XOR (`⊻`) to determine the resulting operator structure.
2.  The function iterates over the qubits where both operators act non-trivially to compute the global phase based on the following multiplication table ($1=X, 2=Z, 3=Y$):

| Type 1 | Type 2 | Phase Update | Math Equivalent |
| :---: | :---: | :--- | :--- |
| 1 ($X$) | 2 ($Z$) | `phase *= -1.0im` | $XZ = -iY$ |
| 2 ($Z$) | 1 ($X$) | `phase *= 1.0im` | $ZX = iY$ |
| 2 ($Z$) | 3 ($Y$) | `phase *= -1.0im` | $ZY = -iX$ |
| 3 ($Y$) | 2 ($Z$) | `phase *= 1.0im` | $YZ = iX$ |
| 3 ($Y$) | 1 ($X$) | `phase *= -1.0im` | $YX = -iZ$ |
| 1 ($X$) | 3 ($Y$) | `phase *= 1.0im` | $XY = iZ$ |

-----

### 2\. Commutator of Pauli Strings

```julia
commutator(p1::PauliString, p2::PauliString) -> (val, p_new)
```

**Description:**
Computes the commutator defined as 

$[P_1, P_2] = P_1 P_2 - P_2 P_1$.

**Returns:**

  * If the strings **commute**: Returns `(0.0im, PauliString(0,0))`.
  * If the strings **anticommute**: Returns `(2 * phase, p_new)`, where `phase` and `p_new` are the result of the multiplication.

**Implementation Details:**
It efficiently checks whether the strings anticommute using symplectic bitwise logic. Two Pauli strings anticommute if the overlap of their X and Z masks has an odd parity:

```julia
isodd(count_ones((p1.x_mask & p2.z_mask) ⊻ (p1.z_mask & p2.x_mask)))
```

  * If **Odd**: They anticommute ($P_1 P_2 = - P_2 P_1$), so the commutator is $2 P_1 P_2$.
  * If **Even**: They commute ($P_1 P_2 = P_2 P_1$), so the commutator is $0$.

-----

### 3\. Operator Addition

```julia
add_ops(A::Operator, B::Operator, scale_B::ComplexF64)
```

**Description:**
Performs a linear combination of two operators, $A$ and $B$, scaling $B$ by a factor:

$ C = A + (\text{scale\_B} \cdot B) $

**Implementation Details:**

1.  Creates a copy of operator $A$ (called $C$).
2.  Iterates over every term in operator $B$.
3.  Adds the coefficient of $B$ (multiplied by `scale_B`) to the corresponding term in $C$.
4.  **Filtering:** If the resulting coefficient absolute value is smaller than `TOLERANCE`, the term is removed from the map to keep the operator sparse and efficient.

**Returns:**

  * `C`: The resulting summed operator.

-----

### 4\. Operator Truncation

```julia
truncate_operator!(O::Operator, max_terms::Int)
```

**Description:**
Reduces the size of operator `O` by keeping only the `max_terms` with the largest coefficients (in absolute value). This is an approximation method used to prevent exponential memory growth during time evolution.

**Implementation Details:**

1.  Converts the operator `O` into a list of `(PauliString, coefficient)` tuples.
2.  Sorts this list by coefficient magnitude (`abs(c)`) in descending order.
3.  Clears the original operator data.
4.  Re-inserts only the top `max_terms`.

**Note:** The `!` at the end of the function name indicates that this is an **in-place operation**, meaning it modifies the input operator `O` directly rather than returning a new one.

-----

## Summary

This module provides the essential mathematical engine for the simulation:

  * **Multiplication:** Handles quantum phase factors correctly.
  * **Commutation:** Determines dynamic evolution direction.
  * **Addition:** Handles linear combinations with noise filtering.
  * **Truncation:** Manages computational resources by approximating complex operators.
