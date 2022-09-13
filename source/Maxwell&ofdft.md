# Things about TD-OF

![1](_static/Maxwell&ofdft.assets/1.png)

|                     | Interacting system               | KS system                                             | Noninteracting boson system     |
| ------------------- | -------------------------------- | ----------------------------------------------------- | ------------------------------- |
| Density             | n($\vec r,t$)                    | n($\vec r,t$)                                         | n($\vec r,t$)                   |
| Wave function       | $\Psi(\vec r_1,...,\vec r_N,t )$ | $\frac{1}{\sqrt{N!}}det[\{\phi_{s,l}(\vec r_l ,t)\}]$ | $\prod^n_l \phi_B(\vec r_l ,t)$ |
| Effective potential | $v(\vec r,t)$                    |                                                       |                                 |
| Hamiltonian         |                                  |                                                       |                                 |

**Pauli potential term** $v_{\mathrm{P}}[n](\vec r , t)$, required to compensate for the neglect of the pauli exclusion principle.

## Time-denpendent Schrodinger-like equation

- **Time dependent Schrodinger-like equation**:

$$
\left [-\frac{\nabla^2}{2} + v_B[n](\vec r,t) \right ]\phi_B(\vec r ,t)=i\frac{\partial}{\partial t}\phi_B (\vec r ,t)
$$

 with initial condition $\phi_B (\vec r ,t) = \phi_B^0(\vec r)$.

- **Time-dependent effective potential**:

$$
v_B (\vec r ,t) = v_P (\vec r ,t) + v_S (\vec r ,t) \label{effective potential}
$$
- **Pauli kinetic energy** :

$$
v_{\mathrm{P}}[n](\vec r ) = \frac{\delta T_p[n]}{\delta n(r)}
$$
where $T_p[n] = T_s[n]-T_s^{vW}[n]$

- **Initial orbital** $\phi_b^0$:
$$
\phi_b^0 = \frac{1}{\sqrt{N}} \sqrt{n_0(r)} 
$$
## Properties of $v_P (\vec r ,t)$

Equation $\ref{effective potential}$. indicates that the role of the time-denpendent Pauli Potential in TD-OFDFT is similar to the role of time-dependent **XC  potential** in TD-DFT because it ensures that the fictitious boson system has the same electronic dynamics as the fermion system.

1. *Functional dependency*
$$
	v_S \equiv v_S[n_0,\Psi_0,\Phi_S^0](\vec r) \\
	v_B \equiv v_B[n_0,\Psi_0,\Phi_B^0](\vec r) \\
	v_P \equiv v_P[n_0,\Phi_S^0,\Phi_B^0](\vec r) 
$$

2. *The zero-force throrem*

	- total momentum of a many-body system:

		$$
		\vec P(t) = \int d\vec {r\ j}(\vec r,t)
		$$

	- the zero-force theorem

		$$
		\int d\vec r n(\vec r , t)\nabla v_p(\vec r,t) = 0
		$$

3. *The one- or two-electron limit*

4. *The relationship between $T_p\ and\ V_p$*
	$$
	E(t) = T_B(t) + T_P(t) + E_(H) + E_XC(t) + \int d \vec r\ v(\vec r,t)n(\vec r,t),  \\
	\text{where}
	\begin{cases}
	  & T_B(t)=-\sqrt{n(\vec r,t)} \frac{\nabla^2}{2} \sqrt{n(\vec r,t)}  \\
	  & T_P(t)=T_S(t)-T_B(t) \\
	  & E_H(t)=\frac{1}{2}\int \ d\vec r \ \int d \vec r' \frac{n(\vec r,t)n(\vec r',t)}{|\vec r -\vec r'|} \\
	  & E_{XC}(t)=T(t)-T_S(t)+W(t)-EH(t) \\
	\end{cases}
	$$
	
	and 
	
	$$
	\frac{d T_{\mathrm{P}}(t)}{d t}=\int d \vec r  \frac{\partial n(\vec r , t)}{\partial t} v_{\mathrm{P}}(\vec r , t)
	$$

5. *Nonadiabaticity and causality*



# Maxwell



