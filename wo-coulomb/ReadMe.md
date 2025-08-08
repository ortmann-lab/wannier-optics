# Calculation of transition matrix elements, Coulomb integrals and local field effects

Publication: [Merkel, K. & Ortmann, F. Linear scaling approach for optical excitations using maximally localized Wannier functions. Journal of Physics: Materials 7, 015001 (2023).]( https://dx.doi.org/10.1088/2515-7639/ad06cd )




## Coulomb integrals

`wo-coulomb.x` calculates the following integrals:
$$ W^{c_1 c_2}_{v_1 v_2}(R_c,R_v, R_D) = \int d^3 x  \int d^3 x' \rho_{c_1 c_2 R_c}(x) W(x-x'-R_D) \rho_{v_1 v_2 R_v}(x'), $$
where the overlap densities for conduction Wannier functions $\rho_{c_1 c_2R_c}$ and valence Wannier functions $\rho_{v_1 v_2R_v}$ are given by
$$\rho_{c_1 c_2R_c}(x) = w_{c_1 0}(x) w_{c_2 R_c}(x), $$
$$ \rho_{v_1 v_2 R_v}(x) = w_{v_1 0}(x) w_{v_2 R_v}(x). $$

Every Coulomb integral is therefore characterized by 4 indexes and 3 unit cell vectors
$(c_1, c_2, v_1, v_2, R_D, R_c, R_v)$.
This convention will be used in all input and output files.


### Connection to exciton Hamiltonian
`wo-optics.x` uses the Coulomb integrals to setup the exciton Hamiltonian. The basis of the relevant exciton subspace (see paper)
is given by three indexes $(c,v,S)$, where $c$ and $v$ are the indexes of conduction and valence Wannier functions as before and $S$
is the electron-hole distance in terms of lattice vectors.

By using the substitution rules $A = R_c$, $S=-R_D$, $S'=-R_D + R_c - R_v$, we can obtain the electron-hole interaction (in the exciton Hamiltonian),
$$ \tilde{H}^\text{SC}_{c_1 v_1 S,\, c_2 v_2 S'} = \sum_{A} \tilde{W}^{c_1 c_2}_{v_1, v_2}(A = R_c,S=-R_D, S'=-R_D+R_c-R_v). $$


## Local field effects
`wo-coulomb.x` can also be used to calculate local field effects,
$$ \tilde{H}^\text{LFE}_{c_1 v_1 S_1,\, c_2 v_2 S_2} = \int dx \int dx' \, w_{c_1 0}(x) w_{v_1,-S_1}(x) \left[\sum_{G\ne 0} \tilde{v}(|G|) e^{iG(x-x')} \right] w_{v_2 0}(x') w_{c_2,-S_2}(x'). $$

Local field effects only have two unit cell vectors. They are therefore described by the indexes $(c_1, c_2, v_1, v_2, S_1, S_2)$.
They can be used directly in the exciton Hamilonian.

To calculate local field effects within the `CoulombIntegral.x` code you only need to specify the
solver with `solver = LocalFieldEffects` in the configuration file. Input and output files
have the same format as for usual Coulomb integrals (see above). Only the naming convention is changed:

$R_D = 0$, $R_c = S_1$, $R_v = S_2$


## Transiton matrix elements (optical transition dipoles)
`wo-coulomb.x` calculates transition dipoles between electron and hole in the form:
$$
\tilde{M}^*_{c_1 v_1 S}(e_q)  = \sqrt{\frac{4\pi e^2}{\epsilon_0\epsilon_r}} e_q \int d^3x \, w_{c_1 0}(x) x w_{v_1,-S}(x)
$$

Please note that the prefactor might not be contained.