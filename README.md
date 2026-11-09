# Stability of patterns in variational systems with continuous symmetries: a computer-assisted approach



Table of contents:


* [Introduction](#introduction)
* [The Problem](#the-problem)
   * [Augmenting the System](#augmenting-the-system)
* [Proof of Periodic Solutions](#proof-of-periodic-solutions)
* [Controlling the spectrum](#controlling-the-spectrum)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](To appear) : "Stability of patterns in variational systems with continuous symmetries: a computer-assisted approach", J.B. van den Berg and D. Blanco [ArXiv Link](To appear).

It provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl). 


# The Problem

The problem we consider is

$$(I_d + \Delta)^2u + \mu u - \gamma u^2 + u^3 = 0$$

The problem is to be solved for $u= u(x), ~x \in (-d_1,d_1) \times (-d_2,d_2)$ with periodic boundary conditions. These conditions are enforced by using Fourier series. We assume that $u$ has a continuous symmetry, and hence, we cannot treat the problem directly with usual computer-assisted techniques.

## Augmenting the system

In Section 2.2 in [[1]](To appear), we write the augmented system

$$ \overline{F}(\beta_1,\beta_2,u) = \begin{bmatrix}
(u,\partial_1 u_0) \\
(u,\partial_2 u_0) \\
F(u) + \beta_1 \partial_1 u_0 + \beta_2 \partial_2 u_0
\end{bmatrix}$$

where $(\beta_1,\beta_2)$ are unfolding parameters. 

# Proof of Periodic Solutions 
We show that the unfolding parameters of $\overline{F}$ are $0$ when we are at a solution, which means the zeros of $\overline{F}$ are the same as those of $F$. Using $\overline{F}$, we can perform computer-assisted proofs using the analysis done in Section 4 of [[1]] (to appear).

Given data of approximate solution, continuous_sym.jl provides the explicit computation of the bounds in Section 4. It also provides a value for $r^{*}$ where the proof is successful. 

# Controlling the spectrum

In Section 4.3 in [[1]](To appear), we derive an additional bound $Z_Q$ which allows us to control the spectrum of $DF(\tilde{u})$. In particular, we want to understand the (in)stability of $\tilde{u}$, but we must do so by studying $D\overline{F}(0,0,\tilde{u})$. To simplify this, we define a different operator $Q(0,0,\tilde{u})$ which is easier to study. The bound $Z_Q$ is essentially the error for studying $Q$ instead of $D\overline{F}$. 

The end of the code continuous_sym.jl provides us with a Gershgorin circle argument to obtain enclosures of the eigenvalues of $D\overline{F}(0,0,\tilde{u})$. With the theory we developed, we can then conclude about the (in)stability. Additionally, we provide an explicit lower bound for the spectral gap in Theorem 3.8. This computation is provided in continuous_sym.jl and is given by 
$$ |\lambda| \geq \frac{1-\delta}{(1-\delta)\|\overline{A}\|_{\ell^2} + \delta \|\overline{A}\|_{\xi}}$$

 
 # Utilisation and References

 The code in continuous_sym.jl can serve to prove other solutions than the ones provided as illustration should one have the numerical candidates. In particular, one can just provide the numerical candidate and attempt a proof.

 To choose which of Theorems 5.1-5.4 to prove, set sol on line 128 of continuous_sym.jl to 1, 2, 3 or 4. The figures of the paper are made by plot_solutions.m (MATLAB) from the same data.
 
 The code is build using the following packages :
 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) (version 0.9)
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [DelimitedFiles](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/)

 The code requires version 0.9 of RadiiPolynomial; it does not work with later versions. To install version 0.9, run using Pkg; Pkg.add(name="RadiiPolynomial", version="0.9") in Julia.

 
 
 # License and Citation
 
This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{continuous_sym.jl,
  author = {Dominic Blanco and J.B. van den Berg},
  title  = {continuous_sym.jl},
  url    = {https://github.com/dominicblanco/continuous_sym.jl},
  note = {\url{https://github.com/dominicblanco/continuous_sym.jl}},
  year   = {2026},
  doi = {10.5281/zenodo.24090635}
}
```
DOI : [10.5281/zenodo.24090635](https://doi.org/10.5281/zenodo.24090635) 


# Contact

You can contact us at :

dominic.blanco@rutgers.edu
janbouwe.vanden.berg@vu.nl
