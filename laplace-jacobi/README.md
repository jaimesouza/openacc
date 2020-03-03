## Examplo 1 - A Laplace Solver Using Simple Jacobi Iteration

This program solves Laplace's equation on a regular grid using simple Jacobi iteration. Laplace's equation in two dimensions is:

![picture](figures/laplace-fig-1.png)

Lets assume the grid A. We can discretize Laplace equation as shown in the figure below.

![picture](figures/laplace-fig-2.png)

### How to compile

´´´ pgcc -fast -ta=tesla -Minfo=all <exemplo>.c ´´´


