# Coupled Oscillator Model parameter optimization with UKF on fMRI data

This project aims to leverage nonlinear system identification techniques, focusing on the Unscented Kalman Filter (UKF), to enhance modeling and understanding of neural dynamics from fMRI data, 
addressing the direct challenges of noise and indirect measurement. This experiment focuses on:

- UKF on different fMRI signals (different regions, different subjects)
- Coupled oscillator model modification - New model
- Optimization methods comparison (Iterative, Simulated Annealing, Nelder-Mead, L-BFGS-B)

## Unscented Kalman Filter - UKF

The Unscented Kalman Filter represents an advancement in filtering technology, adept at handling nonlinearities in data without necessitating linear approximations, making it particularly suited for the complex and noisy data characteristic of fMRI. Retains the exact nonlinearity of F but approximates the a posteriori probability density of the state $xt+∆t$ by a Gaussian.

$$xt+∆t ≈ F(xt) + (∇′Pxx,t∇/2) F(xt),$$

$$Pxx,t+∆t ≈ ∇F(xt) Pxx,t(∇F(xt))′$$

More accurate for nonlinear systems, and computationally much simpler to implement.

## New Coupled Oscillator Model - Two Coupled Pendulums

The equation of motion of the combined system is then given by: 
$$L \ddot{\theta_{1}} =-g \sin \theta_{1}-k L\left(\sin \theta_{1}-\sin \theta_{2}\right),$$
$$L \ddot{\theta_{2}} =-g \sin \theta_{2}+k L\left(\sin \theta_{1}-\sin \theta_{2}\right)$$

There are 3 parameters in the equation:
- g, gravity
- L, length of pendulums ($L_1$ and $L_2$)
- k, spring constant

without the assumption of pendulums having same length, $L_1$ and $L_2$ would replace the L in corresponding equation:

$$L_1 \ddot{\theta_{1}} =-g \sin \theta_{1}-k L_1\left(\sin \theta_{1}-\sin \theta_{2}\right),$$

$$L_2 \ddot{\theta_{2}} =-g \sin \theta_{2}+k L_2\left(\sin \theta_{1}-\sin \theta_{2}\right)$$

![image](https://github.com/skaraoglu/UKF-CoupledOscillatorModel-fMRI/assets/32866050/36c78a07-3e31-4651-8d0e-da8e85da073b)


## fMRI Data

- Voxel A & B data
- Harvard-Oxford atlas dataset
