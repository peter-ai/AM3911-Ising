%{
	Author: Peter Akioyamen
	Student #: 250949002
	Course: AM3911G - Modelling & Simulation
	Professor: Allan MacIsaac
	Assignment #: 3
	Due: March 4, 2020

    A script which plots the average energy per 
    spin of a 2D Ising model (75x75 lattice).
%}

syms AE(T) T B O F k Q1 Q2

boltzmann = 1;
J = 1;

ising_data = readtable("../data/ising_model_2d_75by75_results.csv");
t = ising_data{:, 1};
E = ising_data{:, 2};

N = 50^2;
B = 1 / (boltzmann * T);
k = (2 * (sinh(2 * (B) * J))) / (cosh(2 * (B) * J))^2;
F = 1 / (sqrt(1 - ((k^2) * (sin(O)^2))));
K_1 = int(F, O, 0, pi / 2);
Q1 = -2 * N * J * tanh(2 * B * J);
Q2 = -1 * (N * J * (((sinh(2 * B *J))^(2) - 1) / (sinh(2 * B * J) * cosh(2 * B * J))));
AE(T) = (Q1 + Q2 * ((2/pi) * K_1 - 1)) / N;

scatter(t, E, 24, "b")
hold on 
fplot(AE(T), [1, 4])
legend("MCMC", "Analytical Solution", "interpreter", "latex")
title("2D Ising Model - Exact Solution vs MCMC", "interpreter", "latex")
xlabel("$$\textit{T}$$", "interpreter", "latex")
ylabel("$$<\textit{E}>$$", "interpreter", "latex")
xlim([1, 4])




