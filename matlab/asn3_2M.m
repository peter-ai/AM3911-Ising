%{
	Author: Peter Akioyamen
	Student #: 250949002
	Course: AM3911G - Modelling & Simulation
	Professor: Allan MacIsaac
	Assignment #: 3
	Due: March 4, 2020

    A script which plots the average absolute
    magnetism of a 2D Ising model (75x75 lattice).
%}

boltzmann = 1;
J = 1;
T1 = linspace(1, 2 / log(1+sqrt(2)), 100);
beta = 1 ./ (boltzmann .*  T1);
m1 = (1 - (sinh(2 .* beta .* J)).^(-4)).^(1/8);

T2 = linspace(2 / log(1+sqrt(2)), 4.0, 100);
m2 = zeros(100);

T = zeros(200);
m = zeros(200);
T(1:100) = T1;
T(101:200) = T2;
m(1:100) = m1;
m(101:200) = m2(1:100);

ising_data = readtable("../data/ising_model_2d_75by75_results.csv");
T3 = ising_data{:, 1};
m3 = ising_data{:, 3};

scatter(T3, m3, 24, "b");
%hold on
%plot(T2, m2, "r", "LineWidth", 1.15)
hold on
plot(T, m, "r", "LineWidth", 1.15);

legend("MCMC", "Analytical Solution", "interpreter", "latex")
title("2D Ising Model - Exact Solution vs MCMC", "interpreter", "latex")
xlabel("$$\textit{T}$$", "interpreter", "latex")
ylabel("$$<|\textit{M}|>$$", "interpreter", "latex")
xlim([1, 4])
