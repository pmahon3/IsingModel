%% Assignment 3
% Driver using the MCMC and Lattice objects to simulate a 30x30 Ising model
% at different temperatures and plots the results.
close all;
clear all;

% Simulation Parameters
nRows = 30;
nCols = 30;
nSims = 50;
nChains = 3;
nThin = 10;
nBurnin = 2000;
nIterations = 3000;
J = -1;

% Array initialization for results
Ts = linspace(0.1, 5, nSims);
EResults = zeros(1, nSims);
MResults = zeros(1, nSims);

% Run the simulations
for i = 1:nSims
    lat = Lattice(nRows,nCols, J, Ts(i), false);
    MC = MCMC( lat, 1, nIterations, nThin, nBurnin);
    MC.runChains();
    EResults(i) = MC.E;
    MResults(i) = MC.M;
end

% Plot the results 

% Average energy per spin
plot(Ts, EResults./-(nRows*nCols), 'kx');
xlabel("T");
ylabel("<E>/N");
title("System Energy per Spin");
legend("Simulation");

% Average magnitization per spin
figure;
plot(Ts, MResults./(nRows*nCols), 'kx');

% Compute the Osager-Yang solution
f = @(x) ( 1 - 1./(sinh(2./x)).^2).^(1/8);
hold on;
fplot(f, [0 5]);

% Compute the critical temperature
T_c = 2/(log(1+sqrt(2)));
xline(T_c);
title("Magentization per Spin");
xlabel("T");
ylabel("<|M|>");
legend("Simulation", "Osager-Yang Solution", "T_c");