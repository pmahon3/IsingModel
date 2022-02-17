%% Assignment 3

% Compute the exact free energy per spin 
% Source: en.wikipedia.org/wiki/Square-lattice_ising_model
xDat = linspace(0.01, 5, 100);
yDat = zeros(100);

for i = 1:100
    
    xVal = 1/xDat(i);
    k = 1 / (sinh(2*xVal))^2;
    fInt = @(theta) 1 ./ sqrt( 1 - 4 .* k .* (( 1 + k ).^-2) * sin(theta).^2);
    intVal = integral(fInt, 0, pi/2);
    yDat(i) = -1 * coth(2*xVal)*(1 + (2/pi) * (2 * tanh(2 * xVal)^2 - 1) * intVal );
end

csvwrite("Energy.csv", yDat);
T_c = 2/(log(1+sqrt(2)));
% Compute the Osager-Yang solution
xDat2 = linspace(0.01, T_c, 100);
yDat2 =  (1 - 1./(sinh(2./xDat2)).^4).^(1/8);
csvwrite("Mag.csv", [xDat2', yDat2']);

% Compute the critical temperature

