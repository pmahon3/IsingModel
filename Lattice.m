%%Llattice implements functions which operate on a 2D Ising model lattice

classdef Lattice < handle
     
    % lattice is a matrix representing the lattice object
    % nRows is the number of rows in the lattice
    % nCols is the number of columns in the lattice
    % E is the total Hamiltonian of the lattice
    % M is the total magnetization of the lattice
    % J is the coupling parameter
    % T is the temperature
    properties
        lattice;
        nRows;
        nCols;
        E;
        M;
        J;
        T;
        kb = 1;
    end
    
    methods
        % Lattice is the contructor for the class
        % if initialLattice is false than the spins are randomly
        % initialized. Otherwise a matrix representing the lattice of
        % spins should be passed and set as the lattice. 
        function obj = Lattice( nRows, nCols, J, T, initialLattice)
            
            % Initialize the parameters
            obj.T = T;
            obj.nRows = nRows;
            obj.nCols = nCols;
            obj.J = J;
            obj.lattice = ones(nRows, nCols);
            
            % Randomly initialize the lattice if no initial lattice is
            % passed
            if initialLattice == false
                for i = 1: nRows
                    for j = 1:nCols
                        if rand() < 0.5
                            obj.lattice(i,j) = -1;
                        end
                    end
                end
                obj.totalHamiltonian();
                obj.totalMagnetization();
            % Otherwise set the lattice to the passed lattice
            else 
                obj.lattice = initialLattice;
            end
        end
        
        % totalHamiltonian compute the total Hamiltonian of the lattice,
        %   summing row-wise and then column-wise, with periodic boundaries
        function totalHamiltonian( obj )
            obj.E = 0;
            
            % Perform the sum for all row-wise nearest neighbours
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    
                    if j == obj.nCols
                        obj.E = obj.E + obj.lattice(i,j) * obj.lattice(i,1);
                    else
                        obj.E = obj.E + obj.lattice(i, j) * obj.lattice(i, j+1);   
                    end
                end
            end
            
            % Perform the sum for all column-wise nearest neighbours
            for j = 1:obj.nCols
                for i = 1:obj.nRows
                    if i == obj.nRows
                        obj.E = obj.E + obj.lattice(i,j) * obj.lattice(1, j);
                    else
                        obj.E = obj.E + obj.lattice(i, j) * obj.lattice(i+1, j);
                    end
                end
            end
            
            % Multiply the sum by the coupling constant
            obj.E = (obj.E * obj.J);
        end
        
        % totalMagnetization sums all the spins in the lattice.
        function totalMagnetization( obj )
            obj.M = 0;
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    obj.M = obj.M + obj.lattice(i,j);
                end
            end
        end
        
        % dE computes the energy difference between the current lattice and
        % a lattice with spin [i,j] flipped. This method is modified from:
        % 
        %   MathWorks Physics Team (2020). Ising and Metropolis Algorithm
        %   (https://www.mathworks.com/matlabcentral/fileexchange/62194-
        %    ising-model-and-metropolis-algorithm), MATLAB Central File
        %    Exchange. Retrieved March 9, 2020.
        function dE = dESpin(obj, i,j)
            
            % Find index values of spins above and beside the selected spin
            above = mod( i - 2 , obj.nRows ) + 1;
            below = mod( i, obj.nRows ) + 1; 
            left = mod( j - 2, obj.nCols ) + 1; 
            right = mod( j, obj.nCols ) + 1; 
                 
            % Take the sum of the neighbouring spins
            neighbors = [   obj.lattice( above, j);
                            obj.lattice( i, left );
                            obj.lattice( i, right); 
                            obj.lattice( below, j);
                        ];
                    
            % Compute energy difference if the spin was flipped.
            dE = 2 * obj.J * obj.lattice(i,j) * sum(neighbors);
            
        end
        
        % newLattice selects a spin to flip for this lattice and, if the
        % transition function dictates, proposes this altered lattice as a
        % new lattice for a Monte Carlo step. 
        function newLattice = stateTransition( obj ) 
            
            % Pick random spin
            i = ceil(rand()*obj.nRows);
            j = ceil(rand()*obj.nCols);
            
            % Determine energy change of flipping spin
            dE = obj.dESpin( i, j );
            
            % Calculate probability that this energy change is accepted
            TProb = exp(-dE/(obj.kb*obj.T))/(1+exp(-dE/(obj.kb*obj.T)));
           
            % If the probability is less than X = Uniform(x) than accept
            % the change, otherwise no change. 
            if (TProb < rand() ) 
                newLattice = obj;
                newLattice.lattice(i,j) = -1 * newLattice.lattice(i,j);
                newLattice.E = newLattice.E + dE;
                newLattice.M = newLattice.M + 2 * newLattice.lattice(i,j);
            else
                newLattice = obj;
            end
        end
    end
end