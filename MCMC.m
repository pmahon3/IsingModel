% MC implements monte carlo functions for an ising model simulation

classdef MCMC < handle
    
    % ATTRIBUTES
    % currentLattice is the current state of the Markov chain
    % previousLattice is the previous state of the Markov chain
    % chainsE records the total Hamiltonian for each member of the chain
    % chainsM records the total Magnetization for each member of the chain
    % nChains is the number of chains
    % nIterations is the number of non-burnin steps
    % nThins is interval between succesive MC steps to accept a new state
    % nBurnIn is the number of initial MC steps to reject
    % size is the number of spins in the lattice ( dimensionality of the
    % space)
    % nSamples is the number of accepted MC steps
    % chainsLat stores each lattice accepted to the chains
    % E is the mean total Hamiltonian of all chains
    % M is the mean absolute magnitization of all chains
    properties
     previousLattice;
        currentLattice; 
        chainsE;
        chainsM;
        nChains;
        nIterations;
        nThin; 
        nBurnIn; 
        size;
        nSamples;
        chainsLat;
        
        E;
        M;
    end
    
    methods
        
        % The constructor for the class
        function obj = MCMC( lattice, nChains, nIterations, nThin, nBurnin)
            
            % Initialize parameters
            obj.previousLattice = lattice;
            obj.currentLattice = lattice; 
            obj.nChains = nChains;
            obj.nIterations = nIterations;
            obj.nSamples = nIterations / nThin;
            obj.nThin = nThin;
            obj.nBurnIn = nBurnin;
            obj.size = lattice.nRows * lattice.nCols;
            obj.chainsE = zeros( obj.nSamples, obj.nChains);
            obj.chainsM = zeros( obj.nSamples, obj.nChains);
            obj.M = 0;
            obj.E = 0;
        end
        
        % runChains runs the MCMC simulation for each chain, taking the
        % last lattice of each chain as the initial lattice of the new
        % chain.
        function runChains( obj )
            
            % Burn in phase
            for n = 1:obj.nBurnIn
                for k = 1:obj.size
                         newLattice = obj.currentLattice.stateTransition();
                         obj.currentLattice = newLattice;
                end
            end
            
            obj.previousLattice = obj.currentLattice;
            
            % Sampling phase
            for j = 1:obj.nChains
                obj.chainsE(1, j) = obj.currentLattice.E;
                obj.chainsM(1, j) = obj.currentLattice.M;
                for i = 2:obj.nIterations
                    for k = 1:obj.size
                         newLattice = obj.currentLattice.stateTransition();
                         obj.previousLattice = obj.currentLattice;
                         obj.currentLattice = newLattice;
                    end
                    
                    if ( mod(i, obj.nThin ) == 0 )
                        obj.chainsE(i/obj.nThin, j) = obj.currentLattice.E;
                        obj.chainsM(i/obj.nThin, j) = obj.currentLattice.M;
                        %image(obj.currentLattice.lattice, 'CDataMapping', 'scaled');
                    end
                end
            end
            
            % Compute the mean total Hamiltonian and mean absolute
            % magnetization for all chains
            obj.E = mean(mean(obj.chainsE));
            obj.M = mean(mean(abs(obj.chainsM)));
        end
        
        % Plotting functions for verification
        function plotChainsE(obj)
            for j = 1:obj.nChains
                plot(obj.chainsE(:, j));
                hold on;
            end
        end
        
        function plotChainsM(obj)
            for j = 1:obj.nChains
                plot(abs(obj.chainsM(:, j)));
                hold on;
            end
        end
    end
end
    