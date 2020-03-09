% MC implements monte carlo functions for an ising model simulation

classdef MCMC < handle
    
    % ATTRIBUTES
    % lttc is a lattice object representing the current state of the chain.
    % tmp is the tempurature of the lattice. 
    % kb is the Boltzmann constant
    
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
        
        function obj = MCMC( lattice, nChains, nIterations, nThin, nBurnin)
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
                        image(obj.currentLattice.lattice, 'CDataMapping', 'scaled')
                    end
                end
            end
            
            obj.E = mean(mean(obj.chainsE));
            obj.M = mean(mean(abs(obj.chainsM)));
        end
        
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
    