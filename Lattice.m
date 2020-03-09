%% lattice implements functions which operate on a 2D lattics

classdef Lattice < handle
     
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
        function obj = Lattice( nRows, nCols, J, T, initialLattice)
            obj.T = T;
            obj.nRows = nRows;
            obj.nCols = nCols;
            obj.J = J;
            obj.lattice = ones(nRows, nCols);
            
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
            else 
                obj.lattice = initialLattice;
            end
        end
        
        function totalHamiltonian( obj )
            obj.E = 0;
            
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    
                    if j == obj.nCols
                        obj.E = obj.E + obj.lattice(i,j) * obj.lattice(i,1);
                    else
                        obj.E = obj.E + obj.lattice(i, j) * obj.lattice(i, j+1);   
                    end
                end
            end
            
            for j = 1:obj.nCols
                for i = 1:obj.nRows
                    if i == obj.nRows
                        obj.E = obj.E + obj.lattice(i,j) * obj.lattice(1, j);
                    else
                        obj.E = obj.E + obj.lattice(i, j) * obj.lattice(i+1, j);
                    end
                end
            end
            
            obj.E = (obj.E * obj.J);
        end
        
        function totalMagnetization( obj )
            obj.M = 0;
            for i = 1:obj.nRows
                for j = 1:obj.nCols
                    obj.M = obj.M + obj.lattice(i,j);
                end
            end
        end
        
        function dE = dESpin(obj, i,j)
            
            above = mod( i - 2 , obj.nRows ) + 1;
            below = mod( i, obj.nRows ) + 1; 
            left = mod( j - 2, obj.nCols ) + 1; 
            right = mod( j, obj.nCols ) + 1; 
                 
            neighbors = [   obj.lattice( above, j);
                            obj.lattice( i, left );
                            obj.lattice( i, right); 
                            obj.lattice( below, j);
                        ];
            dE = 2 * obj.J * obj.lattice(i,j) * sum(neighbors);
            
        end
        
        function newLattice = stateTransition( obj ) 
            
            % Pick random spin
            i = ceil(rand()*obj.nRows);
            j = ceil(rand()*obj.nCols);
            
            % Determine energy chain of flipping spin
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