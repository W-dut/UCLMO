classdef UCLMO < ALGORITHM
    % <many> <real> <large/none>
    % optimizer1 ---   1 --- Internal optimisation algorithm. 1 = NSGA-III, 2 = SMPSO
    
    methods
        function main(Algorithm,Problem)
            %Initializing the population
            Population = Problem.Initialization();
			%Initializing the parameter (The maximum number of evaluations performed by the two optimizers)
			[t1,t2,optimiser1] = Algorithm.ParameterSet(500,500,1);
            %Generate reference vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);

            Zmin = min(Population.objs);
            
            while Algorithm.NotTerminated(Population)
                %Optimizer 1, using NSGA-III OR SMPSO as the optimizer, which can be replaced by other optimizers;
				if      optimiser1 == 1
					[Zmin,AllP] = optimiseByNSGAIII(Problem,Population,W,t1,Zmin);   % NSGA-III
				else
					[Zmin,AllP] = optimiseBySMPSO(Problem,Population,t1,Zmin);       % SMPSO
                end
				
                %Optimizer 2, based on a culture-assisted evolutionary algorithm, passed in all individuals in the optimization process AllP;
                [Zmin,Population] = optimiseByculture(Problem,AllP,W,t2,Zmin);
            end
            
        end
    end
end
