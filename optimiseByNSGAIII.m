function [Zmin,P] = optimiseByNSGAIII(Problem,Population,W,evaluations,Zmin)
% ----------------------------------------------------------------------- 
    maximum = Problem.FE + evaluations;
    P = Population;
    
    %% Optimization
    while Problem.FE < maximum
        MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
        Offspring   = OperatorGA(Population(MatingPool));
        P = [P,Offspring];
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        Population = WOF_NSGAIIIEnvironmentalSelection([Population,Offspring],Problem.N,W,Zmin);
    end
end