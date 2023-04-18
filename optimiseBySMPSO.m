function [Zmin,P] = optimiseBySMPSO(Problem,inputPopulation,evaluations,Zmin)

    %% Generate random population
    Population       = inputPopulation;
    Pbest            = Population;
    [Gbest,CrowdDis] = WOF_UpdateGbest(Population,Problem.N);
     P = Population;
    maximum = Problem.FE + evaluations;

    %% Optimization
    while Problem.FE < maximum
        Population       = SMPSO_operator(Problem, [Population,Pbest,Gbest(TournamentSelection(2,Problem.N,-CrowdDis))]);
        P = [P,Population];
        Zmin       = min([Zmin;Population.objs],[],1);
        [Gbest,CrowdDis] = WOF_UpdateGbest([Gbest,Population],Problem.N);
        Pbest            = WOF_UpdatePbest(Pbest,Population);
    end

end

function Pbest = WOF_UpdatePbest(Pbest,Population)
    % Update the local best position of each particle
    replace        = ~all(Population.objs>=Pbest.objs,2);
    Pbest(replace) = Population(replace);
end

function [Gbest,CrowdDis] = WOF_UpdateGbest(Gbest,N)
    % Update the global best set
    Gbest    = Gbest(NDSort(Gbest.objs,1)==1);
    CrowdDis = WOFCrowdingDistance(Gbest.objs);
    [~,rank] = sort(CrowdDis,'descend');
    Gbest    = Gbest(rank(1:min(N,length(Gbest))));
    CrowdDis = CrowdDis(rank(1:min(N,length(Gbest))));
end

function CrowdDis = WOFCrowdingDistance(PopObj)
    % Calculate the crowding distance of each solution in the same front
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end