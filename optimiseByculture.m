function [Zmin,Population] = optimiseByculture(Problem,AllP,W,t2,Zmin)
	alpha = 0.06;
	M = t2*alpha;
    NW = size(W,1);
    pAccept = 0.35;                   
    nAccept = round(pAccept*M);
    culture = [];
    [~,replace]=sort(pdist2(AllP.objs-repmat(Zmin,size(AllP,1),1),W,'cosine'));
    tempP = [];
    for i = 1 : NW
       P =  AllP(replace(1:M,i));

       [~,re] = sort(calPBI(P,Zmin,W(i,:)));
       spop = P(re(1:nAccept));
       cul.min = min(spop.decs,[],1);
       cul.max = max(spop.decs,[],1);
       culture = [culture,cul];
       tempP = [tempP,P(re(1))];
    end
    Population = tempP;    
    eva = Problem.FE;
    while Problem.FE < eva+t2
        for i = 1 : NW          
           dec = Population(i).dec;
           PopDec = unifrnd(culture(i).min,culture(i).max);
           decs = [dec;PopDec];
           Offspring = SOLUTION(OperatorGA(decs));    
           Zmin = min([Zmin;Offspring.objs],[],1);
           AllP = [Population(i),Offspring];
           [~,re] = sort(calPBI(AllP,Zmin,W(i,:)));
           Population(i) = AllP(re(1));
        end
    end

end

function PBI = calPBI(Population,Z,W)
    Z = repmat(Z,size(Population.objs,1),1);
    W = repmat(W,size(Population.objs,1),1);
    normW = sqrt(sum(W.^2,2));
    normP = sqrt(sum((Population.objs-Z).^2,2));
    CosineP = sum((Population.objs-Z).*W,2)./normW./normP;
    PBI = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
end


