function [BioYield,yield,oxygen,EtOH] = getYieldPlot(model,target,GUR,MW)
if nargin<4
    MW = 1;
end
ecM = any(contains(model.rxnNames,'prot_pool'));
oxIndex = find(contains(model.rxnNames,'oxygen exchange (rev'));
EtOHIdx = find(strcmpi(model.rxnNames,'ethanol exchange'));

growthIndex = find(strcmpi(model.rxnNames,'growth'));
BioYield    = [];
yield       = [];
oxygen      = [];
EtOH        = [];
iterations  = 10;
if ~ecM
    glucIndex = find(strcmpi(model.rxnNames,'D-glucose exchange'));
    model = setParam(model,'lb',glucIndex,-1.000001*GUR);
    model = setParam(model,'ub',glucIndex,-0.999999*GUR);
else
    glucIndex = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
    model     = setParam(model,'ub',glucIndex,1.000001*GUR);
    model     = setParam(model,'lb',glucIndex,0.999999*GUR);
end 
model  = setParam(model,'ub',model.rxns(growthIndex),1000);
model  = setParam(model,'lb',model.rxns(growthIndex),0);
model  = setParam(model,'ub',model.rxns(target),1000);
model  = setParam(model,'obj',model.rxns(growthIndex),1);
sol    = solveLP(model);
MiuMax = sol.x(growthIndex);
for i=1:iterations+1
    % Set the objective function to the target exchange reaction
    tempModel = setParam(model,'obj',model.rxns(target),1);
    % Fix dilution rate at every iteration
    Drate     = MiuMax*(i-1)/iterations;
    tempModel = setParam(tempModel,'lb',tempModel.rxns(growthIndex),0.9999*Drate);
    solution  = solveLP(tempModel,1);
    GURsim    = solution.x(glucIndex);
    if ~isempty(solution.f)
        production  = solution.x(target);
        %GUR       = abs(solution.x(glucIndex));
        % Calculate biomass yield (g biomass/g glucose)
        value    = Drate/abs(GURsim*0.180156);
        BioYield = [BioYield; value];
        % Calculate 3HP yield [mmol product/mmol glucose]
        value = production*MW/abs(GURsim*0.180156);
        yield = [yield; value];
        oxygen = [oxygen;solution.x(oxIndex)];
        EtOH   = [EtOH;solution.x(EtOHIdx)];
    end
   
end

end