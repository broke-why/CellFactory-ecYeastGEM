function [bioYield,prod_yield, prod_rate] = calculate_potential(model,growthPos,targetIndex,CS_index,CS_MW)
%Unconstrain growth
tempModel = setParam(model,'lb',growthPos,0);
tempModel = setParam(tempModel,'ub',growthPos,1000);
%fix unit glucose uptake rate
tempModel = setParam(tempModel,'obj',growthPos,1);
tempModel = setParam(tempModel,'ub',CS_index,1);
solution  = solveLP(tempModel);
bioYield  = solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
disp(['The maximum biomass yield is ' num2str(bioYield) '[g biomass/g carbon source]']);

%fix suboptimal growth
tempModel = setParam(tempModel,'lb',growthPos,0.1*solution.x(growthPos));
%Get product yield for a unit glucose uptake rate
tempModel = setParam(tempModel,'obj',targetIndex,1);
solution = solveLP(tempModel);
prod_yield  = solution.x(targetIndex)/(solution.x(CS_index));
disp(['The maximum product yield is ' num2str(prod_yield) '[mol biomass/mol carbon source]']);

%Get product max. rate for a unconstrained glucose uptake rate
tempModel = setParam(tempModel,'ub',CS_index,1000);
solution  = solveLP(tempModel);
prod_rate  = solution.x(targetIndex);
disp(['The maximum productio rate is ' num2str(prod_rate) '[mmol/gDW h]']);
end