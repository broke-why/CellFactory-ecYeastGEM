function candidates = getECCs(candidates,model,target)
model.ub(target) = 1000;
base_sol = solveLP(model,1);
current = pwd;
tolerance = 1E-6;
cd GECKO/geckomat/utilities
candidates.ECCs = zeros(height(candidates),1);
perturbation = 1.1;
for i=1:height(candidates)
    enzyme = candidates.enzymes{i};
    pEusage = candidates.pUsage(i);
    temp_model = model;
    if ~isempty(enzyme)
        option = 2;
        [~,rxnIdx] = getKcat(model,enzyme);
        enzPos   = find(strcmpi(model.mets,['prot_' enzyme]));
        enzUsRxn = find(strcmpi(model.rxns,['draw_prot_' enzyme]));
        temp_model.lb(enzUsRxn) = 0;
        temp_model.ub(enzUsRxn) = 1000;
        kcoeff   = model.S(enzPos,rxnIdx);
        if option == 1
            temp_model = setParam(temp_model,'lb',enzUsRxn,(1-tolerance)*perturbation*pEusage);
            temp_model = setParam(temp_model,'ub',enzUsRxn,(1+tolerance)*perturbation*pEusage);
            new_sol    = solveLP(temp_model);
            factor = 1;
            if isempty(new_sol.f)
                option = 2;
            end
        end
        if option == 2
            temp_model.S(enzPos,rxnIdx) = kcoeff./perturbation;
            factor = perturbation;
        end
        a_i = mean(kcoeff*pEusage);
        new_sol = solveLP(temp_model);
        if ~isempty(new_sol.f)
            a_i_new = kcoeff*factor.*new_sol.x(enzUsRxn);
            delta_V  = (new_sol.x(target) - base_sol.x(target))/base_sol.x(target);
            delta_Ea = (a_i_new-a_i)/a_i;
            %if (a_i_new*delta_Ea*delta_V)~=0 %&
            if delta_Ea>0
                %If solution was feasible then calculate the control coefficient
                objCC = delta_V./delta_Ea;
                if candidates.actions(i)>0
                    objCC = max(objCC);
                else
                    objCC = min(objCC);
                end
                candidates.ECCs(i) = objCC;
            end
        end
    end      
end
cd (current)
end