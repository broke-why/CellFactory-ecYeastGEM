function CCMratio = get_CCM_enzBurden(fluxDist,model)
fluxDist = fluxDist(contains(fluxDist.rxns,'prot_'),:);
%metPathway in CCM
pathways = {'Oxidative phosphorylation' ...
    'Glycolysis / Gluconeogenesis' ...
    'Citrate cycle (TCA cycle)' ...
    'Pentose phosphate pathway' ...
    'Galactose metabolism' ...
    'Pyruvate metabolism'};


enzymes  = model.enzymes;
enzUsage = [];
CCM = [];
for i=1:min(length(enzymes),length(model.pathways))
    index = find(contains(fluxDist.rxns,['prot_' enzymes{i}]));
    enzUsage = [enzUsage;fluxDist.flux(index)*model.MWs(i)];
    presence = false;
    j= 1;
    while presence == false & j <= length(pathways)
        presence = (contains(model.pathways(i),pathways{j}));
        j= j+1;
    end
    if presence > 0
        CCM = [CCM; fluxDist.flux(index)*model.MWs(i)];
    end       
end
enzUsage = sum(enzUsage);
CCM      = sum(CCM);
CCMratio = CCM/enzUsage;
end



