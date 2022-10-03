function newModel = get_Kma_acummulation_model(model,accum_met,rxnStr)
%Add pseudo-accum,ulation reaction reaction:
model.c = zeros(length(model.rxns),1);
if ~strcmpi(model.mets,accum_met)
    disp('Metabolite is not present in the model')
else 
    rxnsToAdd.mets         = {accum_met};
    rxnsToAdd.rxns         = {rxnStr};
    rxnsToAdd.rxnNames     = {rxnStr};
    rxnsToAdd.stoichCoeffs = {-1};
    rxnsToAdd.lb           = 0;
    rxnsToAdd.ub           = 1000;
    rxnsToAdd.c            = 1;
    rxnsToAdd.subSystems   = {''};
    rxnsToAdd.grRules      = {''};
    newModel               = addRxns(model,rxnsToAdd,1,'c',false);
end