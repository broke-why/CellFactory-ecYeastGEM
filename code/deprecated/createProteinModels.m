%createProteinModels
load('../ModelFiles/mat/ecYeastGEM_batch.mat');
proteins = {'insulin' 'filgrastim' 'adalimumab'};
AA_tRNA = readtable('../data/pathway_templates/aminoacids_tRNA.txt','Delimiter','\t');
mets = ecModel_batch.mets(find(startsWith(ecModel_batch.mets,'s_')));
mets = strrep(mets,'s_','');
higherMet = str2double(mets);
higherMet = max(higherMet);
reactions = ecModel_batch.rxns(find(startsWith(ecModel_batch.rxns,'r_')));
reactions = strrep(reactions,'r_','');
higherRxn = str2double(reactions);
higherRxn = max(higherRxn);

for i=1:length(proteins)
    composition = readtable(['../data/pathway_templates/' proteins{i} '.txt'],'Delimiter','\t');
    composition.aminoacid = lower(composition.aminoacid);
    composition.aminoacid = strrep(composition.aminoacid,'aspartic acid','aspartate');
    composition.aminoacid = strrep(composition.aminoacid,'glutamic acid','glutamate');
    synth_rxn   = [];
    prods_synth = [];
    for j=1:height(composition)
        idx = find(strcmpi(AA_tRNA.AA,composition.aminoacid{j}));
        tRNA = AA_tRNA.met{idx};
        free = AA_tRNA.metFree{idx};
        coeff = composition.weight_AA(j);
        if j>1
            synth_rxn   = [synth_rxn ' +'];
            prods_synth = [prods_synth ' +'];
        end
        synth_rxn = [synth_rxn ' ' num2str(coeff) ' ' tRNA];
        prods_synth =  [prods_synth ' ' num2str(coeff) ' ' free];
    end
    newProt_c = ['s_' num2str(higherMet+1)];
    synth_rxn = [strtrim(synth_rxn) ' => ' strtrim(prods_synth) ' + ' newProt_c];
    synth_rxn = strtrim(synth_rxn);
    newProt_e = ['s_' num2str(higherMet+2)];
    trans_rxn = [newProt_c ' => ' newProt_e];
    exch_rxn  = [newProt_e ' => '];
    rxnsToAdd.equations = [{synth_rxn} {trans_rxn} {exch_rxn}]
    % % Metabolites to Add
    metsToAdd.mets          = {newProt_c newProt_e};
    metsToAdd.metNames      = {proteins{i} proteins{i}};
    metsToAdd.compartments  = {'c' 'e'};
    model = addMets(ecModel_batch,metsToAdd);
    % % Rxns to Add
    rxnsToAdd.rxns     = {['r_' num2str(higherRxn+1)] ['r_' num2str(higherRxn+2)] ['r_' num2str(higherRxn+3)]};
    rxnsToAdd.rxnNames = {[proteins{i} ' synthesis'] [proteins{i} ' transport'] [proteins{i} ' exchange']};
    % Define objective and bounds
    rxnsToAdd.c  = [0 0 1];
    rxnsToAdd.lb = [0 0 0];
    rxnsToAdd.ub = [1000 1000 1000];
	rxnsToAdd.grRules = {'' '' ''};
    model.c(:) = 0;
    model = addRxns(model,rxnsToAdd,1);
    [model,pos] = changeMedia_batch(model,'D-glucose exchange (reversible)','YEP');
    sol = solveLP(model);
    if ~isempty(sol.x)
        printFluxes(model,sol.x)
        save(['../ModelFiles/production_ecModels/ec' proteins{i} '.mat'],'model')
    end
end
