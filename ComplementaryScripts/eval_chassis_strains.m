current = pwd;
chemicals_info = readtable('../ComplementaryData/chemicals_info.txt','Delimiter','\t');
strain_conditions = readtable('../ComplementaryData/genetic_background.txt','Delimiter','\t');

comp_classes   = unique(chemicals_info.class);
class_short    = {'alc' 'alk' 'AAs' 'aro' 'bio' 'FAL' 'fla' 'oAc' 'stb' 'ter'};
%load chassis strain files
files = dir('../results');
chassis_mod = readtable('../results/chassis_strains_modifications.txt','Delimiter','\t');
chemicals_info.Name = strtrim(chemicals_info.Name);
for i=5
    cd  (current)
    fileStr = ['../results/chassis_strain_' num2str(i) '_chemicals.txt'];
    if isfile(fileStr)
        results = table();
        chemicals = readtable(fileStr,'Delimiter','\t');
        chemicals.product  = strrep(chemicals.product,'coumaricacid','p-coumaric acid');
        chemicals.product  = strrep(chemicals.product,'taxadienalphaylacetate', 'taxadien_5alpha_yl_acetate');
        chemicals.product  = strrep(chemicals.product,'oleanolate','oleate');
        chemicals.product  = strrep(chemicals.product,'glycyrrhetinicacid','glycyrrhetinic_acid');
        chemicals.product  = strrep(chemicals.product,'betacarotene','_-carotene');
        chemicals.product  = strrep(chemicals.product,'betaionone','_-ionone');
        chemicals.product  = strrep(chemicals.product,'dha','DHA');
        chemicals.product  = strrep(chemicals.product,'nicotianamine','Nicotianamine');
        chemicals.product  = strrep(chemicals.product,'sreticuline','(S)-Reticuline');
        chemicals.product  = strrep(chemicals.product,'epa','EPA');

        %load each model
        for j=1:height(chemicals)
            chem  = chemicals.product(j);
            
            index = find(contains(chemicals_info.Name,chem));
            if ~isempty(index)
                modelStr = chemicals_info.ecModel{index};
                try
                    disp(modelStr)
                    load(['../ModelFiles/production_ecModels/' modelStr])
                catch
                    disp(modelStr)
                    modelStr = strrep(modelStr,'.mat','_WBG.mat');
                    load(['../ModelFiles/production_ecModels/' modelStr])
                end
                targets   = chassis_mod(1:i,:);
                genes = [];
                for k = 1:length(targets.geneTargets)
                    idx   = find(strcmpi(model.geneShortNames,targets.geneTargets{k}));
                    genes = [genes; model.genes(idx)];
                    actions = strcmpi(targets.modifications,'OE')*2;
                    actions = num2cell(actions);
                end
                cd method
                modifications = [genes actions actions];    
                mutant = getMutant(model,modifications,[],false);
                cd ..
                index = find(contains(strain_conditions.ecModel,modelStr));
                if ~isempty(index)
                    media = strain_conditions.media{index};
                    cSource = 'D-glucose exchange (reversible)';
                    
                else
                    media = 'Min';
                    cSource = 'D-glucose exchange (reversible)';
                end
                model  = changeMedia_batch(model,cSource,media);
                mutant = changeMedia_batch(mutant,cSource,media);
                %
                growthPos = find(strcmpi(model.rxnNames,'growth'));
                targetIndex = find(model.c);
                CS_index = find(strcmpi(model.rxnNames,cSource));
                CS_MW = 0.180;
                [bioYield,prod_yield, prod_rate] = calculate_potential(model,growthPos,targetIndex,CS_index,CS_MW);
                [bioYieldM,prod_yieldM, prod_rateM] = calculate_potential(mutant,growthPos,targetIndex,CS_index,CS_MW);
                FC_bY = bioYieldM/bioYield;
                FC_pY = prod_yieldM/prod_yield;
                FC_pR = prod_rateM/prod_rate;
                newRow = [chem, height(chemicals), FC_bY, FC_pY, FC_pR];
                results = [results; newRow];
            else
                disp([num2str(i) ': ' chem])
            end
            disp(' ')
        end
        results.Properties.VariableNames = {'chemicals' 'mod_number' 'FC_bY' 'FC_pY' 'FC_pR'};
        writetable(results,['../results/chassis_strain_' num2str(i) '_prod_capabilities.txt'],'delimiter','\t','QuoteStrings',false)

                
                
                
%         targets   = chassis_mod(1:i,:);
%         cd method
%         
%         genes   = targets.geneTargets;
%         [~,iB]  = ismember(ecModel_batch
%         genes   = 
%         actions = targets.modifications;
%         
%         genes2OE      = candidates.genes(OExs);
%         columns       = cell(length(genes2OE),1);
%         columns(:,1)  = {2};
%         modifications = [genes2OE columns columns];
%         
%         mutantModel = getMutant(model,modifications,enzUsage,message)
    end
    
end

