%Clone GECKO and go to the right branch
clear
current = pwd;
clc

fileNames  = dir('../ModelFiles/production_ecModels');
strain_conditions = readtable('../data/genetic_background.txt','Delimiter','\t');
models = [];
objectives = [];
objIdx = [];
newTable = table();

bioY  = [];


for i=1:length(fileNames)
    cd (current)
    file = fileNames(i).name;
    if contains(file,'.mat')
        modelName = file;%file(1:end-4);
        modelName = strrep(modelName,'_WBG','');
        disp([modelName ':  #' num2str(i)])
        if startsWith(modelName,'ec')
            modelName = lower(modelName);
            load(['../ModelFiles/production_ecModels/' file])
            targetIndex = find(model.c);
            targetRxn   = model.rxns(targetIndex);
            models      = [models;{file}];
            objectives  = [objectives;model.rxnNames(targetIndex)];
            objIdx      = [objIdx;targetIndex];
            CS_MW       = 0.18015;
            
            %find specific conditions and strain background
            idx = find(strcmpi(strain_conditions.ecModel,modelName));
            if ~isempty(idx)
                medium    = strain_conditions.media{idx};
                c_source  = [strain_conditions.carbon_source{idx} ' exchange (reversible)'];
                if strcmpi(c_source,'D-xylose exchange (reversible)')
                    CS_MW     = 0.15013;
                end
                CS_index  = find(strcmpi(model.rxnNames,c_source));
                growthPos = find(strcmpi(model.rxnNames,'growth'));
                
                model = changeMedia_batch(model,c_source,medium);
                WTmodel = model;
                disp('Wild-type')
                [bioYield1,prod_yield1, prod_rate1] = calculate_potential(WTmodel,growthPos,targetIndex,CS_index,CS_MW);
                disp(' ')
                
                try
                    if ~isempty(strain_conditions.deletion(idx))
                        genes2delete  = strsplit(strain_conditions.deletion{idx},',');
                        columns       = cell(length(genes2delete),1);
                        columns(:,1)  = {0};
                        modifications = [genes2delete' columns columns];
                        %Search genes in model
                        [iA,iB] = ismember(modifications(:,1),model.genes);
                        modifications = modifications(iA,:);
                        if ~isempty(modifications)
                            model = getMutant(model,modifications,[],true);
                        end
                    end
                    
                    if ~isempty(strain_conditions.downregulation(idx))
                        genes2dReg    = strsplit(strain_conditions.deletion{idx},',');
                        columns       = cell(length(genes2dReg),2);
                        columns(:,1)  = {2};
                        columns(:,2)  = {0.5};
                        modifications = [genes2dReg' columns columns];
                        %Search genes in model
                        [iA,iB] = ismember(modifications(:,1),model.genes);
                        modifications = modifications(iA,:);
                        if ~isempty(modifications)
                            model = getMutant(model,modifications,[],true);
                        end
                    end
                    
                    if ~isempty(strain_conditions.overexpression(idx))
                        genes2OE  = strsplit(strain_conditions.overexpression{idx},',');
                        columns       = cell(length(genes2OE),1);
                        columns(:,1)  = {2};
                        modifications = [genes2OE' columns columns];
                        %Search genes in model
                        [iA,iB] = ismember(modifications(:,1),model.genes);
                        modifications = modifications(iA,:);
                        if ~isempty(modifications)
                            model = getMutant(model,modifications,[],true);
                        end
                    end
                    disp('Mutant (data)')
                    [bioYield2,prod_yield2, prod_rate2] = calculate_potential(model,growthPos,targetIndex,CS_index,CS_MW);
                    disp(' ')
                catch
                    bioYield2   = 0;
                    prod_yield2 = 0;
                    prod_rate2  = 0;
                end
                %get predicted optimal mutant
                modelName  = modelName(3:(end-4));
                %modelName  = strrep(modelName,'-','_');
                optModel = WTmodel;
                try
                    candidates = readtable(['../results/' modelName '_targets/compatible_genes_results.txt'],'delimiter','\t');
                    dels = find(strcmpi(candidates.actions,'KO'));
                    OExs = find(strcmpi(candidates.actions,'OE'));
                    dReg = find(strcmpi(candidates.actions,'KD'));
                    
                    genes2KO      = candidates.genes(dels);
                    columns       = cell(length(genes2KO),1);
                    columns(:,1)  = {0};
                    modifications = [genes2KO columns columns];
                    %Search genes in model
                    [iA,iB] = ismember(modifications(:,1),model.genes);
                    modifications = modifications(iA,:);
                    if ~isempty(modifications)
                        optModel = getMutant(optModel,modifications,[],true);
                    end
                    
                    genes2OE      = candidates.genes(OExs);
                    columns       = cell(length(genes2OE),1);
                    columns(:,1)  = {2};
                    modifications = [genes2OE columns columns];
                    %Search genes in model
                    [iA,iB] = ismember(modifications(:,1),model.genes);
                    modifications = modifications(iA,:);
                    if ~isempty(modifications)
                        optModel = getMutant(optModel,modifications,[],true);
                    end
                    
                    genes2KD      = candidates.genes(dReg);
                    columns       = cell(length(genes2KD),2);
                    columns(:,1)  = {2};
                    columns(:,2)  = {0.5};
                    modifications = [genes2KD columns columns];
                    %Search genes in model
                    [iA,iB] = ismember(modifications(:,1),model.genes);
                    modifications = modifications(iA,:);
                    if ~isempty(modifications)
                        optModel = getMutant(optModel,modifications,[],true);
                    end
                    disp('Mutant (predictions)')
                    [bioYield3,prod_yield3, prod_rate3] = calculate_potential(optModel,growthPos,targetIndex,CS_index,CS_MW);
                    disp(' ')
                    
                catch
                    bioYield3   = 0;
                    prod_yield3 = 0;
                    prod_rate3  = 0;
                end
                newRow = {modelName bioYield1 bioYield2 bioYield3 prod_yield1 prod_yield2 prod_yield3 prod_rate1 prod_rate2 prod_rate3};
                newTable  = [newTable; newRow];
            end
            
            
            if ~prod_rate3>0
                disp('The model is not able to carry flux through the target reaction with the imposed constraints')
            end
        end
        disp(' ')
        disp(' ')
    end
    %clc
end
varNames = {'chemical' 'bioY1' 'bioY2' 'bioY3' 'proY1' 'proY2' 'proY3' 'proR1' 'proR2' 'proR3'};
newTable.Properties.VariableNames = varNames;
writetable(newTable,'../results/production_capacity_comparison.txt','delimiter','\t','QuoteStrings',false)

