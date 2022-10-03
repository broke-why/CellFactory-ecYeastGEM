%Clone GECKO and go to the right branch
clear
current = pwd;
git clone https://github.com/SysBioChalmers/GECKO.git
%git('clone --depth 1 https://github.com/SysBioChalmers/ecFactory.git')
clc
success    = [];
fail       = [];
models     = [];
objectives = [];
steps      = [];
objIdx     = [];
noFlux     = [];
fileNames  = dir('../ModelFiles/production_ecModels');
strain_conditions = readtable('../data/strain_conditions.txt','Delimiter','\t');
mkdir('../results/ECC')
for i=1:length(fileNames)
    clc
    cd (current)
    file = fileNames(i).name;
    if contains(file,'.mat')
        modelName = file(1:end-4);
        modelName = strrep(modelName,'_WBG','');
        disp([modelName ':  #' num2str(i)])
        if startsWith(modelName,'ec')
            modelName = lower(modelName(3:end));
            load(['../ModelFiles/production_ecModels/' file])
            targetIndex = find(model.c);
            targetRxn   = model.rxns(targetIndex);
            models      = [models;{file}];
            objectives  = [objectives;model.rxnNames(targetIndex)];
            objIdx      = [objIdx;targetIndex];
            disp(' ')
            %Set media conditions
            cd ecFactory/code
            CSname = 'D-glucose exchange (reversible)';
            tempModel = changeMedia_batch(model,CSname,'Min');
            modelParam.rxnTarget = model.rxns(targetIndex);       
            modelParam.CS_MW  = 0.18015;
            CS_index  = find(strcmpi(model.rxnNames,CSname));
            modelParam.CSrxn  = model.rxns{strcmpi(model.rxnNames,CSname)};
            modelParam.growthRxn = model.rxns{strcmpi(model.rxnNames,'biomass pseudoreaction')};
            growthPos = find(strcmpi(tempModel.rxnNames,'growth'));
            %Unconstrain growth
            tempModel = setParam(tempModel,'lb',modelParam.growthRxn,0);
            tempModel = setParam(tempModel,'ub',modelParam.growthRxn,1000);
            %Get biomass yield for a unit glucose uptake rate
            tempModel = setParam(tempModel,'obj',modelParam.growthRxn,1);
            tempModel = setParam(tempModel,'ub',modelParam.CSrxn,1);
            solution  = solveLP(tempModel,1);    
            %Check if model can carry flux for the target rxn
            flux = haveFlux(tempModel,1-12,targetIndex);
            if ~flux | solution.x(growthPos)<=1e-4
                %find specific conditions and strain background
                idx = find(strcmpi(strain_conditions.Chemicals,modelName));
                if ~isempty(idx)
                    medium    = strain_conditions.Medium{idx};
                    c_source  = [strain_conditions.CarbonSource{idx} ' exchange (reversible)'];
                    tempModel = changeMedia_batch(tempModel,c_source,medium);
                    solution  = solveLP(tempModel,1);
                end
            end
            WT_yield  = solution.x(growthPos)/(solution.x(CS_index)*modelParam.CS_MW);
            disp(['The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);
            %reformat chemical name
            modelName = strrep(modelName,' ','_');
            modelName = strrep(modelName,'-','_');
            modelName = strrep(modelName,',','');
            modelName = strrep(modelName,'(','');
            modelName = strrep(modelName,')','');
            modelName = regexprep(modelName,'[0-9]','');
            modelName = strrep(modelName,'__','');
            modelName = lower(modelName);
            while startsWith(modelName,'_')
                modelName = modelName(2:end);
            end
            
            if flux & solution.x(growthPos)>0
                resultsFolder = ['../../../results/production_targets/' modelName '_targets'];
                mkdir(resultsFolder)
                WT_yield = 0.48; %WT yield on glucose minimal media
                expYield = 0.49*WT_yield;
                %try                    
                    [optStrain,candidates,step] = run_ecFactory(tempModel,modelParam,expYield,resultsFolder,false);
                    success = [success; {modelName}];
                    %                 catch
%                     disp('The model is not suitable for robust ecFSEOF')
%                     fail = [fail; {modelName}];
%                 end
                disp(' ')
            else
                noFlux = [noFlux; {modelName}];
                disp('The model is not able to carry flux through the target reaction with the imposed constraints')
            end
        end
    end
end
