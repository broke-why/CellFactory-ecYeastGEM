%Clone GECKO and go to the right branch
clear
current = pwd;
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git checkout feat/add_FSEOF_utilities
cd ..
clc
success = [];
fail = [];
models = [];
objectives = [];
objIdx = [];
fileNames = dir('../ModelFiles/production_ecModels');
for i=1:length(fileNames)
    cd (current)
    file = fileNames(i).name;
    if contains(file,'.mat')
        modelName = file(1:end-4);
        disp([modelName ':  #' num2str(i)])
        if startsWith(modelName,'ec')
            modelName = lower(modelName(3:end));          
            load(['../ModelFiles/production_ecModels/' file])
            model = check_enzyme_fields(model);
            targetIndex = find(model.c);
            models = [models;{file}];
            objectives = [objectives;model.rxnNames(targetIndex)];
            objIdx = [objIdx;targetIndex];
            disp(' ')
            %Set media conditions
            cd GECKO/geckomat/kcat_sensitivity_analysis
            c_source = 'D-glucose exchange (reversible)';
            tempModel = changeMedia_batch(model,c_source);
            CS_MW     = 0.18015;
            CS_index  = find(strcmpi(tempModel.rxnNames,c_source));
            growthPos = find(strcmpi(tempModel.rxnNames,'growth'));
            cd ../../../
            %Unconstrain growth
            tempModel = setParam(tempModel,'lb',growthPos,0);
            tempModel = setParam(tempModel,'ub',growthPos,1000);
            %Get biomass yield for a unit glucose uptake rate
            tempModel = setParam(tempModel,'obj',growthPos,1);
            tempModel = setParam(tempModel,'ub',CS_index,1);
            solution  = solveLP(tempModel,1);
            WT_yield  = solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
            disp(['The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);
            flux = haveFlux(tempModel,1-12,targetIndex);
            if flux
                resultsFolder = ['../results/' modelName '_targets'];
                mkdir(resultsFolder)
                expYield = 0.2;
                try
                    [optStrain,candidates] = robust_ecFSEOF(tempModel,tempModel.rxns(targetIndex),0.122,CS_MW,resultsFolder);
                    success = [success; {modelName}];
                catch
                    disp('The model is not suitable for robust ecFSEOF')
                    fail = [fail; {modelName}];
                end
                disp(' ')
                clc
             end
        end
    end
end
