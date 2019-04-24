% optknock for yeast-GEM

cd ..
cd model
model_new = readCbModel('yeastGEM_3HP.xml')
model1 = model_new


% add revesibility into the model
model1.rev = zeros(length(model1.lb),1)
for i = 1: length(model1.lb)
    if model1.lb(i) == 0
        model1.rev(i) = 0
    else
        model1.rev(i) =1
    end
end

changeCobraSolver('gurobi', 'ALL');

% general analysis
cd ..
cd script

% The test data is from Zongjie Dai
% zjd25
% sbustrate
model1 = changeObjective(model1,'r_2111');
model1=changeRxnBounds(model1,'r_1714', -1.634, 'b'); %% qs glucose
model1=changeRxnBounds(model1,'r_1808', 0.050, 'b'); %% q-glycerol
model1=changeRxnBounds(model1,'r_2033', 0.0166, 'b'); %% q-pyruvate
model1=changeRxnBounds(model1,'r_2056', 0.0, 'b');  %% q-succiate
model1=changeRxnBounds(model1,'r_1634', 1.360, 'b'); %% q-acetate
model1=changeRxnBounds(model1,'r_1761', 0, 'b'); %% q-ethonal
FBAsolution1 = optimizeCbModel(model1,'max','one');
FBAsolution1.f

% change the objective function to test 3HP production
% now the model could predict the 3HP production
model1 = changeObjective(model1,'new5');
solMax = optimizeCbModel(model1,'max');
solMin = optimizeCbModel(model1,'min');
solMax.f

% prepare the optknock step
% find the target to improve the production
biomass = 'r_2111';
%% 
% Define the maximum number of solutions to find (i.e., maximum number of 
% removable reactions that lead to the overproduction of the metabolite of interest)
selectedRxnList = model1.rxns(1:200);


%% 
% Then, calculates the production of metabolites before running optKnock.

% determine 3HP production and growth rate
model1 = changeObjective(model1,'r_2111');
fbaWT = optimizeCbModel(model1);
succFluxWT = fbaWT.x(strcmp(model1.rxns, 'new5'));


% *Aim:* *finding optKnock reactions sets of size 2 for increasing production 
% of 3HP*
%%
fprintf('\n...EXAMPLE 1: Finding optKnock sets of size 2 or less...\n\n')
% Set optKnock options
% The exchange of 3HP will be the objective of the outer problem
options = struct('targetRxn', 'new5', 'numDel', 3);
% We will impose that biomass be at least 50% of the biomass of wild-type
constrOpt = struct('rxnList', {{biomass}},'values', 0.5*fbaWT.f, 'sense', 'G');
% We will try to find 10 optKnock sets of a maximun length of 2
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 1;
threshold = 10;
while nIter < threshold
    
    
    fprintf('...Performing optKnock analysis...\n')
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model1, selectedRxnList, options, constrOpt);
    else
        optKnockSol = OptKnock(model1, selectedRxnList, options, constrOpt, previousSolutions, 1);
    end
    
    % determine 3HP production and growth rate after optimization
    succFluxM1 = optKnockSol.fluxes(strcmp(model1.rxns, 'new5'));
    growthRateM1 = optKnockSol.fluxes(strcmp(model1.rxns, biomass));
    setM1 = optKnockSol.rxnList;
    
    if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions = contPreviousSolutions + 1;
        %printing results
        fprintf('optKnock found a optKnock set of large %d composed by ', length(setM1));
        for j = 1:length(setM1)
            if j == 1
                fprintf('%s', setM1{j});
            elseif j == length(setM1)
                fprintf(' and %s', setM1{j});
            else
                fprintf(', %s', setM1{j});
            end
        end
        fprintf('\n');
        fprintf('The production of 3HP after optimization is %.3f \n', succFluxM1);
        fprintf('The growth rate after optimization is %.3f \n', growthRateM1);
        fprintf('...Performing coupling analysis...\n');
        %[type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model1, setM1, 'new5');
        %fprintf('The solution is of type: %s\n', type);
        %fprintf('The maximun growth rate given the optKnock set is %.3f\n', maxGrowth);
        %fprintf(['The maximun and minimun production of 3HP given the optKnock set is ' ...
        %         '%.3f and %.3f, respectively \n\n'], minProd, maxProd);
        %if strcmp(type, 'growth coupled')
        %    singleProductionEnvelope(model, setM1, 'new5', biomass, 'savePlot', 1, 'showPlot', 1, ...
        %                             'fileName', ['succ_ex1_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
        %end
    else
        if nIter == 1
            fprintf('optKnock was not able to found an optKnock set\n');
        else
            fprintf('optKnock was not able to found additional optKnock sets\n');
        end
        break;
    end
    nIter = nIter + 1;
end










