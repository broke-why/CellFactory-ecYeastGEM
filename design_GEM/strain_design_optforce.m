% optforce for yeast-GEM

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
cd strain_design_cobra

% add constraint
model1 = changeObjective(model1,'r_2111');
model1=changeRxnBounds(model1,'r_1714', -10, 'b'); %% qs glucose
FBAsolution1 = optimizeCbModel(model1,'max','one');
FBAsolution1.f


%% 
% Then, we calculate the maximum specific growth rate and the maximum production 
% rate for 3HP.

growthRate = optimizeCbModel(model1); 
fprintf('The maximum growth rate is %1.2f', growthRate.f);

model1 = changeObjective(model1, 'new5');
maxSucc = optimizeCbModel(model1);
fprintf('The maximum production rate of 3HP is %1.2f', maxSucc.f);

%% 
% *TIP: *The biomass reaction is usually set to 1%-10% of maximum theoretical 
% biomass yield when running the following steps, to prevent solutions without 
% biomass formation.
% 
% # Maximizing product formation
% # Finding MUST sets of second order
% # Finding FORCE sets
%% STEP 2: Define constraints for both wild-type and mutant strain
% *TIMING*: This step should take a few days or weeks, depending on the information 
% available for your species. 
% 
% *CRITICAL STEP*: This is a manual task, so you should search for information 
% in articles or even perform your own experiments. You can also make assumptions 
% for describing the phenotypes of both strains which will make this task a little 
% faster but make sure to have two strains different enough, because you should 
% be able to find differences in reactions ranges. 
% 
% We define constraints for each strain as follows: 
% 
% # The WT strain's biomass function ("R75") is constrained to near the maximum 
% growth rate. 
% # The mutant strain's biomass function is set to zero. 3HP export ('EX_suc') 
% is forced to be the maximum as calculated previously.
%%
constrWT = struct('rxnList', {{'r_2111'}}, 'rxnValues', 0.86, 'rxnBoundType', 'b')
constrMT = struct('rxnList', {{'r_2111', 'new5'}}, 'rxnValues', [0, 17], ...
                  'rxnBoundType', 'bb')
              
              
 %% Step 3: Flux Variability Analysis
% *TIMING*: This task should take from a few seconds to a few hours depending 
% on the size of your reconstruction
% 
% We  run the FVA analysis for both strains
%%
[minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(model1, ...
                                                                     constrWT, constrMT);
disp([minFluxesW, maxFluxesW, minFluxesM, maxFluxesM]);

%% 
% Now, the run the next step of OptForce.
%% Step 4: Find Must Sets
% *TIMING: *This task should take from a few seconds to a few hours depending 
% on the size of your reconstruction
% 
% First, we define an ID for this run. Each time you run the functions associated 
% to the optForce procedure, some folders can be generated to store inputs used 
% in that run. Outputs are stored as well. These folders will be located inside 
% the folder defined by your run ID. Thus, if your runID is ''TestOptForce", the 
% structure of the folders will be the following:
% 
% |├── CurrentFolder|
% 
% ||   ├── TestOptForce|
% 
% ||   |   ├── Inputs|
% 
% ||   |   └── Outputs|
% 
% To avoid the generation of inputs and outputs folders, set |keepInputs 
% = 0|, |printExcel = 0| and |printText = 0|.
% 
% Also, a report of the run is generated each time you run the functions 
% associated to the optForce procedure. So, the idea is to give a different |runID| 
% each time you run the functions, so you will be able to see the report (inputs 
% used, outputs generated, errors in the run) for each run.
% 
% We define then our |runID|.
%%
runID = 'TestOptForceM';
%% 
% Fow now, only functions to find first and second order must sets are supported 
% in this third step. As depicted in *Figure 1*, the first order must sets are 
% MUSTU and MUSTL; and second order must sets are MUSTUU, MUSTLL and MUSTUL.
% 
% *A) Finding first order must sets*
% 
% We define constraints.

constrOpt = struct('rxnList', {{'r_1714', 'r_2111', 'new5'}}, 'values', [-10, 0, 17]');
%% 
% We then run the functions |findMustL| and |findMustU| that will allow 
% us to find |mustU| and |mustL| sets, respectively.
% 
% *i) MustL Set: *

[mustLSet, pos_mustL] = findMustL(model1, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
                                  'outputFileName', 'MustL' , 'printExcel', 1, 'printText', 1, ...
                                  'printReport', 1, 'keepInputs', 1, 'verbose', 0);
%% 
% Note that the folder "TestOptForceM" was created. Inside this folder, 
% two additional folders were created: "InputsMustL" and "OutputsMustL". In the 
% inputs folder you will find all the inputs required to run the the function 
% |findMustL|. Additionally, in the outputs folder you will find the |mustL| set 
% found, which were saved in two files (.xls and .txt). Furthermore, a report 
% which summarize all the inputs and outputs used during your running was generated. 
% The name of the report will be in this format "report-Day-Month-Year-Hour-Minutes". 
% So, you can mantain a chronological order of your experiments. 
% 
% We display the reactions that belongs to the |mustL| set.

disp(mustLSet)
%% 
% *ii) MustU set: *
%%
[mustUSet, pos_mustU] = findMustU(model1, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                  'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
                                  'outputFileName', 'MustU' , 'printExcel', 1, 'printText', 1, ...
                                  'printReport', 1, 'keepInputs', 1, 'verbose', 0);
%% 
% Note that the folders "InputsMustU" and "OutputsFindMustU" were created. 
% These folders contain the inputs and outputs of |findMustU|, respectively. 
% 
% We display the reactions that belongs to the |mustU| set.

disp(mustUSet)
%% 
% *B) Finding second order must sets *
% 
% First, we define the reactions that will be excluded from the analysis. 
% It is suggested to include in this list the reactions found in the previous 
% step as well as exchange reactions.
constrOpt = struct('rxnList', {{'r_1714', 'r_2111', 'new5'}}, 'values', [-10, 0, 17]');
%constrOpt = struct('rxnList', {{'EX_gluc', 'R75', 'EX_suc'}}, 'values', [-100, 0, 155.5]');
exchangeRxns = model1.rxns(cellfun(@isempty, strfind(model1.rxnNames, 'exchange')) == 0);
excludedRxns = unique([mustUSet; mustLSet; exchangeRxns]);
%% 
% Now, we run the functions for finding second order must sets.
% 
% *i) MustUU: *

[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
    findMustUU(model1, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);
%% 
% Note that the folders "InputsMustUU" and "OutputsFindMustUU" were created. 
% These folders contain the inputs and outputs of |findMustUU|, respectively. 
% 
% We display the reactions that belongs to the| mustUU| set

disp(mustUU);

%% 
% *ii) MustLL: *

[mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
    findMustLL(model1, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);
%% 
% Note that the folders "InputsMustLL" and "OutputsFindMustLL" were created. 
% These folders contain the inputs and outputs of |findMustLL|, respectively. 
% 
% We display the reactions that belongs to the |mustLL| set. In this case, 
% |mustLL| is an empty array because no reaction was found in the |mustLL| set.

disp(mustLL);
%% 
% *iii) MustUL: *

[mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = ...
    findMustUL(model1, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUL', 'outputFileName', 'MustUL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);
%% 
% Note that the folders "InputsMustUL" and "OutputsFindMustUL" were created. 
% These folders contain the inputs and outputs of |findMustUL|, respectively. 
% 
% We display the reactions that belongs to the |mustUL| set. In this case, 
% |mustUL| is an empty array because no reaction was found in the |mustUL| set.

disp(mustUL);
%% 
% *TROUBLESHOOTING 1: * "I didn't find any reaction in my must sets"
% 
% *TROUBLESHOOTING 2: * "I got an error when running the |findMustX| functions 
% (X = L or U or LL or UL or UU depending on the case)"
%% Step 5: OptForce
% *TIMING: *This task should take from a few seconds to a few hours depending 
% on the size of your reconstruction
% 
% We define constraints and we define |K| the number of interventions allowed, 
% |nSets| the maximum number of sets to find, and |targetRxn| the reaction producing 
% the metabolite of interest (in this case, 3HP). 
% 
% Additionally, we define the |mustU| set as the union of the reactions that 
% must be upregulated in both first and second order must sets; and |mustL| set 
% as the union of the reactions that must be downregulated in both first and second 
% order must sets .
%%
mustU = unique(union(mustUSet, mustUU));
mustL = unique(union(mustLSet, mustLL));
targetRxn = 'new5';
biomassRxn = 'r_2111';
k = 1;
nSets = 1;
constrOpt = struct('rxnList', {{'r_1714','r_2111'}}, 'values', [-10, 0]);

[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model1, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
             'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
             'runID', runID, 'outputFolder', 'OutputsOptForce', ...
             'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
             'printReport', 1, 'keepInputs', 1, 'verbose', 1);
%% 
% Note that the folders "InputsOptForce" and "OutputsOptForce" were created. 
% These folders contain the inputs and outputs of |optForce|, respectively.
% 
% We display the reactions found by |optForce|

disp(optForceSets)
%% 
% The reaction found was "SUCt", i.e. a transporter for 3HP (a very 
% intuitive solution).
% 
% Next, we will increase |k| and we will exclude "SUCt" from upregulations 
% to find non-intuitive solutions. 
% 
% *TIP: *Sometimes the product is at the end of a long linear pathway. In 
% that case, the recomendation is to also exclude most reactions on the linear 
% pathway. Essential reactions and reactions not associated with any gene (i.e. 
% spontaneous reacitons) should also be excluded. 
% 
% We will only search for the 20 best solutions, but you can try with a higher 
% number.
% 
% We will change the runID to save this second result (K = 2) in a diffetent 
% folder than the previous result (K = 1) 

k = 2;
nSets = 20;
runID = 'TestOptForceM2';
excludedRxns = struct('rxnList', {{'new5'}}, 'typeReg','U');
[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
    optForce(model1, targetRxn, biomassRxn, mustU, mustL, ...
             minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
             'k', k, 'nSets', nSets, 'constrOpt', constrOpt, ...
             'excludedRxns', excludedRxns, ...
             'runID', runID, 'outputFolder', 'OutputsOptForce', ...
             'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
             'printReport', 1, 'keepInputs', 1, 'verbose', 1);

%% 
% Note that the folders "InputsOptForce" and "OutputsOptForce" were created 
% inside TestOptForce2. These folders contain the inputs and outputs of |optForce|, 
% respectively.
% 
% We display the reactions found by |optForce| 
disp(optForceSets)


              
              
              
