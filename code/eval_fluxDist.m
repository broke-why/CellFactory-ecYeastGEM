clear
current = pwd;
%subSystems_GEM = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../data/chemicals_info.txt','Delimiter','\t');
d         = dir('../results');
isub      = [d(:).isdir]; %# returns logical vector
load('../ModelFiles/mat/ecYeastGEM_batch.mat')
load('../ModelFiles/mat/yeastGEM.mat')
GEM  = ravenCobraWrapper(model);
chemicals_info.Name = strtrim(chemicals_info.Name);
prod_capabilities  = table();
rxns = ecModel_batch.rxns;
fluxes = table(rxns);
families = [];
biomass_prod = false;
Prot_cost =[]; 
%Loop through all maxRate folders
mkdir('../results/yieldPlots')
for i=1:height(chemicals_info)
    compound = chemicals_info.Name{i};
    model    = [];
    %try to load GEM
    modelStr = ['../ModelFiles/production_ecModels/' chemicals_info.ecModel{i}];
    try
        modelStr = strrep(modelStr,'.mat','_WBG.mat');
        load(modelStr);
    catch
        modelStr = strrep(modelStr,'_WBG.mat','.mat');
        load(modelStr);
    end
    
    if ~isempty(model)
        ecModel = model;
        %Check presence of objective reaction
        index_ec = find(model.c);
        objRxn   = model.rxns(index_ec);
        indexGEM = find(strcmpi(GEM.rxns,objRxn));
        if ~isempty(index_ec)
            %c_source = 'D-glucose exchange';
            if startsWith(compound,'L-')
                AA = true;
            else
                AA = false;
            end
            ecModel = changeMedia_batch(ecModel,'D-glucose exchange (reversible)','Min',AA);
            %unconstrain growth
            ecModel = setParam(ecModel,'lb',find(strcmpi(ecModel.rxnNames,'growth')),0);
            %Check growth capabilities
            temp = setParam(ecModel,'obj',find(strcmpi(ecModel.rxnNames,'growth')),1);
            sol1 = solveLP(temp);
            gRate = -sol1.f>0;
            %check flux
            sol1    = solveLP(ecModel);
            obj1    = -sol1.f;
            flux_ec = obj1>=1E-6;

            if ~flux_ec | ~gRate
                ecModel = changeMedia_batch(ecModel,'D-glucose exchange (reversible)','YEP',AA);
                %check flux
                sol1    = solveLP(ecModel,1);
                flux_ec = (-sol1.f)>0;
            end
            
            if flux_ec
                disp(['Performing pFBA for: ' compound])
                %Get biomass and product yields for low glucose
                %consumption
                %If objective reaction is also present in GEM, then add GEM
                CS_index  = find(strcmpi(ecModel.rxnNames,'D-glucose exchange (reversible)'));
                growthPos = find(strcmpi(ecModel.rxnNames,'growth'));
                [bioY_ec,proY_ec,rate_ec,fluxDist] = calculate_potential(ecModel,growthPos,index_ec,CS_index,0.180,biomass_prod);
                %Get flux distribution
                cost = fluxDist.flux(strcmpi(fluxDist.rxns,ecModel.rxns(index_ec)));
                [presence,iB]   = ismember(fluxes.rxns,fluxDist.rxns);
                fluxes = fluxes(presence,:);
                fluxDist = fluxDist(iB(presence),:);
                if ~biomass_prod
                    cost = fluxDist.flux(strcmpi(fluxDist.rxns,'prot_pool_exchange'))/cost;
                end
                str = strrep(chemicals_info.ecModel{i},'-','_');
                str = str(3:(end-3));
                str = regexprep(str,'[^a-zA-Z]','');
                str = strrep(str,',','_');
                str = strrep(str,'.mat','');
                str = strrep(str,'(','');
                str = strrep(str,')','');
                str = lower(str);
                              
                eval(['fluxes.' str '=fluxDist.flux;'])
                %now with GEM
                if ~isempty(indexGEM)
                    GEM  = changeMedia_Original(GEM,AA,1,1000);
                    GEM  = setParam(GEM,'obj',indexGEM,1);
                    sol2 = solveLP(GEM);
                    obj2 = -sol2.f;
                    FC   = obj2/obj1;
                    if obj2>0
                        CS_index  = find(strcmpi(GEM.rxnNames,'D-glucose exchange'));
                        growthPos = find(strcmpi(GEM.rxnNames,'growth'));
                        [bioY,proY,rate] = calculate_potential(GEM,growthPos,indexGEM,CS_index,0.180,biomass_prod);
                    end
                else
                    bioY =NaN;proY=NaN;rate=NaN;
                end
                newRow            = [{compound},bioY_ec,proY_ec,rate_ec,bioY,proY,rate];
                prod_capabilities = [prod_capabilities; newRow];
                families = [families;chemicals_info.class(i)];
                Prot_cost =[Prot_cost;cost]; 
            end
        end
    else
        sprintf(['ecModel file for: ' compound ' not found.\n'])
    end
end
%
if biomass_prod
    file1 = '../results/prodCapabilities_allChemicals_wBio.txt';
    file2 = '../results/proteinLimitations_allChemicals_wBio.txt';
    file3 = '../results/fluxDist_distance_allChemicals_wBio.txt';
else
    file1 = '../results/prodCapabilities_allChemicals.txt';
    file2 = '../results/proteinLimitations_allChemicals.txt';
    file3 = '../results/fluxDist_distance_allChemicals.txt';
end
    
prod_capabilities.Properties.VariableNames = {'compound' 'bioYield_ec' 'prodYield_ec' 'prodRate_ec' 'bioYield_gem' 'prodYield_gem' 'prodRate_gem'};
writetable(prod_capabilities,file1,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
%
FC = prod_capabilities.prodYield_ec./prod_capabilities.prodRate_ec;
newTable = table(prod_capabilities.compound,FC,families);
newTable.Properties.VariableNames = {'compound' 'Prod_FC' 'family'};
if ~isempty(Prot_cost)
    newTable.Prot_cost = Prot_cost;
end
writetable(newTable,file2,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
%calculate euclidean distance matrix (flux distributions)
[m,n]   = size(fluxes);
distMat = zeros(n-1,n-1);
for i=2:n
    for j=2:n
        D = norm(table2array(fluxes(:,i)) - table2array(fluxes(:,j)));
        distMat((i-1),(j-1)) = D;
    end
end
%save distance matrix
distMat = array2table(distMat);
distMat.Properties.VariableNames = fluxes.Properties.VariableNames(2:end);
distMat.Properties.RowNames = fluxes.Properties.VariableNames(2:end);
writetable(distMat,file3,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
    