clear
current = pwd;
%subSystems_GEM = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../data/chemicals_info.txt','Delimiter','\t');

d         = dir('../results');
isub      = [d(:).isdir]; %# returns logical vector
load('../ModelFiles/ecYeastGEM_batch.mat')
load('../ModelFiles/yeastGEM.mat')
GEM  = ravenCobraWrapper(model);
chemicals_info.Name = strtrim(chemicals_info.Name);
prod_capabilities  = table();
rxns = ecModel_batch.rxns;
fluxes = table(rxns);
families = [];
biomass_prod = false;
protFree = true;
Prot_cost =[]; 
%precursor mets
precursors = {'D-glucose 6-phosphate' 'D-fructose 6-phosphate' 'ribose-5-phosphate' ...
              'D-erythrose 4-phosphate' 'glyceraldehyde 3-phosphate' '3-phosphonato-D-glycerate(3-)' ...
              'phosphoenolpyruvate' 'pyruvate' 'acetyl-CoA' '2-oxoglutarate' ...
              'succinyl-CoA' 'oxaloacetate'};                  
compartments = [1 1 1 1 1 1 1 1 1 1 9 1];  
prec_mTO = table(precursors');
%Loop through all maxRate folders
mkdir('../results/production_capabilities/yieldPlots')

% metTOsWT = [];
% ecModel = changeMedia_batch(ecModel_batch,'D-glucose exchange (reversible)','Min');
% 
% [~,~,~,fluxDist,~,~] = calculate_potential(ecModel,growthPos,growthPos,CS_index,0.180,true);
% 
% for k=1:length(precursors)
%     midx = find(strcmpi(ecModel_batch.metNames,precursors{k}));
%     midx = midx(find(ecModel_batch.metComps(midx) == compartments(k)));
%     if ~isempty(midx)
%         rxns = find(ecModel_batch.S(midx,:));
%         [iA,rxns2] = ismember(ecModel_batch.rxns(rxns),fluxDist.rxns);
%         rxns = rxns(iA);
%         rxns2 = rxns2(find(rxns2));
%         %compute turnover numbers
%         metTO = 0.5*sum(abs(ecModel_batch.S(midx,rxns))*fluxDist.flux(rxns2));
%     else
%         disp(precursors{k})
%         pause
%     end
%     
%     metTOsWT = [metTOsWT;metTO];
%     
% end


for i=1:height(chemicals_info)
    compound = chemicals_info.Name{i};
    model    = [];
    %try to load GEM
    modelStr = ['../ModelFiles/production_ecModels/' lower(chemicals_info.ecModel{i})];
    try
        modelStr = strrep(modelStr,'.mat','_WBG.mat');
        load(modelStr);
    catch
        modelStr = strrep(modelStr,'_WBG.mat','.mat');
        try
            load(modelStr);
        catch
        end
    end
    
    if ~isempty(model)
        ecModel = model;
        %Check presence of objective reaction
        index_ec = find(model.c);
        indexProt= find(contains(model.rxnNames,'prot_pool_exchange'));
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
                %ecModel   = setParam(ecModel,'ub',indexProt,1000);
                [bioY_ec,proY_ec,rate_ec,fluxDist,cFlux_l,cFlux_h] = calculate_potential(ecModel,growthPos,index_ec,CS_index,0.180,biomass_prod);
                fluxDist_1 = fluxDist;
                enzUsages  = fluxDist(contains(fluxDist.rxns,'prot_'),:);
                %Get burden of CCM enzymes
                CCMratio   = get_CCM_enzBurden(fluxDist,ecModel);
                %Get flux distribution
                cost       = fluxDist.flux(strcmpi(fluxDist.rxns,ecModel.rxns(indexProt)));
                [presence,iB]   = ismember(fluxes.rxns,fluxDist.rxns);
                fluxes = fluxes(presence,:);
                fluxDist = fluxDist(iB(presence),:);
                str = strrep(chemicals_info.ecModel{i},'-','_');
                str = str(3:(end-3));
                str = regexprep(str,'[^a-zA-Z]','');
                str = strrep(str,',','_');
                str = strrep(str,'.mat','');
                str = strrep(str,'(','');
                str = strrep(str,')','');
                str = lower(str);
                eval(['fluxes.' str '=fluxDist.flux;'])    
                %Get metabolic precursors turnover
                metTOs = [];
                for k=1:length(precursors)
                    midx = find(strcmpi(ecModel_batch.metNames,precursors{k}));
                    midx = midx(find(ecModel_batch.metComps(midx) == compartments(k)));
                    if ~isempty(midx)
                        rxns = find(ecModel_batch.S(midx,:));
                        [iA,rxns2] = ismember(ecModel_batch.rxns(rxns),fluxDist_1.rxns);
                        rxns = rxns(iA);
                        rxns2 = rxns2(find(rxns2));
                        %compute turnover numbers
                        metTO = 0.5*sum(abs(ecModel_batch.S(midx,rxns))*fluxDist_1.flux(rxns2));
                    else
                        disp(precursors{k})
                    end 
                    metTOs = [metTOs;metTO];
                end
                eval(['prec_mTO.' str '=metTOs;'])
                %now with GEM
%                 if ~isempty(indexGEM)
%                     GEM  = changeMedia_Original(GEM,AA,1,1000);
%                     GEM  = setParam(GEM,'obj',indexGEM,1);
%                     sol2 = solveLP(GEM);
%                     obj2 = -sol2.f;
%                     FC   = obj2/obj1;
%                     if obj2>0
%                         CS_index  = find(strcmpi(GEM.rxnNames,'D-glucose exchange'));
%                         growthPos = find(strcmpi(GEM.rxnNames,'growth'));
%                         [bioY,proY,rate] = calculate_potential(GEM,growthPos,indexGEM,CS_index,0.180,biomass_prod);
%                     end
%                 else
%                     bioY =NaN;proY=NaN;rate=NaN;
%                 end
                %no prot case
                compound = chemicals_info.internal_ids{i};

                [bioY_np,proY_np,rate_np,~,~,~,protYield] = calculate_potential(ecModel,growthPos,index_ec,CS_index,0.180,biomass_prod,protFree);
                newRow            = [{compound},chemicals_info.Group(i),chemicals_info.class(i),chemicals_info.MW(i),bioY_np,proY_np,rate_np,bioY_ec,proY_ec,rate_ec,cFlux_l,cFlux_h,cost,CCMratio,protYield];
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
mkdir('../results/production_capabilities')
if biomass_prod
    file1 = '../results/production_capabilities/prodCapabilities_allChemicals_wBio.txt';
    file2 = '../results/production_capabilities/proteinLimitations_allChemicals_wBio.txt';
    file3 = '../results/fluxDist_distance_allChemicals_wBio.txt';
elseif protFree
    file1 = '../results/production_capabilities/prodCapabilities_allChemicals.txt';
    file2 = '../results/production_capabilities/proteinLimitations_allChemicals.txt';
    file3 = '../results/fluxDist_distance_allChemicals_noProt.txt';
    file4 = '../results/fluxDistributions_allChemicals.txt';
    file5 = '../results/met_precursors_turnovers_allChemicals.txt';
end
    
prod_capabilities.Properties.VariableNames = {'compound' 'type' 'family' 'MW' 'bioYield_gem' 'prodYield_gem' 'prodRate_gem' 'bioYield_ec' 'prodYield_ec' 'prodRate_ec' 'cFlux_l' 'cFlux_h' 'Pburden' 'CCMratio' 'protScaled'};
writetable(prod_capabilities,file1,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
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
writetable(fluxes,file4,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
writetable(prec_mTO,file5,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')

%writetable(distMat,file3,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
    