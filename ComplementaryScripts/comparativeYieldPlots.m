clear
current = pwd;
%subSystems_model = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../ComplementaryData/chemicals_info.txt','Delimiter','\t');
d         = dir('../results');
isub      = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
load('../ModelFiles/mat/ecYeastGEM_batch.mat')
load('../ModelFiles/mat/yeastGEM.mat')
model = ravenCobraWrapper(model);
subSystems_model = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
models = [];
chemicals_info.Name = strtrim(chemicals_info.Name);
maxYield = table([],[],[],[]);
%Loop through all maxYield folders
mkdir(['../results/yieldPlots'])
for i=1:length(nameFolds)
    folder = nameFolds{i};
    cd (current)
    if contains(folder,'_targets')
        compound  = strrep(folder,'_targets','');
        models    = [models; {['ec' compound]}];
        model_idx = find(strcmpi(chemicals_info.ecModel,['ec' strrep(folder,'_targets','') '.mat']));
        if ~isempty(model_idx)
            if strcmpi(chemicals_info.Group{model_idx},'native')
                met   = chemicals_info.Name{model_idx};
                MWeight = chemicals_info.MW(model_idx)/1000;
                index = find(strcmpi(ecModel_batch.metNames,met));
                %If model refers to a native compound
                if ~isempty(index)
                    %Check presence of exchange reaction
                    str      = [met ' exchange'];
                    index_ec = find(strcmpi(ecModel_batch.rxnNames,str));
                    index    = find(strcmpi(model.rxnNames,str));
                    if ~isempty(index) & ~isempty(index_ec)
                        c_source = 'D-glucose exchange';
                        ecModel  = setParam(ecModel_batch,'obj',index_ec,1);
                        ecModel  = changeMedia_batch(ecModel,[c_source ' (reversible)'],'Min',false,1);
                        index    = find(strcmpi(model.rxnNames,str));
                        model    = setParam(model,'obj',index,1);
                        model    = changeMedia_Original(model,1,1,1);
                        sol1 = solveLP(ecModel);
                        sol2 = solveLP(model);
                        FC = -sol2.f/-sol1.f;
                        maxYield = [maxYield;{met} {-sol1.f} {-sol2.f} {FC}];
                        
                        [BioYield_ec,yield_ec] = getYieldPlot(ecModel,index_ec,1,MWeight);
                        [BioYield,yield] = getYieldPlot(model,index,1,MWeight);
                        figure
                        xStr = 'Biomass yield [g_{biomass}/g_{glucose}]';
                        yStr = 'Product yield [g_{product}/g_{glucose}]';
                        plot(BioYield,yield,BioYield_ec,yield_ec,'LineWidth',5,'Color',[0.8 0.4 0.1])
                        hold on
                        plot(BioYield_ec,yield_ec,'LineWidth',5,'Color',[0.1 0 0.8])
                        xlabel(xStr)
                        ylabel(yStr)
                        set(gca,'FontSize',22)
                        saveas(gcf,['../results/yieldPlots/'  met '_lowGlc_yieldPlot.jpg'])
                        hold off
                        close all
                        
                        [BioYield_ec,yield_ec] = getYieldPlot(ecModel,index_ec,20,MWeight);
                        [BioYield,yield] = getYieldPlot(model,index,20,MWeight);
                        figure
                        set(gca,'FontSize',22)
                        xStr = 'Biomass yield [g_{biomass}/g_{glucose}]';
                        yStr = 'Product yield [g_{product}/g_{glucose}]';
                        plot(BioYield,yield,BioYield_ec,yield_ec,'LineWidth',5,'Color',[0.8 0.4 0.1])
                        hold on
                        plot(BioYield_ec,yield_ec,'LineWidth',5,'Color',[0.1 0 0.8])
                        xlabel(xStr)
                        ylabel(yStr)
                        set(gca,'FontSize',22)
                        saveas(gcf,['../results/yieldPlots/'  met '_highGlc_yieldPlot.jpg'])
                        hold off
                        close all                         
                    end
                end
            end
                
        else
            %disp(folder)
            class = {''};
        end
    end
end



