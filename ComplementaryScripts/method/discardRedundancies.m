function candidates = discardRedundancies(model,candidates)
tempModel = model;
threshold = 1E-12;
%Iterate through each priority level
for priority = [1 2 3]
    disp(['*****Priority level: ' num2str(priority)])
    %Isolate deletion targets for the ith priority level
    targets  = candidates(candidates.priority==priority & ~strcmpi(candidates.actions,'OE'),:);
    if ~isempty(targets)
        %Rank targets by mean(pYield FC, pRate FC)
        %targets.FC = mean([targets.foldChange_yield targets.foldChange_pRate],2);
        %targets    = sortrows(targets,'FC','descend');
        %Remove each of the remaining deletion targets from model
        i=1;
        while i<=height(targets)
            protein = targets.enzymes{i};
            name = targets.shortNames(i);
            %Check if gene still remains in candidates list
            [presence,iB] = ismember(protein,candidates.enzymes);
            if presence
                disp(['Checking protein: ' name{1}])
                %delete gene
                tempModel = removeGenes(model,candidates.genes(iB));
                remaining = candidates(:,[1 2 3]);
                idx = find(strcmpi(remaining.enzymes,protein));
                remaining(idx,:) = [];
                %check if the remaining genes can still carry flux
                for j = 1:height(remaining)
                    enzyme = remaining.enzymes{j};
                    name = remaining.shortNames{j};
                    enzRxn = ['draw_prot_' enzyme];
                    idx = find(strcmpi(tempModel.rxnNames,enzRxn));
                    if ~isempty(idx) & ismember(enzyme,candidates.enzymes)
                        flux = haveFlux(tempModel,threshold,idx);
                        %REmove those enzymes that cannot carry flux
                        %anymore from the deletion targets
                        if ~flux
                            position = find(strcmpi(candidates.enzymes,enzyme));
                            if ~isempty(position) & ~strcmpi(protein,enzyme)
                                disp([' discarding ' name])
                                candidates(position,:) = [];
                                %remaining(j,:) = [];
                            end
                        end
                    end
                end
            end
            i = i+1;
        end
        
        
    end
end
end