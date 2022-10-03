candidates = readtable('../results/production_targets/spermine_targets/compatible_genes_results.txt','delimiter','\t');
%candidates = readtable('../results/spermine_targets/candidates_mech_validated.txt','delimiter','\t');
load('../ModelFiles/ecYeastGEM_batch.mat');

candidates.kcat = cell(height(candidates),1);
cd GECKO/geckomat/utilities    
for i=1:height(candidates)
    enzyme = candidates.enzymes{i};
    [kcat,rxnIdx,rxnName,MW] = getKcat(ecModel_batch,enzyme);
    if length(kcat)>1
        kcat = (unique(kcat));
        str = [];
        for j=1:length(kcat)
            substr = kcat(j);
            substr = num2str(substr);
            if j>1
                str = [str,';',substr];
            else
                str = [str,substr];
            end
        end
    else
        str = num2str(kcat);
    end
    %kcat = strjoin(kcat,';');
    candidates.kcat{i} = str;
end
cd ../../..
results = candidates(:,[1 2 3 4 20 5 6 7 8 12 19 13 16]);
writetable(results,'../results/spermine_targets.txt','delimiter','\t','QuoteStrings',false)