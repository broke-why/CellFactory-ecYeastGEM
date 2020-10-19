cand_1_total  = [];
cand_1_OE  = [];
cand_1_del  = [];

cand_2_total  = [];
cand_2_OE  = [];
cand_2_del  = [];

cand_3_total  = [];
cand_3_OE  = [];
cand_3_del  = [];

del_targets = [];
OE_targets = [];

subSystems_del = [];
subSystems_OE = [];

chem_class_del = [];
chem_class_OE  = [];

models = [];
chemClass = [];
current = pwd;
%fileNames  = dir();
d = dir('../results');
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
load('../ModelFiles/mat/ecYeastGEM_batch.mat')
subSystems_model = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../ComplementaryData/chemicals_info.txt','Delimiter','\t');
for i=1:length(nameFolds)
    folder = nameFolds{i};
    cd (current)
    if contains(folder,'_targets')        
        models    = [models; {['ec' strrep(folder,'_targets','')]}];
        model_idx = find(strcmpi(chemicals_info.ecModel,['ec' strrep(folder,'_targets','') '.mat']));
        if ~isempty(model_idx)
            class = chemicals_info.class(model_idx);
        else
            class = {''};
        end
        chemClass = [chemClass; class];
        try
            cand_1       = readtable(['../results/' folder '/candidates_ecFSEOF.txt'],'Delimiter','\t');
            cand_1_del   = [cand_1_del;sum(cand_1.actions==0)];
            cand_1_OE    = [cand_1_OE;sum(cand_1.actions>0)];
            cand_1_total = [cand_1_total;height(cand_1)];

        catch
            cand_1_total = [cand_1_total;0];
            cand_1_del   = [cand_1_del;0];
            cand_1_OE    = [cand_1_OE;0];
        end
        
        try
            cand_2 = readtable(['../results/' folder '/candidates_mech_validated.txt'],'Delimiter','\t');
            cand_2_total = [cand_2_total;height(cand_2)];
            cand_2_del   = [cand_2_del;sum(cand_2.actions==0)];
            cand_2_OE    = [cand_2_OE;sum(cand_2.actions>0)];
        catch
            cand_2_total = [cand_2_total;0];
            cand_2_del   = [cand_2_del;0];
            cand_2_OE    = [cand_2_OE;0];
        end
              
        try
            cand_3 = readtable(['../results/' folder '/compatible_genes_results.txt'],'Delimiter','\t');
            cand_3_total = [cand_3_total;height(cand_3)];
            
            cand_3_del     = [cand_3_del;sum(contains(cand_3.actions,'deletion'))];
            deletions      = cand_3.enzymes(contains(cand_3.actions,'deletion'));
            subSystems_del = [subSystems_del;mapEnzymeSubSystems(deletions,ecModel_batch)];
            del_targets    = [del_targets;cand_3.shortNames(contains(cand_3.actions,'deletion'))];
            temp           = repelem(class,sum(contains(cand_3.actions,'deletion')));
            chem_class_del = [chem_class_del;temp'];
            
            cand_3_OE     = [cand_3_OE;sum(contains(cand_3.actions,'OE'))];
            OE_targets    = [OE_targets;cand_3.shortNames(contains(cand_3.actions,'OE'))];
            temp          = repelem(class,sum(contains(cand_3.actions,'OE')));
            chem_class_OE = [chem_class_OE;temp'];
            OEs           = cand_3.enzymes(contains(cand_3.actions,'OE'));
            subSystems_OE = [subSystems_OE;mapEnzymeSubSystems(OEs,ecModel_batch)];
        catch
            cand_3_total = [cand_3_total;0];
            cand_3_del   = [cand_3_del;0];
            cand_3_OE    = [cand_3_OE;0];
        end    
    end
    
end
models = strrep(models,'_targets','');
t = table(models,chemClass,cand_1_total,cand_1_OE,cand_1_del,cand_2_total,cand_2_OE,cand_2_del,cand_3_total,cand_3_OE,cand_3_del);
writetable(t,'../results/targets_summary.txt','delimiter','\t');

t = table(del_targets,subSystems_del,chem_class_del);
writetable(t,'../results/all_deletions.txt','delimiter','\t');

t = table(OE_targets,subSystems_OE,chem_class_OE);
writetable(t,'../results/all_OEs.txt','delimiter','\t');
