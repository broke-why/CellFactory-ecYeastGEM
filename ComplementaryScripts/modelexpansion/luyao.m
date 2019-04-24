% The theoretical yield analysis for different pathways of nonnative
% product

% Add all possible reactions to yeastGEM model
addnewpathway_RxnToGEM

% The theoretical yield analysis
% Set constraints
model=changeRxnBounds(model,'r_1714',-5,'l') % maximum glucose uptake rate
model=changeRxnBounds(model,'r_2111',0.1,'l') % minimal growth rate
changeObjective(model,'r_xxxx') % r_xxxx: RxnID of new product's exchange reaction

% Block other added pathways by setting the 'ub' of 0 when analysing a specific
% pathway
FBAsolution=optimizeCbModel(model)




