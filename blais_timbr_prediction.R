% This script is involved in generating biomarker predictions 
% based on gene expression changes using the new algorithm,
% TIMBR (Transcriptionally-Inferred Biomarker Response).
% 
% This script utilizes the COBRA toolbox, rat and human metabolic networks,
% and TIMBR reaction weights which represent reaction-level summaries 
% of gene expression changes. 
% Gene expression fold changes and false discovery-rate adjusted
% q-values (FDR) should be calculated in the R/Bioconductor
% programming environment prior to generating TIMBR reaction weights.

% Most data/code needed to reproduce TIMBR predictions are available 
% as Supplementary Tables/Datasets or on the ratcon GitHub website:
% www.github.com/edikblais/ratcon

% Raw gene expression microarray data can be obtained from Open TG-GATEs:
% Open Toxicogenomics Project-Genomics Assisted Toxicity Evaluation System
% http://www.ncbi.nlm.nih.gov/pubmed/25313160.
% http://toxico.nibio.go.jp/english/index.html
% Raw gene expression microarray data are also available from ArrayExpress:
% E-MTAB-797 - Transcription profiling by array of rat hepatocytes 
% treated with approximately 130 chemicals in vitro
% https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-797/
% E-MTAB-798 - Transcription profiling by array of human hepatocytes 
% treated with approximately 130 chemicals in vitro
% https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-798/

% TIMBR predictions can be reproduced by running the following scripts:
% 1. ncomm_blais_timbr_expression.R
%   Inputs raw expression data
%   Outputs gene expression changes
%   Helper functions: ncomm_helper.R
% 2. ncomm_blais_timbr_weights.R
%   Inputs gene expression changes, rat and human GPR rules
%   Outputs TIMBR reaction weights
%   Helper functions: ncomm_helper.R
% 3. ncomm_blais_timbr_predictions.m
%   Inputs TIMBR reaction weights, rat and human metabolic networks
%   Outputs raw TIMBR predictions
%   Helper functions: timbr.m and ncomm_blais_model2irrev.m
% 4. ncomm_blais_timbr_analysis.R
%   Inputs raw TIMBR predictions
%   Outputs normalized TIMBR production scores, manuscript figures
%   Helper functions: ncomm_helper.R

%% TIMBR software requirements: COBRA toolbox and linear programming solver

% COBRA (COnstraint-Based Reconstruction and Analysis) toolbox:
% We installed the COBRA toolbox version pre2.0.6 from GitHub on 9/23/15:
% https://github.com/opencobra/cobratoolbox

% We installed an academic version of GUROBI Optimizer 6.0 from: 
% http://www.gurobi.com/academia/academia-center

% COBRA toolbox installations include the GLPK solver which could be used
% instead of GUROBI (or IBM's CPLEX).

% Set paths to required software
% Specify path to MATLAB installation folder:
path_cobra_installation = 'C:/Users/edikb/Documents/MATLAB/cobra/';
% Specify path to GUROBI installation folder:
path_gurobi_installation = 'E:/gurobi604/win64/matlab/';
% Specify path to GUROBI license file (gurobi.lic):
path_gurobi_license = 'C:/Users/edikb/';
% Specify path to downlaoded Supplementary Tables/Datasets
path_manuscript_directory = 'E:/sync/Dropbox/2016_01_ncomm_submission/';
% Specify path to COBRA model files in .xls or .sbml format
path_model_directory = 'E:/sync/Dropbox/2016_01_ncomm_submission/';
% Specify output directory for TIMBR predictions
path_timbr_directory = [path_model_directory 'timbr_directory/'];
% Specify path to R/Biocondcutor working directory:
path_rcode_directory = 'E:/sync/Dropbox/2016_01_ncomm_submission/';

addpath(genpath(path_cobra_installation));
addpath(genpath(path_gurobi_installation));
addpath(path_gurobi_license);
addpath(path_model_directory);

% initCobraToolbox % this only needs to be run once after installation
changeCobraSolver('gurobi6'); % or glpk

%% Load human-specific and rat-specific metabolic models
% rno = rattus norvegicus
rno_cobra_file = [path_model_directory 'ncomm_blais_data_rno_cobra.xlsx'];
rno_cobra_load = ncomm_blais_xls2model(rno_cobra_file);

% hsa = homo sapiens
hsa_cobra_file = [path_model_directory 'ncomm_blais_data_hsa_cobra.xlsx'];
hsa_cobra_load = ncomm_blais_xls2model(hsa_cobra_file);

% Each model contains the same sets reactions initially.
% Reactions conisdered "off" have lb and ub equal to 0
off_in_rno = rno_cobra_load.rxns(...
    rno_cobra_load.lb == 0 & rno_cobra_load.ub == 0);
off_in_hsa = hsa_cobra_load.rxns(...
    hsa_cobra_load.lb == 0 & hsa_cobra_load.ub == 0);

% Note that some reactions present in HMR2 were 
% disabled in both models by default:
off_in_both = intersect(off_in_rno,off_in_hsa);

% Remove reactions disabled in each model and set the objective to biomass
rno_cobra = changeObjective(removeRxns(...
    rno_cobra_load,off_in_rno),'RCR99999');
hsa_cobra = changeObjective(removeRxns(...
    hsa_cobra_load,off_in_hsa),'RCR99999');

%% % Save these COBRA models since xls2model takes a while to import a model
save ncomm_blais_data_cobra.mat rno_cobra hsa_cobra;
writeCbToSBML(rno_cobra,[path_model_directory 'ncomm_blais_data_rno_sbml'])
writeCbToSBML(hsa_cobra,[path_model_directory 'ncomm_blais_data_hsa_sbml'])

%% Load COBRA models
load ncomm_blais_data_cobra.mat;

%% Test biomass production with COBRA models
% Because exchange fluxes were formulated as fmol / cell / hour and
% biomass was formulated in units of fmol / cell, 
% the growth rate is measured as growth per hour (1 / hour):
% Note: These values similar to but not equal to values 
% reported in the manuscript. By default, exchange reaction boundary 
% conditions are set to relaxed physiological constraints. 
% Strict physiological constraints are defined in Supplementary Table 4.

rno_growth_rate_cobra = optimizeCbModel(rno_cobra); 
rno_growth_rate = rno_growth_rate_cobra.f;
hsa_growth_rate_cobra = optimizeCbModel(hsa_cobra); 
hsa_growth_rate = hsa_growth_rate_cobra.f;
rno_doubling_time = log(2) / (rno_growth_rate); 
hsa_doubling_time = log(2) / (hsa_growth_rate); 
disp([num2str(rno_growth_rate), ' = growth rate of 0.0510 per hour'])
disp([num2str(hsa_growth_rate), ' = growth rate of 0.0435 per hour']) 
disp([num2str(rno_doubling_time), ' = doubling time of 13.6026 hours'])
disp([num2str(hsa_doubling_time), ' = doubling time of 15.9510 hours'])
% The difference in growth rates is attributed to slight differences in the
% biomass composition for each organism and not likely any underlying
% differences in the network structure.

%% Read SBML models
rno_sbml = readCbModel([path_model_directory ...
    'ncomm_blais_data_rno_sbml.xml']);
hsa_sbml = readCbModel([path_model_directory ...
    'ncomm_blais_data_hsa_sbml.xml']);

%% Test biomass production with SBML models
rno_growth_rate_sbml = optimizeCbModel(rno_sbml); 
hsa_growth_rate_sbml = optimizeCbModel(hsa_sbml); 
disp([num2str(rno_growth_rate_sbml.f), ...
    ' = growth rate of 0.0510 per hour'])
disp([num2str(hsa_growth_rate_sbml.f), ...
    ' = growth rate of 0.0435 per hour']) 
% % Optionally, you can use SBML models instead of COBRA models
% rno_cobra  = rno_sbml
% hsa_cobra  = hsa_sbml

%% Convert COBRA models into irreversible format
% Specify default production requirements for all TIMBR simulations
% Originally, we required that both models produced biomass at a
% growth rate of 1 / week with an ATP maintenance requirement.
% However, we found that it was too difficult to resolve 
% how much TIMBR production scores reflected the global demand for 
% biomass production vs biomarker production.

% obj_value_required = 1 / 7 / 24;
obj_value_required = 0;
% atp_value_required = 5000;
atp_value_required = 0;
% Set unconstrained reaction bounds to +/- 10^6 because the largest
% consumption rate is fairly close to the default 10^3 reaction bound.
default_bound = 1000000; % 10^6
rno_irrev_base = ncomm_blais_model2irrev(rno_cobra);
hsa_irrev_base = ncomm_blais_model2irrev(hsa_cobra);

default_rxns_rno = rno_irrev_base.rxns(rno_irrev_base.ub == 1000);
default_rxns_hsa = hsa_irrev_base.rxns(hsa_irrev_base.ub == 1000);
rno_irrev_default = changeRxnBounds(rno_irrev_base, ...
    default_rxns_rno, default_bound, 'u');
hsa_irrev_default = changeRxnBounds(hsa_irrev_base, ...
    default_rxns_hsa, default_bound, 'u');
rno_irrev_atp = changeRxnBounds(rno_irrev_default,...
    'RCR11017_f',atp_value_required,'l');
hsa_irrev_atp = changeRxnBounds(hsa_irrev_default,...
    'RCR11017_f',atp_value_required,'l');
rno_irrev_biomass = changeRxnBounds(rno_irrev_atp,...
    'RCR99999_f',obj_value_required,'l');
hsa_irrev_biomass = changeRxnBounds(hsa_irrev_atp,...
    'RCR99999_f',obj_value_required,'l');

% Note: as mentioned above, TIMBR predictions were not performed while 
% requiring biomass production.
rno_irrev = rno_irrev_biomass;
hsa_irrev = hsa_irrev_biomass;
%% Determine maxmimum flux through each exchange reaction for iRno
[rno_exchange, rno_demand] = findExcRxns(rno_irrev);
rno_production = rno_irrev.rxns((cellfun(@length, ...
    regexp(rno_irrev.rxns,'_r')) == 0) & rno_exchange);
rno_fva = zeros(length(rno_production),1);
for exchange_index = 1:length(rno_production)
    rno_consumption = setdiff(intersect(...
        regexprep(rno_production(exchange_index),'_f','_r'),...
        rno_irrev.rxns),...
        rno_production(exchange_index));
    if (~isempty(rno_consumption))
        rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_consumption,0,'b'),...
            rno_production(exchange_index)));
        rno_fva(exchange_index,1) = rno_production_fba.f;
    else 
        rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_production(exchange_index)));
        rno_fva(exchange_index,1) = rno_production_fba.f;
    end
    disp(rno_production_fba)
end

%% Determine maxmimum flux through each exchange reaction for iHsa
[hsa_exchange, hsa_demand] = findExcRxns(hsa_irrev);
hsa_production = hsa_irrev.rxns((cellfun(@length, ... 
    regexp(hsa_irrev.rxns,'_r')) == 0) & hsa_exchange);
hsa_fva = zeros(length(hsa_production),1);
for exchange_index = 1:length(hsa_production)
    hsa_consumption = setdiff(intersect(...
        regexprep(hsa_production(exchange_index),'_f','_r'),...
        hsa_irrev.rxns),...
        hsa_production(exchange_index));
    if (~isempty(hsa_consumption))
        hsa_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(changeRxnBounds(...
            hsa_irrev,hsa_production,default_bound,'u'),...
            hsa_consumption,0,'b'),...
            hsa_production(exchange_index)));
        hsa_fva(exchange_index,1) = hsa_production_fba.f;
    else 
        hsa_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(...
            hsa_irrev,hsa_production,default_bound,'u'),...
            hsa_production(exchange_index)));
        hsa_fva(exchange_index,1) = hsa_production_fba.f;
    end
    disp(hsa_production_fba)
end

%% Specify minimimum required flux for each exchange reaction

fba_threshold = 1e-4;

rno_production_ok = rno_fva > fba_threshold;
rno_production_id = rno_production(rno_production_ok,1);
rno_production_requirement = min(0.90 * rno_fva(rno_production_ok,1), 100);
rno_production_count = length(rno_production_id);

hsa_production_ok = hsa_fva > fba_threshold;
hsa_production_id = hsa_production(hsa_production_ok,1);
hsa_production_requirement = min(0.90 * hsa_fva(hsa_production_ok,1), 100);
hsa_production_count = length(hsa_production_id);

%% Load TIMBR reaction weights
% These files are generated by the R script:
% ncomm_blais_timbr_weights.R

rno_timbr_weights_load = readtable([path_rcode_directory, ...
    'ncomm_blais_timbr_weights_rno.txt'],'Delimiter','\t');

hsa_timbr_weights_load = readtable([path_rcode_directory, ...
    'ncomm_blais_timbr_weights_hsa.txt'],'Delimiter','\t');

[rno_model_has_weight, rno_model2weight_index] = ismember(...
    rno_irrev.rxns, rno_timbr_weights_load.rxn_irreversible);
[rno_weight_in_model, rno_weight2model_index] = ismember(...
    rno_timbr_weights_load.rxn_irreversible, rno_irrev.rxns);

if (all(rno_model_has_weight))
    % first 3 columns are organism_id, rxn_id, and rxn_irrerversible
    rno_timbr_weights = rno_timbr_weights_load(rno_model2weight_index,4:end);
    rno_timbr_id = rno_timbr_weights.Properties.VariableNames';
    rno_timbr_count = length(rno_timbr_id);

else
    warning('reaction weights not specified for all reactions in rno')
end

[hsa_model_has_weight, hsa_model2weight_index] = ismember(...
    hsa_irrev.rxns, hsa_timbr_weights_load.rxn_irreversible);
[hsa_weight_in_model, hsa_weight2model_index] = ismember(...
    hsa_timbr_weights_load.rxn_irreversible, hsa_irrev.rxns);

if (all(hsa_model_has_weight))
    % first 3 columns are organism_id, rxn_id, and rxn_irrerversible
    hsa_timbr_weights = hsa_timbr_weights_load(hsa_model2weight_index,4:end);
    hsa_timbr_id = hsa_timbr_weights.Properties.VariableNames';
    hsa_timbr_count = length(hsa_timbr_id);

else
    warning('reaction weights not specified for all reactions in hsa')
end

%% Run TIMBR algorithm and save results as .txt files
% Estimate the global network demand of producing each exchange metabolite
% under treatment and control conditions for various compounds. 
% Saved outputs will be analyzed in the R/Bioconductor script:
% ncomm_blais_timbr_analysis.R

for rno_timbr_index = 1:rno_timbr_count
    
    rno_timbr_network_demand = zeros(rno_timbr_count, 1);
    disp(rno_timbr_id{hsa_timbr_index});
    for rno_production_index = 1:rno_production_count
        disp(rno_production_id(rno_production_index,1));
        rno_timbr_network_demand(rno_production_index,1) = ...
            timbr(rno_irrev, ...
            rno_production_id(rno_production_index,1), ...
            rno_production_requirement(rno_production_index,1), ...
            rno_timbr_weights(:,rno_timbr_index));
    end
    rno_timbr_file_name = [path_timbr_directory ...
        'ncomm_blais_timbr_' rno_timbr_id{hsa_timbr_index} '.txt'];
    
    rno_timbr_table =  table(...
        rno_production_id,...
        rno_timbr_network_demand,...
        'VariableNames',{'timbr_id' 'timbr_value'});
    writetable(rno_timbr_table,rno_timbr_file_name,...
        'Delimiter','\t');
end

for hsa_timbr_index = 1:hsa_timbr_count
    hsa_timbr_network_demand = zeros(hsa_timbr_count, 1);
    disp(hsa_timbr_id{hsa_timbr_index});
    for hsa_production_index = 1:hsa_production_count
        disp(hsa_production_id(hsa_production_index,1));
        hsa_timbr_network_demand(hsa_production_index,1) = ...
            timbr(hsa_irrev, ...
            hsa_production_id(hsa_production_index,1), ...
            hsa_production_requirement(hsa_production_index,1), ...
            hsa_timbr_weights(:,hsa_timbr_index));
    end
    hsa_timbr_file_name = [path_timbr_directory ...
        'ncomm_blais_timbr_' hsa_timbr_id{hsa_timbr_index} '.txt'];
    
    hsa_timbr_table =  table(...
        hsa_production_id,...
        hsa_timbr_network_demand,...
        'VariableNames',{'timbr_id' 'timbr_value'});
    writetable(hsa_timbr_table,hsa_timbr_file_name,...
        'Delimiter','\t');
end
