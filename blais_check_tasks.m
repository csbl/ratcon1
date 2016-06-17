% This MATLAB script imports iRno and iHsa into the RAVEN toolbox and 
% performs the checkTasks function for 327 metabolic tasks described 
% in Supplementary Table 4. Most tasks capture metabolic functions that 
% are common to both rats and humans. Some tasks are meant to fail 
% in both organisms while others are meant to fail in only one organism, 
% consistent with known species-specific differences between 
% rat and human metabolism.

% Software requirements: RAVEN toolbox and a linear programming solver

% RAVEN (Reconstruction, Analysis, and Visualization of Metabolic Networks)
% We installed the RAVEN toolbox version 1.08 from the BIOMET website:
% http://biomet-toolbox.org/index.php?page=downtools-raven

% We installed an academic version of Mosek from: 
% http://mosek.com/resources/download/
% and obtained an academic license from:
% http://license.mosek.com/cgi-bin/student.py

% RAVEN toolbox installation does not include a pre-packaged solver.

%% Set paths to required software
% Specify path to Mosek installation folder:
path_mosek_installation = 'E:/mosek/7/toolbox/r2013a';
% Specify path to GUROBI license file (gurobi.lic):
path_mosek_license = 'C:/Users/edikb/';
% Specify path to MATLAB installation folder:
path_raven_installation = 'C:/Users/edikb/Documents/MATLAB/RAVEN/';
% Specify path to R/Biocondcutor working directory:
path_rcode_directory = 'E:/sync/Dropbox/2016_01_ncomm_submission/';
path_manuscript_directory = 'E:/sync/Dropbox/2016_01_ncomm_submission/';
path_model_directory = 'E:/sync/Dropbox/2016_01_ncomm_submission/';


%% RAVEN toolbox and set solver to Mosek
addpath(genpath(path_mosek_installation));
addpath(path_mosek_license);
addpath(genpath(path_raven_installation));
addpath(path_manuscript_directory);
addpath(path_model_directory);
%% Load human-specific and rat-specific metabolic models
% rno = rattus norvegicus
rno_raven_file = [path_model_directory 'ncomm_blais_data_rno_raven.xlsx'];
rno_raven_load = importExcelModel(rno_raven_file,false);
% hsa = homo sapiens
hsa_raven_file = [path_model_directory 'ncomm_blais_data_hsa_raven.xlsx'];
hsa_raven_load = importExcelModel(hsa_raven_file,false);

% note: expect warnings related to metabolites compositions, genes not associated
% with reactions, and metabolites never being. they do not affect the
% outcome of the metabolic tasks.

%% Run checkTasks with Supplementary Table 4
check_tasks_file = [path_manuscript_directory 'ncomm_blais_supplementary_table4.xls'];

disp('The following 12 rat-specific tasks should report: PASS (should fail):');
disp('[278,279,280,281,282,283,284,313,316,319,321,327]');
check_tasks_report_rno = checkTasks(rno_raven_load, check_tasks_file,true,true);

disp('The following 2 human-specific tasks should report: PASS (should fail):');
disp('[314,325]');
check_tasks_report_hsa = checkTasks(hsa_raven_load, check_tasks_file,true,true);

% No other tasks should be printed to the command line other than those
% explicitly mentioned.

%% End
