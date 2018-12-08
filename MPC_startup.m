%% Setup the MATLAB path to contain the proper folders for the Model Predictive Control Toolbox

% Find out where the toolbox is
filename = [mfilename, '.m'];
scriptDir = which(filename);
scriptDir = strrep(scriptDir, filename, '');

% Add all subfolders to the path
%addpath( genpath(scriptDir) );

%% Add the Toolbox folder to the path
addpath( scriptDir );

%% Add the packages to the path
addpath( [scriptDir, filesep, 'matrices']);

addpath( genpath([scriptDir, filesep, 'optimization']) );
addpath( genpath([scriptDir, filesep, 'spectralProperties']) );

% Remove the testbench directory from the path
rmpath( genpath([scriptDir, filesep, 'testbenchs']) );

% Remove the variables used
clear filename scriptDir
