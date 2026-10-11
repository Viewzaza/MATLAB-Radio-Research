% startup.m - runs when MATLAB starts (this folder is the MATLAB user folder).
% Added for the FieldFox suite: the ITU-R P.618 digital maps used by the
% Satellite Communications Toolbox (p618PropagationLosses) for rain, gas,
% cloud and scintillation losses (downloaded from MathWorks, 337 MB).
p618Maps = fullfile(fileparts(mfilename('fullpath')), 'support', 'p618');
if isfolder(p618Maps), addpath(p618Maps); end
clear p618Maps
