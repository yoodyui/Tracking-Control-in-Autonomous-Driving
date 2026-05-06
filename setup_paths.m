% SETUP_PATHS Add shared helper functions for the textbook examples.
%
% Run this script from the codes folder, or let each chapter example call it
% automatically. It keeps reusable model, geometry, and waypoint utilities in
% one shared location.

code_root = fileparts(mfilename('fullpath'));
addpath(fullfile(code_root, 'helpers'));
