function startup
%STARTUP Configure the MATLAB path for this project.
%
%   This file is run automatically by MATLAB at start-up whenever MATLAB is
%   launched in (or with this folder as) the project root. It adds the
%   source folders that hold the reusable functions so that scripts in the
%   project root can call them without any per-script addpath directives.
%
%   To set this up as a persistent project, launch MATLAB from this folder,
%   or set this folder as the MATLAB start-up folder. The path additions are
%   not saved permanently (savepath is intentionally not called), so the
%   change is local to each MATLAB session.

    projectRoot = fileparts(mfilename('fullpath'));

    addpath(fullfile(projectRoot, 'src', 'dst'));
    addpath(fullfile(projectRoot, 'src', 'domains'));

    % Local user's chebfun
    addpath(fullfile(userpath, 'chebfun'));

end
