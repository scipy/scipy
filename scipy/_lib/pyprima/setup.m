function setup(varargin)
%SETUP sets the package up for MATLAB.
%
%   This script can be called in the following ways.
%
%   setup  % Compile all the solvers with the default options
%   setup(options)  % Compile all the solvers with `options`
%   setup path  % Add the paths needed to use the package
%   setup clean  % Remove all the compiled MEX files
%   setup uninstall  % Uninstall the package
%
%   In addition, the script can also be used as follows to compile a solver
%   specified by a string `solver_name`. Note, however, this is only intended
%   for developers. End users should NEVER do it. After doing this, any solver
%   other than the one being compiled WILL BECOME UNUSABLE.
%
%   setup(solver_name)  % Compile a solver with the default options
%   setup(solver_name, options)  % Compile a solver with `options`
%
%   Possible options for compilation:
%   - half: whether to compile the half precision of the Fortran solvers (default: false)
%   - single: whether to compile the single precision of the Fortran solvers (default: true)
%   - quadruple: whether to compile the quadruple precision of the Fortran solvers (default: false)
%   - classical: whether to compile the classical variant of the Fortran solvers (default: true)
%   - debug: whether to compile the debugging version of the Fortran solvers (default: false)
%   - debug_only: whether to compile only the debugging version (default: false);
%     `debug_only` prevails if both `debug` and `debug_only` are present
%   - verbose: whether to be verbose during the compilation (default: false)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   TODO: When the MATLAB implementation is available, we should support
%   setting options.fortran = false and skipping the compilation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REMARKS:
%
%   1. Since MEX is the standard way of calling Fortran code in MATLAB, you
%   need to have MEX properly configured for compile Fortran before using
%   the package. It is out of the scope of this package to help the users
%   to configure MEX.
%
%   2. To run this script, you need to have write access to the directory that
%   contains this script and its subdirectories.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University.
%
%   Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
%
%   Started in July 2020.
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup starts

% Name of the package. It will be used as a stamp to be included in the path_string. Needed only
% if `savepath` fails.
package_name = 'prima';

% Check the version of MATLAB.
if verLessThan('matlab', '9.4')   % MATLAB R2018a = MATLAB 9.4
    fprintf('\nSorry, this package does not support MATLAB R2017b or earlier releases.\n\n');
    return
end

% The full paths to several directories needed for the setup.
cpwd = pwd();  % Current directory, which may not be the directory containing this setup script
setup_dir = fileparts(mfilename('fullpath')); % The directory containing this setup script
fortd = fullfile(setup_dir, 'fortran'); % Directory of the Fortran source code
matd = fullfile(setup_dir, 'matlab'); % Matlab directory
gateways = fullfile(matd, 'mex_gateways'); % Directory of the MEX gateway files
interfaces = fullfile(matd, 'interfaces'); % Directory of the interfaces
mexdir = fullfile(interfaces, 'private'); % The private subdirectory of the interfaces
tests = fullfile(matd, 'tests'); % Directory containing some tests
tools = fullfile(matd, 'setup_tools'); % Directory containing some tools, e.g., interform.m
fortd_interform = fullfile(fortd, '.interform'); % Directory of the intersection-form Fortran code
gateways_interform = fullfile(gateways, '.interform'); % Directory of the intersection-form MEX gateways

% We need write access to `setup_dir` (and its subdirectories). Return if we do not have it.
% N.B.: This checking is NOT perfect because of the following --- but it is better than nothing.
% 1. `fileattrib` may not reflect the attributes correctly, particularly on Windows. See
% https://www.mathworks.com/matlabcentral/answers/296657-how-can-i-check-if-i-have-read-or-write-access-to-a-directory
% 2. Even if we have write access to `setup_dir`, we may not have the same access to its subdirectories.
[~, attribute] = fileattrib(setup_dir);
if ~attribute.UserWrite
    fprintf('\nSorry, we cannot continue because we do not have write access to\n\n%s\n\n', setup_dir);
    return
end

% `tools` contains some functions needed in the sequel.
addpath(tools);

% Parse the input.
[solver_list, options, action, wrong_input] = parse_input(varargin);

% Exit if wrong input detected. Error messages have been printed during the parsing.
if wrong_input
    rmpath(tools);
    error('prima:InvalidInput', 'setup: The input is invalid.');
end

% Remove the compiled MEX files (and the files generated along side) if requested.
if strcmp(action, 'clean')
    clean_mex(mexdir, 'verbose');
    clean_generated_files(fortd_interform, gateways_interform, tools, mexdir);
    rmpath(tools);
    return
end

% Uninstall the package if requested.
if strcmp(action, 'uninstall')
    uninstall_prima(package_name);
    % `rmpath(tools)` is not needed, as `uninstall_prima` already does it.
    return
end

% Add the path and return if requested.
if strcmp(action, 'path')
    add_save_path(interfaces, package_name);
    % Create `all_precisions.m` and `all_variants.m` under `tools` according to the content of
    % `mexdir`. They reflect the precisions ('half', 'single', 'double', 'quadruple') and variants
    % ('modern', 'classical') of the solvers available under `mexdir`.
    create_all_precisions(mexdir);
    create_all_variants(mexdir);
    % Some tools are shared between the setup and the runtime. They are all ready now. Copy them
    % from `tools` (the setup directory) to `mexdir` (the runtime directory).
    copy_shared_tools(tools, mexdir);
    rmpath(tools);
    fprintf('\nPath added.\n\n')
    return
end


%%%%%%%%%%%%%%% If we arrive here, then the user requests us to compile the solvers. %%%%%%%%%%%%%%%


% All previously compiled solvers are removed by the following lines.
clean_mex(mexdir);
clean_generated_files(fortd_interform, gateways_interform, tools, mexdir);


% Create `all_precisions.m` and `all_variants.m` under `tools` according to `options`.
% They reflect the precisions ('half', 'single', 'double', 'quadruple') and variants ('modern','classical')
% of the solvers available after the compilation.
create_all_precisions(options);
create_all_variants(options);
% Some tools are shared between the setup and the runtime. They are all ready now. Copy them from
% from `tools` (the setup directory) to `mexdir` (the runtime directory).
copy_shared_tools(tools, mexdir);

% Check whether MEX is properly configured.
fprintf('\nVerifying the setup of MEX ... \n');
language = 'Fortran'; % Language to compile
mex_well_conf = try_mex_setup(language);
if mex_well_conf == 0
    fprintf('\nMATLAB needs you to set MEX up for Fortran.');
    fprintf('\nTry ''help mex'' or ask a MATLAB expert about MEX setup.');
    fprintf('\nNote: MEX is a product of MathWorks. Its configuration is not part of this package.\n\n');
    return

elseif mex_well_conf == -1
    fprintf('\nmex(''-setup'', ''%s'') runs successfully but we cannot verify that MEX works properly.', language);
    fprintf('\nWe will try to continue.\n\n');
else
    fprintf('\nMEX is correctly set up.\n\n');
end

% Generate the intersection-form Fortran source code
% We need to do this because MEX accepts only the (obsolescent) fixed-form Fortran code on Windows.
% Intersection-form Fortran code can be compiled both as free form and as fixed form.
fprintf('Refactoring the Fortran code ... ');
interform(fortd);
interform(gateways);
fprintf('Done.\n\n');

% Compilation starts
fprintf('Compilation starts. It may take some time ...\n');
% Change directory to mexdir. All the intermediate files produced by the compilation (e.g., .mod)
% will be dumped to this directory. They will be removed when the compilation finishes.
cd(mexdir);
exception = [];
try
    compile(solver_list, mexdir, fortd_interform, gateways_interform, options);
catch exception
    % Do nothing for the moment.
end

%% Remove the intersection-form Fortran files unless we are debugging.
%if ~debug_flag
%    rmdir(fortd_interform, 's');
%    rmdir(gateways_interform, 's');
%end

cd(cpwd); % Change directory back to cpwd

if ~isempty(exception)
    rethrow(exception);  % Rethrow any exception caught during the compilation.
end

% Compilation ends successfully if we arrive here.
fprintf('Package compiled successfully!\n');

% Check whether all the solvers are compiled. It may not be the case if the user has specified
% a particular solver to compile. This must be done before removing `tools` from the MATLAB path.
all_solvers_compiled = isempty(setdiff(all_solvers(), solver_list));

% Remove `tools` from the MATLAB path. This should be done before calling `add_save_path`.
rmpath(tools);

% Add `interfaces` to the MATLAB path, and then try saving the path.
path_saved = add_save_path(interfaces, package_name);

fprintf('\nThe package is ready to use.\n');
fprintf('\nYou may now try ''help prima'' for information on the usage of the package.\n');

if all_solvers_compiled
    addpath(tests);
    fprintf('\nYou may also run ''testprima'' to test the package on a few examples.\n');
end

if ~path_saved  % `add_save_path` failed to save the path.
    add_path_string = sprintf('addpath(''%s'');', interfaces);
    fprintf('\n***** To use the package in other MATLAB sessions, do ONE of the following. *****\n');
    fprintf('\n- Append the following line to your startup script');
    fprintf('\n  (see https://www.mathworks.com/help/matlab/ref/startup.html for information):\n');
    fprintf('\n    %s\n', add_path_string);
    fprintf('\n- OR come to the current directory and run ''setup path'' when you need the package.\n');
end

fprintf('\n');

% setup ends
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function path_saved = add_save_path(path_string, path_string_stamp)
%ADD_SAVE_PATH adds the path indicated by PATH_STRING to the MATLAB path and then tries saving path.
% PATH_STRING_STAMP is a stamp used when writing PATH_STRING to the user's startup.m file, which is
% needed only if `savepath` fails.
% N.B.: Why not putting this function as an individual file in the `tools` directory? Because we
% need it after `tools` is removed from the path.

if nargin < 2
    path_string_stamp = sprintf('Added by %s', mfilename);
end

if ~exist(path_string, 'dir')
    error('prima:PathNotExist', 'The string %s does not correspond to an existing directory.', path_string);
end

addpath(path_string);

% Try saving the path in the system path-defining file at sys_pathdef. If the user does not have
% writing permission for this file, then the path will not saved.
% N.B. Do not save the path to the pathdef.m file under userpath. This file is not loaded by default
% at startup. See
% https://www.mathworks.com/matlabcentral/answers/269482-is-userdata-pathdef-m-for-local-path-additions-supported-on-linux
orig_warning_state = warning;
warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path
sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
path_saved = (savepath(sys_pathdef) == 0);
warning(orig_warning_state); % Restore the behavior of displaying warnings

% If path not saved, try editing the startup.m of this user. Do this only if userpath is nonempty.
% Otherwise, we will only get a startup.m in the current directory, which will not be executed
% when MATLAB starts from other directories. On Linux, the default value of userpath is
% ~/Documents/MATLAB, but it will be '' if this directory does not exist. We refrain from creating
% this directory in that case.
if ~path_saved && numel(userpath) > 0
    user_startup = fullfile(userpath, 'startup.m');
    add_path_string = sprintf('addpath(''%s'');', path_string);
    full_add_path_string = sprintf('%s  %s %s', add_path_string, '%', path_string_stamp);

    % First, check whether full_add_path_string already exists in user_startup or not.
    if exist(user_startup, 'file')
        startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
        path_saved = any(strcmp(startup_text_cells, full_add_path_string));
    end

    if ~path_saved
        % We first check whether the last line of the user startup script is an empty line (or the
        % file is empty or even does not exist at all). If yes, we do not need to put a line break
        % before the path adding string.
        if exist(user_startup, 'file')
            startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
            last_line_empty = isempty(startup_text_cells) || (isempty(startup_text_cells{end}) && ...
                isempty(startup_text_cells{max(1, end-1)}));
        else
            last_line_empty = true;
        end

        file_id = fopen(user_startup, 'a');  % Open/create file for writing. Append data to the end.
        if file_id ~= -1 % If FOPEN cannot open the file, it returns -1; We keep silent if it fails.
            if ~last_line_empty  % The last line of user_startup is not empty
                fprintf(file_id, '\n');  % Add a new empty line
            end
            fprintf(file_id, '%s', full_add_path_string);
            fclose(file_id);
            % Check that full_add_path_string is indeed added to user_startup.
            if exist(user_startup, 'file')
                startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                path_saved = any(strcmp(startup_text_cells, full_add_path_string));
            end
        end
    end
end

% add_save_path ends
return
