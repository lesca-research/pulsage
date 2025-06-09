function psa_install(mode, extra)
%PSA_INSTALL Installation of Pulsage plugin for Brainstorm (linux and windows).
%   PSA_INSTALL(MODE) install Pulsage processes to brainstorm user folder:
%     - linux: $HOME/.brainstorm/process
%     - windows: C:\Documents and Settings\username\.brainstorm
%   It backups existing scripts in brainstorm user folder.
%
%   IMPORTANT: assume brainstorm has been installed and its functions are
%   available in matlab's path.
%
%   If MODE == 'copy' then all processes are copied to the installation
%   folder. This is recommended for an end-user installation (no frequent
%   code updates). 
%   If MODE == 'link' (linux only) then symbolink links pointing to scripts 
%   in the  source folder are created in the installation folder. 
%   This is recommeded when source code needs to be updated often. 
%   Any modification in the source folder will be available in brainstorm.
%   IMPORTANT: for *new scripts*, installation has to be run again to create
%   new symbolic links.
%   PSA_INSTALL(MODE, EXTRA)
%   Install extra scripts that override some brainstorm functions.
%   WARNING: these are mostly in-dev features so they may be unstable.
%            They are strongly dependent on the current version of
%            brainstorm, so it's better to have the most up-to-date
%            version.
%            DO NOT install these scripts unless you know what you're doing ;)
%   EXTRA can be a string or cell of string corresponding to extra
%   installation scenarios. Files specified in MANIFEST.<extra> will then
%   be installed.
%
%  To cleanly uninstall pulsage, run psa_uninstall() (see psa_unsintall.m)
%
if nargin < 1           
   mode = 'copy'; 
end

if nargin < 2
    extra = {};
elseif ~iscellstr(extra)
    if ~ischar(extra)
        error('Argument "extra" must be a string');
    end
    extra = {extra};
end

if nargin < 4
    dry = 0;
end

%% Check Brainstorm installation
try
    bst_folder = bst_get('BrainstormUserDir');
catch
    msg = ['Could not find Brainstorm installation. '...
           'Check that matlab path contains Brainstorm folders'];
    throw(MException('Pulsage:Installation', msg));
end

bst_process_folder = fullfile(bst_folder, 'process');
if ~exist(bst_process_folder, 'dir')
    display(['Could not find Brainstorm process folder "' ...
             bst_process_folder '". Check brainstorm installation']);
    return;
end
addpath(fullfile(pwd, 'dist_tools'));
install_package('pulsage', 'bst_plugin', bst_process_folder, mode, extra, dry);
