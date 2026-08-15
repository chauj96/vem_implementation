function [cell_struct, face_struct, V3, cells3D] = readVTU(filename)
% read a VTU mesh using Python/PyVista

    pythonExe = findPythonWithPyVista();

    rootDir = fileparts(fileparts(mfilename('fullpath')));
    pythonScript = fullfile(rootDir, 'meshes', 'mesh_reader.py');
    matfile = [tempname, '.mat'];

    command = sprintf('"%s" "%s" "%s" "%s"', pythonExe, pythonScript, filename, matfile);

    status = system(command);

    if status ~= 0
        error('Python mesh reader failed.');
    end

    data = load(matfile);

    cell_struct = data.cell_struct;
    face_struct = data.face_struct;
    V3 = data.V3;
    cells3D = data.cells3D;

    delete(matfile);

end

%% helper
function pythonExe = findPythonWithPyVista()
% find a Python executable with the required mesh-reader dependencies
%
% MATLAB started from the Finder/Dock inherits a minimal PATH that usually
% excludes Conda, so interpreters are located on the file system rather than
% through the shell. set VEM_PYTHON to bypass the search

    persistent cachedPython

    if ~isempty(cachedPython) && (isfile(cachedPython) || ~ismember(filesep, cachedPython))

        if hasDependencies(cachedPython)
            pythonExe = cachedPython;
            return;
        end

        cachedPython = '';

    end

    candidates = {};

    % explicit user override
    override = getenv('VEM_PYTHON');

    if ~isempty(override)
        candidates{end+1} = override;
    end

    homeDir = char(java.lang.System.getProperty('user.home'));

    % Conda / Mamba installation roots
    rootPatterns = { ...
        fullfile(homeDir, 'miniconda*'), ...
        fullfile(homeDir, 'anaconda*'), ...
        fullfile(homeDir, 'miniforge*'), ...
        fullfile(homeDir, 'mambaforge*'), ...
        fullfile(homeDir, 'opt', 'miniconda*'), ...
        fullfile(homeDir, 'opt', 'anaconda*'), ...
        '/opt/miniconda*', ...
        '/opt/anaconda*', ...
        '/opt/homebrew/Caskroom/miniconda/base'};

    installRoots = {};

    for k = 1:numel(rootPatterns)

        matches = dir(rootPatterns{k});

        for j = 1:numel(matches)

            if matches(j).isdir
                installRoots{end+1} = fullfile(matches(j).folder, matches(j).name); %#ok<AGROW>
            end

        end

    end

    % search order: vem env, base, PATH, other envs
    envCandidates = {};
    baseCandidates = {};
    otherCandidates = {};

    for k = 1:numel(installRoots)

        root = installRoots{k};

        baseCandidates{end+1} = interpreterIn(root); %#ok<AGROW>

        envs = dir(fullfile(root, 'envs', '*'));

        for j = 1:numel(envs)

            if ~envs(j).isdir || startsWith(envs(j).name, '.')
                continue;
            end

            python = interpreterIn(fullfile(envs(j).folder, envs(j).name));

            if strcmp(envs(j).name, 'vem')
                envCandidates{end+1} = python; %#ok<AGROW>
            else
                otherCandidates{end+1} = python; %#ok<AGROW>
            end

        end

    end

    if ispc
        pathCandidates = {'python', 'py'};
    else
        pathCandidates = {'python3', 'python'};
    end

    candidates = [candidates, envCandidates, baseCandidates, pathCandidates, otherCandidates];

    for k = 1:numel(candidates)

        python = candidates{k};

        if isempty(python)
            continue;
        end

        if ismember(filesep, python) && ~isfile(python)
            continue;
        end

        if hasDependencies(python)
            cachedPython = python;
            pythonExe = python;
            return;
        end

    end

    error('%s', sprintf(['Could not find a Python environment with PyVista, NumPy and SciPy installed.\n' ...
           'Create one with:  conda env create -f environment.yml\n' ...
           'Or point MATLAB at an existing interpreter with:\n' ...
           '  setenv(''VEM_PYTHON'', ''/full/path/to/python'')']));

end

function python = interpreterIn(envPath)
% interpreter inside a Conda environment directory

    if ispc
        python = fullfile(envPath, 'python.exe');
    else
        python = fullfile(envPath, 'bin', 'python');
    end

end

function tf = hasDependencies(python)
% check that the interpreter can import the mesh-reader dependencies

    if ispc
        nullOutput = '> NUL 2>&1';
    else
        nullOutput = '> /dev/null 2>&1';
    end

    command = sprintf('"%s" -c "import pyvista, numpy, scipy" %s', python, nullOutput);

    tf = system(command) == 0;

end
