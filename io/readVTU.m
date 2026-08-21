function [cell_struct, face_struct, V3, cells3D] = readVTU(filename)
% read a VTU mesh using Python/PyVista
%
% The reader rebuilds the global face structure with per-cell, per-face Python
% loops, which costs ~17 s on a 1e4-cell mesh and dominates a run. The result
% depends only on the file, so it is cached beside it as <name>.vemcache.mat,
% keyed on the file's size and modification time. Set VEM_NOCACHE=1 to bypass.

    cacheFile = [filename '.vemcache.mat'];
    useCache = isempty(getenv('VEM_NOCACHE'));

    key = cacheKey(filename);

    if useCache && isfile(cacheFile)

        try

            cached = load(cacheFile);

            if isfield(cached, 'key') && isequal(cached.key, key)

                cell_struct = cached.cell_struct;
                face_struct = cached.face_struct;
                V3 = cached.V3;
                cells3D = cached.cells3D;

                return;

            end

        catch
            % unreadable or stale cache: fall through and rebuild
        end

    end

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

    if useCache

        try
            save(cacheFile, 'key', 'cell_struct', 'face_struct', 'V3', 'cells3D', '-v7.3');
        catch
            % a read-only mesh directory is not a reason to fail the run
        end

    end

end

function key = cacheKey(filename)
% identity of the source file: size and modification time

    d = dir(filename);

    assert(~isempty(d), 'Mesh file not found: %s', filename);

    key = struct('bytes', d.bytes, 'datenum', d.datenum);

end

%% helper
function pythonExe = findPythonWithPyVista()
% find a Python executable with the required mesh-reader dependencies

    if ispc
        nullOutput = '> NUL 2>&1';
        nullError = '2> NUL';
    else
        nullOutput = '> /dev/null 2>&1';
        nullError = '2> /dev/null';
    end
    
    % find Conda
    homeDir = char(java.lang.System.getProperty('user.home'));
    
    if ispc
        condaCandidates = { ...
            'conda', ...
            fullfile(homeDir, 'miniconda3', 'Scripts', 'conda.exe'), ...
            fullfile(homeDir, 'anaconda3', 'Scripts', 'conda.exe'), ...
            fullfile(homeDir, 'miniconda', 'Scripts', 'conda.exe'), ...
            fullfile(homeDir, 'anaconda', 'Scripts', 'conda.exe')};
    else
        condaCandidates = { ...
            'conda', ...
            fullfile(homeDir, 'miniconda3', 'bin', 'conda'), ...
            fullfile(homeDir, 'anaconda3', 'bin', 'conda'), ...
            fullfile(homeDir, 'miniconda', 'bin', 'conda'), ...
            fullfile(homeDir, 'anaconda', 'bin', 'conda')};
    end
    
    condaInfo = '';
    
    for k = 1:numel(condaCandidates)
    
        conda = condaCandidates{k};
    
        command = sprintf('"%s" env list %s', conda, nullError);
        [status, out] = system(command);
    
        if status == 0
            condaInfo = out;
            break;
        end
    
    end
    
    % check VEM Conda environment first
    if ~isempty(condaInfo)
    
        lines = splitlines(condaInfo);
    
        for k = 1:numel(lines)
    
            line = strtrim(lines{k});
    
            if isempty(line) || startsWith(line, '#')
                continue;
            end
    
            tokens = regexp(line, '\s+', 'split');
    
            if ~strcmp(tokens{1}, 'vem')
                continue;
            end
    
            envPath = tokens{end};
    
            if ispc
                python = fullfile(envPath, 'python.exe');
            else
                python = fullfile(envPath, 'bin', 'python');
            end
    
            if isfile(python)
    
                command = sprintf('"%s" -c "import pyvista, numpy, scipy" %s', python, nullOutput);
                status = system(command);
    
                if status == 0
                    pythonExe = python;
                    return;
                end
    
            end
    
        end
    
    end
    
    % check Python available on PATH
    if ispc
        candidates = {'python', 'py'};
    else
        candidates = {'python3', 'python'};
    end
    
    for k = 1:numel(candidates)
    
        python = candidates{k};
    
        command = sprintf('"%s" -c "import pyvista, numpy, scipy" %s', python, nullOutput);
        status = system(command);
    
        if status == 0
            pythonExe = python;
            return;
        end
    
    end
    
    % check other Conda environments
    if ~isempty(condaInfo)
    
        lines = splitlines(condaInfo);
    
        for k = 1:numel(lines)
    
            line = strtrim(lines{k});
    
            if isempty(line) || startsWith(line, '#')
                continue;
            end
    
            tokens = regexp(line, '\s+', 'split');
            envPath = tokens{end};
    
            if ispc
                python = fullfile(envPath, 'python.exe');
            else
                python = fullfile(envPath, 'bin', 'python');
            end
    
            if ~isfile(python)
                continue;
            end
    
            command = sprintf('"%s" -c "import pyvista, numpy, scipy" %s', python, nullOutput);
    
            status = system(command);
    
            if status == 0
                pythonExe = python;
                return;
            end
    
        end
    
    end
    
    error('Could not find a Python environment with PyVista, NumPy and SciPy installed');

end