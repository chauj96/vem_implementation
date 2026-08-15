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

    if ispc
        candidates = {'python', 'py'};
        nullOutput = '> NUL 2>&1';
        nullError = '2> NUL';
    else
        candidates = {'python3', 'python'};
        nullOutput = '> /dev/null 2>&1';
        nullError = '2> /dev/null';
    end
    
    % check Python available on PATH
    for k = 1:numel(candidates)
    
        python = candidates{k};
    
        command = sprintf('"%s" -c "import pyvista, numpy, scipy" %s', python, nullOutput);
    
        status = system(command);
    
        if status == 0
            pythonExe = python;
            return;
        end
    
    end
    
    % find Conda
    homeDir = char(java.lang.System.getProperty('user.home'));
    
    if ispc
        condaCandidates = { ...
            'conda', ...
            fullfile(homeDir, 'miniconda3', 'Scripts', 'conda.exe'), ...
            fullfile(homeDir, 'anaconda3', 'Scripts', 'conda.exe')};
    else
        condaCandidates = { ...
            'conda', ...
            fullfile(homeDir, 'miniconda3', 'bin', 'conda'), ...
            fullfile(homeDir, 'anaconda3', 'bin', 'conda')};
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
    
    % check Conda environments
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