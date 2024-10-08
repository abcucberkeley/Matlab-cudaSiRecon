function siReconDataWrapper(dataPaths,ChannelPatterns, configFiles,varargin)

% Remember to carefully check your settings here and in your config files

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @iscell);
ip.addRequired('ChannelPatterns', @iscell);
ip.addRequired('configFiles', @iscell);

% PSF files or OTF files must have the same size as channel patterns
ip.addParameter('psfFiles', {''}, @iscell);
ip.addParameter('otfFiles', {''}, @iscell);

ip.addParameter('DeskewData', false, @islogical);
ip.addParameter('DeskewPsfs', false, @islogical);
ip.addParameter('DeskewCpusPerTask', 8, @isnumeric);
ip.addParameter('np', 5, @isnumeric); % number of phases
ip.addParameter('ndir', 3, @isnumeric); % number of phases
ip.addParameter('angle', 32.45, @isnumeric);
ip.addParameter('dz', 0.4, @isnumeric);
ip.addParameter('xyPixelSize', .098, @isnumeric);

% makeOTF Parameters
ip.addParameter('moAngle', 1.57, @isnumeric);
ip.addParameter('background', 100, @isnumeric);
ip.addParameter('zres', .2146, @isnumeric);

ip.addParameter('outFol', 'GPUsirecon', @ischar);
ip.addParameter('OverwriteSIrecon', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('overlap', 128, @isnumeric);
ip.addParameter('chunkSize', [512,0,0], @isnumeric);
ip.addParameter('SlurmParam', '', @ischar); % Slurm Parameter. Default below.
ip.addParameter('cpusPerTask', 5, @isnumeric); % Number of cores to request from Slurm.
ip.addParameter('memNeeded', 120, @isnumeric); % Memory to request from slurm. In GB.

ip.parse(dataPaths, ChannelPatterns, configFiles, varargin{:});

pr = ip.Results;

psfFiles = pr.psfFiles;
otfFiles = pr.otfFiles;

DeskewData = pr.DeskewData;
DeskewPsfs = pr.DeskewPsfs;
DeskewCpusPerTask = pr.DeskewCpusPerTask;
np = pr.np;
ndir = pr.ndir;
angle = pr.angle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;

moAngle = pr.moAngle;
background = pr.background;
zres = pr.zres;

outFol = pr.outFol;
parseCluster = pr.parseCluster;
overlap = pr.overlap;
chunkSize = pr.chunkSize;
SlurmParam = pr.SlurmParam;
memNeeded = pr.memNeeded;
cpusPerTask = pr.cpusPerTask;
OverwriteSIrecon = pr.OverwriteSIrecon;

% Number of config files must be the same as the number of channel patterns
if(numel(configFiles) ~= numel(ChannelPatterns))
    disp('Number of Config Files must be the same as the number of Channel Patterns.');
    return;
end

% PSF Files or OTF files need to be provided
if(numel(psfFiles{1}) == 0 && numel(otfFiles{1}) == 0)
    disp('This function requires either otf files or psf files.');
    return;
end

genOTFs = false;

% Check if we need to generate otf files or not
if(numel(psfFiles{1}) > 0)
    if(numel(psfFiles) ~= numel(ChannelPatterns))
        disp('Number of PSF Files must be the same as the number of Channel Patterns.');
        return;
    end
    genOTFs = true;
else
    if(numel(otfFiles) ~= numel(ChannelPatterns))
        disp('Number of OTF Files must be the same as the number of Channel Patterns.');
        return;
    end
end

% Default SlurmParam
if isempty(SlurmParam)
    SlurmParam = ['-p abc_a100 --qos abc_high -n1 --mem=' sprintf('%d',memNeeded) 'G --gres=gpu:a100:1'];
end

%dataPaths = {'/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/DS/'};
%{
dataPaths = {'/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/PSFs/488_NA0p4_sig0p1_highSN/DS/', ...
'/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/PSFs/488_NA0p46_sig0p1_512xPix_wait4ms_1p65msInt/DS/', ...
'/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/PSFs/488_NA0p46_sig0p1_512xPix_noWait_2p5msInt/DS/', ...
'/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/PSFs/488_NA0p46_sig0p1_320xPix_wait4ms_2p6Int/DS', ...
'/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/PSFs/488_NA0p46_sig0p1_320xPix_noWait_3p4msInt/DS/'};
%}
%ChannelPatterns = {'CamB'};
%parseCluster = true;
%outFol = 'GPUsirecon';



% if Deskew is needed then Deskew images first
if DeskewData || DeskewPsfs
    nFiles = 0;

    if DeskewData
        for i = 1:numel(dataPaths)
            for cPatt = 1:numel(ChannelPatterns)
                fnames = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);
                fnames = {fnames.name};
                nFiles = nFiles + numel(fnames);
            end
        end
    end

    % Check to see if any channel patterns are actually found
    if nFiles == 0
        ME = MException('dataInput:noFiles', ...
            'No files matching any of the given channel patterns were found');
        throw(ME);
    end


    if DeskewPsfs
        nFiles = nFiles + numel(psfFiles);
    end

    inputFullpaths = cell(nFiles, 1);
    outputFullpaths = cell(nFiles, 1);
    funcStrs = cell(nFiles, 1);
    curr = 1;

    if DeskewData
        for i = 1:numel(dataPaths)
            jobLogDir = dataPaths{i};
            dsFol = [dataPaths{i} '/DS/'];
            if  ~exist(dsFol, 'dir')
                mkdir(dsFol);
                fileattrib(dsFol, '+w', 'g');
            end
            for cPatt = 1:numel(ChannelPatterns)
                fnames = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);
                fnames = {fnames.name};
                if parseCluster
                    if  ~exist(jobLogDir, 'dir')
                        warning('The job log directory does not exist, use %s/job_logs as job log directory.', dataPaths{i})
                        jobLogDir = sprintf('%s/job_logs', dataPaths{i});
                        mkdir(jobLogDir);
                        fileattrib(jobLogDir, '+w', 'g');
                    end
                    job_log_fname = [jobLogDir, '/job_%A_%a.out'];
                    job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
                end
                for j = 1: numel(fnames)
                    [pathstr, fsname, ext] = fileparts(fnames{j});
                    dataFullpath = [dataPaths{i} filesep fnames{j}];
                    dataReconFullpath = [dataPaths{i} '/DS/' fsname ext];
                    inputFullpaths{curr} = dataFullpath;
                    outputFullpaths{curr} = dataReconFullpath;

                    funcStrs{curr} =  sprintf('cd /clusterfs/nvme/matthewmueller/Matlab-cudaSiRecon/src/cudaSirecon/;deskewPhasesFrame(''%s'',%.5f,%.5f,''SkewAngle'',%.5f,''nphases'',%.5f)', dataFullpath,xyPixelSize,dz,angle,np);
                    curr = curr+1;
                end
            end
        end
    end

    if DeskewPsfs
        for i = 1:numel(psfFiles)
            [pathstr, fsname, ext] = fileparts(psfFiles{i});
            jobLogDir = pathstr;
            dsFol = [pathstr '/DS/'];
            if  ~exist(dsFol, 'dir')
                mkdir(dsFol);
                fileattrib(dsFol, '+w', 'g');
            end
            if parseCluster
                if  ~exist(jobLogDir, 'dir')
                    warning('The job log directory does not exist, use %s/job_logs as job log directory.', dataPaths{i})
                    jobLogDir = sprintf('%s/job_logs', pathstr);
                    mkdir(jobLogDir);
                    fileattrib(jobLogDir, '+w', 'g');
                end
                job_log_fname = [jobLogDir, '/job_%A_%a.out'];
                job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
            end

            dataFullpath = [pathstr filesep fsname ext];
            dataReconFullpath = [pathstr '/DS/' fsname ext];
            inputFullpaths{curr} = dataFullpath;
            outputFullpaths{curr} = dataReconFullpath;

            funcStrs{curr} =  sprintf('cd /clusterfs/nvme/matthewmueller/Matlab-cudaSiRecon/src/cudaSirecon/;deskewPhasesFrame(''%s'',%.5f,%.5f,''SkewAngle'',%.5f,''nphases'',%.5f)', dataFullpath,xyPixelSize,dz,angle,np);
            curr = curr+1;
        end
    end

    maxJobNum = inf;
    taskBatchNum = 1;
    slurmDeskewParam = '-p abc --qos abc_high -n1 --mem-per-cpu=21000';

    is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', DeskewCpusPerTask, 'SlurmParam', slurmDeskewParam, ...
        'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'parseCluster', parseCluster, ...
        'MatlabLaunchStr', 'module load matlab/r2023a; matlab -nodisplay -nosplash -nodesktop -r', ...
        'masterCompute', false);

    if ~all(is_done_flag)
        slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', DeskewCpusPerTask, 'SlurmParam', slurmDeskewParam, ...
            'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'parseCluster', parseCluster, ...
            'MatlabLaunchStr', 'module load matlab/r2023a; matlab -nodisplay -nosplash -nodesktop -r', ...
            'masterCompute', false);
    end
    if DeskewData
        for i = 1:numel(dataPaths)
            dataPaths{i} = [dataPaths{i} '/DS/'];
        end
    end
    if DeskewPsfs
        for i = 1:numel(psfFiles)
            [pathstr, fsname, ext]  = psfFiles{i};
            psfFiles{i} = [pathstr '/DS/' fsname ext];
        end
    end

end

% Generate OTFs if needed
if genOTFs
    nFiles = numel(psfFiles);
    inputFullpaths = cell(nFiles, 1);
    outputFullpaths = cell(nFiles, 1);
    funcStrs = cell(nFiles, 1);
    curr = 1;

    for i = 1:numel(psfFiles)
        [pathstr, fsname, ext] = fileparts(psfFiles{i});
        jobLogDir = pathstr;
        if parseCluster
            if  ~exist(jobLogDir, 'dir')
                warning('The job log directory does not exist, use %s/job_logs as job log directory.', dataPaths{i})
                jobLogDir = sprintf('%s/job_logs', pathstr);
                mkdir(jobLogDir);
                fileattrib(jobLogDir, '+w', 'g');
            end
            job_log_fname = [jobLogDir, '/job_%A_%a.out'];
            job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
        end

        dataFullpath = [pathstr filesep fsname ext];
        dataReconFullpath = [pathstr '/OTFs/' fsname '_otf' ext];
        inputFullpaths{curr} = dataFullpath;
        outputFullpaths{curr} = dataReconFullpath;

        funcStrs{curr} =  sprintf('cd /clusterfs/nvme/matthewmueller/Matlab-cudaSiRecon/src/cudaSirecon/;makeOTF(''%s'',%d,%.5f,%d,%.5f,%.5f)', dataFullpath,np,moAngle,background,xyPixelSize,zres);
        curr = curr+1;
    end


    maxJobNum = inf;
    taskBatchNum = 1;
    slurmDeskewParam = '-p abc --qos abc_high -n1 --mem-per-cpu=21000';

    is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', DeskewCpusPerTask, 'SlurmParam', slurmDeskewParam, ...
        'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'parseCluster', parseCluster, ...
        'MatlabLaunchStr', 'module load matlab/r2023a; module load cudasirecon/1.0.0; matlab -nodisplay -nosplash -nodesktop -r', ...
        'masterCompute', false);

    if ~all(is_done_flag)
        slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', DeskewCpusPerTask, 'SlurmParam', slurmDeskewParam, ...
            'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'parseCluster', parseCluster, ...
            'MatlabLaunchStr', 'module load matlab/r2023a; module load cudasirecon/1.0.0; matlab -nodisplay -nosplash -nodesktop -r', ...
            'masterCompute', false);
    end

    for i = 1:numel(psfFiles)
        [pathstr, fsname, ext]  = psfFiles{i};
        otfFiles{i} = [pathstr '/OTFs/' fsname '_otf' ext];
        otfFol = [pathstr '/OTFs/'];
        if  ~exist(otfFol, 'dir')
            mkdir(otfFol);
            fileattrib(otfFol, '+w', 'g');
        end
    end

end


% Run sim recon
nFiles = 0;
for i = 1:numel(dataPaths)
    for cPatt = 1:numel(ChannelPatterns)
        fnames = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);
        fnames = {fnames.name};
        nFiles = nFiles + numel(fnames);
    end
end

% Check to see if any channel patterns are actually found
if nFiles == 0
    ME = MException('dataInput:noFiles', ...
        'No files matching any of the given channel patterns were found');
    throw(ME);
end

inputFullpaths = cell(nFiles, 1);
outputFullpaths = cell(nFiles, 1);
funcStrs = cell(nFiles, 1);
curr = 1;

for i = 1:numel(dataPaths)
    jobLogDir = dataPaths{i};
    oFol = [dataPaths{i} filesep outFol filesep];
    if  ~exist(oFol, 'dir')
        mkdir(oFol);
        fileattrib(oFol, '+w', 'g');
    end
    for cPatt = 1:numel(ChannelPatterns)
        fnames = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);
        fnames = {fnames.name};
        if parseCluster
            if  ~exist(jobLogDir, 'dir')
                warning('The job log directory does not exist, use %s/job_logs as job log directory.', dataPaths{i})
                jobLogDir = sprintf('%s/job_logs', dataPaths{i});
                mkdir(jobLogDir);
                fileattrib(jobLogDir, '+w', 'g');
            end
            job_log_fname = [jobLogDir, '/job_%A_%a.out'];
            job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
        end
        for j = 1: numel(fnames)
            [pathstr, fsname, ext] = fileparts(fnames{j});
            dataFullpath = [dataPaths{i} filesep fnames{j}];
            dataReconFullpath = [dataPaths{i} filesep outFol filesep fsname ext];
            inputFullpaths{curr} = dataFullpath;
            outputFullpaths{curr} = dataReconFullpath;
            idx(curr) = ~exist([dataPaths{i} filesep outFol filesep fsname '_recon' ext],'file');
            funcStrs{curr} =  sprintf(['cd /clusterfs/nvme/matthewmueller/Matlab-cudaSiRecon/src/cudaSirecon/;siReconWrapper(''%s'',''%s'',''%s'',''%s'',''%s'',[%d,%d,%d],%d, %d, %d)'], dataPaths{i},fsname,outFol,otfFiles{cPatt},configFiles{cPatt},chunkSize(1),chunkSize(2),chunkSize(3),overlap, background, ndir);
            curr = curr+1;
        end
    end
end

if OverwriteSIrecon
    idx = ones(size(idx));
end

maxJobNum = inf;
taskBatchNum = 1;

if nnz(idx)>0
    is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths(idx), outputFullpaths(idx), ...
        funcStrs(idx), 'cpusPerTask', cpusPerTask, 'SlurmParam', SlurmParam, ...
        'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'parseCluster', parseCluster, ...
        'MatlabLaunchStr', 'module load matlab/r2023a; module load cudasirecon/1.0.0; matlab -nodisplay -nosplash -nodesktop -r', ...
        'masterCompute', false, 'maxTrialNum', 1);

%     if ~all(is_done_flag)
%         slurm_cluster_generic_computing_wrapper(inputFullpaths(idx), outputFullpaths(idx), ...
%             funcStrs(idx), 'cpusPerTask', cpusPerTask, 'SlurmParam', SlurmParam, ...
%             'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'parseCluster', parseCluster, ...
%             'MatlabLaunchStr', 'module load matlab/r2023a; module load cudasirecon/1.0.0; matlab -nodisplay -nosplash -nodesktop -r', ...
%             'masterCompute', false, 'maxTrialNum', 1);
%     end
end

end
