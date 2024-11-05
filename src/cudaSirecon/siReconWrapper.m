function siReconWrapper(inFol, inN, outFol, otfF, configF, chunkSize, overlap, background, ndirs, nphases, occupancyRatio)
addpath(genpath(fileparts(mfilename('fullpath'))));
cd '/clusterfs/nvme/matthewmueller/PetaKit5D';
setup;
gpuDevice;

outFullPath = [inFol sprintf('/%s/%s_recon.tif',outFol,inN)];
mkdir([inFol '/' outFol]);
if ~exist(outFullPath,'file')
    % GL: this cudaSiReconChunk.m currently (2023/11/20) messes with z and
    % is only suitable for LLS-SIM
    out = cudaSireconChunk(inFol, inN, otfF, configF, outFullPath, 'chunkSize', chunkSize, ...
                           'overlap', overlap, 'background', background, 'ndirs', ndirs, ...
                           'nphases', nphases, 'occupancyRatio', occupancyRatio);
    writetiff(out, outFullPath);
else
    disp('LLS-SIM reconstructed file exists!!!')
end

end