function siReconWrapper(inFol,inN,outFol,otfF,configF,chunkSize,overlap, bak, ndir)
if nargin < 8
    bak = 0;
end
addpath(genpath('/clusterfs/nvme/matthewmueller/cudasireconMex/src/cudaSirecon/'));
cd '/clusterfs/nvme/ABCcode/XR_Repository'
setup([],true);
gpuDevice
%inFol = '/clusterfs/nvme/Gokul/LatticeSIMtest/DS';
%inN = 'test3';
%otfF = '/clusterfs/nvme/Gokul/LatticeSIMtest/OTF_camA_v3_px.tif';
%configF = [inFol '/mattConfig'];
%outFullPath = [inFol '/GPUsirecon/out3_zf1.5.tif'];
%inFol = '/clusterfs/nvme/Data/20211005_latticeSIM_wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV3/DS/';
%inN = 'RAW_exp01_CamA_ch0_CAM1_stack0149_488nm_0845113msec_0020289375msecAbs_000x_000y_000z_0000t';
%inN = 'test3_bk';
% otfF = '/clusterfs/nvme/Gokul/LatticeSIMtest/OTF_camB_v3_px.tif';
%otfF = '/clusterfs/nvme/Data/20220616_SIM_tiledVSsingle_CD9_Mito/psf/otf_560_single_C2.tif';
%configF = '/clusterfs/nvme/Gokul/LatticeSIMtest/configDS.txt';
% cd /clusterfs/nvme/matthewmueller/cudasireconMex/src/cudaSirecon/;
%fileList = dir('/clusterfs/nvme/Gokul/LatticeSIMtest/DS/*.txt');
%%
%fileList = dir([inFol '*.txt']);
%cs = [512,512,25];
%ol = 128;
%outFol = 'GPUsirecon_512c_128o_eq1_apom1_560bead_w0p005';
% for ii = 1:size(cs,1)
%for i = 1:size(fileList,1)
%outFullPath = [inFol sprintf('/GPUsirecon/%s/%s_recon_cs%d_ol%d.tif',outFol,inN,cs(ii,1),ol(ii))];
outFullPath = [inFol sprintf('/%s/%s_recon.tif',outFol,inN)];
mkdir([inFol filesep outFol]);
%         if isfile(outFullPath)
%             continue;
%         end
%,'chunkSize',[308,404,0],'overlap',0, testing if this is read only
if ~exist(outFullPath,'file')
    % GL: this cudaSiReconChunk.m currently (2023/11/20) messes with z and
    % is only suitable for LLS-SIM
    % out = cudaSireconChunk(inFol,inN,otfF,configF,outFullPath,'chunkSize',chunkSize,'overlap',overlap, 'bak', bak);
    out = cudaSireconChunk(inFol,inN,otfF,configF,outFullPath,'chunkSize',chunkSize,'overlap',overlap, 'bak', bak,'ndirs', ndir);
    %parallelWriteTiff(outFullPath,out,'w');
    %         bfsave(out,outFullPath,'Compression','Uncompressed');
    writetiff(out,outFullPath);
else
    disp('LLS-SIM reconstructed file exists!!!')
end
%writetiff(out(1:end/2,:,:),[outFullPath(1:end-4) '_firstHalfY.tif']);
%writetiff(out((end/2)+1:end,:,:),[outFullPath(1:end-4) '_secondHalfY.tif']);
%end
% end

%outFullPath = ['/clusterfs/nvme/matthewmueller/cudasireconMex/src/cudaSirecon/' 'test00.tif'];
%cudaSireconChunk(inFol,inN,otfF,configF,outFullPath,'chunk',true);
end