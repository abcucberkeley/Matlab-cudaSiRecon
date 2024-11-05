function [recIm] = cudaSireconChunk(inFol, inN, otfF, configF, outFullPath, varargin)

%export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib
%deskew=-32.45

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inFol', @ischar);
ip.addRequired('inN', @ischar);
ip.addRequired('otfF', @ischar);
ip.addRequired('configF', @ischar);
ip.addRequired('outFullPath', @ischar);
ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('ndirs', 1, @isnumeric);
ip.addParameter('overlap', 128, @isnumeric);
ip.addParameter('zoomfactor', 2, @isnumeric); % for xy
ip.addParameter('chunkSize', [320,320,0], @isnumeric);
ip.addParameter('background', 99, @isnumeric);
ip.addParameter('padpx', 0, @isnumeric);
ip.addParameter('occupancyRatio', 0.0, @isnumeric);
ip.parse(inFol, inN, otfF, configF, outFullPath, varargin{:});

pr = ip.Results;

nphases = pr.nphases;
ndirs = pr.ndirs;
ol = pr.overlap;
chunkSize = pr.chunkSize;
chunkSize(3) = chunkSize(3) * nphases*ndirs; % user defined in actual z planes 
background = pr.background;
padpx = pr.padpx;
occupancyRatio = pr.occupancyRatio;
zf = pr.zoomfactor;

im = readtiff([inFol '/' inN '.tif']);
[sy, sx, sz] = size(im);

% mask and interpolate
mask = im(:,:,1:nphases:end);
mask = max(mask, 0);
mask = logical(mask);
imask = imresize3(mask,[ceil(sy*zf),ceil(sx*zf), sz/(nphases*ndirs)], 'nearest');
clear mask;

im = single(im) - background;
im = max(im, 0);
otf = readtiff(otfF);

cd /clusterfs/nvme/matthewmueller/Matlab-cudaSiRecon/src/cudaSirecon;

if(chunkSize(1) > sx)
    chunkSize(1) = sx;
end
if(chunkSize(2) > sy)
    chunkSize(2) = sy;
end
if(chunkSize(3) > sz)
    chunkSize(3) = sz;
end
if(chunkSize(1) == 0)
    chunkSize(1) = sx;
end
if(chunkSize(2) == 0)
    chunkSize(2) = sy;
end
if(chunkSize(3) == 0)
    chunkSize(3) = sz;
end

nyc = floor(sy/(chunkSize(2)-ol));
nxc = floor(sx/(chunkSize(1)-ol));
nzc = floor(sz/(chunkSize(3)-(ol*nphases*ndirs)));

ymin = zeros(nyc,1);
ymax = zeros(nyc,1);

ymin_out = zeros(nyc,1);
ymax_out = zeros(nyc,1);

xmin = zeros(nxc,1);
xmax = zeros(nxc,1);

xmin_out = zeros(nxc,1);
xmax_out = zeros(nxc,1);

zmin = zeros(nzc,1);
zmax = zeros(nzc,1);

zmin_out = zeros(nzc,1);
zmax_out = zeros(nzc,1);

edgey = xmin;
edgex = edgey;
edgez = edgex;

nn = 0;
for h = 1:nxc
    for j = 1:nyc
        for k = 1:nzc
            nn = nn + 1;
            xmin(nn) = (1+((h-1)*(chunkSize(1)-ol)));
            ymin(nn) = (1+((j-1)*(chunkSize(2)-ol)));
            zmin(nn) = (1+((k-1)*(chunkSize(3)-(ol*nphases*ndirs))));
            if j==1
                ymin(nn) = 1;
                ymin_out(nn) = 1;
                edgey(nn) = 1;
            end
    
            if h==1
                xmin(nn) = 1;
                xmin_out(nn) = 1;
                edgex(nn) = 1;
            end
    
            if k==1
                zmin(nn) = 1;
                zmin_out(nn) = 1;
                edgez(nn) = 1;
            end
    
            if j < nyc
                %             ymax(nn) = (s(2)*j-(ol*(j)));
                %             ymax_out(nn) = (s(2)*j-(ol*(j)))*zf;
                ymax(nn) = (chunkSize(2)*j-(ol*(j-1)));
                ymax_out(nn) = (chunkSize(2)*j-(ol*(j-1)))*zf;
            else
                ymax(nn) = sy;
                ymax_out(nn) = sy*zf;
                edgey(nn) = 1;
            end
    
            if h < nxc
                %             xmax(nn) = (s(1)*h-(ol*(h)));
                %             xmax_out(nn) = (s(1)*h-(ol*(h)))*zf;
                xmax(nn) = (chunkSize(1)*h-(ol*(h-1)));
                xmax_out(nn) = (chunkSize(1)*h-(ol*(h-1)))*zf;
            else
                xmax(nn) = sx;
                xmax_out(nn) = sx*zf;
                edgex(nn) = 1;
            end
            
            if k < nzc
                zmax(nn) = (chunkSize(3)*k-((ol*nphases*ndirs)*(k-1)));
                zmax_out(nn) = (chunkSize(3)*k-((ol*nphases*ndirs)*(k-1)))/(nphases*ndirs);
            else
                zmax(nn) = sz;
                zmax_out(nn) = sz/(nphases*ndirs);
                edgez(nn) = 1;
            end
        end
    end
end
a = 1:nn;

recIm = zeros(ceil(sy*zf),ceil(sx*zf), sz/(nphases*ndirs),'single');

% Don't use a dimension's edges if we are not chunking in that dimension
if(chunkSize(1) == sx)
    edgex(:) = 0;
end
if(chunkSize(2) == sy)
    edgey(:) = 0;
end
if(chunkSize(3) == sz)
    edgez(:) = 0;
end

if nn > 1
    for rr = a(logical(edgey+edgex+edgez))
        subImage = im(ymin(rr):ymax(rr), xmin(rr):xmax(rr), zmin(rr):zmax(rr));
        [csy, csx, csz] = size(subImage);
        ck_im = zeros(csy+padpx*2, csx+padpx*2, csz, 'single');
        ck_im(padpx+1:csy+padpx, padpx+1:csx+padpx, :) = subImage;
        OR = nnz(subImage) / numel(subImage);
        
        if OR>=occupancyRatio
            %         outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
            %             configF, im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:), otf);
            outImg = cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
                configF, ck_im, otf);
            outImg = outImg(2*padpx+1:csy*2+2*padpx, 2*padpx+1:csx*2+2*padpx, :);

            zboundsL = (zmin(rr)-1)/(nphases*ndirs)+1-(ol*nphases*ndirs)/2:zmax(rr)/(nphases*ndirs)-(ol*nphases*ndirs)/2;
            if zmin(rr) == 1
                zboundsL = (zmin(rr)-1)/(nphases*ndirs)+1:zmax(rr)/(nphases*ndirs)-ol/2;
            elseif zmax(rr) == sz
                zboundsL = (zmin(rr)-1)/(nphases*ndirs)+1+ol/2:zmax(rr)/(nphases*ndirs);
            end

            zboundsR = ':';
            if zmin(rr) == 1
                [~,~,outZ] = size(outImg);
                zboundsR = 1:outZ-ol/2;
            elseif zmax(rr) == sz
                [~,~,outZ] = size(outImg);
                zboundsR = 1+ol/2:outZ;
            end


            if ymin(rr) ==1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) < sx % top left edge
                recIm(ymin(rr)*zf-1:ymax(rr)*zf-ol/2, xmin(rr)*zf-1:xmax(rr)*zf-ol/2, zboundsL) = outImg(1:end-ol/2, 1:end-ol/2, zboundsR);
            elseif ymin(rr) ==1 && xmin(rr) >1 && ymax(rr) < sy && xmax(rr) == sx % top right edge
                recIm(ymin(rr)*zf-1:ymax(rr)*zf-ol/2, xmin(rr)*zf-1+ol/2:xmax(rr)*zf, zboundsL) = outImg(1:end-ol/2, ol/2+1:end, zboundsR);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) == sy && xmax(rr) < sx % bottom left edge
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf, xmin(rr)*zf-1:xmax(rr)*zf-ol/2, zboundsL) = outImg(ol/2+1:end, 1:end-ol/2, zboundsR);
            elseif ymin(rr) >1 && xmin(rr) >1 && ymax(rr) == sy && xmax(rr) == sx % bottom right edge
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf, xmin(rr)*zf-1+ol/2:xmax(rr)*zf, zboundsL) = outImg(ol/2+1:end, ol/2+1:end, zboundsR);
            elseif ymin(rr) ==1 && xmin(rr) >1 && ymax(rr) < sy && xmax(rr) < sx % top middle
                recIm(ymin(rr)*zf-1:ymax(rr)*zf-ol/2, xmin(rr)*zf-1+ol/2:xmax(rr)*zf-ol/2, zboundsL) = outImg(1:end-ol/2, ol/2+1:end-ol/2, zboundsR);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) < sx % left middle
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf-ol/2, xmin(rr)*zf-1:xmax(rr)*zf-ol/2, zboundsL) = outImg(ol/2+1:end-ol/2, 1:end-ol/2, zboundsR);
            elseif ymin(rr) >1 && xmin(rr) >1 && ymax(rr) < sy && xmax(rr) == sx % right middle
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf-ol/2, xmin(rr)*zf-1+ol/2:xmax(rr)*zf, zboundsL) = outImg(ol/2+1:end-ol/2, ol/2+1:end, zboundsR);
            elseif ymin(rr) >1 && xmin(rr) >1 && ymax(rr) ==sy && xmax(rr) < sx % bottom middle
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf, xmin(rr)*zf-1+ol/2:xmax(rr)*zf-ol/2, zboundsL) = outImg(ol/2+1:end, ol/2+1:end-ol/2, zboundsR);
            elseif ymin(rr) ==1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) == sx % top center edge
                recIm(ymin(rr)*zf-1:ymax(rr)*zf-ol/2, xmin(rr)*zf-1:xmax(rr)*2,zboundsL) = outImg(1:end-ol/2, 1:end,zboundsR);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) == sx % middle center edge
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf-ol/2, xmin(rr)*zf-1:xmax(rr)*zf,zboundsL) = outImg(ol/2+1:end-ol/2, 1:end,zboundsR);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) == sy && xmax(rr) == sx % bottom center edge
                recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf, xmin(rr)*zf-1:xmax(rr)*zf,zboundsL) = outImg(1+ol/2:end, 1:end,zboundsR);
            end
        end
    end

    for rr = a(~logical(edgey+edgex+edgez))
        subImage = im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),zmin(rr):zmax(rr));
        [csy, csx, ~] = size(subImage);
        OR = nnz(subImage)/numel(subImage);
        if OR>=occupancyRatio
            %         outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
            %             configF, im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:), otf);

            outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
                configF, subImage, otf);
            outImg = outImg(2*padpx+1:csy*2+2*padpx, 2*padpx+1:csx*2+2*padpx, :);

            recIm(ymin(rr)*zf-1+ol/2:ymax(rr)*zf-ol/2, xmin(rr)*zf-1+ol/2:xmax(rr)*zf-ol/2, (zmin(rr)-1)/(nphases*ndirs)+1+ol/2:zmax(rr)/(nphases*ndirs)-ol/2) = outImg(ol/2+1:end-ol/2, ol/2+1:end-ol/2, ol/2+1:end-ol/2);
        end
    end
else
    recIm= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
        configF, im, otf);
end

% mask to zero
recIm = recIm .* imask;
end