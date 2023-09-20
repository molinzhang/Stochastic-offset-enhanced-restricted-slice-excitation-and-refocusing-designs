function pAD = refocus()
clear
IniVar = matfile('./IniVars.mat');
IniVar.Properties.Writable = true;
IniVar1 = matfile('./IniVars.mat');
IniVar1.Properties.Writable = true;
IniVar2 = matfile('./IniVars.mat');
IniVar2.Properties.Writable = true;

%%
cube = IniVar.cube;
cube1 = IniVar1.cube;
cube2 = IniVar2.cube;

%% for the first cube
% load your own binary mask. It should have -TAR: target ROI. -MASK: target ROI + background. Alias de-selecting method is recommended for cartesian sampling.
MASK1 = load('./mask.mat'); 
MASK2 = MASK1.Mask;
MASK = MASK2; %logical(ones(size(MASK1)));

CUBE = cube;

CUBE.mask = logical(MASK);
CUBE.nM = sum(sum(sum(MASK)));

CUBE.dim = size(MASK);
CUBE.fov = [36,36,24];

b0Map = nan(CUBE.dim);
b0Map(MASK == 1) = 0;
CUBE.b0Map_ = b0Map(MASK == 1); 
CUBE.b0Map = b0Map;
[Xv, Yv, Zv] = meshgrid(-18:CUBE.res(1):17.99, -18:CUBE.res(2):17.99, -12:CUBE.res(3):11.99);

CUBE.loc = nan([CUBE.dim 3]);

inter = CUBE.loc(:,:,:,1);
inter(MASK == 1) = Xv(MASK ==1);
CUBE.loc(:,:,:,1) = inter;
inter = CUBE.loc(:,:,:,2);
inter(MASK == 1) = Yv(MASK ==1);
CUBE.loc(:,:,:,2) = inter;
inter = CUBE.loc(:,:,:,3);
inter(MASK == 1) = Zv(MASK ==1);
CUBE.loc(:,:,:,3) = inter;
CUBE.loc_ =  [Xv(MASK ==1), Yv(MASK== 1), Zv(MASK ==1)];

Mag = zeros(([CUBE.dim 3]));
mag = Mag(:,:,:,3);
mag(MASK == 1) = 1;
Mag(:,:,:,1) = mag;
CUBE.M = Mag;
mag1= Mag(:,:,:,1);
mag2 = Mag(:,:,:,2);
mag3 = Mag(:,:,:,3);

CUBE.M_ = [mag1(MASK ==1), mag2(MASK ==1), mag3(MASK == 1)];
%%% assume T1, T2, gamma is all the same;
T1 = CUBE.T1_(1); T2 = CUBE.T2_(1); gamma = CUBE.gam_(1);
CUBE.T1_ = zeros(CUBE.nM,1);
CUBE.T1_(:,1) = T1;
CUBE.T2_ = zeros(CUBE.nM,1);
CUBE.T2_(:,1) = T2;
CUBE.gam_ = zeros(CUBE.nM,1);
CUBE.gam_(:,1) = gamma;

%gamma : Hz/Gauss

cube = CUBE;

%% for the second cube
MASK2 = MASK1.Mask;

MASK = MASK2; %logical(ones(size(MASK1)))

CUBE1 = cube1;
CUBE1.mask = logical(MASK);
CUBE1.nM = sum(sum(sum(MASK)));

CUBE1.dim = size(MASK);
CUBE1.fov = [36, 36, 24];

b0Map = nan(CUBE1.dim);
b0Map(MASK == 1) =0;
CUBE1.b0Map_ = b0Map(MASK == 1); 
CUBE1.b0Map = b0Map;
[Xv, Yv, Zv] = meshgrid(-18:CUBE1.res(1):17.99, -18:CUBE1.res(2):17.99, -12:CUBE1.res(3):11.99);
CUBE1.loc = nan([CUBE1.dim 3]);

inter = CUBE1.loc(:,:,:,1);
inter(MASK == 1) = Xv(MASK ==1);
CUBE1.loc(:,:,:,1) = inter;
inter = CUBE1.loc(:,:,:,2);
inter(MASK == 1) = Yv(MASK ==1);
CUBE1.loc(:,:,:,2) = inter;
inter = CUBE1.loc(:,:,:,3);
inter(MASK == 1) = Zv(MASK ==1);
CUBE1.loc(:,:,:,3) = inter;
CUBE1.loc_ =  [Xv(MASK ==1), Yv(MASK== 1), Zv(MASK ==1)];

Mag = zeros(([CUBE1.dim 3]));
mag = Mag(:,:,:,3);
mag(MASK == 1) =1 ;
Mag(:,:,:,2) = mag;
CUBE1.M = Mag;
mag1= Mag(:,:,:,1);
mag2 = Mag(:,:,:,2);
mag3 = Mag(:,:,:,3);

CUBE1.M_ = [mag1(MASK ==1), mag2(MASK ==1), mag3(MASK == 1)];
%%% assume T1, T2, gamma is all the same;
T1 = CUBE1.T1_(1); T2 = CUBE1.T2_(1); gamma = CUBE1.gam_(1);
CUBE1.T1_ = zeros(CUBE1.nM,1);
CUBE1.T1_(:,1) = T1;
CUBE1.T2_ = zeros(CUBE1.nM,1);
CUBE1.T2_(:,1) = T2;
CUBE1.gam_ = zeros(CUBE1.nM,1);
CUBE1.gam_(:,1) = gamma;

%gamma : Hz/Gauss

cube1 = CUBE1;

%% for the third cube, excitation now
MASK2 = MASK1.Mask;
MASK = logical(MASK2); %logical(ones(size(MASK1)))
CUBE2 = cube2;


CUBE2.mask = logical(MASK);
CUBE2.nM = sum(sum(sum(MASK)));

CUBE2.dim = size(MASK);
CUBE2.fov = [36, 36, 24];

b0Map = nan(CUBE2.dim);
b0Map(MASK == 1) = 0; %-100; %-460 * 8 -60; %-80%- 4 * 1.2 * 0.48933 * 42.58 * 10 * 1.8; %-48 * 4 * 0.7242 / 4  *9.5 * 42.58 ;%; %-1035 - 650; %+ rand(CUBE.nM,1)* 10 ; %300 for ismrm

CUBE2.b0Map_ = b0Map(MASK == 1); 
CUBE2.b0Map = b0Map;
[Xv, Yv, Zv] = meshgrid(-18:CUBE2.res(1):17.99, -18:CUBE2.res(2):17.99, -12:CUBE2.res(3):11.99);
CUBE2.loc = nan([CUBE2.dim 3]);

inter = CUBE2.loc(:,:,:,1);
inter(MASK == 1) = Xv(MASK ==1);
CUBE2.loc(:,:,:,1) = inter;
inter = CUBE2.loc(:,:,:,2);
inter(MASK == 1) = Yv(MASK ==1);
CUBE2.loc(:,:,:,2) = inter;
inter = CUBE2.loc(:,:,:,3);
inter(MASK == 1) = Zv(MASK ==1);
CUBE2.loc(:,:,:,3) = inter;
CUBE2.loc_ =  [Xv(MASK ==1), Yv(MASK== 1), Zv(MASK ==1)];

Mag = zeros(([CUBE2.dim 3]));
mag = Mag(:,:,:,3);
mag(MASK == 1) =1 ;
Mag(:,:,:,3) = mag;
CUBE2.M = Mag;
mag1= Mag(:,:,:,1);
mag2 = Mag(:,:,:,2);
mag3 = Mag(:,:,:,3);

CUBE2.M_ = [mag1(MASK ==1), mag2(MASK ==1), mag3(MASK == 1)];
%%% assume T1, T2, gamma is all the same;
T1 = CUBE2.T1_(1); T2 = CUBE2.T2_(1); gamma = CUBE2.gam_(1);
CUBE2.T1_ = zeros(CUBE2.nM,1);
CUBE2.T1_(:,1) = T1;
CUBE2.T2_ = zeros(CUBE2.nM,1);
CUBE2.T2_(:,1) = T2;
CUBE2.gam_ = zeros(CUBE2.nM,1);
CUBE2.gam_(:,1) = gamma;

%gamma : Hz/Gauss

cube2 = CUBE2;

%%
IniVar.cube = CUBE;

%%modify target and pulse
 
Target = IniVar.target_OV90;
TAR = double(MASK1.Target);
TAR_all = double(MASK1.Mask);


%TAR_refo = TAR;
%TAR_refo(TAR == 1) = -1;

TAR_refo = TAR_all;
TAR_refo(TAR == 1) = -1;

Target.d = cat(4, TAR_all, zeros(size(TAR)), 1 - TAR_all); %Target for Mx: (1,0,0) -> (1,0,0)
Target.d1 = cat(4, zeros(size(TAR)), TAR_refo, 1 - TAR_all); %Target for My: (0,1,0) -> (0,-1,0)
%Target.d2 = cat(4, zeros(size(TAR)), zeros(size(TAR)), -ones(size(TAR))); %Target for Mz: (0,0,1) -> (0,0,+1/-1)
Target.d2 = cat(4, zeros(size(TAR)), TAR, 1 - TAR);


% %% Assign different weights for different locations.
% Mask_weight = double(MASK);
% Slice_of_interest = squeeze(Mask_weight(45,:,:)); % 19, 44
% Slice_of_interest(Slice_of_interest == 1) = 1.1;%1.1; %2;
% Mask_weight(45,:,:) = Slice_of_interest;

% Slice_of_interest = squeeze(Mask_weight(46,:,:)); % 19, 44
% Slice_of_interest(Slice_of_interest == 1) = 1;%0.1; %2;
% Mask_weight(46,:,:) = Slice_of_interest;
% Slice_of_interest = squeeze(Mask_weight(44,:,:)); % 19, 44
% Slice_of_interest(Slice_of_interest == 1) = 1;%0.1; %2;
% Mask_weight(44,:,:) = Slice_of_interest;

Mask_weight = double(MASK);
Slice_of_interest = squeeze(Mask_weight(:,:,31)); % 19, 44
Slice_of_interest(Slice_of_interest == 1) = 1.1; %1.1; % 4
Mask_weight(:,:,31) = Slice_of_interest;

Mask_weight(Mask_weight == 1) =0.1; %0.1;%0.0064  %0.1;
Mask_weight(TAR == 1) = 0;

W_back = double(MASK);
W_back(TAR == 1) = 0; 

Target.weight = double(TAR) * 4;%double(MASK);%MASK1.Mask);%double(TAR);%Mask_weight;
Target.weight1 = Mask_weight;%Mask_weight; %double(MASK);%MASK1.Mask);%Mask_weight;%double(MASK);%double(TAR);%Mask_weight;

Ex_mask = double(Mask_weight);
Ex_mask(TAR == 1) = 16;

Target.weight2 = Ex_mask;%double(MASK);

IniVar.target_OV90 = Target;



%% % Load initialization. Gauss pulse.
RF1 = load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/rf_gauss_false.mat');
RF = RF1.rf_gau;
non_linear = zeros(64,1);
non_linear = repmat(non_linear, 1, size(RF,2));

%% % Load optimized RF
% load('/home/molin/shimm_nick/op_refoc/mrm_results/sagital/pulse1.mat')
% RF = rf;
% non_linear = gr(4:end, :);


kk = IniVar.pIni_OV90;

kk.rf = RF;
kk.nr = non_linear;
kk.lr = zeros(3,size(RF,2));
kk.lr(3,:) = 0.2; % add linear gradient
kk.gr = [kk.lr; kk.nr];
kk.dt = 8e-6; %4e-6;
kk.smax = 1.0;
kk.rfmax = 1.0;
IniVar.pIni_OV90 = kk;

pAD = OV90(cube, cube1, cube2, IniVar.target_OV90, IniVar.pIni_OV90);

end

function pAD = OV90(cube, cube1, cube2, target, pIni)
  pAD = adpulses.opt.arctanAD(target, cube, cube1, cube2, pIni, 'err_meth', 'l2xy' ...
                           , 'doClean',false, 'gpuID',0);

end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   