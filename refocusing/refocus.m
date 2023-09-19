function pAD = refocus()
clear
IniVar = matfile('/home/molin/shimm_nick/decompose_refocus/demo/IniVars.mat');
IniVar.Properties.Writable = true;
IniVar1 = matfile('/home/molin/shimm_nick/decompose_refocus/demo/IniVars.mat');
IniVar1.Properties.Writable = true;
IniVar2 = matfile('/home/molin/shimm_nick/decompose_refocus/demo/IniVars.mat');
IniVar2.Properties.Writable = true;
%IniVar3 = matfile('/home/molin/shimm_nick/decompose_refocus/demo/IniVars.mat');
%IniVar3.Properties.Writable = true;


%%
%%% this is a much complicated way to rebuild the CUBE. The class SpinCube
%%% and SpinArray could do this automatically by the property and function
%%% bewtween public/dependent and compact properties. I modify them into
%%% all public property. But some of them still use the function of
%%% automatically cacultation from M_, mask to M. Modify this after ISMRM.
cube = IniVar.cube;
cube1 = IniVar1.cube;
cube2 = IniVar2.cube;
cube3 = IniVar3.cube;

%% for the first cube
%MASK1 = load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/Mask_sagital_pregnancy_small.mat'); %r7mask.mat;%Mask_checkslice.mat');%MASK_12mmROI_slice.mat');%Mask_pregnancy_ismrm2022.mat'); %; % %% MASK
%MASK1 = load('/home/mmatlaolin/shimm_nick/decompose_refocus/demo/ISMRM/compare/phantom_3mm_mask_x_noshift.mat');
MASK1 = load('/home/molin/shimm_nick/op_rewind/MASK_z.mat');
MASK2 = MASK1.Mask; %Mask;%.Target;
%MASK2 = ones(size(MASK2));
%MASK = zeros(size(MASK2));

%MASK(19, 15, 10) = MASK2(19, 15 ,10);

%load('/home/molin/shimm_nick/op_rewind/fullmask.mat')
%MASK2 = D1;
MASK = MASK2; %logical(ones(size(MASK1)));


%MASK = ones(size(MASK));
CUBE = cube;

CUBE.mask = logical(MASK);
CUBE.nM = sum(sum(sum(MASK)));

CUBE.dim = size(MASK);
CUBE.fov = [0.4091*3, 0.4091*3, 24];%[36,36,24]; %[21, 21, 15];

b0Map = nan(CUBE.dim);
b0Map(MASK == 1) = 0;%-100; %-80; %- 460.0 * 8 - 60; 
%b0Map(MASK == 1) = b0Map(MASK == 1) +50;
CUBE.b0Map_ = b0Map(MASK == 1); 
CUBE.b0Map = b0Map;
%[Xv, Yv, Zv] = meshgrid(-18:CUBE.res(1):17.99, -18:CUBE.res(2):17.99, -12:CUBE.res(3):11.99);
[Xv, Yv, Zv] = meshgrid(-0.8178:CUBE.res(1):0.0004, -0.8178:CUBE.res(2):0.0004, -12:CUBE.res(3):11.99);
%%% need to think about it more. Currently, Mask(:,45,31) = 1. But this
%%% corresponds to the change in Yv.
%[Xv, Yv, Zv] = meshgrid(-10.5:CUBE.res(1):10.49, -10.5:CUBE.res(2):10.49, -7.5:CUBE.res(3):7.49);
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

%mag1(19,:,:) = 1;
%mag3(19,:,:) = 0;

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
MASK2 = MASK1.Mask;%Mask; %.Target;
%MASK2 = ones(size(MASK2));
%MASK = zeros(size(MASK2));
%MASK(19,15,10) = MASK2(19,15,10);
%load('/home/molin/shimm_nick/op_rewind/fullmask.mat')
%MASK2 = D1;
MASK = MASK2; %logical(ones(size(MASK1)))
%MASK = ones(size(MASK));
CUBE1 = cube1;
CUBE1.mask = logical(MASK);
CUBE1.nM = sum(sum(sum(MASK)));

CUBE1.dim = size(MASK);
CUBE1.fov = [0.4091*3, 0.4091*3, 24];%[36, 36, 24];%[21, 21, 15]; % 36, 36, 24

b0Map = nan(CUBE1.dim);
b0Map(MASK == 1) =0; %-100; %-460 * 8 -60; %-80; %; 
%b0Map(MASK == 1) = b0Map(MASK == 1) +50;
CUBE1.b0Map_ = b0Map(MASK == 1); 
CUBE1.b0Map = b0Map;
%[Xv, Yv, Zv] = meshgrid(-18:CUBE1.res(1):17.99, -18:CUBE1.res(2):17.99, -12:CUBE1.res(3):11.99);
[Xv, Yv, Zv] = meshgrid(-0.8178:CUBE.res(1):0.0004, -0.8178:CUBE.res(2):0.0004, -12:CUBE.res(3):11.99);
%[Xv, Yv, Zv] = meshgrid(-10.5:CUBE1.res(1):10.49, -10.5:CUBE1.res(2):10.49, -7.5:CUBE1.res(3):7.49);
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

%mag2(19,:,:) = 1;
%mag3(19,:,:) = 0;

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
%MASK2 = ones(size(MASK2));
%MASK = zeros(size(MASK2));
%MASK(MASK2 == 1) = 1;
%sliceS = MASK(19,:,:);
%sliceS(sliceS == 1) = 0;
%MASK(19,:,:) = sliceS;
MASK = logical(MASK2); %logical(ones(size(MASK1)))
%MASK(19,:,:) = 0;
CUBE2 = cube2;
%t = MASK1.Target;
%MASK(t == 1) = 0;

CUBE2.mask = logical(MASK);
CUBE2.nM = sum(sum(sum(MASK)));

CUBE2.dim = size(MASK);
CUBE2.fov = [0.4091*3, 0.4091*3, 24];%[36, 36, 24];%[21, 21, 15];%[36, 36, 24];

b0Map = nan(CUBE2.dim);
b0Map(MASK == 1) = 0; %-100; %-460 * 8 -60; %-80%- 4 * 1.2 * 0.48933 * 42.58 * 10 * 1.8; %-48 * 4 * 0.7242 / 4  *9.5 * 42.58 ;%; %-1035 - 650; %+ rand(CUBE.nM,1)* 10 ; %300 for ismrm
%b0Map(MASK == 1) = b0Map(MASK == 1) +50;
CUBE2.b0Map_ = b0Map(MASK == 1); 
CUBE2.b0Map = b0Map;
%[Xv, Yv, Zv] = meshgrid(-18:CUBE2.res(1):17.99, -18:CUBE2.res(2):17.99, -12:CUBE2.res(3):11.99);
[Xv, Yv, Zv] = meshgrid(-0.8178:CUBE.res(1):0.0004, -0.8178:CUBE.res(2):0.0004, -12:CUBE.res(3):11.99);
%[Xv, Yv, Zv] = meshgrid(-10.5:CUBE2.res(1):10.49, -10.5:CUBE2.res(2):10.49, -7.5:CUBE2.res(3):7.49); %meshgrid(-18:CUBE2.res(1):17.99, -18:CUBE2.res(2):17.99, -12:CUBE2.res(3):11.99);
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
%%% for test
%MASK1 = load('/home/molin/shimm_nick/Dynamic_resemble/demo/ISMRM/compare/masks_slice.mat'); %% MASK
%MASK2 = MASK1.Target;
%%MASK2 = ones(size(MASK2));
%MASK = MASK2; %logical(ones(size(MASK1)));
%CUBE3 = cube3;

%CUBE3.mask = logical(MASK);
%CUBE3.nM = sum(sum(sum(MASK)));

%CUBE3.dim = size(MASK);
%CUBE3.fov = [36, 36, 24];

%b0Map = nan(CUBE3.dim);
%b0Map(MASK == 1) = - 460.0 * 10; %- 4 * 1.2 * 0.48933 * 42.58 * 10 * 1.8; %-48 * 4 * 0.7242 / 4  *9.5 * 42.58 ;%; %-1035 - 650; %+ rand(CUBE.nM,1)* 10 ; %300 for ismrm
%b0Map(MASK == 1) = b0Map(MASK == 1) +50;
%CUBE3.b0Map_ = b0Map(MASK == 1); 
%CUBE3.b0Map = b0Map;
%[Xv, Yv, Zv] = meshgrid(-18:CUBE3.res(1):17.99, -18:CUBE3.res(2):17.99, -12:CUBE3.res(3):11.99);
%CUBE3.loc = nan([CUBE3.dim 3]);

%inter = CUBE3.loc(:,:,:,1);
%inter(MASK == 1) = Xv(MASK ==1);
%CUBE3.loc(:,:,:,1) = inter;
%inter = CUBE3.loc(:,:,:,2);
%inter(MASK == 1) = Yv(MASK ==1);
%CUBE3.loc(:,:,:,2) = inter;
%inter = CUBE3.loc(:,:,:,3);
%inter(MASK == 1) = Zv(MASK ==1);
%CUBE3.loc(:,:,:,3) = inter;
%CUBE3.loc_ =  [Xv(MASK ==1), Yv(MASK== 1), Zv(MASK ==1)];

%Mag = zeros(([CUBE3.dim 3]));
%mag = Mag(:,:,:,3);
%mag(MASK == 1) = 1;
%Mag(:,:,:,1) = mag;
%CUBE3.M = Mag;
%mag1= Mag(:,:,:,1);
%mag2 = Mag(:,:,:,2);
%mag3 = Mag(:,:,:,3);


%mag1(19,:,:) = -1 + 2*rand(30,20);
%mag2(19,:,:) = (-1 + 2* randi([0 1], 30, 20)) .* squeeze(sqrt(1 - mag1(19,:,:).^2));
%CUBE3.M_ = [mag1(MASK ==1), mag2(MASK ==1), mag3(MASK == 1)];
%%%% assume T1, T2, gamma is all the same;
%T1 = CUBE3.T1_(1); T2 = CUBE3.T2_(1); gamma = CUBE3.gam_(1);
%CUBE3.T1_ = zeros(CUBE3.nM,1);
%CUBE3.T1_(:,1) = T1;
%CUBE3.T2_ = zeros(CUBE3.nM,1);
%CUBE3.T2_(:,1) = T2;
%CUBE3.gam_ = zeros(CUBE3.nM,1);
%CUBE3.gam_(:,1) = gamma;

%%gamma : Hz/Gauss

%cube3 = CUBE3;


%%
IniVar.cube = CUBE;

%%modify target and pulse
 
Target = IniVar.target_OV90;
%load('r_4balltarget.mat');
TAR = double(MASK1.Target);

TAR_all = double(MASK1.Mask);


%TAR_refo = TAR;
%TAR_refo(TAR == 1) = -1;

TAR_refo = TAR_all;
TAR_refo(TAR == 1) = -1;
%TAR_refo(16, :, :) = -1;

Target.d = cat(4, TAR_all, zeros(size(TAR)), 1 - TAR_all); %Target for Mx: (1,0,0) -> (1,0,0)
Target.d1 = cat(4, zeros(size(TAR)), TAR_refo, 1 - TAR_all); %Target for My: (0,1,0) -> (0,-1,0)
%Target.d2 = cat(4, zeros(size(TAR)), zeros(size(TAR)), -ones(size(TAR))); %Target for Mz: (0,0,1) -> (0,0,+1/-1)
Target.d2 = cat(4, zeros(size(TAR)), TAR, 1 - TAR);

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

load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/refocus_test/Check_nosense/slice_refoc.mat');
%load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/refocus_test/Check_nosense/pulse2.mat');
load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/refocus_test/selective/field_correct/pulse1.mat');

load('/home/molin/shimm_nick/op_refoc/mrm_results/sagital/pulse1.mat')

%RF1 = load('../../rf_185.mat'); % T = 6ms, B = 1.5kHz, bd = 300Hz
%RF = RF1.rf2 * 7; %refocs * 7; % 50 % 7
%RF1 = load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/refocus_test/Check_nosense/SLR_refocus.mat');
RF1 = load('/home/molin/shimm_nick/decompose_refocus/demo/ISMRM/compare/rf_gauss_false.mat');
RF = RF1.rf_gau;
%RF = RF1.refocs * 5;
%RF = zeros(1, 500);

%RF(1,1001:1500) = RF1;
RF = rf;

%non_linear = load('/home/molin/shimm_nick/Dynamic_resemble/demo/ISMRM/compare/current_onlybrain_soft_64coil.mat');
%non_linear = non_linear.amps;matla
%non_linear(:) = 0;
non_linear = zeros(64,1);
non_linear = repmat(non_linear, 1, size(RF,2));
non_linear = gr(4:end, :);
%non_linear = zeros(64, size(RF,2));

kk = IniVar.pIni_OV90;

kk.rf = RF;
kk.nr = non_linear;
kk.lr = zeros(3,size(RF,2)); %gr(1:3, :); %
kk.lr(3,:) = 0.2;%0.7; %2 * 4.8933;%7242 / 160; %%%%%% lr(2,:) or lr(3,:)
%kk.lr(2,1001:1500) = 0.001;
%kk.lr = gr(1:3,:);
kk.gr = [kk.lr; kk.nr]; %gr + 5 * rand(size(gr));%
kk.dt = 8e-6; %4e-6;
kk.smax = 10000000000;
kk.rfmax = 1.0;
IniVar.pIni_OV90 = kk;



pAD = OV90(cube, cube1, cube2, IniVar.target_OV90, IniVar.pIni_OV90);

%IV180(cube, IniVmask_out(mask_out > 0.000001) = 1;ar.target_IV180, IniVar.pIni_IV180);

%IV180M(cube, IniVar.target_IV180M, IniVar.pIni_IV180M);

end

function pAD = OV90(cube, cube1, cube2, target, pIni)
  pAD = adpulses.opt.arctanAD(target, cube, cube1, cube2, pIni, 'err_meth', 'l2xy' ...
                           , 'doClean',false, 'gpuID',0);

  figure 
  plot_res(pIni, pAD, cube, target, 'xy');
  plot_res(pIni, pAD, cube1, target, 'xy');
  %plot_res(pIni, pAD, cube2, target, 'xy');
  %suptitle('OV90');
end

function IV180(cube, target, pIni)
  %pAD = adpulses.opt.arctanAD(target, cube, pIni, 'err_meth', 'l2z' ...
  %                            , 'doClean',false, 'gpuID',0);
   
  %figure
  plot_res(pIni, PIni, cube, target, 'z');
  suptitle('IV180');
end

function IV180M(cube, target, pIni)
  pAD = adpulses.opt.arctanAD(target, cube, pIni, 'err_meth', 'l2z' ...
                              , 'doClean',false, 'gpuID',0);
 
  figure
  plot_res(pIni, pAD, cube, target, 'z');
  suptitle('IV180M');
end

%% Utils
function plot_res(p1, p2, cube, target, xy_z) 
  fn_sim1 = @(p)cube.applypulse(p, 'doCim',true, 'doEmbed',true, 'b1Map_', 0.9 + 0.05 * rand(cube.nM,1));
  fn_sim = @(p)cube.applypulse(p, 'doCim',true, 'doEmbed',true);
  [MT_1, MT_2] = deal(fn_sim(p1), fn_sim(p2));
  
  if strcmpi(xy_z, 'xy')
    % xy component and transversal MLS NRMSE
    fn_res = @(MT)MT(:,:,:,1) + 1i*MT(:,:,:,2);
    fn_tile = @(P)abs(tile3dto2d(P));
    cl_res = [0, 1];
  elseif strcmpi(xy_z, 'z')
    % z component and longitudinal LS NRMSE
    fn_res = @(MT)MT(:,:,:,3);
    fn_tile = @(P)tile3dto2d(P);
    cl_res = [-1, 1];
  else, error('Unknown xy_z type');
  end
  
  d = fn_res(target.d);
  d(cube.mask == 0) = nan;
  %[MT_1, MT_2] = deal(fn_res(MT_1), fn_res(MT_2));
  MT_1 = MT_2(:,:,:,1);
  MT_2 = MT_2(:,:,:,2);
  cl_res = [-1,1];
  %MT_2 = MT_2 .* target.d(:,:,:,1);
  %d(isnan(d)) = -1;
  %MT_1(isnan(MT_1)) = -1;
  
  
  AAA = abs(MT_2((cube.mask == 1) & (target.d(:,:,:,1) == 1)));
  BBB = abs(MT_2((cube.mask == 1) & (target.d(:,:,:,1) == 0)));
  
  min_AAA = min(min(AAA));
  max_BBB = max(max(BBB));
  
  %loss seperately
  TAR = target.d(:,:,:,1);
  loss_i = norm(BBB - TAR((cube.mask == 1)&(TAR == 0))) ^2;
  
  % reshapes
  [d, MT_1, MT_2] = deal(fn_tile(d), fn_tile(MT_1), fn_tile(MT_2));
  
  %subplot(131);
  %imagesc(d); colorbar;
  %caxis(cl_res); axis('equal'); pbaspect([1,1,1]); title('target');
  
  subplot(121);
  imagesc(MT_1); colorbar;
  caxis(cl_res); axis('equal'); pbaspect([1,1,1]); title('Uniform B1+');

  subplot(122);
  imagesc(MT_2); colorbar;
  caxis(cl_res); axis('equal'); pbaspect([1,1,1]); title('Random B1+');
end

function P = tile3dto2d(P)
  [nx, ny, nz] = size(P);
  nz_rt = sqrt(nz);
  [nc, nr] = deal(ceil(nz_rt), floor(nz_rt));
  nc = 10; nr = 6;
  P = cat(3, P, zeros([nx, ny, nr*nc-nz])); % -> (nx, ny, nc*nr)
  P = reshape(P, nx, ny*nc, []); % -> (nx, nc*ny, nr)
  P = reshape(permute(P, [2, 1, 3]), ny*nc, nx*nr).'; % -> (nr*nx, nc*ny);
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   