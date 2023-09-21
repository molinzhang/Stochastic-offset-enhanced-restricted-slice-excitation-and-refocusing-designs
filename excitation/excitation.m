function pAD = excitation()
clear
IniVar = matfile('IniVars.mat');
IniVar.Properties.Writable = true;

%%
%%% this is a much complicated way to rebuild the CUBE. The class SpinCube
%%% and SpinArray could do this automatically by the property and function
%%% bewtween public/dependent and compact properties. I modify them into
%%% all public property. But some of them still use the function of
%%% automatically cacultation from M_, mask to M. Modify this after ISMRM.
cube = IniVar.cube;
MASK1 = load('/home/molin/shimm_nick/slice_excite/demo/ISMRM/compare/mask_true_noalias.mat'); %% MASK
MASK2 = MASK1.Mask;
%MASK2 = ones(size(MASK2));

MASK = MASK2; %logical(ones(size(MASK1)));
CUBE = cube;

CUBE.mask = logical(MASK);
CUBE.nM = sum(sum(sum(MASK)));

CUBE.dim = size(MASK);
CUBE.fov = [36, 36, 24];

b0Map = nan(CUBE.dim);
b0Map(MASK == 1) =-48 * 4 * 0.7242 / 4  *9.5 * 42.58;%-1200.12; %192.5 for slice shim;%for 64 coils only brain, -923.12; %- 4 * 1.2 * 0.7242 * 42.58 * 10 * 1.8; %-48 * 4 * 0.7242 / 4  *9.5 * 42.58 ;%; %-1035 - 650; %+ rand(CUBE.nM,1)* 10 ; %300 for ismrm
%b0Map(MASK == 1) = b0Map(MASK == 1) +50;
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
mag(MASK == 1) =1 ;
Mag(:,:,:,3) = mag;
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
IniVar.cube = CUBE;

%%modify target and pulse
 
Target = IniVar.target_OV90;
%load('r_4balltarget.mat');
TAR = double(MASK1.Target);
Target.d = cat(4, zeros(size(TAR)), TAR, 1 - TAR);



Target.weight = double(MASK2);
%Target.weight((Target.weight == 1) & (TAR == 1)) = 1;%1 - sum(sum(sum(TAR))) / CUBE.nM;%0.932;
%Target.weight((Target.weight == 1) & (TAR == 0)) = 0.2;%sum(sum(sum(TAR))) / CUBE.nM; %0.068;
A1 = squeeze(Target.weight(56, :, :)); % 19 % 56
B1 = squeeze(TAR(56, :, :));
A1((A1 == 1) & (B1 == 0)) = 2;

%Target.weight((Target.weight(19, :, :) == 1) & (TAR(19, :, :) == 0)) = 20;
%A1 = squeeze(Target.weight(16:22, :, :)); % 19 % 56
%B1 = squeeze(TAR(16:22, :, :));
%A1((A1 == 1) & (B1 == 0)) = 1.2;
%Target.weight(16:22,:,:) = A1;
Target.weight((Target.weight == 1) & (TAR == 1)) = 1;
Target.weight((Target.weight == 1) & (TAR == 0)) = 0.1;
Target.weight(56,:,:) = A1;

IniVar.target_OV90 = Target;

load('/home/molin/shimm_nick/slice_excite/demo/ISMRM/compare/New_consider/pulse_2_1_0.2.mat');

RF = load('../../rf_185.mat'); % rf_500: T = 6ms, B = 1.5kHz, bd = 300Hz % rf_185: T = 2ms, B = 925Hz
RF = RF.rf2 * 50; % 50 % 7 % pulse_rf
%RF = rf;

non_linear = load('/home/molin/shimm_nick/Dynamic_resemble/demo/ISMRM/compare/current_onlybrain_soft_64coil.mat');
non_linear = non_linear.amps;
non_linear(:) = 0;
%non_linear = zeros(64,1);
non_linear = repmat(non_linear, 1, size(RF,2));

%non_linear = gr(4:end, :);

kk = IniVar.pIni_OV90; 

kk.rf = RF;
kk.nr = non_linear;
kk.lr = zeros(3,size(RF,2)); %gr(1:3, :); %
kk.lr(2,:) =  0.7242; %8;%15;%7242 / 160;
%kk.lr = gr(1:3,:);
kk.gr = [kk.lr; kk.nr]; %gr + 5 * rand(size(gr));%

%kk.nr = non_linear(:,1:250);

kk.dt = 4e-6; 
kk.smax = 10000000000;
kk.gmax = 51;
IniVar.pIni_OV90 = kk;



pAD = OV90(cube, IniVar.target_OV90, IniVar.pIni_OV90);

%IV180(cube, IniVar.target_IV180, IniVar.pIni_IV180);

%IV180M(cube, IniVar.target_IV180M, IniVar.pIni_IV180M);

end

function pAD = OV90(cube, target, pIni)
  pAD = adpulses.opt.arctanAD(target, cube, pIni, 'err_meth', 'ml2xy' ...
                            , 'doClean',false, 'gpuID',1);

  figure 
  plot_res(pIni, pAD, cube, target, 'xy');
  %suptitle('OV90');
end

function IV180(cube, target, pIni)
  pAD = adpulses.opt.arctanAD(target, cube, pIni, 'err_meth', 'l2z' ...
                              , 'doClean',false, 'gpuID',0);
   
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
  [MT_1, MT_2] = deal(fn_res(MT_1), fn_res(MT_2));
  %MT_2 = MT_2 .* target.d(:,:,:,1);
  %d(isnan(d)) = -1;
  %MT_1(isnan(MT_1)) = -1;
  
  
  
  AAA = abs(MT_2((cube.mask == 1) & (target.d(:,:,:,1) == 1)));
  BBB = abs(MT_2((cube.mask == 1) & (target.d(:,:,:,1) == 0)));
  
  min_AAA = min(min(AAA));
  max_BBB = max(max(BBB));
  
  %loss seperately
  TAR = target.d(:,:,:,1);
  loss_o = norm(AAA - TAR((cube.mask == 1)&(TAR == 1))) ^2;
  loss_i = norm(BBB - TAR((cube.mask == 1)&(TAR == 0))) ^2;
  
  % reshapes
  [d, MT_1, MT_2] = deal(fn_tile(d), fn_tile(MT_1), fn_tile(MT_2));
  
  subplot(131);
  imagesc(d,'AlphaData',~isnan(d)); colorbar;
  caxis(cl_res); axis('equal'); pbaspect([1,1,1]); title('target');
  
  subplot(132);
  imagesc(MT_1,'AlphaData',~isnan(d)); colorbar;
  caxis(cl_res); axis('equal'); pbaspect([1,1,1]); title('Uniform B1+');

  subplot(133);
  imagesc(MT_2,'AlphaData',~isnan(d)); colorbar;
  caxis(cl_res); axis('equal'); pbaspect([1,1,1]); title('Random B1+');
end

function P = tile3dto2d(P)
  [nx, ny, nz] = size(P);
  nz_rt = sqrt(nz);
  [nc, nr] = deal(ceil(nz_rt), floor(nz_rt));
  nc = 6; nr = 10;
  P = cat(3, P, zeros([nx, ny, nr*nc-nz])); % -> (nx, ny, nc*nr)
  P = reshape(P, nx, ny*nc, []); % -> (nx, nc*ny, nr)
  P = reshape(permute(P, [2, 1, 3]), ny*nc, nx*nr).'; % -> (nr*nx, nc*ny);
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   