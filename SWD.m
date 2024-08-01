% surface wave dispersion curve inversion with multi-offset data and fundamental mode
% Zhang Chang zhangchang23@mails.jlu.edu.cn; JIng Li, KAUST, jing.li@kaust.edu.sa 
% Improved efficiency and readability of the code
% C++ Eigen accelerate 2DFD; Modification for inversion logic
%% 
clc;
clear
close all;
addpath('./core/');
model_dir='C_model';
% load(['./model/',model_dir,'/low_vel.mat']);
load('model8_2.mat')
vs = model8_2;
% vs = imresize(vs,[40,300],'nearest');
[nz,nx]=size(vs);
vs=single(vs);
vp=vs*1.732;
vs_d = vs;
vsmin=min(vs(:));vsmax=max(vs(:));
fr=30;
dx=1;%dx=(min(vs_d(:))/fr/12);
dt=dx/max(vp(:))*0.5;
dtx=dt/dx;
for i=1:nx
%     vs(:,i)=linspace((vsmin+vsmax)/2,(vsmin+vsmax)/2,nz);
    vs(:,i)=linspace(min(vs(:)),max(vs(:)),nz);
end
load('vs3_0300_ADWD.mat', 'vs3_0300_init_WD')
vs = vs3_0300_init_WD;
vp=vs*1.732;
pickMethod=1;  %1==FDC 2==argmax
nt=floor((nx*dx/min(vs(:))/dt+1000)/10)*10-500;  % time step
[s,nw]=ricker(fr,dt,nt); s =single(s); % source wavelet
nbc=40;   % boundary layer
%% define acquisition geometry
ds=6; sx=single(1:ds:nx); sz=zeros(size(sx),'single')+1;[~,ns]=size(sx);
dg=2;gx=single(1:dg:nx);  gz=zeros(size(gx),'single')+1;  ng=numel(gx);
M=ds/dg;refsx=floor(sx/dg)+1;
dg=dg*dx;
%% define wavefield record
dt_wf=dt*5; nt_wf=floor((nt-1)*dt/dt_wf)+1;
pur=-0.05;
offset=floor(1.5*nz*dx);   % multi-offset &the maximum offset
offset=12;
offmin=3;offmax=floor(1.5*nz*dx);
NN=(offset-offmin)/3;   % offset change number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_type=0;
fd_order=22;fsz=0;source_type='w';
isfs=1;

parallel_init(ns);


%% Set Radon Transform papramter to cal. phase velocity 
vmin=floor(min(vs_d(:))/2);
vmax=floor(max(vs_d(:))*1.2);
np=vmax-vmin+1;
df1=0.1;
df=1/nt/dt;
N=df/df1;
df=df1;
fre1=(df*(0:floor(nt)-1));
fmin=10;
fmax=80;
fmin=fre1(find(fre1==fmin));
fmax=fre1(find(fre1==fmax));
m=0;FK=0;w=3;  % dot not need to change,w=3 is the min. traces for your multi-offset dispersion. You can increase the number
err=0.01;
parfor is=1:ns
  [~,seismo_v_d(:,:,is)]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs_d,isfs,fsz,fd_order,source_type,parameter_type);  
end
imagesc(seismo_v_d(:,:,1))
[a,b]=size(seismo_v_d(:,:,1));  
iteration=190;

%% ++++++++Start iteration 
%%--------------------------------------------------------------
k =1;kk=1;
vs_all=zeros(nz,nx,100);
dk_vs_all=zeros(nz,nx,100);
ind=find(fre1==fmin);
freq=(fmin:df:fmax);   
[~,npair]=size(freq);
win=npair+2*ind;
SoftArgNorm = 1e+4;
offsets = zeros(iteration,1);

for i =1:npair
    space_M(:,i) = linspace(1,np,np);
end
[s,nw]=ricker(fr,dt,nt); s =single(s); % source wavelet

while (k<=iteration)  
    tic;
    display(['Elastic_LSM, k=',num2str(k),' iteration=',num2str(iteration)]);

fre=2*pi*linspace(fmin,fmax,(fmax-fmin)/df+1);   

ml = zeros(np,npair,ns);
ml1 = zeros(np,npair,ns);
%% Cal. the obsdata dispersion curve two sides(Left and right)
line_mask_r = ones(npair,ns);
line_mask_l = ones(npair,ns);
cr_0 = 1.*ones(npair,ns);
cr_0l = 1.*ones(npair,ns);
parfor is=1:ns-floor(m/M)-round(w/M)
    [ml(:,:,is),dataLen]=RTr(seismo_v_d(:,:,is),is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M);
    [line_mask_r(:,is),~,cr_0(:,is)] = LHDispPick(ml(:,:,is),npair,vmin,cr_0(:,is),line_mask_r(:,is),dataLen,pickMethod);
end
parfor is=round(w/M)+floor(m/M)+1:ns
    [ml1(:,:,is),dataLen]=RTl(seismo_v_d(:,:,is),is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M);
    [line_mask_l(:,is),~,cr_0l(:,is)] = LHDispPick(ml1(:,:,is),npair,vmin,cr_0l(:,is),line_mask_l(:,is),dataLen,pickMethod);
end

cr_pre_r = 1.*ones(npair,ns);
cr_pre_l = 1.*ones(npair,ns);
   g_cl=zeros(nz,nx);
   g_cm=zeros(nz,nx);
   g_illum=zeros(nz,nx);
   %% 
   
parfor is=1:ns
   [~,seismo_v,wavefield_gradient]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf);
   pre_mask_r = zeros(npair,ns);
   pre_mask_l = zeros(npair,ns);
   saveForBackwardr = 0;
   saveForBackwardl = 0;
        if is<=ns-floor(m/M)-round(w/M)
            [mlr,dataLen,saveForBackwardr]=RTrAD(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M); % Cal. the predicetd data dispersion curve for two sides 
            [pre_mask_r,res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,npair,vmin,cr_0(:,is),line_mask_r(:,is),dataLen,pickMethod);
            saveForBackwardr.cr_r = cr_pre_r(:,is)-vmin;
        end

        if is>=round(w/M)+floor(m/M)+1
            [mll,dataLen,saveForBackwardl]=RTlAD(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M); % Cal. the predicetd data dispersion curve for two sides 
            [pre_mask_l,res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,npair,vmin,cr_0l(:,is),line_mask_l(:,is),dataLen,pickMethod);
            saveForBackwardl.cr_l = cr_pre_l(:,is)-vmin;
        end
        grad_outputr = cr_pre_r(:,is)-cr_0(:,is);
        grad_outputl = cr_pre_l(:,is)-cr_0l(:,is);

%     [seismo_v_d1]=weight_data_muti3(seismo_v_d(:,:,is),seismo_v,is,dt,df,offset,dg,0,w,M,m,ns,refsx,win,np,vmin,vmax,fmin,fmax, ...
%     cr_0(:,is),cr_0l(:,is),cr_pre_r(:,is),cr_pre_l(:,is),FK,ind,pre_mask_r,pre_mask_l); %The key step to calculate the backprohated data

%     [seismo_v_d1]=ADWDgrad_norm(nt,ng,ns,npair,is,w,m,M,SoftArgNorm,grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl);
% mask = envelope(norm_trace(seismo_v),300,'analytic');
% mask(mask>0.1)=1; 
% mask(mask~=1)=0; 
% [seismo_v_d1]=ADWDgrad_norm(nt,ng,ns,npair,is,w,m,M,SoftArgNorm,grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl);
[seismo_v_d1]=ADWDgrad_1(nt,ng,ns,npair,is,w,m,M,SoftArgNorm,grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl);
% [seismo_v_d1]=ADWDgrad3(nt,ng,ns,npair,is,w,m,M,SoftArgNorm,grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl);

% seismo_v_d1 = seismo_v_d1 .* mask;  
% [seismo_v_d1]=ADWDgradlimband(nt,ng,ns,npair,is,w,m,M,SoftArgNorm,grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl);
        [cl_img,cm_img,illum_div]=e2drtm_eigen(wavefield_gradient,single(seismo_v_d1),is,nbc,nt,dtx,dx,dt,gx,gz,s,vp,vs,isfs,fsz,fd_order,parameter_type,dt_wf);
        g_cl = g_cl+cl_img;g_cm = g_cm+cm_img;g_illum = g_illum+illum_div;
end

    % Use the penality method to build the objective function
  residual(k)=mean(mean(res_r));%+mean(mean(res_l));
    display(['residual = ',num2str( residual(k) ),' k=',num2str(k)]);
    res0=residual(k); 
    g_cl=g_cl./g_illum;g_cm=g_cm./g_illum;
    dk_vs = -4*vs.*g_cl+2*vs.*g_cm;
    dk_vs=single(smooth2a(double(dk_vs),1,6));  % Smooth the Vs gradient

    if k==1
        f1=0.5;
    end
    v_mean=(sum(vs(:).*vs(:)))^0.5;
    g_mean=(sum(dk_vs(:).*dk_vs(:)))^0.5;
    alpha=v_mean/g_mean*pur;  
    display(['v_mean=',num2str(v_mean),' g_mean=',num2str(g_mean),' alpha=',num2str(alpha)]);
    vs1=vs+alpha*f1*dk_vs;
    vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax; 
    %vs1(vp./vs1<1.3)=vp(vp./vs1<1.3)/1.3;
    %% 

parfor is=1:ns
    [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf);
    if is<=ns-floor(m/M)-round(w/M)
        [mlr,dataLen]=RTr(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M); % Cal. the predicetd data dispersion curve for two sides 
        [pre_mask_r,res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,npair,vmin,cr_0(:,is),line_mask_r(:,is),dataLen,pickMethod);

    end
    if is>=round(w/M)+floor(m/M)+1
        [mll,dataLen]=RTl(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M);
        [pre_mask_l,res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,npair,vmin,cr_0l(:,is),line_mask_l(:,is),dataLen,pickMethod);
    end


end
   res1=mean(mean(res_r));%+mean(mean(res_l));
    display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
 clear seismo_v cr_1 deta_c deta_c3
 
    if res1>res0
        while res1>res0 && f1>0.0001
            f2=f1; res2=res1;
            f1=f1*0.5;
            vs1=vs+alpha*f1*dk_vs;
            vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax; 
                %vs1(vp./vs1<1.3)=vp(vp./vs1<1.3)/1.3;

       
parfor is=1:ns
   [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf);
    if is<=ns-floor(m/M)-round(w/M)
        [mlr,dataLen]=RTr(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M); % Cal. the predicetd data dispersion curve for two sides 
        [pre_mask_r,res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,npair,vmin,cr_0(:,is),line_mask_r(:,is),dataLen,pickMethod);

    end
    if is>=round(w/M)+floor(m/M)+1
        [mll,dataLen]=RTl(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M);
        [pre_mask_l,res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,npair,vmin,cr_0l(:,is),line_mask_l(:,is),dataLen,pickMethod);
    end

end
   res1=mean(mean(res_r));%+mean(mean(res_l));
    display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
 clear seismo_v cr_1 deta_c deta_c3
        end
    else
        f2=f1*2;
        vs1=vs+alpha*f2*dk_vs; 
        vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax; 
            %vs1(vp./vs1<1.3)=vp(vp./vs1<1.3)/1.3;

        %vp=vs1*1.7;
        %% 
        
parfor is=1:ns
    [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf);
    if is<=ns-floor(m/M)-round(w/M)
        [mlr,dataLen]=RTr(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M); % Cal. the predicetd data dispersion curve for two sides 
        [pre_mask_r,res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,npair,vmin,cr_0(:,is),line_mask_r(:,is),dataLen,pickMethod);

    end
    if is>=round(w/M)+floor(m/M)+1
        [mll,dataLen]=RTl(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M);
        [pre_mask_l,res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,npair,vmin,cr_0l(:,is),line_mask_l(:,is),dataLen,pickMethod);
    end

end
   res2=mean(mean(res_r));%+mean(mean(res_l));
 clear seismo_v deta_c cr_1 deta_c3
        display(['f2= ',num2str(f2),' res2= ',num2str(res2)]);

end
    gama=(f1^2*(res0-res2)+f2^2*(res1-res0))/(2*res0*(f1-f2)+2*res1*f2-2*res2*f1);

    if isinf(gama)
        gama=0;
    end
    display(['gama= ',num2str(gama),' numerical step_length= ',num2str(gama*alpha)]);
    vs1=vs+alpha*gama*dk_vs;
    vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax; 
    %vs1(vp./vs1<1.3)=vp(vp./vs1<1.3)/1.3;

     %vp=vs1*1.7;
     %% 

parfor is=1:ns
   [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf);
    if is<=ns-floor(m/M)-round(w/M)
        [mlr,dataLen]=RTr(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M); % Cal. the predicetd data dispersion curve for two sides 
        [pre_mask_r,res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,npair,vmin,cr_0(:,is),line_mask_r(:,is),dataLen,pickMethod);

    end
    if is>=round(w/M)+floor(m/M)+1
        [mll,dataLen]=RTl(seismo_v,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,dg,offset,m,FK,M);
        [pre_mask_l,res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,npair,vmin,cr_0l(:,is),line_mask_l(:,is),dataLen,pickMethod);
    end

end

   res3=mean(mean(res_r));%+mean(mean(res_l));
 clear seismo_v deta_c cr_1 deta_c3
    display(['res3= ',num2str(res3)]);
    if (res3>res1 || res3>res2)
        if res1>res2
            res0=res2;
            gama=f2; 
        else
            res0=res1;
            gama=f1; 
        end
        vs1=vs+alpha*gama*dk_vs;      
       vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax; 
    %vs1(vp./vs1<1.3)=vp(vp./vs1<1.3)/1.3;

     %vp=vs1*1.7;
    else
        res0=res3;
    end
    vs=vs1;
   vs_all(:,:,k) = vs1; 
    dk_vs_all(:,:,k) = dk_vs; 
    offsets(k) = offset;

  if (k>1) && ((residual(k) - res0)/residual(k)) < err
    offset = offset - 3  
    f1=0.5;
  end
       k = k +1;
       kk = kk + 1;
  
   if offset==3 ||offset==51
       break
   end

  clear cr_0 cr_0l

end

 toc; 

dfs = mean(mean((vs-vs_d).^2))
% for i=1:nx
%     vs3_udlm_init_WD(:,i)=linspace(vsmin,vsmax,nz);
% end
% residual3_udlm_WD = residual;
% vs3_udlm_WD = vs;
% vs3_all_udlm_WD = vs_all;
% dk_vs3_all_udlm_WD = dk_vs_all;
% vs3_d_udlm_WD = vs_d;
% offsets3_udlm_WD = offsets;
% save('vs3_udlm_ADWD.mat','vs3_udlm_init_WD','vs3_udlm_WD',...
%     'vs3_all_udlm_WD','dk_vs3_all_udlm_WD','vs3_d_udlm_WD',...
%     'offsets3_udlm_WD','residual3_udlm_WD','-append')