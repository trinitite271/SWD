function [seismo_v_d1,deta_c]=weight_data_muti4(seismo_v_d,seismo_v,is,dt,df,offset,dg,F1,w,M,m,ns,refsx,win,freq,npair,cr_0,cr_0l,cr_pre_r,cr_pre_l,FK,ind,line_mask_r,line_mask_l)
 [nt,b]=size(seismo_v_d(:,:,1));
fre=2*pi*freq; 
 a=nt;
 %% Calcuate the weight source
 ntpad=floor(1/dt/df);
 %spad=zeros(ntpad);  
 % Caluate the weight residual data
maxoffset=b;


if  is<=ns-floor(m/M)-round(w/M)
deta_cc2=zeros(win,1);
deta_cc2(ind:ind+npair-1) = fre'.*(cr_0-cr_pre_r)./(cr_0.*cr_pre_r).*line_mask_r;
    nend=m+offset-1+refsx(is);
    if nend>maxoffset
      nend=maxoffset;
    end
    newoffset=nend-(refsx(is)+m)+1;
    tempdata=zeros(ntpad,newoffset);
    y(1:newoffset)=(m+linspace(0,newoffset-1,newoffset))*dg;% Define the offset (m)
    if F1==1
    tempdata(1:nt,:)=seismo_v_d(:,round(m+refsx(is)):nend);
    else
    tempdata(1:nt,:)=seismo_v(:,round(m+refsx(is)):nend);    
    end
    deta_seismo_v1(:,1:newoffset)=fft(tempdata); %TRansform shot gather to FK domain
    y1(:)=deta_cc2(:); % Define the delta k
    
for ii=1:newoffset
      if  F1==1
     deta_seismo_v4(:,ii)=(((-1i*y1(:)*((y(ii))).*exp(-1i*y1(:)*((y(ii))))/pi/2)));   %% observed data
    else
     deta_seismo_v4(:,ii)=(-1i*y1(:)*y(ii))/pi/2;   %% predicted data
      end
end


clear y y1 tempdata
%for omega=1:floor(ntpad/2)
    for omega=1:win
%     if omega>=nfmin && omega<=nfmax
        deta_seismo_v1(omega+1,:)=deta_seismo_v1(omega+1,:).*deta_seismo_v4(omega,:);
%     else
%         deta_seismo_v11(omega+1,:,:)=0;
%     end
   deta_seismo_v1(ntpad+1-omega,:)=conj(deta_seismo_v1(omega+1,:));
end
end

   seismo_v_d1=zeros(nt,maxoffset); % Define the backprogate data
%  seismo_v_d1=seismo_v_d;
if is<=ns-floor(m/M)-round(w/M)
    nend=m+offset-1+refsx(is);
    if nend>maxoffset
      nend=maxoffset;
    end
    newoffset=nend-(m+refsx(is))+1;
   deta_seismo_vv(:,:)=(ifft(deta_seismo_v1(:,:)));
   seismo_v_d1(:,round(m+refsx(is)):nend)=deta_seismo_vv(1:nt,1:newoffset);
end
    
clear deta_seismo_vv


%%% Left side part

if is>round(w/M)+floor(m/M)+1
deta_cc22 = zeros(win,1);
deta_cc22(ind:ind+npair-1) = fre'.*(cr_0l-cr_pre_l)./(cr_0l.*cr_pre_l).*line_mask_l;
   nend=refsx(is)-offset-m;
  if nend<1
      nend=0;
  end
    newoffset=refsx(is)-m-nend;
    tempdata=zeros(ntpad,newoffset);
    y(1:newoffset)=(m+linspace(newoffset-1,0,newoffset))*dg;% Define the offset (m)
    if F1==1
    tempdata(1:nt,:)=seismo_v_d(:,nend+1:refsx(is)-m);
    else
       tempdata(1:nt,:)=seismo_v(:,nend+1:refsx(is)-m);  
    end
    deta_seismo_v11(:,1:newoffset)=fft(tempdata); %TRansform shot gather to FK domain
    y1(:)=deta_cc22(:); % Define the delta k
    
for ii=1:newoffset
    if  F1==1
     deta_seismo_v44(:,ii)=(((-1i*y1(:)*((y(ii))).*exp(-1i*y1(:)*((y(ii))))/pi/2)));   %% observed data
    else
     deta_seismo_v44(:,ii)=(-1i*y1(:)*y(ii))/pi/2;   %% predicted data
     end
end



%for omega=1:floor(ntpad/2)
    for omega=1:win
%     if omega>=nfmin && omega<=nfmax
        deta_seismo_v11(omega+1,:)=deta_seismo_v11(omega+1,:).*deta_seismo_v44(omega,:);
%     else
%         deta_seismo_v11(omega+1,:,:)=0;
%     end
   deta_seismo_v11(ntpad+1-omega,:)=conj(deta_seismo_v11(omega+1,:));
end
%end


%%% Left side part
%if is>=round(w/M)+floor(m/M)+1
    
    nend=refsx(is)-offset-m;
  if nend<1
      nend=0;
  end
  newoffset=refsx(is)-m-nend;
   deta_seismo_vv(:,:)=(ifft(deta_seismo_v11(:,:)));
   seismo_v_d1(:,nend+1:refsx(is)-m)=deta_seismo_vv(1:nt,1:newoffset);
end
%clear deta_seismo_vv deta_seismo_v1 deta_seismo_v11 deta_seismo_v4 deta_seismo_v44 deta_cc deta_cc1 y y1 tempdata
end