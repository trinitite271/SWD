clc
clear 
close all
load topo
load topo_1
[nz,nx]=size(vs);
vs2=719*ones(35,nx);
%  for i=1:nx
%      vs1(1:nz-fsz(i)+1,i)=vs(fsz(i):end,i);
%  end

 for i=1:nx
      vs2(fsz(i)+1:nz,i)=vs(1:end-fsz(i),i);
 end
vs2(:,101:end)=vs2(:,1:20);
vs2(19:26,:)=flipud(vs2(19:26,:));
% load topo_result_NEW
% vs2(:,22:32)=vs(:,22:32);
figure
imagesc(smooth2a(vs2,1,3))
colormap(jet)
% figure
% %vs=(vs1(1:35,:));
% imagesc(vs)
% colormap(jet)