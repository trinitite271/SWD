clc
clear 
close all
load topo8
%fsz=fsz(1:160,:);
load tt
 load ca1
%load fig8
%load topob
%load rot_ini
%load rot_result
%fsz=13-fsz;
[nz,nx]=size(vs);
vs1=3411*ones(40,nx);
for i=1:nx
    vs1(1:nz-fsz(i)+1,i)=vs(fsz(i):end,i);
end

% for i=161:nx
%     vs1(:,i)=vs1(:,160);
% end
% for i=1:nx
%     vs1(fsz(i):end,i)=vs(1:nz-fsz(i)+1,i);
% end
vs2=(smooth2a(vs1(2:end,:),2,2));
figure
imagesc(vs2(1:end,:))
colormap(jet)

% vs=(vs1(1:35,:));imagesc(vs)
% colormap(jet)