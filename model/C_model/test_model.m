clc
clear
load model
vp=vp(1:2:68,1:2:end);
vs=vp;
for i=1:34
    for j=1:180
        if vs(i,j)==666
            vs(i,j)=450;
        end
    end
end
% for i=1:34
%     for j=1:180
%         if vs(i,j)==736
%             vs(i,j)=500;
%         end
%     end
% end
vs=smoothn(vs,'robust');
vp=vs*2;

save low_vel1.mat vp vs