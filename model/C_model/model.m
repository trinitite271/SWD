clc
clear
mi=300;
ma=1500;
z=5;z1=21;
nx=120;
for i=1:nx
      vs(1:z,i)=linspace(mi,ma,z);
end
for i=1:nx
      vs(z+1:z1+z,i)=linspace(ma,ma+200,z1);
end
%vs(12:15,:)=500;
save Hi.mat vs