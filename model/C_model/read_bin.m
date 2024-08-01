%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.

%  file:       read_bin
%  purpose:    read in binary file 
%  assumption: None
%  input:      filename   ------ filename for reading
%              n1-n5      ------ dimensions from first to fifth
%                                (Note: doon't spt dimensions higher than 5
%  output:     data       ------ data read after5 sorting

%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Oct 21, 2011

function data=read_bin(filename,n1,n2,n3,n4,n5)

fid=fopen(filename,'rb');

if nargin == 2 % readin a vector
   temp = fread(fid,n1,'float');fclose(fid);
   data = temp;
elseif nargin == 3 % readin a 2D matrix
   temp = fread(fid,n1*n2,'float');fclose(fid);
   data = reshape(temp,n1,n2);
elseif nargin == 4  % readin a 3D matrix
   temp = fread(fid,n1*n2*n3,'float');fclose(fid);
   data = reshape(temp,n1,n2,n3);
elseif nargin == 5 % readin 2 4D matrix
   temp = fread(fid,n1*n2*n3*n4);fclose(fid);
   data = reshape(temp,n1,n2,n3,n4);
elseif nargin == 6 % readin a 5D matrix
   temp = fread(fid,n1*n2*n3*n4*n5);fclose(fid);
   data = reshape(temp,n1,n2,n3,n4,n5);
elseif nargin == 1 % wrong
   disp('You have to define the dimension of your data!');
else % readin higher dimension matrix
   disp('Higher than 6 dimension matrix readin is not supported by this code!');
end

end

    
    
