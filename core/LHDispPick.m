
function [pre_mask,res,cr_pre] = LHDispPick(ml,npair,vmin,cr_0,line_mask,offset,method)
mode=1;
if offset>100
    [~,line_start_stop,line_part]=pickdisp_xl(ml);
    [line_start_stop,line_part]=CEDC(ml);
else
    [~,p_max] = max(ml);
    llp = length(p_max);
    line_part{1}=[p_max;1:llp];
    line_start_stop=[1,llp];
end

if mode==1
    if method==1
    [idx,~,~]=pickFDC1(ml,10,2,4);
    cr_pre = idx' + vmin;
    pre_mask = ones(npair,1);
    res = mean(((cr_pre - cr_0).^2).*line_mask);
    else
    cr_pre = 1.*ones(npair,1);
    cr_pre(line_start_stop(mode,1):line_start_stop(mode,2),1) = line_part{mode}(1,:)+vmin;
    res = mean(((cr_pre - cr_0).^2).*line_mask);
    pre_mask = ones(npair,1);
    pre_mask = pre_mask .* line_mask;
    end
else 
    if length(line_part)>=mode
            cr_pre = 1.*ones(npair,1);
            cr_pre(line_start_stop(mode,1):line_start_stop(mode,2),1) = line_part{mode}(1,:)+vmin;
            res = mean(((cr_pre - cr_0).^2).*line_mask);
            pre_mask = ones(npair,1);
            pre_mask = pre_mask .* line_mask;
    else
        cr_pre = 1.*ones(npair,1);
        res = mean(((cr_pre - cr_0).^2).*line_mask);
        pre_mask = ones(npair,1) .* line_mask;
    end
end
    
end