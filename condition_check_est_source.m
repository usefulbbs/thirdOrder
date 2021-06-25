function res = condition_check_est_source(Xf,ep)
res=1;
for i = 1:size(Xf,2)-1
    if(Xf(i)==0)
        res =0;
        break;
    end
    for j = i+1:size(Xf,2)
        if(abs(angle(Xf(j)/Xf(i))) > ep & abs(angle(Xf(j)/Xf(i)) < pi-ep))
            res=0;
            break;
        end
    end
end
end