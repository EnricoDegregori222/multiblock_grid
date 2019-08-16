function [lh2] = strio(fileTEC, lh2, outputs)
    for i = 1:length(outputs)
        lh1=lh2;        lh2=lh1+1;
        fwrite(fileTEC, int32(unicode2native(outputs(i))), 'int');
    end
    lh1=lh2;        lh2=lh1+1;
    fwrite(fileTEC, 0, 'int');
end