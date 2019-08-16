function [lh2] = techead_write(fileTEC, mb, mbk, mq, lh2, lxi, let, lze, output, var, time, gridmino, gridmaxo)
if mb==0
    fwrite(fileTEC, '#', 'char*1');
    fwrite(fileTEC, '!', 'char*1');
    fwrite(fileTEC, 'T', 'char*1');
    fwrite(fileTEC, 'D', 'char*1');
    fwrite(fileTEC, 'V', 'char*1');
    fwrite(fileTEC, '1', 'char*1');
    fwrite(fileTEC, '1', 'char*1');
    fwrite(fileTEC, '2', 'char*1');
    lh1 = lh2; lh2 = lh1+1;
    fwrite(fileTEC, 1, 'int');
    lh1=lh2; lh2=lh1+1;
    fwrite(fileTEC, 1, 'int');
    lh2 = strio(fileTEC, lh2, output);
    lh1=lh2;        lh2=lh1+1;
    fwrite(fileTEC, mq, 'int');
    lh2 = strio(fileTEC, lh2, 'x');
    lh2 = strio(fileTEC, lh2, 'y');
    lh2 = strio(fileTEC, lh2, 'z');
%     for i=1:length(var)
%         cinput = var(i); lh2 = strio(fileTEC, lh2, cinput);
%     end
    for k=0:mbk
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, 299.0, 'single');
        czonet = "z00";
        czonet = strcat(czonet, int2str(k));
        lh2 = strio(fileTEC, lh2, czonet);
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, -1, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, k+1, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, time, 'double');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, -1, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, 0, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, 0, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, 0, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, 0, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, lxi(k+1)+1, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, let(k+1)+1, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, lze(k+1)+1, 'int');
        lh1=lh2;    lh2=lh1+1;
        fwrite(fileTEC, 0, 'int');
    end
    lh1=lh2;    lh2=lh1+1;
    fwrite(fileTEC, 357.0, 'single');
end
lh1=lh2;    lh2=lh1+1;
fwrite(fileTEC, 299.0, 'single');
for k=1:mq
    lh1=lh2;    lh2=lh1+1;
    fwrite(fileTEC, 1, 'int');
end
lh1=lh2;    lh2=lh1+1;
fwrite(fileTEC, 0, 'int');
lh1=lh2;    lh2=lh1+1;
fwrite(fileTEC, 0, 'int');
lh1=lh2;    lh2=lh1+1;
fwrite(fileTEC, -1, 'int');
for k=1:(mq)
    lh1=lh2;    lh2=lh1+1;
    fwrite(fileTEC, gridmino(k), 'double');
    lh1=lh2;    lh2=lh1+1;
    fwrite(fileTEC, gridmaxo(k), 'double');
end
end