function [A, B, tc] = tgauss(n)
% 代码生成兼容版本的 tgauss(n)
% n = 2 到 80 时，分别调用 tgauss2, tgauss3, ..., tgauss80
if n == 2
    [A, B, tc] = tgauss2();
elseif n == 3
    [A, B, tc] = tgauss3();
elseif n == 4
    [A, B, tc] = tgauss4();
elseif n == 5
    [A, B, tc] = tgauss5();
elseif n == 6
    [A, B, tc] = tgauss6();
elseif n == 7
    [A, B, tc] = tgauss7();
elseif n == 8
    [A, B, tc] = tgauss8();
elseif n == 9
    [A, B, tc] = tgauss9();
elseif n == 10
    [A, B, tc] = tgauss10();
else
    error('Unsupported n value: %d', n);
end
end
%%
% function [A,B,tc]=tgauss(n)
% %A=zeros(n+1,n+1);
% eval(strcat('tgauss',num2str(n)));
% end