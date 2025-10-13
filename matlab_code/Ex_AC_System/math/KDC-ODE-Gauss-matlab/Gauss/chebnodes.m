function [A,B,tc]=chebnodes(n)
which = 0;

if (which == 1)
  th=-pi:pi/n:0; tc=(1 + cos(th'))/2;
  A =tgalo(n);
else
  [A,B,tc]=tgauss(n);
end
% th=[0:-pi*2/(2*n-1):-pi -pi]; tc=(1 + cos(th(end:-1:1)'))/2;
% A=tgalo(n);
return