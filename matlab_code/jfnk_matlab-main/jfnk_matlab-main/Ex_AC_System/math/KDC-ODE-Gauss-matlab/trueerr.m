function err=trueerr(y,t)

%t is a scalar time.
%y is a row vector
%evaluate the true err at t if you know the eact solution at t,
%otherwise let err=10 (a dummy value) is ok.
global iformulation;
switch iformulation
  case {0,1}
    err=max(abs(y-analsolu(t)));
  otherwise
end;

return;