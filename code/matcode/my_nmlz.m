function u = my_nmlz(u)
% normalize
maxu = max(max(u));
minu = min(min(u));
if maxu < -minu
    u = u / minu;
else
    u = u / maxu;
end
