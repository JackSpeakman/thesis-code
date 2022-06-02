function c = cShift(c,s)
% shifts a colour c by s
s = s(:);

c = (1+s).*(s<=0).*c+(1-s).*(s>0).*c+s.*(s>0);
end
