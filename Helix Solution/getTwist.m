function twist = getTwist()
global x ne refTwist refLen
theta_f = x(4:4:end);
theta_e = [0; theta_f(1:end-1)];
deltam = theta_f - theta_e;
getUndeformedTwist = zeros(ne,1);
twist = zeros(ne, 1);
for c=1:ne
    twist(c) = (deltam(c) + refTwist(c) - getUndeformedTwist(c))/refLen(c);
end
