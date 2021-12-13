function stretch = getStretch()
global x ne refLen

edgeLen = zeros(ne, 1);
for c=1:ne
    dx = x(4*c+1:4*c+3) - x(4*(c-1)+1:4*(c-1)+3);
    edgeLen(c) = norm(dx);
end

stretch = edgeLen./refLen - 1;
