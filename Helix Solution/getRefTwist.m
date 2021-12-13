function refTwist = getRefTwist(d1, tangent, refTwist)

ne = length(refTwist);
% Twist = zeros(ne-1, 1);

for c=2:ne
    u0 = d1(c-1,:);
    u1 = d1(c,:);
    t = tangent(c,:);
    t0 = tangent(c-1,:);
    
    ut = parallel_transport(u0, t0, t);
    ut = rotateAxisAngle(ut, t, refTwist(c) );
    
    % get signed angle
    sgnAngle = signedAngle(ut, u1, t);
    refTwist(c) = refTwist(c) + sgnAngle;
    
%     %get twist
%     Twist(c) = x(4*(c+1)) - x(4*c) + refTwist(c);
    
    
end

end
