function kappa = getKappa()
global x nv ne m1 m2 refLen

tangentL = zeros(ne,3); % tangent
for c=1:ne
    dx = x(4*c+1:4*c+3) - x( 4*(c-1) + 1: 4*(c-1)+3);
    tangentL(c,1:3) = dx / norm(dx);
end

kb = zeros(nv, 3);
for c=1:nv
    if c==1 || c==nv
        kb(c,:) = [0 0 0];
    else
        t0 = tangentL(c-1,:);
        t1 = tangentL(c,:);
        kb(c,:) = 2.0 * cross(t0, t1) / (1.0 + dot(t0, t1));
    end
end

% Compute Kappa (based on vertex)
kappaVertex = zeros(nv, 2);
for c=2:ne
    
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);    
    kappa1 = 0.5 * dot( kb(c,:), m2e + m2f);
    kappa2 = -0.5 * dot( kb(c,:), m1e + m1f);
    kappaVertex(c,1) = kappa1;
    kappaVertex(c,2) = kappa2;
end

kappa = zeros(ne, 2); % Edge based
for c=1:ne
    if c==1
        kappa(c,:) = 0.5 * kappaVertex(c,:) / refLen(c);
    else
        kappa(c,:) = 0.5 * (kappaVertex(c,:)+kappaVertex(c+1,:))/refLen(c);        
    end
end

end
