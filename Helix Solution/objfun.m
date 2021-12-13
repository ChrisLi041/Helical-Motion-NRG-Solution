function [xUncons, error] = objfun(xUncons)

%% Global variables
global m dt x0 uUncons xCons
global unconsInd consInd
global Fb Fs Ft Fg tol maximum_iter
global d1 refTwist ScaleSolver

%%
x0Uncons = x0(unconsInd); % previous position
mUncons = m(unconsInd); % mass vector sans constrained dofs
mMat = diag(mUncons); % convert mass vector to mass matrix

%%
iter = 0; % number of iterations
normf = tol*ScaleSolver*10; % norm of function value (initialized to a value higher 
% than tolerance) 
error = 1; % Start with a 'good' simulation (error=1 means no error)
while (normf>tol*ScaleSolver)
    xCurrentIterate = x0;
    xCurrentIterate(consInd) = xCons; 
    xCurrentIterate(unconsInd) = xUncons; % dofs at this iteration
    
    % Update directors & get reference twist
    tangent = computeTangent(xCurrentIterate);
    [d1Iterate, d2Iterate] = computeTimeParallel(d1, x0, xCurrentIterate);
    refTwistIterate = getRefTwist(d1Iterate, tangent, refTwist); % Compute reference twist
    theta = xCurrentIterate(4:4:end);
    [m1, m2] = computeMaterialDirectors(d1Iterate, d2Iterate, theta);
    
    % Get forces
    [Fb, Jb] = getFb(xCurrentIterate, m1, m2);
    [Fg, Jg] = getFg(xCurrentIterate);
    [Fs, Js] = getFs(xCurrentIterate);
    [Ft, Jt] = getFt(xCurrentIterate, refTwistIterate);

    Forces = (Fb + Fs + Ft + Fg);
    Forces = Forces(unconsInd);
    
    % Equation of motion
    f = mUncons .* (xUncons - x0Uncons)/dt^2 - mUncons.*uUncons/dt - Forces;
    %% Force
    
    % Manipulate the Jacobians
    Jforces = Jb + Js + Jt + Jg;
    Jforces = Jforces(unconsInd, unconsInd);
%     J= mMat/dt^2 - mMat/dt - Jforces;
    J= mMat/dt^2 - Jforces;
    
    % Sir Newton's update
    xUncons = xUncons -  J\f;
    
    % Get the norm
    normfNew = norm(f);
    
    % Update iteration number
    iter = iter+1;
    fprintf('Iter=%d, error=%e\n', iter-1, normfNew);
    normf = normfNew;
    
    if (iter>maximum_iter)
        error = -1; % return with an error signal
        return
    end
end

