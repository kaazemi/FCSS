function [X Y] = calc_Calcium_IP(coeffs)
num_coeffs = length(coeffs);
if num_coeffs < 5; error('There should be 5 coefficients for Calcium parameter estimation'); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to be careful about the order
cx2 = coeffs(1); if cx2 <= 0; error('Non-Convex Problem'); end
cy2 = coeffs(2); if cy2 <= 0; error('Non-Convex Problem'); end
cxy = coeffs(3);
cx  = coeffs(4);
cy  = coeffs(5);

%%%%%%%%%%%%%%%%%%%%%%%%%
A = [2*cx2 cxy; cxy 2*cy2]\[-cx; -cy];
X = A(1); Y = A(2);
%%%%%%%%%%%%%%%% If not satisfying constraints go and check the edges of the region
if check_constraints(X,Y)
    
    r(1) = abs(X+sqrt(X^2+4*Y))/2;
    r(2) = abs(X-sqrt(X^2+4*Y))/2;
    X = X*.99/max(r);
    Y = Y*.99^2/max(r)^2;
end

end

function c = check_constraints(X,Y)
c = abs(X+sqrt(X^2+4*Y))/2 > 0.999 || abs(X-sqrt(X^2+4*Y))/2 > 0.99 || X^2+4*Y < 0;
end

% function a = f(X,Y,coeffs)
% 
% if ~check_constraints(X,Y); a = coeffs(1)*X^2+coeffs(2)*Y^2+coeffs(3)*X*Y+coeffs(4)*X+coeffs(5)*Y;
% else a = Inf;   end
%     
% end







