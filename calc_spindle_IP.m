function [X Y] = calc_spindle_IP(coeffs,fs,min_freq,max_freq)
num_coeffs = length(coeffs);
if num_coeffs < 5; error('There should be 5 coefficients for spindle parameter estimation'); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3; min_freq = 12; max_freq = 14; end
if nargin < 2; fs = 200; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to be careful about the order
cx2 = coeffs(1); if cx2 <= 0; error('Non-Convex Problem'); end
cy2 = coeffs(2); if cy2 <= 0; error('Non-Convex Problem'); end
cxy = coeffs(3);
cx  = coeffs(4);
cy  = coeffs(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ly = -.99^2;
uy = -.95^2;
% lx2 = -4 * ly * cos(2*pi*min_freq/fs);
% ux2 = -4 * uy * cos(2*pi*max_freq/fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = - a^2; X = 2a cos(2pi f/fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization Problem : (Non-Convex)
% minimize       cx2   X^2 +  cy2  Y^2  + cxy XY + cx X + cy Y
% subject to    lambda3:  -4Y cos(2*pi max_freq/fs)^2    <=  X^2 <= -4Y cos(2*pi  min_freq/fs)^2 : lambda4                  
%               lambda1:                            ly  <=  Y <=   uy                        : lambda2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 cy2 Y + cxy X + cy + lambda2 -  lambda1 + 4 lambda4 cos(2*pi  min_freq/fs)^2 - 4 lambda3 cos(2*pi  max_freq/fs)^2 = 0
% 2 cx2 X + cxy Y + cx + 2 lambda4 X - 2 lambda3 X = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complementary Slackness
% (2 cy2 Y + cxy X + cy)(Y-ub)(Y-lb)(X^2 + 4Y cos(2*pi  min_freq/fs) )^2( X^2 + 4Y cos(2*pi  max_freq/fs) )^2 = 0
% (2 cx2 X + cxy Y + cx)(X^2 + 4Y cos(2*pi  min_freq/fs) )^2( X^2 + 4Y cos(2*pi  max_freq/fs) )^2 = 0


%%%%%%%%%%%%%%%%%%%%%%%%%
A = [2*cx2 cxy; cxy 2*cy2]\[-cx; -cy];
X = A(1); Y = A(2);
%%%%%%%%%%%%%%%% If not satisfying constraints go and check the edges of the region
if check_constraints(X,Y,uy,ly,min_freq,max_freq,fs)
    
    % case 1
    Y(1) = ly;
    X(1) = (-cx-cxy*Y(1))/(2*cx2);
    a(1) = f(X(1),Y(1),uy,ly,min_freq,max_freq,fs,coeffs);
    
    Y(2) = ly;
    X(2) = sqrt(-4*Y(2))*cos(2*pi*  min_freq/fs);
    a(2) = f(X(2),Y(2),uy,ly,min_freq,max_freq,fs,coeffs);
    
    Y(3) = ly;
    X(3) = sqrt(-4*Y(3))*cos(2*pi*  max_freq/fs);
    a(3) = f(X(3),Y(3),uy,ly,min_freq,max_freq,fs,coeffs);
    
    % case 2
    Y(4) = uy;
    X(4) = (-cx-cxy*Y(4))/(2*cx2);
    a(4) = f(X(4),Y(4),uy,ly,min_freq,max_freq,fs,coeffs);
     
    Y(5) = uy;
    X(5) = sqrt(-4*Y(5))*cos(2*pi*  min_freq/fs);
    a(5) = f(X(5),Y(5),uy,ly,min_freq,max_freq,fs,coeffs);
    
    Y(6) = uy;
    X(6) = sqrt(-4*Y(6))*cos(2*pi*  max_freq/fs);
    a(6) = f(X(6),Y(6),uy,ly,min_freq,max_freq,fs,coeffs);
    
    % case 3
    syms z;
    S = solve(2*cx2*z - 3/4*cxy*z^2/cos(2*pi*min_freq/fs)^2 + cy2*z^3/4/cos(2*pi*min_freq/fs)^4+cx-cy*z/2/cos(2*pi*min_freq/fs)^2 ...
        ,z>0,z<2*(cos(2*pi*min_freq/fs)));
    
    for i=1:length(S)
    X(end+1) = S(i);    
    Y(end+1) = -X(end).^2/4/cos(2*pi*min_freq/fs)^2;
    a(end+1) = f(X(end),Y(end),uy,ly,min_freq,max_freq,fs,coeffs);
    end
    
    S = solve(2*cx2*z - 3/4*cxy*z^2/cos(2*pi*max_freq/fs)^2 + cy2*z^3/4/cos(2*pi*max_freq/fs)^4+cx-cy*z/2/cos(2*pi*max_freq/fs)^2 ...
        ,z>0,z<2*(cos(2*pi*min_freq/fs)));
    
        for i=1:length(S)
    X(end+1) = S(i);    
    Y(end+1) = -X(end).^2/4/cos(2*pi*max_freq/fs)^2;
    a(end+1) = f(X(end),Y(end),uy,ly,min_freq,max_freq,fs,coeffs);
        end

        [m, I] = min(a);
        X = X(I);
        Y = Y(I);
%         check_constraints(X,Y,uy,ly,min_freq,max_freq,fs)
%         a = sqrt(-Y)
%         f = fs*acos(X/2/a)/2/pi
end

end

function c = check_constraints(X,Y,uy,ly,min_freq,max_freq,fs)
c = (Y > uy) || (Y < ly) || (X^2 > -4*Y* cos(2*pi*  min_freq/fs)^2) || (X^2 < -4*Y* cos(2*pi*  max_freq/fs)^2);
end

function a = f(X,Y,uy,ly,min_freq,max_freq,fs,coeffs)

if ~check_constraints(X,Y,uy,ly,min_freq,max_freq,fs); a = coeffs(1)*X^2+coeffs(2)*Y^2+coeffs(3)*X*Y+coeffs(4)*X+coeffs(5)*Y;
else a = Inf;   end
    

end







