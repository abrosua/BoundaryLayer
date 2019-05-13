function [x, y, gama, v, cp] = panelMethod(xb,yb,m,mp1,alpha,bc)
  % NOTE:
  % this subroutine solves a set of algebraic equations c(i,j)*x(j) = a(i),
  % i=1,2,...,N. it is taken from p.114 of chow (1979)

  % Geometry parameterization
  % coordinates (x,y) of control point and panel length s are computed
  %     for each of the vortex panels.
  % rhs represents the right-hand side of eq.(15.47).
  
  for i = 1:m;
    ip1       = i+1; % index of next panel
    x(i)      = 0.5*(xb(i)+xb(ip1)); % x coordinate of panel's midpoint
    y(i)      = 0.5*(yb(i)+yb(ip1)); % y coordinate of panel's midpoint
    s(i)      = sqrt( (xb(ip1)-xb(i))^2 + (yb(ip1)-yb(i))^2); % length
    theta(i)  = atan2((yb(ip1)-yb(i)),(xb(ip1)-xb(i)));
    sine(i)   = sin(theta(i)); % sin(tetha)
    cosine(i) = cos(theta(i)); % cos(tetha)
    rhs(i)    = sin(theta(i)-alpha);
  end

  for i = 1:m
    for j = 1:m
      if i == j
        % representing Kutta condition
        cn1(i,j) = -1.0;
        cn2(i,j) =  1.0;
        ct1(i,j) = pi/2;
        ct2(i,j) = pi/2;
      else
        a = -(x(i)-xb(j))*cosine(j) - (y(i)-yb(j))*sine(j);
        b = (x(i)-xb(j))^2 + (y(i)-yb(j))^2;
        c = sin(theta(i)-theta(j));
        d = cos(theta(i)-theta(j));
        e = (x(i)-xb(j))*sine(j) - (y(i)-yb(j))*cosine(j);
        f = log (1+s(j)*(s(j)+2*a)/b);
        g = atan2(e*s(j),b+a*s(j)); 
        p = (x(i)-xb(j))*sin(theta(i)-2*theta(j)) +...
            (y(i)-yb(j))*cos(theta(i)-2*theta(j));            
        q = (x(i)-xb(j))*cos(theta(i)-2*theta(j)) -...
            (y(i)-yb(j))*sin(theta(i)-2*theta(j));
        cn2(i,j) = d+0.5*q*f/s(j) - (a*c+d*e)*g/s(j);
        cn1(i,j) = 0.5*d*f +c*g -cn2(i,j);
        ct2(i,j) = c+0.5*p*f/s(j) + (a*d-c*e)*g/s(j);
        ct1(i,j) = 0.5*c*f-d*g-ct2(i,j);
      end
    end
  end
  
  % compute influence coefficients in eq (5.47) and (5.49).
  for i = 1:m;
    an(i,1) = cn1(i,1);
    an(i,mp1) = cn2(i,m);
    at(i,1) = ct1(i,1);
    at(i,mp1) = ct2(i,m);
    for j = 2:m
      an(i,j) = cn1(i,j) + cn2(i,(j-1));
      at(i,j) = ct1(i,j) + ct2(i,(j-1));
    end
  end
  
  an(mp1,1) = 1.0;
  an(mp1,mp1) = 1.0;
  for j = 2:m
    an(mp1,j) = 0.0;
  end
  rhs(mp1) = 0.0;
  % solve eq.(5.47) for dimensionless strength gama using cramer's rule.
  % compute and print dimensionless velocity and pressure coefficient at
  %     every control point.
  gama=an\(rhs' + bc);
  
  for i = 1:m
    v(i) = cos(theta(i)-alpha);
    for j = 1:mp1
        v(i) = v(i) +at(i,j)*gama(j) ;
        cp(i) = (1- v(i)^2); % flipping the plot 
    end
  end

end