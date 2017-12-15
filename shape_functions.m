function [ phi, dphidx, dphidy ] = shape_functions (t)

  area = t(1,1) * ( t(2,2) - t(2,3) ) ...
       + t(1,2) * ( t(2,3) - t(2,1) ) ...
       + t(1,3) * ( t(2,1) - t(2,2) );

  phi(1) =     (  ( t(1,3) - t(1,2) ) * ( p(2,1:n) - t(2,2) )     ...
                    - ( t(2,3) - t(2,2) ) * ( p(1,1:n) - t(1,2) ) );
  dphidx(1) =   - ( t(2,3) - t(2,2) );
  dphidy(1) =     ( t(1,3) - t(1,2) );

  phi(2) =     (  ( t(1,1) - t(1,3) ) * ( p(2,1:n) - t(2,3) )     ...
                    - ( t(2,1) - t(2,3) ) * ( p(1,1:n) - t(1,3) ) );
  dphidx(2) =   - ( t(2,1) - t(2,3) );
  dphidy(2) =     ( t(1,1) - t(1,3) );

  phi(3) =     (  ( t(1,2) - t(1,1) ) * ( p(2,1:n) - t(2,1) )     ... % 3. Spalte 1. Zeile +
                    - ( t(2,2) - t(2,1) ) * ( p(1,1:n) - t(1,1) ) ); % 3. Spalte 2. Zeile
  dphidx(3) =   - ( t(2,2) - t(2,1) ); % y_12 3. Spalte 2. Zeile
  dphidy(3) =     ( t(1,2) - t(1,1) ); % x_21 3. Spalte 1. Zeile
%
%  Normalize.
%
  phi(1:3)    = phi(1:3) / area;
  dphidx(1:3) = dphidx(1:3) / area;
  dphidy(1:3) = dphidy(1:3) / area;

  return
end