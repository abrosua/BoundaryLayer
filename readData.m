function [xb,yb,m,mp1,Uinf,alpha] = readData(filename)
  % This function reads data from given airfoil coordinate data and,
  % reads a user-based input for freestream velocity and angle of attack.
  
  % Specifying a coordinates (xb,yb) of boundary points on airfoil surface.
  % The last point coincides with the first.
  fileID = fopen(filename, 'r');
  formatSpec = '%f';
  sizeA = [2 inf];
  A = fscanf(fileID, formatSpec, sizeA);
  fclose(fileID);
  
  xy = A';
  % Defining column-wise coordinate
  xb0=xy(:,1)'; yb0=xy(:,2)';

  % Calculating total panel point 
  m   = length(xb0)-1; % TE & LE is treated as 1 point 
  mp1 = m+1;
  
  % Reverse point index, clock-wise direction
  for i=1:mp1
      xb(i)=xb0(mp1+1-i);
      yb(i)=yb0(mp1+1-i);
  end

  fprintf('Input parameters:\n');
  % setting freestream velocity through user input by a control condition 
  Uinf = input('Freestream velocity (1-100m/s) : ');
  while (Uinf<=0) || (Uinf>100)
    Uinf = input('Incorrect input! Freestream velocity (1-100 m/s): ');
  end
  
  % setting up the angle of attack(AoA)
  alphaDeg = input('Angle of attack (-5 to 10 deg): ');
  while (alphaDeg>10) || (alphaDeg<-4)
    alphaDeg = input('Incorrect input! Angle of attack (-5 to 10deg): ');
  end
  alpha = deg2rad(alphaDeg); % convert into radians

end
