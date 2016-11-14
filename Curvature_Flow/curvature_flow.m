function main()
    close all; clear all;
    
    %User configurable parameters
    animated    = 1;    %A true value cause the circle to shrink to point
    nPoints     = 2000; %Number of points in a curve boundary
    
    %Developer configurable parameters
    if(animated)
        gapPlot     = 50;   %Num iterations between successive curve-plot
        nCurves     = 1;    %Use a single curve
        timeStep    = 3;
        linewidth   = 1;
    else
        gapPlot     = 5000; %Num iterations between successive curve-plot
        nCurves     = 3;    %Number of different curves to be tested
        nPlots      = 3;    %Determines how many plots will be displayed per curve
        timeStep    = 2;
        linewidth   = 1.5;  %Width of line in plot
    end
    
    %Compute number of iterations based on configured parameters
    if(~animated)
        nItr        = gapPlot * nPlots;
    end
    
    %Loop over all curves
    for curveID=1:nCurves
        %Generate a random curve
        [x, y, thetaFine] = get_rand_curve(nPoints); 
        x = [x(end - 1); x ; x(2)];
        y = [y(end - 1); y ; y(2)];
        thetaFine = [thetaFine(end - 1); thetaFine ; thetaFine(2)];
    
        %Draw initial curve
        figure, 
        pltID = 1;
        %plot(x, y, 'blue');
        plot(x(2:end-1), y(2:end-1), 'blue', 'LineWidth',linewidth);
        hold on;
        
        %Loop over iterations for a curve
        itrID   = 1; 
        halt    = 0;
        while(~halt)
            itrID       = itrID + 1;
            
            % To make the curve appear circular we create a cycle by making 
            % first and last two points same
            x(end)      = x(3);
            x(end-1)    = x(2);
            x(1)        = x(end-2);
            
            y(end)      = y(3);
            y(end-1)    = y(2);
            y(1)        = y(end-2);            
            
            %Compute curvature K
            dx      = gradient(x);
            ddx     = gradient(dx);
            dy      = gradient(y);
            ddy     = gradient(dy);
            ds      = sqrt(dx.^2 + dy.^2);
            K       = (ds .* (dx .* ddy - ddx .* dy))./power(ds, 3);
            K(K<0)  = 0;
            
            % Update coordinates
            x       = (x - timeStep * K .* (dy));
            y       = (y - timeStep * K .* (-dx));            
            rFine   = sqrt(x.^2 + y.^2);
            [x, y]  = get_updated_curve(rFine, thetaFine);
            
            %Display each curve nPlots times on screen
            if (0 == mod(itrID, gapPlot))
                plot(x(2:end-1,1), y(2:end-1, 1), 'red','LineWidth',linewidth);
                pltID = pltID + 1;
                pause(0.1);
            end

            if(~animated)
                halt = (itrID >= nItr);
            end
        end
        
        if(~animated)
             if(nPlots == 1)
                 legend('Init', '1');
             elseif(nPlots == 2)
                 legend('Init', '1', '2');
             elseif(nPlots == 3)
                 legend('Init', '1', '2', '3');
            end
        end

        hold off; %Hold off for possible next curve
    end
end %Main function

%{
To avoid generating a curve whose boundary has large fluctuations, we first
create a curve with small number of boundary points (Coarse) and then interpolate
the points within the boundary (Fine) with the help of a cubic spline.
%}
function [x, y, thetaFine] = get_rand_curve(N)
    numCoarsePts  = 20; %Chosen by hit and trial
    thetaCoarse   = linspace(0, 2 * pi, numCoarsePts);
    thetaFine     = linspace(0, 2 * pi, N);

    % 1 is added to avoid self-intersecting curves.
    % 5 is just for a bigger scale. Chosen arbitrarily.
    rCoarse     = 5 * (1 + rand(size(thetaCoarse)));
    rCoarse(end)= rCoarse(1);
    rFine       = interp1(thetaCoarse, rCoarse, thetaFine, 'spline');
    
    x           = rFine .* cos(thetaFine);
    y           = rFine .* sin(thetaFine);
    
    %Convert to column vectors
    x           = x(:); 
    y           = y(:);
    thetaFine   = thetaFine(:);
end

%{
Get updated curve from rFine and thetaFine values
%}
function [x, y] = get_updated_curve(rFine, thetaFine)
    x           = rFine .* cos(thetaFine);
    y           = rFine .* sin(thetaFine);
end