function Prob1_04()
    close all; clear all;
    
    %User configurable parameters
    nCurves     = 3;    %Number of different curves to be tested
    nPoints     = 2000; %Number of points in a curve boundary
    nPlots      = 3;    %Determines how many plots will be displayed per curve
    timeStep    = 1;
    
    %Developer configurable parameters
    gapPlot     = 100;  %Num iterations between successive curve-plot
    
    %Compute number of iterations based on configured parameters
    nItr        = gapPlot * nPlots;
    
    %Loop over all curves
    for curveID=1:nCurves
        %Generate a random curve
        [x, y, rFine, thetaFine] = get_rand_curve(nPoints); 

        %Draw initial curve
        figure, 
        pltID = 1;
        plot(x, y);
        hold on;
        
        %Loop over iterations for a curve
        for itrID=1:nItr
            %Compute curvature K
            dx      = gradient(x);
            ddx     = gradient(dx);
            dy      = gradient(y);
            ddy     = gradient(dy);
            ds      = sqrt(dx.^2 + dy.^2);
            K       = (dx .* ddy - ddx .* dy)./power((dx.^2 + dy.^2), 1.5);

            %Compute Normal to the curve
            N       = [-dy./ds, dx./ds];
            
            %Compute the amount by which coordinates has to change in next iteration
            delta   = -[K .* N(:, 2), K .* N(:, 1)];

            % Update coordinates
            x       = (x - timeStep * delta(:, 1));
            y       = (y - timeStep * delta(:, 2));            
            rFine   = sqrt(x.^2 + y.^2);
            [x y]   = get_updated_curve(rFine, thetaFine);
            
            %Display each curve nPlots times on screen
            if (0 == mod(itrID, gapPlot))
                plot(x, y);
                pltID = pltID + 1;
            end            
        end
        
        if(nPlots == 1)
            legend('Init', '1');
        elseif(nPlots == 2)
            legend('Init', '1', '2');
        elseif(nPlots == 3)
            legend('Init', '1', '2', '3');
        end

        hold off;
    end
end %Main function

%{
To avoid generating a curve whose boundary has large fluctuations, we first
create a curve with small number of boundary points (Coarse) and then interpolate
the points within the boundary (Fine) with the help of a cubic spline.
%}
function [x, y, rFine, thetaFine] = get_rand_curve(N)
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
end

%{
Get updated curve from rFine and thetaFine values
%}
function [x, y] = get_updated_curve(rFine, thetaFine)
    x           = rFine .* cos(thetaFine);
    y           = rFine .* sin(thetaFine);
end