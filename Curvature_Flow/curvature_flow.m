function Prob1_00()
    close all; clear all;
    
    xlimit = 20; 
    ylimit = 20;
    
    nCurves     = 1;
    nPoints     = 200;
    nTimeSteps  = 3;
    
    gapTimeStep = 5;
    
    nItr        = gapTimeStep * nTimeSteps; %Total number of iterations
    
    for curveID=1:nCurves      
        [x, y] = gen_rand_curve(nPoints); 
        x = [x(end-1) x x(2)];
        y = [y(end-1) y y(2)];

        %x = smooth(x);
        %y = smooth(y);
        
        %disp([x y]);    
        
        figure, 
        pltID = 1;
        axis([-xlimit xlimit -ylimit ylimit]);
        subplot(2, 2, pltID);
        plot(x, y);
        title(pltID);
        pltID = pltID + 1;
        %plot(x(3:end-2), y(3:end-2));
        
        for itrID=1:nItr
            %fprintf('Iteration: %d\n', itrID);
            %Compute the curvature
            dx      = gradient(x);
            ddx     = gradient(dx);
            dy      = gradient(y);
            ddy     = gradient(dy);
            ds      = sqrt(dx.*dx + dy.*dy);
            num     = dx .* ddy - ddx .* dy;
            denom   = dx .* dx + dy .* dy;
            denom   = sqrt(denom);
            denom   = denom .* denom .* denom;
            K            = num ./ denom;
            %K(denom < 0) = NaN;
            %K(isnan(K))  = 0;

            N       = [-dy./ds, dx./ds];
            
            %for i=1:size(N,1)
            %    N(i, :) = N(i, :)/sqrt(N(i,1)*N(i,1) + N(i,2)*N(i,2));
            %end
            
            K= abs(K);
            
            %delta   = -[K .* N(:, 1), K .* N(:, 2)];%original
            delta   = [K .* N(:, 2), K .* N(:, 1)];
            delta(isnan(delta)) = 0;

            x = x + delta(:, 1);
            y = y + delta(:, 2);
            
            %Display each curve nTimeSteps times on screen
            if (0 == mod(itrID, gapTimeStep))
                subplot(2, 2, pltID), plot(x(2:end-1), y(2:end-1));
                title(pltID);
                axis([-xlimit xlimit -ylimit ylimit])
                pltID = pltID + 1;
                %fprintf('Size K: (%d, %d), Size N: (%d, %d)\n', size(K,1), size(K,2), size(N,1), size(N,2));
            end            
        end
        
        
        %{
        % In each iteration find the normal and curvature at each point and also 
        % disp as Curv*Normal
        itr = 300;
        for i = 1:itr
            dx = gradient(x);
            ddx = gradient(dx);
            dy = gradient(y);
            ddy = gradient(dy);

            temp = dx .* ddy - ddx .* dy;
            ds = sqrt(dx.*dx + dy.*dy);

            K = temp./(ds.*ds.*ds);
            K(isnan(K)) = 0;

            Nor = [-dy./ds,dx./ds];
            Disp = -[K.*Nor(:,1),K.*Nor(:,2)];
            Disp(isnan(Disp)) = 0;

            %x = smooth(x + Disp(:,1));
            %y = smooth(y + Disp(:,2));
            
            x = x + Disp(:,1);
            y = y + Disp(:,2);

            %Plot the curve at 3 time steps of 100 itr.
            if(i == 100 || i == 200 || i == 300)
                figure,
                plot(y,x,'r')
                %axis([-50 300 -50 300])
            end   
        end
        %}
        
    end
end %Main function




function [x, y] = gen_rand_curve(N)
    n           = 15;
    t           = linspace(0, 2 * pi, n);
    tt          = linspace(0, 2 * pi, N);

    r           = 1 + 5 * rand(size(t));
    r(end)      = r(1);
    rr          = spline(t, [0 r 0], tt);
    x           = rr .* cos(tt);
    y           = rr .* sin(tt);

    %% display
    %figure, plot(x,y);
end