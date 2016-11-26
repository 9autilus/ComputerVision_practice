function main()
    close all; clear all;

    %User configurable parameters
    %inputImageFileName = 'rose.jpg';
    inputImageFileName = 'pacman.jpg';
    
    %Developer configurable parameters
    nPoints     = 40;
    maxItr      = 120;  %max iteration to be performed
    gapPlot     = 40;   %How many iteration between successive plots
    gaussSigma  = 3;    %std_dev for guassian smoothing of image for external forces
    alpha       = 0.1;%0.5;  %Value based on hit-and-trial. Values are specific to this image.
    beta        = 0.01; %1   %Value based on hit-and-trial. Values are specific to this image.
    delta       = 5000; %Value based on hit-and-trial. Values are specific to this image.
    time_step   = .02;  %Value based on hit-and-trial. Values are specific to this image.
    linewidth   = 1.5;
    
    %Get a grayscale image
    imgColor = double(imread(inputImageFileName))/255.0;
    img      = rgb2gray(imgColor);
    [x, y]   = get_points(inputImageFileName, nPoints);
    x = x(:); y = y(:); %Convert to column vector
    x(end) = x(1); y(end) = y(1); %Make curve closed
    
    % Display input image and starting points of curve
    %figure, imshow(imgColor), hold on, plot(x, y, 'red', 'LineWidth',linewidth), hold off;
    [u, v] = get_grad_vector_field(img);
 
    return;
    
    %Get the gradient term of Energy Equation
    smoothImg   = imgaussfilt(img, gaussSigma);
    [gradImgx, gradImgy] = gradient(smoothImg);
    squaredNormGrad = (gradImgx .^2 + gradImgy .^2);
    [gradSquaredNormGradx, gradSquaredNormGrady] = gradient(squaredNormGrad);
    
    D2 = get_D2(nPoints); %Second Order discrete derivation Matrix D2
    D4 = get_D4(nPoints); % Fourth Order discrete derivation Matrix D4
    
    A       = eye(nPoints) + time_step * (-alpha * D2 + beta * D4);
    A_inv   = inv(A);
    
    hold on;
    for itrID = 1:maxItr
        interp_Px = interp2(gradSquaredNormGradx, x, y);
        interp_Py = interp2(gradSquaredNormGrady, x, y);
        interp_Px = interp_Px(:); interp_Py = interp_Py(:);  %Convert to column vector
    
        x = A_inv * (x + time_step * delta * interp_Px);
        y = A_inv * (y + time_step * delta * interp_Py);
        
        x(end) = x(1); y(end) = y(1); %Make curve closed by joining first and last points
        
        if(mod(itrID, gapPlot) == 0)
            plot(x, y, 'LineWidth',linewidth);
        end
    end
    
    if(3 == int32(maxItr/gapPlot))
        legend('Init', '1', '2', '3');
    end
    
    hold off;
end

function [x, y] = get_points(inputImageFileName, nPoints)
    if(1)   %Use hard-coded point values
        assert((nPoints == 40), 'Error: nPoints must be 40 with hard-coded values.');
        if strcmp(inputImageFileName, 'pacman.jpg')
            x = [120.3873  135.3249  150.2626  165.4577  175.7596  184.2586  172.1539 161.0795  150.5201  142.0211  135.3249  142.5362  151.2928  157.7314 165.9728  176.2746  181.9406  172.6690  157.7314  145.1117  123.4779 108.7978   94.3753   76.8622   65.5302   56.2586   47.5020   44.4115  44.4115   48.5322   56.0010   64.7575   75.3169   82.7857   89.9970 94.1177   99.5262  102.8742  106.9950  112.9185];
            y = [42.1057   43.9245   48.9263   55.2922   64.3863   73.9352   88.9405  98.4893  108.9476  117.1323  124.4076  132.1377  139.4130  147.5977 155.3277  166.6954  172.1519  183.9742  193.9778  199.8890  202.6172 202.1625  198.0702  189.4307  178.5178  166.6954  148.5071  127.1359 105.3099   89.8499   79.3917   67.1146   57.1110   53.0187   48.9263 48.4716   46.1980   46.1980   44.8339    43.9245];
        elseif strcmp(inputImageFileName, 'rose.jpg')
            x = [119.6147  123.2203  129.9165  144.0815  153.8682  196.8783  234.2223  253.7958  217.7394  166.4879  170.0936 144.8541  128.3712  157.9889  168.2907  178.8501  174.2143  149.7475  126.3109  103.3893   83.0433   70.4235 58.5765   62.9547   77.6348   68.1056   82.0131  102.6167  112.6610  111.6308  110.8581  102.3592   83.3008  60.1217   29.9889   11.7032   41.8360   57.0312   74.0292  110.0855];
            y = [328.4432  311.1359  297.4094  263.3917  261.0044  239.5195  231.1643  229.9707  209.6794  205.5018  179.8393 189.9849  167.9032  135.0790   92.1092   80.1732   48.5426   31.8321   25.8641   21.0897   22.8801   33.6226 59.2851   76.5924  105.8357  124.9334  139.8535  146.4183  147.0151  173.2744  196.5497  212.6634  206.0986 232.3579  253.8428  272.3437  300.3934  265.1821  299.1998  302.1838];
        else
            assert(0, 'Error: Hard-coded values only work with pacman.jpg or rose.jpg currently.');
        end
    else    %Take points mannually
        fprintf('Taking %d points mannually from user.\n', nPoints);
        xy2d = Get2DPoints(inputImageFileName, nPoints);
        x = xy2d(1, :);
        x = xy2d(2, :);
        disp(xy2d);
    end
end

function D2 = get_D2(nPoints)
    V   = -2 * ones(nPoints, 1);
    W   = 1 * ones(nPoints-1, 1);
    D2  = diag(V)+diag(W, 1) + diag(W, -1);
    
    D2(1, nPoints) = 1;
    D2(nPoints, 1) = 1;
end

function D4 = get_D4(nPoints)
    V   = 6 * ones(nPoints, 1);
    W   = -4 * ones(nPoints - 1, 1);
    U   = ones(nPoints - 2, 1);
    D4  = diag(V) + diag(W, 1) + diag(W, -1) + diag(U, 2) + diag(U, -2);
    
    D4(1, nPoints)      = -4;
    D4(nPoints, 1)      = -4;
    D4(1, nPoints-1)    = 1;
    D4(nPoints - 1, 1)  = 1;
    D4(2, nPoints)      = 1;
    D4(nPoints, 2)      = 1;
end

function [u, v] = get_grad_vector_field(img)
    time_step = .1;
    mu = 1;
    maxItr = 100;
    gaussSigma  = 1;    %std_dev for guassian smoothing of image for external forces
    
    ht = size(img, 1);
    wd = size(img, 2);
    
    u = zeros(ht, wd);
    v = zeros(ht, wd);
    
    smoothImg   = imgaussfilt(img, gaussSigma);
    
    [fx, fy]    = gradient(smoothImg);
    b           = fx.^2 + fy.^2;
    c1          = b .* fx;
    c2          = b .* fy;

    r = mu * time_step;
    
    for itrID = 1:maxItr
        u = (1 - time_step * b) .* u ...
            + r * ([u(2:end,:); zeros(1,wd)] + [zeros(1,wd); u(1:end-1, :)] + [u(:, 2:end), zeros(ht,1)] + [zeros(ht,1), u(:, 1:end-1)] - 4 * u) ...
            + time_step * c1;
        v = (1 - time_step * b) .* v ...
            + r * ([v(2:end,:); zeros(1,wd)] + [zeros(1,wd); v(1:end-1, :)] + [v(:, 2:end), zeros(ht,1)] + [zeros(ht,1), v(:, 1:end-1)] - 4 * v) ...
            + time_step * c2;
    end
    
    figure,
    subplot(1,2,1), imshow(smoothImg);
    hold on;
    quiver(fx,fy);
    hold off;
    
    subplot(1,2,2), imshow(smoothImg);
    hold on;
    quiver(u,v);
    hold off;
    
end