function upscaledArray = interpolate(shortArray,longArray, plotResult)
    %inerpolate like this: https://se.mathworks.com/matlabcentral/answers/784166-interpolate-data-in-array-to-match-specified-length-of-data-in-another-vector
    m = size(shortArray,1);
    n = size(longArray,1);
    upscaledArray = zeros(n,size(shortArray,2));
    xi = (1:n)'; %x indexes
    x = linspace(1,n,m)'; 
    for i = 1:size(shortArray,2)
        array_i = interp1(x,shortArray(:,i),xi); 
        upscaledArray(:,i) = array_i; 
    end

    %option of plotting to verify the result
    %knowledge that OptiTrack operates at 120Hz and Franka robot at 1kHz used for time vectors
    if plotResult == true 
        for i = 1:size(shortArray,2)
            subplot(size(shortArray,2),1,i)
            hold on
            plot(0:1/120:1/120*(length(shortArray)-1), shortArray(:,i), 'color','r', 'LineWidth',2.5, 'DisplayName', 'shortArray')
            plot(0:1/1000:1/1000*(length(upscaledArray)-1), upscaledArray(:,i), 'color','b', 'DisplayName', 'upscaledArray')
            legend
            hold off
        end
    end

end