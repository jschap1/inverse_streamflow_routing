% Creates a custom colormap for ISR paper
%
% Goes from green (large) to black (zero) to red (negative)

function cmap = custom_colormap()

% Define custom colormap
negativeColor = [1, 0, 0];  % Red for smaller numbers
zeroColor = [0, 0, 0];      % Black for zero
positiveColor = [0, 1, 0];  % Green for larger numbers

% Number of colormap steps (adjust as needed)
numSteps = 256;

% Create custom colormap matrix
customColormap = zeros(numSteps, 3);

% Interpolate colors for the colormap
for i = 1:numSteps
    if i <= numSteps / 2
        % Interpolate from red to black for smaller numbers
        customColormap(i, :) = (numSteps / 2 - i) / (numSteps / 2) * negativeColor + i / (numSteps / 2) * zeroColor;
    else
        % Interpolate from black to green for larger numbers
        customColormap(i, :) = (i - numSteps / 2) / (numSteps / 2) * zeroColor + (numSteps - i) / (numSteps / 2) * positiveColor;
    end
end

% Make sure values greater than or equal to zero are set to green
%customColormap(numSteps/2+1:end, :) = positiveColor;

% Create your image data (replace this with your actual data)
data = randn(100, 100);  % Example random data

% Create the imagesc plot with the custom colormap
imagesc(data);
colormap(customColormap);

% Add colorbar for reference
colorbar;
cmap = customColormap;

return