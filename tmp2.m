% Create data for animation
nFrames = 50;  % Number of frames
nRows = 50;    % Number of rows
nCols = 50;    % Number of columns

% Create a figure
figure;

% Loop through frames
for k = 1:nFrames
    % Generate some data (you can replace this with your actual data)
    data = sin((1:nRows)' * (1:nCols)) + k/10;  % Example of a dynamic matrix
    
    % Create the heatmap
    imagesc(data);
    
    % Adjust the color scale (optional)
    colorbar;
    
    % Set axis properties
    axis tight;
    axis off;
    
    % Update the plot
    drawnow;  % This forces MATLAB to update the figure window
    
    % Optional: Pause for a short period to control frame rate
    pause(0.05);  % Pause for 50 milliseconds between frames
end
