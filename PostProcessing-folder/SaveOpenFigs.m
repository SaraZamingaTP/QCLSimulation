function SaveOpenFigs()

    figHandles = findall(0, 'Type', 'figure');
    
    % Specify the directory where you want to save the figures
    saveDir = 'D:\PhD\MyPapers\NUSOD-Extended\images2';
    
    % Loop through each figure and save it
    for i = 1:length(figHandles)
        fig = figHandles(i);
        figName = sprintf('figure_%d.png', i); % You can use a different format or extension
        fullFilePath = fullfile(saveDir, figName);
        saveas(fig, fullFilePath);
    end

end