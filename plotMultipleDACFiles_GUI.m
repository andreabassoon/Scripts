function plotMultipleDACFiles_GUI
    % Create the GUI window with title including date/time and user info
    dateStr = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    userName = getenv('USERNAME');
    
    % Create main figure with normalized units
    hFig = figure('Name', sprintf('Analisi Spettrale - %s - User: %s', dateStr, userName), ...
                  'NumberTitle', 'off', ...
                  'Units', 'normalized', ...
                  'Position', [0.1, 0.1, 0.8, 0.8], ...
                  'ResizeFcn', @resizeGUI);
    
    % Create main panel to contain all controls
    mainPanel = uipanel('Parent', hFig, ...
                       'Units', 'normalized', ...
                       'Position', [0, 0, 1, 1], ...
                       'BorderType', 'none');
    
    % Create status bar at the bottom first
    statusBar = uicontrol(mainPanel, 'Style', 'text', ...
                         'Units', 'normalized', ...
                         'Position', [0.01, 0, 0.98, 0.03], ...
                         'String', 'Pronto', ...
                         'HorizontalAlignment', 'left');
    
    % Define updateStatus function early
    function updateStatus(message)
        if isvalid(statusBar)  % Check if statusBar still exists
            set(statusBar, 'String', [datestr(now, 'HH:MM:SS') ' - ' message]);
            drawnow;
        end
    end
    
    % Control panel for buttons at the top
    controlPanel = uipanel('Parent', mainPanel, ...
                          'Units', 'normalized', ...
                          'Position', [0, 0.9, 1, 0.1], ...
                          'BorderType', 'none');
    
    % Create buttons with normalized units
    uicontrol(controlPanel, 'Style', 'pushbutton', ...
              'String', 'Carica File', ...
              'Units', 'normalized', ...
              'Position', [0.01, 0.2, 0.08, 0.6], ...
              'BackgroundColor', [0.8 0.8 1.0], ...
              'Callback', @caricaFileCallback, ...
              'TooltipString', 'Carica i file DAC per l''analisi');
    
    % Analysis type selection menu
    analysisMenu = uicontrol(controlPanel, 'Style', 'popupmenu', ...
        'String', {'Seleziona Analisi', ...
                  'Analisi dello spettro a intervalli temporali', ...
                  'Analisi dell''intensità della PL nel tempo in un range di lunghezze d''onda', ...
                  'Visualizzazione contour plot delle matrici caricate'}, ...
        'Units', 'normalized', ...
        'Position', [0.1, 0.2, 0.25, 0.6], ...
        'Callback', @updateUI, ...
        'TooltipString', 'Seleziona il tipo di analisi da eseguire');
    
    % Execute analysis button
    uicontrol(controlPanel, 'Style', 'pushbutton', ...
              'String', 'Esegui Analisi', ...
              'Units', 'normalized', ...
              'Position', [0.36, 0.2, 0.08, 0.6], ...
              'BackgroundColor', [0.8 1.0 0.8], ...
              'Callback', @eseguiAnalisiCallback, ...
              'TooltipString', 'Esegui l''analisi selezionata');
    
    % Create plotting panel
    plotPanel = uipanel('Parent', mainPanel, ...
                       'Units', 'normalized', ...
                       'Position', [0, 0.05, 1, 0.85], ...
                       'BorderType', 'none');
    
    % Create axes with normalized units
    hAxes1 = axes('Parent', plotPanel, ...
                 'Units', 'normalized', ...
                 'Position', [0.05, 0.5, 0.9, 0.45], ...
                 'Visible', 'off');
    
    hAxes2 = axes('Parent', plotPanel, ...
                 'Units', 'normalized', ...
                 'Position', [0.05, 0.02, 0.9, 0.45], ...
                 'Visible', 'off');
    
    % Parameter controls
    uicontrol(controlPanel, 'Style', 'text', ...
              'Units', 'normalized', ...
              'Position', [0.45, 0.5, 0.08, 0.3], ...
              'String', 'Range WL:', ...
              'TooltipString', 'Intervallo di lunghezze d''onda');
    
    wlRangeBox = uicontrol(controlPanel, 'Style', 'edit', ...
                          'Units', 'normalized', ...
                          'Position', [0.54, 0.2, 0.08, 0.6], ...
                          'String', '[650 760]', ...
                          'Callback', @updatePlot, ...
                          'TooltipString', 'Inserisci l''intervallo di lunghezze d''onda [min max]');
    
    % Time range controls
    uicontrol(controlPanel, 'Style', 'text', ...
              'Units', 'normalized', ...
              'Position', [0.63, 0.5, 0.08, 0.3], ...
              'String', 'Range Tempo:', ...
              'TooltipString', 'Intervallo di tempo');
    
    timeRangeBox = uicontrol(controlPanel, 'Style', 'edit', ...
                            'Units', 'normalized', ...
                            'Position', [0.72, 0.2, 0.08, 0.6], ...
                            'String', '[0 1900]', ...
                            'Callback', @updatePlot, ...
                            'TooltipString', 'Inserisci l''intervallo di tempo [min max]');
    
    % Loading text
    loadingText = uicontrol(mainPanel, 'Style', 'text', ...
                           'Units', 'normalized', ...
                           'Position', [0.01, 0.85, 0.25, 0.03], ...
                           'String', '', ...
                           'Visible', 'off');
    
    % Global variables
    global S filenames filepath backgroundPoints;
    S = [];
    filenames = {};
    filepath = '';
    
    % Resize function
    function resizeGUI(~,~)
        % This function will be called when the figure is resized
        drawnow;
    end
    

    function caricaFileCallback(~, ~)
        try
            set(loadingText, 'String', 'Caricamento file...', 'Visible', 'on');
            drawnow;
            WLR = [400, 700];
            [Load, filenames, filepath, S] = BrowseFiles(WLR);
            if Load
                updateStatus('File caricati con successo!');
                msgbox('File caricati con successo!', 'Successo');
            else
                updateStatus('Nessun file caricato.');
                msgbox('Nessun file caricato.', 'Errore', 'error');
            end
        catch err
            updateStatus(['Errore: ' err.message]);
            errordlg(['Errore nel caricamento: ' err.message], 'Errore');
        end
        set(loadingText, 'Visible', 'off');
    end

    function updateUI(~, ~)
        analysisType = get(analysisMenu, 'Value');
        switch analysisType
            case 2
                set(wlRangeBox, 'Visible', 'off');
                set(timeRangeBox, 'Visible', 'on');
            case 3
                set(wlRangeBox, 'Visible', 'on');
                set(timeRangeBox, 'Visible', 'off');
            otherwise
                set(wlRangeBox, 'Visible', 'off');
                set(timeRangeBox, 'Visible', 'off');
        end
        updateStatus('Interfaccia aggiornata');
    end

    function updatePlot(~, ~)
        analysisType = get(analysisMenu, 'Value');
        if analysisType == 2 || analysisType == 3
            eseguiAnalisiCallback();
        end
    end

    function eseguiAnalisiCallback(~, ~)
        try
            analysisType = get(analysisMenu, 'Value');
            updateStatus('Esecuzione analisi...');
            switch analysisType
                case 2
                    interval = str2num(get(timeRangeBox, 'String'));
                    plotSpettroIntervalliTemporali(interval);
                case 3
                    wlRange = str2num(get(wlRangeBox, 'String'));
                    plotIntensitaTempo(wlRange);
                case 4
                    plotContour();
            end
            updateStatus('Analisi completata');
        catch err
            updateStatus(['Errore: ' err.message]);
            errordlg(['Errore nell''analisi: ' err.message], 'Errore');
        end
    end

function plotSpettroIntervalliTemporali(interval)
    if isempty(S)
        updateStatus('Nessun dato da analizzare');
        return;
    end
    
    % Store data in persistent variables so it's available for all sub-functions
    persistent spectraData currentInterval
    currentInterval = interval;
    spectraData = cell(1, numel(S));
    for i = 1:numel(S)
        timeIdx = (S{i}.TIME >= currentInterval(1) & S{i}.TIME <= currentInterval(2));
        spectraData{i}.x = S{i}.WL;
        spectraData{i}.y = mean(S{i}.Data(timeIdx, :), 1);
        spectraData{i}.name = filenames{i};
    end
    
    % Plot the initial data
    plotSpectralData(spectraData, currentInterval, false); % false = not normalized
    
    % Show analysis options menu
    showSpectralOptions(spectraData, currentInterval);
end

function plotSpectralData(spectraData, interval, isNormalized)
    axes(hAxes1);
    set(hAxes1, 'Visible', 'on');
    cla(hAxes1); % Clear previous plot
    hold(hAxes1, 'on');
    
    for i = 1:numel(spectraData)
        yData = spectraData{i}.y;
        if isNormalized
            yData = yData / max(yData);
        end
        plot(hAxes1, spectraData{i}.x, yData, 'DisplayName', spectraData{i}.name);
    end
    
    hold(hAxes1, 'off');
    xlabel(hAxes1, 'Lunghezza d''onda (nm)');
    if isNormalized
        ylabel(hAxes1, 'Intensità normalizzata (u.a.)');
    else
        ylabel(hAxes1, 'Intensità media (u.a.)');
    end
    title(hAxes1, sprintf('Spettro a intervalli temporali [%.1f, %.1f] s', interval(1), interval(2)));
    legend(hAxes1, 'show', 'Location', 'best');
    grid(hAxes1, 'on');
end

function showSpectralOptions(spectraData, interval)
    options = {'Normalizzare gli spettri', ...
              'Fit con gaussiane multiple', ...
              'Cambia intervallo temporale', ...
              'Salva spettri come txt', ...  % New option
              'Fine'};
    
    while true
        [option, ok] = listdlg('ListString', options, ...
                             'SelectionMode', 'single', ...
                             'Name', 'Opzioni di analisi spettrale', ...
                             'PromptString', 'Seleziona un''opzione:', ...
                             'ListSize', [300 150]);
        if ~ok
            break;
        end
        
        switch option
            case 1
                plotSpectralData(spectraData, interval, true);
            case 2
                fitGaussians(spectraData);
            case 3
                changeTimeInterval(interval);
            case 4
                saveSpectraToTxt(spectraData);  % New function
            case 5
                break;
        end
    end
end

function saveSpectraToTxt(spectraData)
    % Get save location from user
    [fileName, pathName] = uiputfile('*.txt', 'Salva spettri come txt');
    if fileName == 0
        return;
    end
    
    % Prepare data for saving
    % First column will be wavelength, subsequent columns will be intensities
    x = spectraData{1}.x(:);  % Assuming all spectra share same x-axis
    data = zeros(length(x), numel(spectraData) + 1);
    data(:,1) = x;
    
    % Create header
    header = 'Wavelength(nm)';
    
    % Add data from each spectrum
    for i = 1:numel(spectraData)
        data(:,i+1) = spectraData{i}.y(:);
        header = [header, sprintf('\t%s', spectraData{i}.name)];
    end
    
    % Save to file
    fullPath = fullfile(pathName, fileName);
    fid = fopen(fullPath, 'w');
    fprintf(fid, '%s\n', header);
    dlmwrite(fullPath, data, '-append', 'delimiter', '\t', 'precision', '%.6f');
    fclose(fid);
end

function fitGaussians(spectraData)
    hold(hAxes1, 'on');
    
    % Ask for number of peaks to fit
    answer = inputdlg('Numero di picchi da fittare:', 'Parametri di fit', [1 35]);
    if isempty(answer)
        return;
    end
    numPeaks = str2double(answer{1});
    
    % Get initial peak positions from first spectrum
    i = 1;  % Use first dataset for initial peak selection
    yData = double(spectraData{i}.y(:)');
    xData = double(spectraData{i}.x(:)');
    yData_norm = yData / max(yData);
    
    % Initialize peaks manually
    initialPeaks = zeros(1, numPeaks);
    initialAmps = zeros(1, numPeaks);
    
    % Single figure for all peak selections
    f = figure('Name', 'Selezione Picchi');
    for p = 1:numPeaks
        clf(f);
        plot(xData, yData_norm);
        xlabel('Lunghezza d''onda (nm)');
        ylabel('Intensità normalizzata (u.a.)');
        title(['Seleziona il picco ' num2str(p) ' di ' num2str(numPeaks) ' (stima iniziale)']);
        grid on;
        
        [peakX, peakY] = ginput(1);
        initialPeaks(p) = peakX;
        initialAmps(p) = peakY;
    end
    close(f);
    
    % Create initial parameters
    initialParams = [];
    for p = 1:numPeaks
        initialParams = [initialParams, initialAmps(p), initialPeaks(p), 30];
    end
    
    % Initialize storage for fit parameters and R-squared values
    fitParamsArray = cell(1, numel(spectraData));
    allRsquared = zeros(1, numel(spectraData));
    
    % Initialize arrays for peak data and uncertainties
    peakCenters = zeros(numel(spectraData), numPeaks);
    peakAmplitudes = zeros(numel(spectraData), numPeaks);
    centerErrors = zeros(numel(spectraData), numPeaks);
    amplitudeErrors = zeros(numel(spectraData), numPeaks);
    sampleNames = cell(numel(spectraData), 1);
    
    % Fit all spectra
    for i = 1:numel(spectraData)
        try
            % Get current spectrum data
            yData = double(spectraData{i}.y(:)');
            xData = double(spectraData{i}.x(:)');
            yData_norm = yData / max(yData);
            
            % Perform optimization
            options = optimset('Display', 'off', 'MaxFunEvals', 1000*length(initialParams));
            [finalParams, fval] = fminsearch(@(params) fitError(params, xData, yData_norm, numPeaks), initialParams, options);
            
            % Store parameters
            fitParamsArray{i} = finalParams;
            
            % Calculate fit and R-squared
            yFit = calculateMultiGaussian(finalParams, xData, numPeaks);
            SSres = sum((yData_norm - yFit).^2);
            SStot = sum((yData_norm - mean(yData_norm)).^2);
            allRsquared(i) = 1 - SSres/SStot;
            
            % Calculate uncertainties based on residuals
            residuals = yData_norm - yFit;
            chi_squared = sum(residuals.^2);
            dof = length(yData_norm) - length(finalParams); % degrees of freedom
            sigma = sqrt(chi_squared / dof);
            
            % Store peak parameters and estimated uncertainties
            sampleNames{i} = spectraData{i}.name;
            for p = 1:numPeaks
                paramIdx = (p-1)*3 + 1;
                peakAmplitudes(i, p) = finalParams(paramIdx);
                peakCenters(i, p) = finalParams(paramIdx + 1);
                
                % Estimate uncertainties using fit quality and parameter values
                amplitudeErrors(i, p) = sigma * abs(finalParams(paramIdx)) * 0.05; % 5% of amplitude
                centerErrors(i, p) = sigma * 0.5; % 0.5 nm uncertainty in position
            end
            
            % Plot results for current spectrum
            if i == 1
                axes(hAxes1);
                cla(hAxes1);
            end
            
            colorOrder = get(gca, 'ColorOrder');
            currentColor = colorOrder(mod(i-1, size(colorOrder,1))+1, :);
            
            plot(hAxes1, xData, yData_norm, '-', 'Color', currentColor, ...
                 'DisplayName', [spectraData{i}.name ' - Data']);
            plot(hAxes1, xData, yFit, '--', 'Color', currentColor, ...
                 'DisplayName', [spectraData{i}.name ' - Fit']);
            
            % Plot individual peaks
            for p = 1:numPeaks
                paramIdx = (p-1)*3 + 1;
                amplitude = finalParams(paramIdx);
                center = finalParams(paramIdx + 1);
                width = abs(finalParams(paramIdx + 2));
                
                yIndividual = amplitude * exp(-((xData-center)/width).^2);
                
                plot(hAxes1, xData, yIndividual, ':', 'Color', currentColor * 0.8, ...
                     'DisplayName', sprintf('%s - Peak %d (%.1f nm)', spectraData{i}.name, p, center));
            end
            
        catch err
            errordlg(['Errore nel fit per spettro ' num2str(i) ': ' err.message], 'Errore');
            disp(err.message);
        end
    end
    
    % Update plot aesthetics
    xlabel(hAxes1, 'Lunghezza d''onda (nm)');
    ylabel(hAxes1, 'Intensità normalizzata (u.a.)');
    ylim(hAxes1, [0 1.1]);
    grid(hAxes1, 'on');
    legend(hAxes1, 'show', 'Location', 'best');
    
    % Create scatter plot for each peak position with error bars
    for p = 1:numPeaks
        figure('Name', sprintf('Peak %d Positions', p), ...
               'NumberTitle', 'off', ...
               'Position', [200+30*p, 200+30*p, 800, 400]);
        
        x_positions = 1:numel(spectraData);
        currentPeakPositions = peakCenters(:, p);
        currentErrors = centerErrors(:, p);
        
        % Create error bar plot
        errorbar(x_positions, currentPeakPositions, currentErrors, 'o', ...
                'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
                'MarkerEdgeColor', [0 0.4470 0.7410], ...
                'LineWidth', 1.5, ...
                'MarkerSize', 8);
        
        grid on;
        set(gca, 'XTick', x_positions);
        set(gca, 'XTickLabel', sampleNames);
        
        if numel(sampleNames) > 5
            xtickangle(45);
        end
        
        xlabel('Nome Campione');
        ylabel('Posizione del picco (nm)');
        title(sprintf('Posizione del picco %d per ciascun campione\n%s - %s', ...
              p, '2025-03-06 16:00:55', 'andreabassoon'));
        
        xlim([0.5, numel(spectraData) + 0.5]);
        y_range = max(currentPeakPositions + currentErrors) - min(currentPeakPositions - currentErrors);
        if y_range == 0
            y_range = max(currentPeakPositions) * 0.1;
        end
        ylim([min(currentPeakPositions - currentErrors) - 0.1*y_range, ...
              max(currentPeakPositions + currentErrors) + 0.1*y_range]);
    end
    
    % Add plot for peak intensity ratio if exactly 2 peaks were fitted
    if numPeaks == 2
        figure('Name', 'Rapporto delle Intensità dei Picchi', ...
               'NumberTitle', 'off', ...
               'Position', [260, 260, 800, 400]);
        
        % Calculate ratio and its uncertainty using error propagation
        peakRatio = peakAmplitudes(:, 1) ./ peakAmplitudes(:, 2);
        % Error propagation formula for ratio A/B: δ(A/B) = |A/B| * sqrt((δA/A)^2 + (δB/B)^2)
        ratioErrors = abs(peakRatio) .* sqrt((amplitudeErrors(:, 1)./peakAmplitudes(:, 1)).^2 + ...
                                            (amplitudeErrors(:, 2)./peakAmplitudes(:, 2)).^2);
        
        % Create error bar plot for ratio
        errorbar(x_positions, peakRatio, ratioErrors, 'o', ...
                'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
                'MarkerEdgeColor', [0.6350 0.0780 0.1840], ...
                'LineWidth', 1.5, ...
                'MarkerSize', 8);
        
        grid on;
        set(gca, 'XTick', x_positions);
        set(gca, 'XTickLabel', sampleNames);
        
        if numel(sampleNames) > 5
            xtickangle(45);
        end
        
        xlabel('Nome Campione');
        ylabel('Rapporto delle Intensità (Picco 1/Picco 2)');
        title(sprintf('Rapporto delle Intensità dei Picchi\n%s - %s', ...
              '2025-03-06 16:00:55', 'andreabassoon'));
        
        xlim([0.5, numel(spectraData) + 0.5]);
        y_range = max(peakRatio + ratioErrors) - min(peakRatio - ratioErrors);
        if y_range == 0
            y_range = max(peakRatio) * 0.1;
        end
        ylim([min(peakRatio - ratioErrors) - 0.1*y_range, ...
              max(peakRatio + ratioErrors) + 0.1*y_range]);
    end
    
    % Prompt to save parameters
    saveParams = questdlg('Vuoi salvare i parametri del fit?', ...
                         'Salva Parametri', ...
                         'Si', 'No', 'Si');
    if strcmp(saveParams, 'Si')
        [fileName, pathName] = uiputfile('*.txt', 'Salva parametri di fit');
        if fileName ~= 0
            saveAllFitParameters(fullfile(pathName, fileName), spectraData, fitParamsArray, allRsquared, numPeaks);
        end
    end
    
    hold(hAxes1, 'off');
end

function saveAllFitParameters(fullPath, spectraData, allFitParams, allRsquared, numPeaks)
    % Open file for writing
    fid = fopen(fullPath, 'w');
    
    % Write timestamp and user info
    fprintf(fid, 'Data Analysis Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, 'User: %s\n\n', getenv('USERNAME'));
    
    % Write header
    fprintf(fid, 'Spectrum Name\tR-squared');
    for p = 1:numPeaks
        fprintf(fid, '\tPeak%d_Amplitude\tPeak%d_Center\tPeak%d_Width', p, p, p);
    end
    fprintf(fid, '\n');
    
    % Write parameters for each spectrum
    for i = 1:numel(spectraData)
        % Write spectrum name and R-squared
        fprintf(fid, '%s\t%.6f', spectraData{i}.name, allRsquared(i));
        
        % Write parameters for each peak
        params = allFitParams{i};
        for p = 1:numPeaks
            paramIdx = (p-1)*3 + 1;
            amplitude = params(paramIdx);
            center = params(paramIdx + 1);
            width = abs(params(paramIdx + 2));
            fprintf(fid, '\t%.6f\t%.6f\t%.6f', amplitude, center, width);
        end
        fprintf(fid, '\n');
    end
    
    fclose(fid);
    
    % Show confirmation message
    msgbox(sprintf('Parametri salvati in:\n%s', fullPath), 'Salvataggio Completato', 'help');
end

function err = fitError(params, x, y, numPeaks)
    % Calculate error between data and fit
    yFit = calculateMultiGaussian(params, x, numPeaks);
    err = sum((y - yFit).^2);  % Sum of squared errors
end

function y = calculateMultiGaussian(params, x, numPeaks)
    % Initialize output array
    y = zeros(size(x));
    
    % Add contribution from each Gaussian
    for p = 1:numPeaks
        paramIdx = (p-1)*3 + 1;
        amplitude = params(paramIdx);
        center = params(paramIdx + 1);
        width = abs(params(paramIdx + 2));  % Ensure width is positive
        
        % Add this Gaussian to the total
        y = y + amplitude * exp(-((x-center)/width).^2);
    end
end











    function plotIntensitaTempo(wlRange)
        if isempty(S)
            updateStatus('Nessun dato da analizzare');
            return;
        end
        
        intensities = cell(1, numel(S));
        axes(hAxes2);
        set(hAxes2, 'Visible', 'on');
        hold(hAxes2, 'on');
        for i = 1:numel(S)
            wlIdx = (S{i}.WL >= wlRange(1) & S{i}.WL <= wlRange(2));
            intensities{i} = mean(S{i}.Data(:, wlIdx), 2);
            plot(hAxes2, S{i}.TIME, intensities{i}, 'DisplayName', filenames{i});
        end
        hold(hAxes2, 'off');
        xlabel(hAxes2, 'Tempo (s)');
        ylabel(hAxes2, 'Intensità media (u.a.)');
        title(hAxes2, 'Intensità della PL nel tempo');
        legend(hAxes2, 'show', 'Location', 'best');
        grid(hAxes2, 'on');
        showMenuOptions(intensities);
    end

function showMenuOptions(intensities)
    options = {'Sottrarre il background', ...
              'Normalizzare tutte le curve', ...
              'Visualizzare in scala logaritmica', ...
              'Visualizzare in scala lineare', ...  % New option
              'Fit Triesponenziale', ...
              'Salva le curve finali in un file txt', ...
              'Fine'};
    while true
        [option, ok] = listdlg('ListString', options, ...
                             'SelectionMode', 'single', ...
                             'Name', 'Opzioni di analisi', ...
                             'PromptString', 'Seleziona un''opzione:', ...
                             'ListSize', [300 150]);
        if ~ok
            break;
        end
        
        switch option
            case 1
                x = inputdlg('Inserisci il numero di punti per calcolare il background:', ...
                           'Background', [1 50]);
                if ~isempty(x)
                    backgroundPoints = str2double(x{1});
                    sottraiBackground(intensities, backgroundPoints);
                end
            case 2
                normalizzaCurve(intensities);
            case 3
                visualizzaScalaLogaritmica(intensities);
            case 4
                visualizzaScalaLineare(intensities);  % New function call
            case 5
                % Exit the menu loop before starting the fit
                delete(findobj('Tag', 'OptionDialog'));
                fitTriesponenziale(intensities);
                return; % Exit the function after fitting
            case 6
                salvaCurveFinali(intensities);
            case 7
                break;
        end
    end
end

function visualizzaScalaLineare(intensities)
    % Get the current axes
    hAxes = gca;
    
    % Set y-axis to linear scale
    set(hAxes, 'YScale', 'linear');
    
    % Refresh the plot
    drawnow;
    
    % Update title if needed
    title(hAxes, 'Curve di decadimento - Scala lineare');
end

    function sottraiBackground(intensities, backgroundPoints)
        if isempty(backgroundPoints) || backgroundPoints <= 0
            updateStatus('Numero di punti non valido per il background');
            return;
        end
        
        axes(hAxes2);
        cla(hAxes2);
        hold(hAxes2, 'on');
        for i = 1:numel(intensities)
            background = mean(intensities{i}(1:backgroundPoints));
            intensities{i} = intensities{i} - background;
            plot(hAxes2, S{i}.TIME, intensities{i}, 'DisplayName', filenames{i});
        end
        hold(hAxes2, 'off');
        xlabel(hAxes2, 'Tempo (s)');
        ylabel(hAxes2, 'Intensità media (background sottratto)');
        title(hAxes2, 'Intensità della PL nel tempo (background sottratto)');
        legend(hAxes2, 'show', 'Location', 'best');
        grid(hAxes2, 'on');
        updateStatus('Background sottratto con successo');
    end

    function normalizzaCurve(intensities)
        axes(hAxes2);
        cla(hAxes2);
        hold(hAxes2, 'on');
        for i = 1:numel(intensities)
            intensities{i} = intensities{i} / max(intensities{i});
            plot(hAxes2, S{i}.TIME, intensities{i}, 'DisplayName', filenames{i});
        end
        hold(hAxes2, 'off');
        xlabel(hAxes2, 'Tempo (s)');
        ylabel(hAxes2, 'Intensità normalizzata');
        title(hAxes2, 'Intensità della PL nel tempo (normalizzata)');
        legend(hAxes2, 'show', 'Location', 'best');
        grid(hAxes2, 'on');
        updateStatus('Curve normalizzate');
    end
    function visualizzaScalaLogaritmica(intensities)
        axes(hAxes2);
        cla(hAxes2);
        hold(hAxes2, 'on');
        for i = 1:numel(intensities)
            plot(hAxes2, S{i}.TIME, intensities{i}, 'DisplayName', filenames{i});
        end
        set(hAxes2, 'YScale', 'log');
        hold(hAxes2, 'off');
        xlabel(hAxes2, 'Tempo (s)');
        ylabel(hAxes2, 'Intensità media (scala logaritmica)');
        title(hAxes2, 'Intensità della PL nel tempo (scala logaritmica)');
        legend(hAxes2, 'show', 'Location', 'best');
        grid(hAxes2, 'on');
        updateStatus('Visualizzazione in scala logaritmica');
    end

function fitTriesponenziale(intensities)
    % Create the dialog window with increased size for bounds
    d = figure('Name', 'Parametri di Fit', ...
              'NumberTitle', 'off', ...
              'Position', [300, 300, 700, 700], ... % Increased height for new parameters
              'MenuBar', 'none', ...
              'ToolBar', 'none');
    
    % Parameters panel with increased size
    uipanel('Title', 'Parametri', 'Position', [0.05, 0.25, 0.9, 0.7]);
    
    % Headers for columns
    columnHeaders = {'Parametro', 'Valore Iniziale', 'Limite Inferiore', 'Limite Superiore', 'Fissa'};
    columnPositions = [0.02, 0.2, 0.4, 0.6, 0.8];
    for i = 1:length(columnHeaders)
        uicontrol('Style', 'text', ...
                 'String', columnHeaders{i}, ...
                 'Position', [columnPositions(i)*700, 650, 120, 20], ... % Adjusted position
                 'HorizontalAlignment', 'center');
    end
    
    % Default values structure - added A3 and t3
    defaults = struct(...
        'A1', struct('val', '0.5', 'lb', '0', 'ub', '1'), ...
        't1', struct('val', '250', 'lb', '200', 'ub', '400'), ...
        'A2', struct('val', '0.3', 'lb', '0', 'ub', '1'), ...
        't2', struct('val', '2500', 'lb', '2000', 'ub', '3500'), ...
        'A3', struct('val', '0', 'lb', '0', 'ub', '0'), ...
        't3', struct('val', '0', 'lb', '2000', 'ub', '5000'), ...
        'y0', struct('val', '0.01', 'lb', '0.001', 'ub', '0.05'));
    
    % Create parameter fields with bounds - added A3 and t3
    paramNames = {'A1', 't1', 'A2', 't2', 'A3', 't3', 'y0'};
    handles = struct();
    for i = 1:length(paramNames)
        [handles.(paramNames{i}).val, ...
         handles.(paramNames{i}).lb, ...
         handles.(paramNames{i}).ub, ...
         handles.(paramNames{i}).fix] = createParamFieldWithBounds(...
            paramNames{i}, ...
            0.95 - 0.12*i, ... % Adjusted spacing
            defaults.(paramNames{i}).val, ...
            defaults.(paramNames{i}).lb, ...
            defaults.(paramNames{i}).ub);
    end
    
    % OK and Cancel buttons
    uicontrol('Style', 'pushbutton', ...
             'String', 'OK', ...
             'Position', [230, 20, 100, 30], ...
             'Callback', @onOK);
    
    uicontrol('Style', 'pushbutton', ...
             'String', 'Annulla', ...
             'Position', [370, 20, 100, 30], ...
             'Callback', @(~,~) close(d));
    
    % Nested function to create parameter fields with bounds
    function [h_val, h_lb, h_ub, h_fix] = createParamFieldWithBounds(label, y, defaultVal, defaultLB, defaultUB)
        % Create label
        uicontrol('Style', 'text', ...
                 'String', label, ...
                 'Position', [20, y*500, 80, 25], ...
                 'HorizontalAlignment', 'right');
        
        % Create input fields
        h_val = uicontrol('Style', 'edit', ...
                         'String', defaultVal, ...
                         'Position', [140, y*500, 120, 25], ...
                         'BackgroundColor', 'white');
        
        h_lb = uicontrol('Style', 'edit', ...
                        'String', defaultLB, ...
                        'Position', [280, y*500, 120, 25], ...
                        'BackgroundColor', 'white');
        
        h_ub = uicontrol('Style', 'edit', ...
                        'String', defaultUB, ...
                        'Position', [420, y*500, 120, 25], ...
                        'BackgroundColor', 'white');
        
        % Create checkbox
        h_fix = uicontrol('Style', 'checkbox', ...
                         'Position', [560, y*500, 20, 25]);
    end
    
    function [val, lb, ub] = getParamAndBounds(valHandle, lbHandle, ubHandle)
        val = str2double(get(valHandle, 'String'));
        lb = str2double(get(lbHandle, 'String'));
        ub = str2double(get(ubHandle, 'String'));
    end
    
    function validateParams(params, bounds)
        % Check for valid numbers - added A3 and t3
        fields = {'A1', 't1', 'A2', 't2', 'A3', 't3', 'y0'};
        for i = 1:length(fields)
            f = fields{i};
            if isnan(params.(f)) || isnan(bounds.([f '_lb'])) || isnan(bounds.([f '_ub']))
                error(['Il parametro ' f ' o i suoi limiti non sono numeri validi']);
            end
            
            % Check bounds consistency
            if bounds.([f '_lb']) >= bounds.([f '_ub'])
                error(['Il limite inferiore di ' f ' deve essere minore del limite superiore']);
            end
            if params.(f) < bounds.([f '_lb']) || params.(f) > bounds.([f '_ub'])
                error(['Il valore iniziale di ' f ' deve essere compreso tra i limiti']);
            end
        end
    end
    
    function onOK(~, ~)
        try
            % Get parameter values, bounds, and fixed states
            params = struct();
            bounds = struct();
            fixedParams = zeros(1,7); % Changed to 7 parameters
            
            % Get values for each parameter
            for i = 1:length(paramNames)
                pName = paramNames{i};
                [params.(pName), bounds.([pName '_lb']), bounds.([pName '_ub'])] = ...
                    getParamAndBounds(handles.(pName).val, handles.(pName).lb, handles.(pName).ub);
                fixedParams(i) = get(handles.(pName).fix, 'Value');
            end
            
            % Validate parameters
            validateParams(params, bounds);
            
            % Close dialog and perform fit
            delete(d);
            fitResults = performFit(intensities, params, bounds, fixedParams);
            saveFitParameters(fitResults);
            plotTauAvg(fitResults);
            updateStatus('Fit triesponenziale completato con successo');
            
        catch err
            errordlg(['Errore nel fit: ' err.message], 'Errore');
            updateStatus(['Errore nel fit: ' err.message]);
        end
    end
end
% Here's Part 3 - The updated performFit function with error calculation:

function fitResults = performFit(intensities, params, bounds, fixedParams)
    fitResults = struct('filenames', {filenames}, ...
                      'fitParams', {cell(1,numel(intensities))}, ...
                      'tau_avg_values', zeros(1,numel(intensities)), ...
                      'tau_avg_errors', zeros(1,numel(intensities)));
    
    axes(hAxes2);
    cla(hAxes2);
    hold(hAxes2, 'on');
    
    for i = 1:numel(intensities)
        % Select only data with positive time
        timeIdx = find(S{i}.TIME > 0, 1):length(S{i}.TIME);
        xData = S{i}.TIME(timeIdx);
        yData = intensities{i}(timeIdx);
        
        % Normalize data
        yData = yData / max(yData);
        
        % Initial parameters and bounds - added A3 and t3
        startPoint = [params.A1, params.t1, params.A2, params.t2, params.A3, params.t3, params.y0];
        lb = [bounds.A1_lb, bounds.t1_lb, bounds.A2_lb, bounds.t2_lb, bounds.A3_lb, bounds.t3_lb, bounds.y0_lb];
        ub = [bounds.A1_ub, bounds.t1_ub, bounds.A2_ub, bounds.t2_ub, bounds.A3_ub, bounds.t3_ub, bounds.y0_ub];
        
        % Apply fixed parameters
        for j = 1:7 % Changed to 7 parameters
            if fixedParams(j)
                lb(j) = startPoint(j);
                ub(j) = startPoint(j);
            end
        end
        
        % Fit function - added third exponential term
        fitFunc = @(p,x) p(1)*exp(-x/p(2)) + p(3)*exp(-x/p(4)) + p(5)*exp(-x/p(6)) + p(7);
        
        % Perform fit
        options = optimoptions('lsqcurvefit', 'Display', 'off', ...
                             'Algorithm', 'levenberg-marquardt', ...
                             'SpecifyObjectiveGradient', false);
        [fitParams, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(fitFunc, startPoint, xData, yData, lb, ub, options);
        
        % Store results
        fitResults.fitParams{i} = fitParams;
        
        % Calculate tau_avg using standard weighted average - added third term
        A1 = fitParams(1); t1 = fitParams(2);
        A2 = fitParams(3); t2 = fitParams(4);
        A3 = fitParams(5); t3 = fitParams(6);
        tau_avg = (A1*t1 + A2*t2 + A3*t3)/(A1 + A2 + A3);
        fitResults.tau_avg_values(i) = tau_avg;
        
        % Calculate uncertainties with improved robustness
        yFit = fitFunc(fitParams, xData);
        residuals = yData - yFit;
        chi_squared = sum(residuals.^2);
        dof = length(yData) - length(fitParams);
        sigma = sqrt(chi_squared / dof);
        
        % Calculate parameter uncertainties with condition number check
        try
            JtJ = jacobian'*jacobian;
            cond_num = cond(JtJ);
            
            if cond_num > 1e15
                lambda = 1e-6 * trace(JtJ)/size(JtJ,1);
                covariance = sigma^2 * inv(JtJ + lambda*eye(size(JtJ)));
            else
                covariance = sigma^2 * inv(JtJ);
            end
            
            param_errors = sqrt(abs(diag(covariance)));
            
            max_relative_error = 0.5;
            for p = 1:length(param_errors)
                if param_errors(p) > abs(fitParams(p)) * max_relative_error
                    param_errors(p) = abs(fitParams(p)) * max_relative_error;
                end
            end
            
        catch
            param_errors = abs(fitParams) * 0.1;
        end
        
        % Error propagation for tau_avg with improved stability - added third term
        denominator = (A1 + A2 + A3);
        
        if abs(denominator) < 1e-10
            tau_error = max(fitResults.tau_avg_values) * 0.1;
        else
            % Calculate partial derivatives for three-component weighted average
            dTau_dA1 = (t1*denominator - (A1*t1 + A2*t2 + A3*t3))/(denominator^2);
            dTau_dt1 = A1/denominator;
            dTau_dA2 = (t2*denominator - (A1*t1 + A2*t2 + A3*t3))/(denominator^2);
            dTau_dt2 = A2/denominator;
            dTau_dA3 = (t3*denominator - (A1*t1 + A2*t2 + A3*t3))/(denominator^2);
            dTau_dt3 = A3/denominator;
            
            % Calculate tau_avg error including all components
            tau_error = sqrt((dTau_dA1*param_errors(1))^2 + ...
                           (dTau_dt1*param_errors(2))^2 + ...
                           (dTau_dA2*param_errors(3))^2 + ...
                           (dTau_dt2*param_errors(4))^2 + ...
                           (dTau_dA3*param_errors(5))^2 + ...
                           (dTau_dt3*param_errors(6))^2);
            
            max_tau_error = tau_avg * 0.5;
            if tau_error > max_tau_error
                tau_error = max_tau_error;
            end
        end
        
        fitResults.tau_avg_errors(i) = tau_error;
        
        % Plot results
        plot(hAxes2, xData, yData, 'DisplayName', [filenames{i} ' (dati)']);
        plot(hAxes2, xData, yFit, '--', ...
             'DisplayName', [filenames{i} ' (fit)']);
    end
    
    xlabel(hAxes2, 'Tempo (s)');
    ylabel(hAxes2, 'Intensità normalizzata');
    title(hAxes2, sprintf('Fit triesponenziale\\n%s - %s', ...
          '2025-03-10 13:21:54', 'andreabassoon'));
    legend(hAxes2, 'show', 'Location', 'best');
    grid(hAxes2, 'on');
end
function plotTauAvg(fitResults)
    figure('Name', 'Tau Average Values', ...
           'NumberTitle', 'off', ...
           'Position', [200, 200, 800, 400]);
    
    % Create x-axis numeric positions for the scatter plot
    x_positions = 1:numel(fitResults.filenames);
    
    % Create error bar plot
    errorbar(x_positions, fitResults.tau_avg_values, fitResults.tau_avg_errors, 'o', ...
            'MarkerFaceColor', [0.3010 0.7450 0.9330], ...    % Light blue fill
            'MarkerEdgeColor', [0 0.4470 0.7410], ...         % Darker blue edge
            'LineWidth', 1.5, ...                             % Error bar line width
            'MarkerSize', 8);                                 % Marker size
    
    % Add grid
    grid on;
    
    % Set x-axis ticks and labels
    set(gca, 'XTick', x_positions);
    set(gca, 'XTickLabel', fitResults.filenames);
    
    % Rotate labels if needed
    if numel(fitResults.filenames) > 5
        xtickangle(45);
    end
    
    % Add labels and title
    xlabel('Nome Campione');
    ylabel('tau_{avg} (s)');
    title(sprintf('Valori di tau_{avg} per ciascun campione\n%s - %s', ...
          '2025-03-07 13:14:54', 'andreabassoon'));
    
    % Adjust axes limits to add some padding
    xlim([0.5, numel(fitResults.filenames) + 0.5]);
    y_range = max(fitResults.tau_avg_values + fitResults.tau_avg_errors) - ...
              min(fitResults.tau_avg_values - fitResults.tau_avg_errors);
    if y_range == 0
        y_range = max(fitResults.tau_avg_values) * 0.1;
    end
    ylim([min(fitResults.tau_avg_values - fitResults.tau_avg_errors) - 0.1*y_range, ...
          max(fitResults.tau_avg_values + fitResults.tau_avg_errors) + 0.1*y_range]);
end

function saveFitParameters(fitResults)
    [saveFile, savePath] = uiputfile('*.txt', 'Salva i parametri di fit');
    if ischar(saveFile)
        try
            fullPath = fullfile(savePath, saveFile);
            fid = fopen(fullPath, 'w');
            fprintf(fid, 'Data analisi: %s\n', '2025-03-10 13:21:54');
            fprintf(fid, 'Utente: %s\n\n', 'andreabassoon');
            fprintf(fid, 'Nome Campione\tA1\tt1\tA2\tt2\tA3\tt3\ty0\ttau_avg\ttau_avg_error\n');
            
            for i = 1:numel(fitResults.filenames)
                params = fitResults.fitParams{i};
                fprintf(fid, '%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
                    fitResults.filenames{i}, params(1), params(2), params(3), ...
                    params(4), params(5), params(6), params(7), ...
                    fitResults.tau_avg_values(i), fitResults.tau_avg_errors(i));
            end
            
            fclose(fid);
            updateStatus(['Parametri di fit salvati in ' fullPath]);
            msgbox('Parametri di fit salvati con successo!', 'Salvataggio completato');
        catch err
            errordlg(['Errore nel salvataggio: ' err.message], 'Errore');
            updateStatus(['Errore nel salvataggio: ' err.message]);
        end
    end
end

function salvaCurveFinali(intensities)
    [saveFile, savePath] = uiputfile('*.txt', 'Salva le curve finali');
    if ischar(saveFile)
        try
            fullPath = fullfile(savePath, saveFile);
            fid = fopen(fullPath, 'w');
            
            % Write header with metadata
            fprintf(fid, 'Data analisi: %s\n', '2025-03-07 13:14:54');
            fprintf(fid, 'Utente: %s\n\n', 'andreabassoon');
            
            % Get current plot data
            h = findobj(hAxes2, 'Type', 'line');
            
            % Write header for columns
            fprintf(fid, 'Tempo');
            for i = numel(h):-1:1  % Reverse order to match the original data order
                fprintf(fid, '\t%s', get(h(i), 'DisplayName'));
            end
            fprintf(fid, '\n');
            
            % Get the length of time vectors
            maxLength = length(S{1}.TIME);
            timeVector = S{1}.TIME;
            
            % Write data row by row
            for j = 1:maxLength
                % Write time value
                fprintf(fid, '%f', timeVector(j));
                
                % Write processed intensity values from the plot
                for i = numel(h):-1:1  % Reverse order to match the original data order
                    yData = get(h(i), 'YData');
                    if j <= length(yData)
                        fprintf(fid, '\t%f', yData(j));
                    else
                        fprintf(fid, '\t');
                    end
                end
                fprintf(fid, '\n');
            end
            
            fclose(fid);
            updateStatus(['Curve salvate in ' fullPath]);
            msgbox(['Curve salvate in ' fullPath], 'Salvataggio completato');
        catch err
            errordlg(['Errore nel salvataggio: ' err.message], 'Errore');
            updateStatus(['Errore nel salvataggio: ' err.message]);
        end
    end
end
    function plotContour()
        if isempty(S)
            updateStatus('Nessun dato da visualizzare');
            return;
        end
        
        numFiles = numel(S);
        numCols = min(3, numFiles);
        numRows = ceil(numFiles/numCols);
        
        figure('Name', 'Contour Plot', ...
               'NumberTitle', 'off', ...
               'Position', [100, 100, 300*numCols, 300*numRows]);
        
        for i = 1:numFiles
            subplot(numRows, numCols, i);
            [T, W] = meshgrid(S{i}.TIME, S{i}.WL);
            contourf(T, W, S{i}.Data', 20, 'LineColor', 'none');
            colormap(jet);
            colorbar;
            xlabel('Tempo (s)');
            ylabel('Lunghezza d''onda (nm)');
            title(filenames{i});
        end
        
        sgtitle(sprintf('Visualizzazione contour plot delle matrici caricate\n%s - %s', ...
               '2025-03-03 00:45:51', 'andreabassoon'));
        updateStatus('Contour plot completato');
    end
end

function [Load, filenames, filepath, S] = BrowseFiles(WLR)
    [filenames, filepath, indx] = uigetfile({'*.dac'}, 'File Selector', 'MultiSelect', 'on');
    if isequal(filenames, 0) % no file loaded
        Load = 0;
        S = [];
    else % files loaded
        Load = 1;
        if ischar(filenames) % single file selected
            filenames = {filenames};
        end
        S = cell(1, numel(filenames));
        for i = 1:numel(filenames)
            try
                s = importdata(fullfile(filepath, filenames{i}));
                % different file formats are possible. This section is able to
                % distinguish among the different types.
                if size(fieldnames(s), 1) == 2
                    if size(s.data, 1) == 1
                        [token, remain] = strtok(s.textdata{1, end});
                        s.textdata{1, end} = token;
                        S{i}.WL = str2double(strrep(s.textdata(1, 2:end), ',', '.'))';
                        time0 = str2double(strrep([remain; s.textdata(2:end, 1)], ',', '.'));
                        S{i}.Data = [s.data; str2double(strrep(s.textdata(2:end, 2:end), ',', '.'))];
                    elseif size(s.textdata, 1) == 2
                        S{i}.Data = s.data;
                        s.textdata = strrep(s.textdata, ',', '.');
                        time0 = str2double(s.textdata(2:end, 1));
                        S{i}.WL = str2double(s.textdata(1, 2:end))';
                    elseif size(s.textdata, 1) == 1 + size(s.data, 1) && size(s.textdata, 2) == 1 + size(s.data, 2)
                        S{i}.Data = s.data;
                        s.textdata = strrep(s.textdata, ',', '.');
                        time0 = str2double(s.textdata(2:end, 1));
                        S{i}.WL = str2double(s.textdata(1, 2:end))';
                    end
                elseif size(fieldnames(s), 1) == 3
                    if size(s.data, 1) > 10
                        S{i}.Data = s.data(2:end, :);
                        S{i}.WL = s.data(1, :)';
                        time0 = str2double(s.textdata(2:end, 1));
                    end
                end
                % find the zero
                [~, pos] = max(sum(S{i}.Data(:, (S{i}.WL > WLR(1) & S{i}.WL < WLR(2))), 2));
                S{i}.TIME = time0 - time0(pos);
            catch err
                warndlg(['Errore nel caricamento del file: ' filenames{i}], 'Errore');
                disp(['Error loading file ' filenames{i} ': ' err.message]);
            end
        end
    end
end
    function plotIntensitaTempo(wlRange)
        if isempty(S)
            updateStatus('Nessun dato da analizzare');
            return;
        end
        
        intensities = cell(1, numel(S));
        axes(hAxes2);
        set(hAxes2, 'Visible', 'on');
        hold(hAxes2, 'on');
        for i = 1:numel(S)
            wlIdx = (S{i}.WL >= wlRange(1) & S{i}.WL <= wlRange(2));
            intensities{i} = mean(S{i}.Data(:, wlIdx), 2);
            plot(hAxes2, S{i}.TIME, intensities{i}, 'DisplayName', filenames{i});
        end
        hold(hAxes2, 'off');
        xlabel(hAxes2, 'Tempo (s)');
        ylabel(hAxes2, 'Intensità media (u.a.)');
        title(hAxes2, 'Intensità della PL nel tempo');
        legend(hAxes2, 'show', 'Location', 'best');
        grid(hAxes2, 'on');
        showMenuOptions(intensities);
    end

    function showMenuOptions(intensities)
        persistent isProcessing
        if isempty(isProcessing)
            isProcessing = false;
        end
        
        if isProcessing
            return;
        end
        
        options = {'Sottrarre il background', ...
                  'Normalizzare tutte le curve', ...
                  'Visualizzare in scala logaritmica', ...
                  'Fit biesponenziale', ...
                  'Salva le curve finali in un file txt', ...
                  'Fine'};
        
        while true
            if isProcessing
                break;
            end
            
            [option, ok] = listdlg('ListString', options, ...
                                 'SelectionMode', 'single', ...
                                 'Name', 'Opzioni di analisi', ...
                                 'PromptString', 'Seleziona un''opzione:', ...
                                 'ListSize', [300 150]);
            if ~ok
                break;
            end
            
            switch option
                case 1
                    x = inputdlg('Inserisci il numero di punti per calcolare il background:', ...
                               'Background', [1 50]);
                    if ~isempty(x)
                        backgroundPoints = str2double(x{1});
                        sottraiBackground(intensities, backgroundPoints);
                    end
                case 2
                    normalizzaCurve(intensities);
                case 3
                    visualizzaScalaLogaritmica(intensities);
                case 4
                    isProcessing = true;
                    % Close the options dialog
                    delete(findobj('Tag', 'OptionDialog'));
                    % Start the fitting process
                    fitBiesponenziale(intensities);
                    isProcessing = false;
                    return;
                case 5
                    salvaCurveFinali(intensities);
                case 6
                    break;
            end
        end
    end

    function fitBiesponenziale(intensities)
        % Create the dialog window
        d = figure('Name', 'Parametri di Fit', ...
                  'NumberTitle', 'off', ...
                  'Position', [300, 300, 400, 500], ...
                  'MenuBar', 'none', ...
                  'ToolBar', 'none');
        
        % Parameters panel
        uipanel('Title', 'Parametri', 'Position', [0.05, 0.3, 0.9, 0.65]);
        
        % Create parameter input fields with labels
        edtA1 = createParamField('A1:', 0.1, 0.85, '1');
        edtT1 = createParamField('t1:', 0.1, 0.75, '1');
        edtA2 = createParamField('A2:', 0.1, 0.65, '1');
        edtT2 = createParamField('t2:', 0.1, 0.55, '1');
        edtY0 = createParamField('y0:', 0.1, 0.45, '0');
        
        % Create checkboxes for fixing parameters
        chkA1 = createCheckbox('Fissa A1', 0.1, 0.35);
        chkT1 = createCheckbox('Fissa t1', 0.1, 0.30);
        chkA2 = createCheckbox('Fissa A2', 0.1, 0.25);
        chkT2 = createCheckbox('Fissa t2', 0.1, 0.20);
        chkY0 = createCheckbox('Fissa y0', 0.1, 0.15);
        
        % OK and Cancel buttons
        uicontrol('Style', 'pushbutton', ...
                 'String', 'OK', ...
                 'Position', [80, 20, 100, 30], ...
                 'Callback', @onOK);
        
        uicontrol('Style', 'pushbutton', ...
                 'String', 'Annulla', ...
                 'Position', [220, 20, 100, 30], ...
                 'Callback', @(~,~) delete(d));
        
        % Wait for the dialog to complete
        uiwait(d);
        
        function h = createParamField(label, x, y, defaultVal)
            uicontrol('Style', 'text', ...
                     'String', label, ...
                     'Position', [x*400, y*500, 100, 25]);
            h = uicontrol('Style', 'edit', ...
                         'String', defaultVal, ...
                         'Position', [(x+0.25)*400, y*500, 100, 25]);
        end
        
        function h = createCheckbox(label, x, y)
            h = uicontrol('Style', 'checkbox', ...
                         'String', label, ...
                         'Position', [x*400, y*500, 150, 25]);
        end
        
        function onOK(~, ~)
            try
                % Get parameter values directly from the handles
                params = struct();
                params.A1 = str2double(get(edtA1, 'String'));
                params.t1 = str2double(get(edtT1, 'String'));
                params.A2 = str2double(get(edtA2, 'String'));
                params.t2 = str2double(get(edtT2, 'String'));
                params.y0 = str2double(get(edtY0, 'String'));
                
                % Get checkbox states
                fixedParams = [
                    get(chkA1, 'Value'),
                    get(chkT1, 'Value'),
                    get(chkA2, 'Value'),
                    get(chkT2, 'Value'),
                    get(chkY0, 'Value')
                ];
                
                % Validate parameters
                if any(isnan([params.A1, params.t1, params.A2, params.t2, params.y0]))
                    error('Tutti i parametri devono essere numeri validi');
                end
                
                % Close dialog and perform fit
                delete(d);
                fitResults = performFit(intensities, params, fixedParams);
                saveFitParameters(fitResults);
                plotTauAvg(fitResults);
                updateStatus('Fit biesponenziale completato con successo');
                
                % Show the analysis options menu again after fitting is complete
                showMenuOptions(intensities);
                
            catch err
                errordlg(['Errore nel fit: ' err.message], 'Errore');
                updateStatus(['Errore nel fit: ' err.message]);
            end
        end
    end