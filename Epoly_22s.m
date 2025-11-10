function Epoly_22s()
% EPOLY_22S - Enhanced reflectivity calculator with TMM and Berreman optimization
%
% This version supports both TMM (isotropic) and Berreman 4x4 (anisotropic) methods
% for calculating reflectivity and optimizing optical parameters.
%
% Features:
%   • TMM Method - Standard isotropic Transfer Matrix Method
%   • Berreman Method - FULL 4x4 matrix method with proper eigenvalue decomposition
%     (Translated from berreman.py - includes Layer, HalfSpace, Structure classes)
%   • Particle Swarm (All) - Optimize all parameters including thickness
%   • Particle Swarm (Fixed) - Only optimize n, k, and angle shift
%   • Global Search - Multiple lsqnonlin starting points
%   • Advanced (lsqnonlin) - Fast single optimization
%   • Separate optimization for TE and TM modes (different angle shifts allowed)
%   • Original experimental data files (eploy_22s_Spol.xlsx, eploy_22s_Ppol.xlsx)
%
% New Berreman Implementation Files (required):
%   - BerremanLayer.m - Layer class with eigenvalue decomposition
%   - BerremanHalfSpace.m - Semi-infinite isotropic media
%   - BerremanStructure.m - Scattering matrix assembly
%   - BerremanStackModel.m - High-level interface
%   - berreman_4x4_full.m - Wrapper function matching Python behavior
%
% -----------------------------------------------------------------------------

    %----------------------------
    % App State
    %----------------------------
    app = struct();
    app.Defaults.ito = 15;       % nm
    app.Defaults.pedot = 50.15;   % nm
    app.Defaults.wavelength = 561.0; % nm
    app.Defaults.voltage = 0.0;      % mV
    app.Defaults.n_real = 1.41;
    app.Defaults.k_manual = 0.05;
    app.Defaults.delta_k = -0.01;
    app.Defaults.n_real_te = 1.41;   % TE (in-plane, ordinary)
    app.Defaults.n_real_tm = 1.269;  % TM (out-of-plane, extraordinary)
    app.Defaults.k_te = 0.05;
    app.Defaults.k_tm = 0.05;
    app.Defaults.shift_te = 0.0;
    app.Defaults.shift_tm = 0.0;
    app.Defaults.calculation_method = 'TMM'; % 'TMM' or 'Berreman'
    app.OptimizedShift = 0.0; % global (single-mode) optimized shift
    app.LastOptimizationResults = struct(); % Store last optimization results for export
    app.AllOptimizationResults = struct(); % Store all optimization results for detailed export

    % Default Excel filename as requested
    app.DefaultExcel = '300rpm extinction coefficient.xlsx';
    app.PEDOT = load_pedot_data_improved(app.DefaultExcel); % improved loading

    %----------------------------
    % Main UI - Concise Interface
    %----------------------------
    % Create main window - small and simple
    fig = uifigure('Name','Epoly_22s - Optimization Tool',...
        'Position',[100 100 400 200]);
    
    % Optimization buttons panel
    pnlMain = uipanel(fig, 'Title', 'Optimization Controls', ...
        'Units', 'pixels', 'Position', [10, 10, 380, 180]);
    
    % TE Mode optimization button
    btnOptTE = uicontrol('Parent', pnlMain, 'Style', 'pushbutton', ...
        'String', 'Optimize TE Mode', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', 'Position', [20, 130, 160, 40], ...
        'Callback', @(s,e) on_optimize_te());
    
    % TM Mode optimization button
    btnOptTM = uicontrol('Parent', pnlMain, 'Style', 'pushbutton', ...
        'String', 'Optimize TM Mode', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', 'Position', [200, 130, 160, 40], ...
        'Callback', @(s,e) on_optimize_tm());
    
    % Settings button - opens pop-up settings window
    btnSettings = uicontrol('Parent', pnlMain, 'Style', 'pushbutton', ...
        'String', 'Settings...', ...
        'Units', 'pixels', 'Position', [20, 80, 340, 35], ...
        'Callback', @(s,e) open_settings_window());
    
    % Export button
    btnExport = uicontrol('Parent', pnlMain, 'Style', 'pushbutton', ...
        'String', 'Export Results', ...
        'FontSize', 11, ...
        'Units', 'pixels', 'Position', [20, 35, 340, 35], ...
        'Callback', @(s,e) on_export_results());
    
    % Settings window (created when needed, stored in app structure)
    settingsFig = [];
    
    % Initialize variables for settings window (will be created on first open)
    edITO = []; edPEDOT = []; edWavelength = []; edVoltage = [];
    btnLoadExcel = []; lblVoltInfo = []; lblAvailV = [];
    rbTMM = []; rbBerreman = []; lblMethodInfo = [];
    chkManual = []; edNreal = []; edK = []; lblRIDisplay = [];
    edNte = []; edNtm = []; edKte = []; edKtm = [];
    edShiftTE = []; edShiftTM = [];
    rbTE = []; rbTM = [];
    
    % Initialize display
    refresh_voltage_info();
    update_ri_display();

    %====================
    % Callbacks / helpers (same structure as Python)
    %====================
    function open_settings_window()
        % Create settings pop-up window with all detailed controls
        if ~isempty(settingsFig) && isvalid(settingsFig)
            % If window already exists, bring it to front
            figure(settingsFig);
            return;
        end
        
        % Create new settings window
        settingsFig = uifigure('Name', 'Settings - Epoly_22s', ...
            'Position', [150 50 540 820]);
        
        % Create scrollable panel
        scrollContainer = uipanel(settingsFig, 'Title', 'Settings', ...
            'Position', [0 0 540 820]);
        scrollPanel = uipanel(scrollContainer);
        scrollPanel.Position = [0 0 520 1600];
        scrollPanel.Scrollable = 'on';
        
        y_pos = 1500;
        panel_width = 500;
        left_margin = 10;
        
        %--- Layer thickness panel
        pnlThick = uipanel('Parent', scrollPanel, ...
            'Title', 'Layer Thicknesses (nm)', ...
            'Units', 'pixels', ...
            'Position', [left_margin, y_pos-110, panel_width, 100]);
        
        uicontrol('Parent', pnlThick, 'Style', 'text', ...
            'String', 'ITO:', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 60, 80, 20]);
        edITO = uicontrol('Parent', pnlThick, 'Style', 'edit', ...
            'String', num2str(app.Defaults.ito), ...
            'Units', 'pixels', 'Position', [100, 60, 80, 25]);
        
        uicontrol('Parent', pnlThick, 'Style', 'text', ...
            'String', 'PEDOT:PSS:', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 30, 80, 20]);
        edPEDOT = uicontrol('Parent', pnlThick, 'Style', 'edit', ...
            'String', num2str(app.Defaults.pedot), ...
            'Units', 'pixels', 'Position', [100, 30, 80, 25]);
        
        uicontrol('Parent', pnlThick, 'Style', 'text', ...
            'String', 'Glass | ITO | PEDOT:PSS | Water', 'FontAngle', 'italic', ...
            'Units', 'pixels', 'Position', [200, 60, 250, 20]);
        
        uicontrol('Parent', pnlThick, 'Style', 'pushbutton', ...
            'String', 'Reset to Defaults', ...
            'Units', 'pixels', 'Position', [300, 25, 150, 30], ...
            'Callback', @(src,evt) reset_thickness());
        
        y_pos = y_pos - 130;
        
        %--- Wavelength settings
        pnlWl = uipanel('Parent', scrollPanel, ...
            'Title', 'Wavelength Settings', ...
            'Units', 'pixels', ...
            'Position', [left_margin, y_pos-80, panel_width, 70]);
        
        uicontrol('Parent', pnlWl, 'Style', 'text', ...
            'String', 'Wavelength (nm):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 35, 120, 20]);
        edWavelength = uicontrol('Parent', pnlWl, 'Style', 'edit', ...
            'String', num2str(app.Defaults.wavelength), ...
            'Units', 'pixels', 'Position', [140, 35, 80, 25], ...
            'Callback', @(s,e) update_ri_display());
        
        btnLoadExcel = uicontrol('Parent', pnlWl, 'Style', 'pushbutton', ...
            'String', 'Load PEDOT Excel...', ...
            'Units', 'pixels', 'Position', [240, 35, 200, 25], ...
            'Callback', @(src,evt) on_load_excel());
        
        lblVoltInfo = uicontrol('Parent', pnlWl, 'Style', 'text', ...
            'String', '', 'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [10, 5, 450, 20]);
        
        y_pos = y_pos - 100;
        
        %--- Voltage settings
        pnlV = uipanel('Parent', scrollPanel, ...
            'Title', 'Voltage Settings', ...
            'Units', 'pixels', ...
            'Position', [left_margin, y_pos-60, panel_width, 50]);
        
        uicontrol('Parent', pnlV, 'Style', 'text', ...
            'String', 'Voltage (mV):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 15, 120, 20]);
        edVoltage = uicontrol('Parent', pnlV, 'Style', 'edit', ...
            'String', num2str(app.Defaults.voltage), ...
            'Units', 'pixels', 'Position', [140, 15, 80, 25], ...
            'Callback', @(s,e) update_ri_display());
        
        lblAvailV = uicontrol('Parent', pnlV, 'Style', 'text', ...
            'String', '', 'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [240, 15, 220, 20]);
        
        y_pos = y_pos - 80;
        
        %--- Calculation Method Selection
        pnlMethod = uipanel('Parent', scrollPanel, ...
            'Title', 'Calculation Method', ...
            'Units', 'pixels', ...
            'Position', [left_margin, y_pos-90, panel_width, 80]);
        
        bgMethod = uibuttongroup('Parent', pnlMethod, ...
            'Units', 'pixels', ...
            'Position', [10, 35, 480, 30], ...
            'BorderType', 'none', ...
            'SelectionChangedFcn', @(s,e) on_method_change());
        
        rbTMM = uicontrol('Parent', bgMethod, 'Style', 'radiobutton', ...
            'String', 'TMM (Isotropic)', ...
            'Units', 'pixels', 'Position', [10, 5, 150, 20], ...
            'Value', 1);
        rbBerreman = uicontrol('Parent', bgMethod, 'Style', 'radiobutton', ...
            'String', 'Berreman 4x4 (Anisotropic)', ...
            'Units', 'pixels', 'Position', [200, 5, 200, 20]);
        
        lblMethodInfo = uicontrol('Parent', pnlMethod, 'Style', 'text', ...
            'String', 'TMM: Single n,k for all directions | Berreman: Separate n_te, n_tm (anisotropic)', ...
            'FontAngle', 'italic', 'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [10, 5, 480, 25]);
        
        y_pos = y_pos - 110;
        
        %--- PEDOT refractive index panel
        pnlRI = uipanel('Parent', scrollPanel, ...
            'Title', 'PEDOT:PSS Refractive Index', ...
            'Units', 'pixels', ...
            'Position', [left_margin, y_pos-380, panel_width, 370]);
        
        % Manual input toggle
        chkManual = uicontrol('Parent', pnlRI, 'Style', 'checkbox', ...
            'String', 'Manual input', 'Value', 0, ...
            'Units', 'pixels', 'Position', [10, 340, 150, 20], ...
            'Callback', @(s,e) on_manual_toggle());
        
        % Single mode section
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'Single Mode:', 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [10, 310, 200, 20]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'Real part (n):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 280, 90, 20]);
        edNreal = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.n_real), ...
            'Units', 'pixels', 'Position', [110, 280, 70, 25], ...
            'Callback', @(s,e) update_ri_display());
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'Imaginary part (k):', 'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [200, 280, 110, 20]);
        edK = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.k_manual), 'Enable', 'off', ...
            'Units', 'pixels', 'Position', [320, 280, 70, 25], ...
            'Callback', @(s,e) update_ri_display());
        
        lblRIDisplay = uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', '', 'ForegroundColor', [0 0 1], ...
            'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [10, 250, 450, 25]);
        
        % Separator
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'For Combined TE/TM Plot:', 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [10, 220, 200, 20]);
        
        % TE and TM mode sections
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'TE Mode:', 'FontAngle', 'italic', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [10, 190, 80, 20]);
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'TM Mode:', 'FontAngle', 'italic', ...
            'HorizontalAlignment', 'left', ...
            'Units', 'pixels', 'Position', [250, 190, 80, 20]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'n (TE):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 160, 60, 20]);
        edNte = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.n_real_te), ...
            'Units', 'pixels', 'Position', [80, 160, 70, 25]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'n (TM):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [250, 160, 60, 20]);
        edNtm = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.n_real_tm), ...
            'Units', 'pixels', 'Position', [320, 160, 70, 25]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'k (TE):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 130, 60, 20]);
        edKte = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.k_te), 'Enable', 'off', ...
            'Units', 'pixels', 'Position', [80, 130, 70, 25]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'k (TM):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [250, 130, 60, 20]);
        edKtm = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.k_tm), 'Enable', 'off', ...
            'Units', 'pixels', 'Position', [320, 130, 70, 25]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'Shift (TE):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [10, 100, 60, 20]);
        edShiftTE = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.shift_te), ...
            'Units', 'pixels', 'Position', [80, 100, 70, 25]);
        
        uicontrol('Parent', pnlRI, 'Style', 'text', ...
            'String', 'Shift (TM):', 'HorizontalAlignment', 'right', ...
            'Units', 'pixels', 'Position', [250, 100, 60, 20]);
        edShiftTM = uicontrol('Parent', pnlRI, 'Style', 'edit', ...
            'String', num2str(app.Defaults.shift_tm), ...
            'Units', 'pixels', 'Position', [320, 100, 70, 25]);
        
        y_pos = y_pos - 400;
        
        %--- Polarization settings
        pnlPol = uipanel('Parent', scrollPanel, ...
            'Title', 'Polarization', ...
            'Units', 'pixels', ...
            'Position', [left_margin, y_pos-60, panel_width, 50]);
        
        bgPol = uibuttongroup('Parent', pnlPol, ...
            'Units', 'pixels', ...
            'Position', [10, 10, 480, 30], ...
            'BorderType', 'none', ...
            'SelectionChangedFcn', @(s,e) update_ri_display());
        
        rbTE = uicontrol('Parent', bgPol, 'Style', 'radiobutton', ...
            'String', 'TE (s-pol)', ...
            'Units', 'pixels', 'Position', [10, 5, 100, 20], ...
            'Value', 1);
        rbTM = uicontrol('Parent', bgPol, 'Style', 'radiobutton', ...
            'String', 'TM (p-pol)', ...
            'Units', 'pixels', 'Position', [150, 5, 100, 20]);
        
        % Initialize display in settings window
        refresh_voltage_info();
        update_ri_display();
    end
    
    function [oldStates] = disable_longrun_controls()
        ctrls_all = [btnOptTE, btnOptTM];
        if ~isempty(btnLoadExcel) && isvalid(btnLoadExcel)
            ctrls_all = [ctrls_all, btnLoadExcel];
        end
        oldStates = cell(1, numel(ctrls_all));
        for i=1:numel(ctrls_all)
            oldStates{i} = get(ctrls_all(i), 'Enable');
        end
        set(ctrls_all, 'Enable','off');
    end
    
    function restore_controls(oldStates)
        ctrls_all = [btnOptTE, btnOptTM];
        if ~isempty(btnLoadExcel) && isvalid(btnLoadExcel)
            ctrls_all = [ctrls_all, btnLoadExcel];
        end
        for i=1:numel(ctrls_all)
            if i <= length(oldStates)
                set(ctrls_all(i), 'Enable', oldStates{i});
            end
        end
    end
 
    function reset_thickness()
        if ~isempty(edITO) && isvalid(edITO)
            set(edITO, 'String', num2str(app.Defaults.ito));
        end
        if ~isempty(edPEDOT) && isvalid(edPEDOT)
            set(edPEDOT, 'String', num2str(app.Defaults.pedot));
        end
    end

    function on_load_excel()
        [f,p] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, 'Select PEDOT extinction coefficient Excel');
        if isequal(f,0); return; end
        full = fullfile(p,f);
        try
            oldStates = disable_longrun_controls();
            d = uiprogressdlg(fig,'Title','Loading Excel','Indeterminate','on',...
                'Message','Reading PEDOT:PSS extinction data...');
            app.PEDOT = load_pedot_data_improved(full);
            close(d);
            restore_controls(oldStates);
            refresh_voltage_info();
            update_ri_display();
        catch ME
            if exist('d','var') && isvalid(d), close(d); end
            if exist('oldStates','var'), restore_controls(oldStates); end
            uialert(fig,ME.message,'Error loading Excel');
        end
    end

    function refresh_voltage_info()
        if isempty(lblAvailV) || ~isvalid(lblAvailV)
            return;
        end
        v = app.PEDOT.voltages;
        v = v(~isnan(v));
        if ~isempty(v)
            set(lblAvailV, 'String', sprintf('Available: %.0f to %.0f mV',min(v),max(v)));
        else
            set(lblAvailV, 'String', 'No voltage data loaded');
        end
        if ~isempty(lblVoltInfo) && isvalid(lblVoltInfo)
            if ~isempty(app.PEDOT.wavelengths)
                set(lblVoltInfo, 'String', sprintf('Wavelengths: %.0f–%.0f nm  |  Columns: %d voltages',...
                    min(app.PEDOT.wavelengths), max(app.PEDOT.wavelengths), numel(app.PEDOT.voltages)));
            else
                set(lblVoltInfo, 'String', '');
            end
        end
    end
    
    % Helper function to safely get UI control values (uses defaults if controls don't exist)
    function val = get_ui_value(control, default_val)
        if ~isempty(control) && isvalid(control)
            val = str2double(get(control, 'String'));
            if isnan(val)
                val = default_val;
            end
        else
            val = default_val;
        end
    end

    function on_manual_toggle()
        isMan = get(chkManual, 'Value');
        if isMan
            set(edK, 'Enable', 'on');
            set(edKte, 'Enable', 'on');
            set(edKtm, 'Enable', 'on');
        else
            set(edK, 'Enable', 'off');
            set(edKte, 'Enable', 'off');
            set(edKtm, 'Enable', 'off');
        end
        update_ri_display();
    end

    % on_use_dk_toggle removed - delta k plotting feature removed

    function on_method_change()
        if isempty(rbTMM) || ~isvalid(rbTMM)
            return;
        end
        if get(rbTMM, 'Value')
            app.Defaults.calculation_method = 'TMM';
            if ~isempty(lblMethodInfo) && isvalid(lblMethodInfo)
                set(lblMethodInfo, 'String', 'TMM: Single n,k for all directions (isotropic)');
            end
        else
            app.Defaults.calculation_method = 'Berreman';
            if ~isempty(lblMethodInfo) && isvalid(lblMethodInfo)
                set(lblMethodInfo, 'String', 'Berreman: Separate n_te, n_tm (anisotropic)');
            end
        end
    end

    function update_ri_display(~,~)
        if isempty(edWavelength) || ~isvalid(edWavelength) || ...
           isempty(edVoltage) || ~isvalid(edVoltage) || ...
           isempty(edNreal) || ~isvalid(edNreal) || ...
           isempty(chkManual) || ~isvalid(chkManual) || ...
           isempty(lblRIDisplay) || ~isvalid(lblRIDisplay)
            return;
        end
        try
            wl = get_ui_value(edWavelength, app.Defaults.wavelength);
            V = get_ui_value(edVoltage, app.Defaults.voltage);
            nreal = str2double(get(edNreal, 'String'));
            if get(chkManual, 'Value')
                if isempty(edK) || ~isvalid(edK)
                    return;
                end
                k = str2double(get(edK, 'String'));
                set(lblRIDisplay, 'String', sprintf('Manual input: n = %.3f + %.6fi', nreal, k));
            else
                k = get_k_improved(app.PEDOT, wl, V);
                set(lblRIDisplay, 'String', sprintf('At %.1f nm, %.1f mV: n = %.3f + %.6fi', wl, V, nreal, k));
            end
        catch
            if ~isempty(lblRIDisplay) && isvalid(lblRIDisplay)
                set(lblRIDisplay, 'String', '');
            end
        end
    end

    %====================
    % Optimization Functions (same as Python)
    %====================

    function on_optimize_te()
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Optimizing','Indeterminate','on','Message','Optimizing TE...');
        try
            [best, ok] = optimize_mode('s', d, 'Finishing...');
            if ~ok
                uialert(fig,'TE experimental data (eploy_22s_Spol.xlsx) not found or invalid.','No Data','warn');
                return;
            end
            
            % Update UI fields based on calculation method (if settings window is open)
            if ~isempty(edNte) && isvalid(edNte)
                if isfield(best, 'calculation_method') && strcmp(best.calculation_method, 'Berreman')
                    % Berreman TE: Show n_te, k_te only
                    set(edNte, 'String', num2str(best.n_te));
                    if ~isempty(edKte) && isvalid(edKte), set(edKte, 'String', num2str(best.k_te)); end
                    if ~isempty(edShiftTE) && isvalid(edShiftTE), set(edShiftTE, 'String', num2str(best.shift)); end
                    if ~isempty(edNreal) && isvalid(edNreal), set(edNreal, 'String', num2str(best.n_te)); end
                    if ~isempty(edK) && isvalid(edK), set(edK, 'String', num2str(best.k_te)); end
                else
                    % TMM: Show single n, k
                    set(edNte, 'String', num2str(best.n_real));
                    if ~isempty(edKte) && isvalid(edKte), set(edKte, 'String', num2str(best.kappa)); end
                    if ~isempty(edShiftTE) && isvalid(edShiftTE), set(edShiftTE, 'String', num2str(best.shift)); end
                    if ~isempty(edNreal) && isvalid(edNreal), set(edNreal, 'String', num2str(best.n_real)); end
                    if ~isempty(edK) && isvalid(edK), set(edK, 'String', num2str(best.kappa)); end
                end
                
                if ~isempty(chkManual) && isvalid(chkManual)
                    set(chkManual, 'Value', 1);
                    on_manual_toggle();
                end
                if ~isempty(edITO) && isvalid(edITO), set(edITO, 'String', num2str(best.ito_thickness)); end
                if ~isempty(edPEDOT) && isvalid(edPEDOT), set(edPEDOT, 'String', num2str(best.pedot_thickness)); end
            end
            
            app.OptimizedShift = best.shift;
            app.LastOptimizationResults.TE = best;
            uialert(fig,format_best_message('TE',best),'TE Optimization Complete');
        catch ME
            uialert(fig,ME.message,'TE Optimization Error');
        end
        if isvalid(d), close(d); end
        restore_controls(oldStates);
    end

    function on_optimize_tm()
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Optimizing','Indeterminate','on','Message','Optimizing TM...');
        try
            [best, ok] = optimize_mode('p', d, 'Finishing...');
            if ~ok
                uialert(fig,'TM experimental data (eploy_22s_Ppol.xlsx) not found or invalid.','No Data','warn');
                return;
            end
            
            % Update UI fields based on calculation method (if settings window is open)
            if ~isempty(edNtm) && isvalid(edNtm)
                if isfield(best, 'calculation_method') && strcmp(best.calculation_method, 'Berreman')
                    % Berreman TM: Show both n_te, k_te AND n_tm, k_tm
                    if ~isempty(edNte) && isvalid(edNte), set(edNte, 'String', num2str(best.n_te)); end
                    if ~isempty(edKte) && isvalid(edKte), set(edKte, 'String', num2str(best.k_te)); end
                    set(edNtm, 'String', num2str(best.n_tm));
                    if ~isempty(edKtm) && isvalid(edKtm), set(edKtm, 'String', num2str(best.k_tm)); end
                    if ~isempty(edShiftTM) && isvalid(edShiftTM), set(edShiftTM, 'String', num2str(best.shift)); end
                    if ~isempty(edNreal) && isvalid(edNreal), set(edNreal, 'String', num2str(best.n_tm)); end
                    if ~isempty(edK) && isvalid(edK), set(edK, 'String', num2str(best.k_tm)); end
                else
                    % TMM: Show single n, k
                    set(edNtm, 'String', num2str(best.n_real));
                    if ~isempty(edKtm) && isvalid(edKtm), set(edKtm, 'String', num2str(best.kappa)); end
                    if ~isempty(edShiftTM) && isvalid(edShiftTM), set(edShiftTM, 'String', num2str(best.shift)); end
                    if ~isempty(edNreal) && isvalid(edNreal), set(edNreal, 'String', num2str(best.n_real)); end
                    if ~isempty(edK) && isvalid(edK), set(edK, 'String', num2str(best.kappa)); end
                end
                
                if ~isempty(chkManual) && isvalid(chkManual)
                    set(chkManual, 'Value', 1);
                    on_manual_toggle();
                end
                if ~isempty(edITO) && isvalid(edITO), set(edITO, 'String', num2str(best.ito_thickness)); end
                if ~isempty(edPEDOT) && isvalid(edPEDOT), set(edPEDOT, 'String', num2str(best.pedot_thickness)); end
            end
            
            app.OptimizedShift = best.shift;
            app.LastOptimizationResults.TM = best;
            uialert(fig,format_best_message('TM',best),'TM Optimization Complete');
        catch ME
            uialert(fig,ME.message,'TM Optimization Error');
        end
        if isvalid(d), close(d); end
        restore_controls(oldStates);
    end
    

    function on_export_results()
        % Check if we have optimization results with local minima
        has_local_minima = false;
        if isfield(app.LastOptimizationResults, 'TE') && isfield(app.LastOptimizationResults.TE, 'all_local_minima')
            has_local_minima = true;
        end
        if isfield(app.LastOptimizationResults, 'TM') && isfield(app.LastOptimizationResults.TM, 'all_local_minima')
            has_local_minima = true;
        end
        
        if has_local_minima
            % Ask user what to export
            choice = questdlg(['What would you like to export?' newline newline ...
                '1. All Local Minima - Export all local minima found during optimization' newline ...
                '2. Best Results Only - Export only the best optimization results'], ...
                'Export Options', ...
                'All Local Minima', 'Best Results Only', 'Cancel', ...
                'All Local Minima');
            
            if strcmp(choice, 'Cancel')
                return;
            elseif strcmp(choice, 'All Local Minima')
                export_all_local_minima();
                return;
            end
        end
        
        % Check if we have optimization results to export
        hasResults = false;
        if isfield(app, 'LastOptimizationResults')
            if isfield(app.LastOptimizationResults, 'TE') || isfield(app.LastOptimizationResults, 'TM')
                hasResults = true;
            end
        end
        
        if ~hasResults
            msgbox('No optimization results to export. Please run optimization first.','No Data','warn');
            return;
        end
        
        % For single results, no need to warn about large datasets
        totalResults = 1; % Single result
        
        % Warn about large datasets
        if totalResults > 1000000
            choice = uiconfirm(fig, sprintf(['Large dataset detected: %d results\n\n' ...
                'This may cause memory issues or exceed Excel limits.\n' ...
                'Options:\n' ...
                '1. Export only top 100,000 results (recommended)\n' ...
                '2. Export all results (may fail)\n' ...
                '3. Cancel'], totalResults), ...
                'Large Dataset Warning', 'Options', {'Top 100k', 'All Results', 'Cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 3);
            
            if strcmp(choice, 'Cancel')
                return;
            elseif strcmp(choice, 'Top 100k')
                exportTopResults = true;
            else
                exportTopResults = false;
            end
        else
            exportTopResults = false;
        end
        
        [filename, pathname] = uiputfile({'*.xlsx','Excel files (*.xlsx)'}, 'Save optimization results as...', 'optimization_results.xlsx');
        if isequal(filename,0)
            return; % User cancelled
        end
        
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Exporting Results','Indeterminate','on',...
            'Message','Preparing data for export...');
        
        try
            fullpath = fullfile(pathname, filename);
            
            % Prepare data for export (single results)
            headers = {'Mode', 'ITO_Thickness_nm', 'PEDOT_Thickness_nm', 'n_real', 'k_imaginary', 'Angle_Shift_deg', 'Error', 'Valid_Points', 'Method', ...
                'Mean_Error_Percent', 'Std_Error_Percent', 'RMSE_Percent', 'MAE_Percent', 'Max_Error_Percent', 'Min_Error_Percent'};
            allData = {};
            rowCount = 0;
            
            % Export TE result if available
            if isfield(app.LastOptimizationResults, 'TE')
                te = app.LastOptimizationResults.TE;
                rowCount = rowCount + 1;
                % Get percentage error values
                if isfield(te, 'percent_error_stats') && isfield(te.percent_error_stats, 'mean_error_percent')
                    pe = te.percent_error_stats;
                    pe_vals = {pe.mean_error_percent, pe.std_error_percent, pe.rmse_percent, ...
                        pe.mae_percent, pe.max_error_percent, pe.min_error_percent};
                else
                    pe_vals = {NaN, NaN, NaN, NaN, NaN, NaN};
                end
                allData{rowCount,1} = 'TE';
                allData{rowCount,2} = te.ito_thickness;
                allData{rowCount,3} = te.pedot_thickness;
                allData{rowCount,4} = te.n_real;
                allData{rowCount,5} = te.kappa;
                allData{rowCount,6} = te.shift;
                allData{rowCount,7} = te.err;
                allData{rowCount,8} = te.num_valid_points;
                allData{rowCount,9} = te.optimization_method;
                for col = 10:15
                    allData{rowCount, col} = pe_vals{col-9};
                end
            end
            
            % Export TM result if available
            if isfield(app.LastOptimizationResults, 'TM')
                tm = app.LastOptimizationResults.TM;
                rowCount = rowCount + 1;
                % Get percentage error values
                if isfield(tm, 'percent_error_stats') && isfield(tm.percent_error_stats, 'mean_error_percent')
                    pe = tm.percent_error_stats;
                    pe_vals = {pe.mean_error_percent, pe.std_error_percent, pe.rmse_percent, ...
                        pe.mae_percent, pe.max_error_percent, pe.min_error_percent};
                else
                    pe_vals = {NaN, NaN, NaN, NaN, NaN, NaN};
                end
                allData{rowCount,1} = 'TM';
                allData{rowCount,2} = tm.ito_thickness;
                allData{rowCount,3} = tm.pedot_thickness;
                allData{rowCount,4} = tm.n_real;
                allData{rowCount,5} = tm.kappa;
                allData{rowCount,6} = tm.shift;
                allData{rowCount,7} = tm.err;
                allData{rowCount,8} = tm.num_valid_points;
                allData{rowCount,9} = tm.optimization_method;
                for col = 10:15
                    allData{rowCount, col} = pe_vals{col-9};
                end
            end
            
            % Create table and export
            T = array2table(allData, 'VariableNames', headers);
            writetable(T, fullpath);
            
            uialert(fig, sprintf('Optimization results exported successfully!\nFile: %s\nTotal results: %d', fullpath, rowCount), 'Export Complete');
            
        catch ME
            uialert(fig, ME.message, 'Export Error');
        end
        
        if isvalid(d), close(d); end
        restore_controls(oldStates);
    end
    
    function export_top_results_only(fullpath, totalResults)
        % Export only the top 100,000 results to avoid memory issues
        headers = {'Mode', 'ITO_Thickness_nm', 'PEDOT_Thickness_nm', 'n_real', 'k_imaginary', 'Angle_Shift_deg', 'Error', 'Valid_Points'};
        exportData = cell(100000, 8); % Pre-allocate for 100k results
        rowCount = 0;
        
        % Process TE results (top 50k)
        if isfield(app.AllOptimizationResults, 'TE')
            teResults = app.AllOptimizationResults.TE;
            [~, sortedIdx] = sort([teResults.error]); % Sort by error (ascending)
            maxTE = min(50000, length(teResults));
            
            for i = 1:maxTE
                idx = sortedIdx(i);
                rowCount = rowCount + 1;
                exportData{rowCount,1} = 'TE';
                exportData{rowCount,2} = teResults(idx).ito_thickness;
                exportData{rowCount,3} = teResults(idx).pedot_thickness;
                exportData{rowCount,4} = teResults(idx).n_real;
                exportData{rowCount,5} = teResults(idx).kappa;
                exportData{rowCount,6} = teResults(idx).shift;
                exportData{rowCount,7} = teResults(idx).error;
                exportData{rowCount,8} = teResults(idx).valid_points;
            end
        end
        
        % Process TM results (top 50k)
        if isfield(app.AllOptimizationResults, 'TM')
            tmResults = app.AllOptimizationResults.TM;
            [~, sortedIdx] = sort([tmResults.error]); % Sort by error (ascending)
            maxTM = min(50000, length(tmResults));
            
            for i = 1:maxTM
                idx = sortedIdx(i);
                rowCount = rowCount + 1;
                exportData{rowCount,1} = 'TM';
                exportData{rowCount,2} = tmResults(idx).ito_thickness;
                exportData{rowCount,3} = tmResults(idx).pedot_thickness;
                exportData{rowCount,4} = tmResults(idx).n_real;
                exportData{rowCount,5} = tmResults(idx).kappa;
                exportData{rowCount,6} = tmResults(idx).shift;
                exportData{rowCount,7} = tmResults(idx).error;
                exportData{rowCount,8} = tmResults(idx).valid_points;
            end
        end
        
        % Trim unused rows and create table
        exportData = exportData(1:rowCount, :);
        T = array2table(exportData, 'VariableNames', headers);
        writetable(T, fullpath);
    end
    
    function export_all_results_chunked(fullpath, totalResults)
        % Export all results using chunking to handle large datasets
        headers = {'Mode', 'ITO_Thickness_nm', 'PEDOT_Thickness_nm', 'n_real', 'k_imaginary', 'Angle_Shift_deg', 'Error', 'Valid_Points'};
        chunkSize = 50000; % Process 50k rows at a time
        allData = {};
        
        % Process TE results in chunks
        if isfield(app.AllOptimizationResults, 'TE')
            teResults = app.AllOptimizationResults.TE;
            numTE = length(teResults);
            numChunks = ceil(numTE / chunkSize);
            
            for chunk = 1:numChunks
                startIdx = (chunk-1) * chunkSize + 1;
                endIdx = min(chunk * chunkSize, numTE);
                chunkData = cell(endIdx - startIdx + 1, 8);
                
                for i = startIdx:endIdx
                    rowIdx = i - startIdx + 1;
                    chunkData{rowIdx,1} = 'TE';
                    chunkData{rowIdx,2} = teResults(i).ito_thickness;
                    chunkData{rowIdx,3} = teResults(i).pedot_thickness;
                    chunkData{rowIdx,4} = teResults(i).n_real;
                    chunkData{rowIdx,5} = teResults(i).kappa;
                    chunkData{rowIdx,6} = teResults(i).shift;
                    chunkData{rowIdx,7} = teResults(i).error;
                    chunkData{rowIdx,8} = teResults(i).valid_points;
                end
                
                allData = [allData; chunkData];
                
                % Update progress
                if isvalid(d)
                    d.Message = sprintf('Processing TE results... (%d/%d chunks)', chunk, numChunks);
                end
            end
        end
        
        % Process TM results in chunks
        if isfield(app.AllOptimizationResults, 'TM')
            tmResults = app.AllOptimizationResults.TM;
            numTM = length(tmResults);
            numChunks = ceil(numTM / chunkSize);
            
            for chunk = 1:numChunks
                startIdx = (chunk-1) * chunkSize + 1;
                endIdx = min(chunk * chunkSize, numTM);
                chunkData = cell(endIdx - startIdx + 1, 8);
                
                for i = startIdx:endIdx
                    rowIdx = i - startIdx + 1;
                    chunkData{rowIdx,1} = 'TM';
                    chunkData{rowIdx,2} = tmResults(i).ito_thickness;
                    chunkData{rowIdx,3} = tmResults(i).pedot_thickness;
                    chunkData{rowIdx,4} = tmResults(i).n_real;
                    chunkData{rowIdx,5} = tmResults(i).kappa;
                    chunkData{rowIdx,6} = tmResults(i).shift;
                    chunkData{rowIdx,7} = tmResults(i).error;
                    chunkData{rowIdx,8} = tmResults(i).valid_points;
                end
                
                allData = [allData; chunkData];
                
                % Update progress
                if isvalid(d)
                    d.Message = sprintf('Processing TM results... (%d/%d chunks)', chunk, numChunks);
                end
            end
        end
        
        % Write to Excel
        if isvalid(d)
            d.Message = 'Writing to Excel file...';
        end
        
        T = array2table(allData, 'VariableNames', headers);
        writetable(T, fullpath);
    end
    
    function export_all_local_minima()
        % Export all local minima found during Particle Swarm optimization
        [filename, pathname] = uiputfile({'*.xlsx','Excel files (*.xlsx)'}, 'Save all local minima as...', 'all_local_minima.xlsx');
        if isequal(filename,0)
            return; % User cancelled
        end
        
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Exporting Local Minima','Indeterminate','on',...
            'Message','Preparing local minima data for export...');
        
        try
            fullpath = fullfile(pathname, filename);
            
            % Prepare data for export
            allData = {};
            rowCount = 0;
            
            % Process TE local minima
            if isfield(app.LastOptimizationResults, 'TE') && isfield(app.LastOptimizationResults.TE, 'all_local_minima')
                te_minima = app.LastOptimizationResults.TE.all_local_minima;
                for i = 1:length(te_minima)
                    rowCount = rowCount + 1;
                    params = te_minima(i).parameters;
                    if length(params) == 3 % Fixed thickness mode
                        allData{rowCount,1} = 'TE';
                        allData{rowCount,2} = app.LastOptimizationResults.TE.ito_thickness; % Fixed
                        allData{rowCount,3} = app.LastOptimizationResults.TE.pedot_thickness; % Fixed
                        allData{rowCount,4} = params(1); % n_real
                        allData{rowCount,5} = params(2); % kappa
                        allData{rowCount,6} = params(3); % angle_shift
                        allData{rowCount,7} = te_minima(i).fval; % error
                        allData{rowCount,8} = te_minima(i).iteration; % iteration found
                        allData{rowCount,9} = te_minima(i).particle_id; % particle ID
                    else % All parameters mode
                        allData{rowCount,1} = 'TE';
                        allData{rowCount,2} = params(1); % ITO thickness
                        allData{rowCount,3} = params(2); % PEDOT thickness
                        allData{rowCount,4} = params(3); % n_real
                        allData{rowCount,5} = params(4); % kappa
                        allData{rowCount,6} = params(5); % angle_shift
                        allData{rowCount,7} = te_minima(i).fval; % error
                        allData{rowCount,8} = te_minima(i).iteration; % iteration found
                        allData{rowCount,9} = te_minima(i).particle_id; % particle ID
                    end
                end
            end
            
            % Process TM local minima
            if isfield(app.LastOptimizationResults, 'TM') && isfield(app.LastOptimizationResults.TM, 'all_local_minima')
                tm_minima = app.LastOptimizationResults.TM.all_local_minima;
                for i = 1:length(tm_minima)
                    rowCount = rowCount + 1;
                    params = tm_minima(i).parameters;
                    if length(params) == 3 % Fixed thickness mode
                        allData{rowCount,1} = 'TM';
                        allData{rowCount,2} = app.LastOptimizationResults.TM.ito_thickness; % Fixed
                        allData{rowCount,3} = app.LastOptimizationResults.TM.pedot_thickness; % Fixed
                        allData{rowCount,4} = params(1); % n_real
                        allData{rowCount,5} = params(2); % kappa
                        allData{rowCount,6} = params(3); % angle_shift
                        allData{rowCount,7} = tm_minima(i).fval; % error
                        allData{rowCount,8} = tm_minima(i).iteration; % iteration found
                        allData{rowCount,9} = tm_minima(i).particle_id; % particle ID
                    else % All parameters mode
                        allData{rowCount,1} = 'TM';
                        allData{rowCount,2} = params(1); % ITO thickness
                        allData{rowCount,3} = params(2); % PEDOT thickness
                        allData{rowCount,4} = params(3); % n_real
                        allData{rowCount,5} = params(4); % kappa
                        allData{rowCount,6} = params(5); % angle_shift
                        allData{rowCount,7} = tm_minima(i).fval; % error
                        allData{rowCount,8} = tm_minima(i).iteration; % iteration found
                        allData{rowCount,9} = tm_minima(i).particle_id; % particle ID
                    end
                end
            end
            
            if rowCount == 0
                uialert(fig,'No local minima data found to export.','No Data');
                return;
            end
            
            % Create table and export
            headers = {'Mode', 'ITO_Thickness_nm', 'PEDOT_Thickness_nm', 'n_real', 'k_imaginary', 'Angle_Shift_deg', 'Error', 'Iteration_Found', 'Particle_ID'};
            T = array2table(allData, 'VariableNames', headers);
            writetable(T, fullpath);
            
            uialert(fig, sprintf('All local minima exported successfully!\nFile: %s\nTotal solutions: %d', fullpath, rowCount), 'Export Complete');
            
        catch ME
            uialert(fig, ME.message, 'Export Error');
        end
        
        if isvalid(d), close(d); end
        restore_controls(oldStates);
    end

    function msg = compose_results_message(bestTE,hasTE, bestTM,hasTM)
        parts = {};
        if hasTE
            parts{end+1} = format_best_message('TE',bestTE);
        else
            parts{end+1} = 'TE: No experimental data found';
        end
        if hasTM
            parts{end+1} = format_best_message('TM',bestTM);
        else
            parts{end+1} = 'TM: No experimental data found';
        end
        msg = strjoin(parts, sprintf('\n\n'));
    end

    function msg = format_best_message(mode, best)
        % Check if Berreman method was used
        is_berreman = isfield(best, 'calculation_method') && strcmp(best.calculation_method, 'Berreman');
        
        % Check if percentage error stats are available
        has_percent_errors = isfield(best, 'percent_error_stats') && ...
            isfield(best.percent_error_stats, 'mean_error_percent') && ...
            ~isnan(best.percent_error_stats.mean_error_percent);
        
        % Helper function to append percentage error info
        function err_str = get_percent_error_string()
            if ~has_percent_errors
                err_str = '';
                return;
            end
            pe = best.percent_error_stats;
            err_str = sprintf(['\n  Percentage Error Statistics:\n' ...
                '    Mean Error: %.2f%%\n' ...
                '    Std Dev: %.2f%%\n' ...
                '    RMSE: %.2f%%\n' ...
                '    MAE: %.2f%%\n' ...
                '    Max Error: %.2f%%\n' ...
                '    Min Error: %.2f%%'], ...
                pe.mean_error_percent, pe.std_error_percent, pe.rmse_percent, ...
                pe.mae_percent, pe.max_error_percent, pe.min_error_percent);
        end
        
        if is_berreman
            % Format message for Berreman results - depends on polarization
            if strcmp(mode, 'TE')
                % TE mode: Only show in-plane indices (n_te, k_te)
            msg = sprintf(['%s Mode (Berreman):\n' ...
                    '  n_te = %.4f + %.4fi  (in-plane, ordinary)\n' ...
                    '  Angle shift = %.2f°\n' ...
                    '  ITO thickness = %.1f nm\n' ...
                    '  PEDOT thickness = %.1f nm\n' ...
                    '  Error = %.3e%s'], ...
                    mode, best.n_te, best.k_te, best.shift, ...
                    best.ito_thickness, best.pedot_thickness, best.err, get_percent_error_string());
            else
                % TM mode: Show both in-plane and out-of-plane indices
                msg = sprintf(['%s Mode (Berreman):\n' ...
                    '  n_te = %.4f + %.4fi  (in-plane, ordinary)\n' ...
                    '  n_tm = %.4f + %.4fi  (out-of-plane, extraordinary)\n' ...
                '  Angle shift = %.2f°\n' ...
                '  ITO thickness = %.1f nm\n' ...
                '  PEDOT thickness = %.1f nm\n' ...
                '  Error = %.3e%s'], ...
                mode, best.n_te, best.k_te, best.n_tm, best.k_tm, best.shift, ...
                best.ito_thickness, best.pedot_thickness, best.err, get_percent_error_string());
            end
        elseif isfield(best, 'optimization_method')
            % Format message for TMM results - show separate labels for TE and TM
            if strcmp(mode, 'TE')
                % TE mode: Show n_te, k_te
                if contains(best.optimization_method, 'lsqnonlin')
                    if contains(best.optimization_method, 'fixed_thickness')
                        opt_type = 'Fixed Thickness';
                        thickness_info = sprintf('(ITO=%.1f nm, PEDOT=%.1f nm - FIXED)', best.ito_thickness, best.pedot_thickness);
                    else
                        opt_type = 'All Parameters';
                        thickness_info = sprintf('ITO=%.1f nm, PEDOT=%.1f nm', best.ito_thickness, best.pedot_thickness);
                    end
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  Angle shift = %.2f°\n  %s\n  Error = %.3e\n  Iterations: %d%s', ...
                        mode, opt_type, best.n_real, best.kappa, best.shift, thickness_info, best.err, best.iterations, get_percent_error_string());
                elseif strcmp(best.optimization_method, 'global_search')
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Starting points: %d, Iterations: %d%s', ...
                        mode, 'Global Search', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.num_starting_points, best.iterations, get_percent_error_string());
                elseif contains(best.optimization_method, 'particleswarm')
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Iterations: %d%s', ...
                        mode, 'Particle Swarm', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.iterations, get_percent_error_string());
                else
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Combinations tested: %d%s', ...
                        mode, 'Grid Search', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.total_combinations, get_percent_error_string());
                end
            else
                % TM mode: Show n_te, k_te, n_tm, k_tm (TMM uses same values for both)
                if contains(best.optimization_method, 'lsqnonlin')
                    if contains(best.optimization_method, 'fixed_thickness')
                        opt_type = 'Fixed Thickness';
                        thickness_info = sprintf('(ITO=%.1f nm, PEDOT=%.1f nm - FIXED)', best.ito_thickness, best.pedot_thickness);
                    else
                        opt_type = 'All Parameters';
                        thickness_info = sprintf('ITO=%.1f nm, PEDOT=%.1f nm', best.ito_thickness, best.pedot_thickness);
                    end
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  n_tm = %.4f + %.4fi\n  Angle shift = %.2f°\n  %s\n  Error = %.3e\n  Iterations: %d%s', ...
                        mode, opt_type, best.n_real, best.kappa, best.n_real, best.kappa, best.shift, thickness_info, best.err, best.iterations, get_percent_error_string());
                elseif strcmp(best.optimization_method, 'global_search')
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  n_tm = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Starting points: %d, Iterations: %d%s', ...
                        mode, 'Global Search', best.n_real, best.kappa, best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.num_starting_points, best.iterations, get_percent_error_string());
                elseif contains(best.optimization_method, 'particleswarm')
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  n_tm = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Iterations: %d%s', ...
                        mode, 'Particle Swarm', best.n_real, best.kappa, best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.iterations, get_percent_error_string());
                else
                    msg = sprintf('%s Mode (%s - TMM):\n  n_te = %.4f + %.4fi\n  n_tm = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Combinations tested: %d%s', ...
                        mode, 'Grid Search', best.n_real, best.kappa, best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.total_combinations, get_percent_error_string());
                end
            end
        else
            % Fallback for old format - show based on mode
            if strcmp(mode, 'TE')
                msg = sprintf('%s Mode:\n  n_te = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e%s', ...
                    mode, best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, get_percent_error_string());
            else
                msg = sprintf('%s Mode:\n  n_te = %.4f + %.4fi\n  n_tm = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e%s', ...
                    mode, best.n_real, best.kappa, best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, get_percent_error_string());
            end
        end
    end

    function [best, ok] = optimize_mode(pol, dlg, nextMessage)
        [expData, ok] = load_experimental_data(pol);
        if ~ok
            best = struct(); return;
        end
        
        % Get the current calculation method
        calc_method = app.Defaults.calculation_method;
        
        % Choose optimization method
        choice = uiconfirm(fig, sprintf(['Choose optimization method for %s mode (%s):\n\n' ...
            '1. Particle Swarm (All) - Optimize all parameters including thickness\n' ...
            '2. Particle Swarm (Fixed) - Only optimize n, k, and angle shift\n' ...
            '3. Global Search - Multiple lsqnonlin starting points\n' ...
            '4. Advanced (lsqnonlin) - Fast single optimization'], pol, calc_method), ...
            'Optimization Method', 'Options', {'Particle Swarm (All)', 'Particle Swarm (Fixed)', 'Global Search', 'Advanced (All)'}, ...
            'DefaultOption', 1);
        
        if strcmp(choice, 'Particle Swarm (All)')
            best = optimize_mode_particleswarm(pol, expData, dlg, calc_method);
            ok = true;
        elseif strcmp(choice, 'Particle Swarm (Fixed)')
            best = optimize_mode_particleswarm_fixed_thickness(pol, expData, dlg, calc_method);
            ok = true;
        elseif strcmp(choice, 'Global Search')
            best = optimize_mode_global(pol, expData, dlg, calc_method);
            ok = true;
        elseif strcmp(choice, 'Advanced (All)')
            best = optimize_mode_lsqnonlin(pol, expData, dlg, 'all', calc_method);
            ok = true;
        else % Cancelled or closed
            best = struct();
            ok = false;
            return;
        end
        
        if isvalid(dlg), dlg.Message = nextMessage; end
    end
    
    function best = optimize_mode_lsqnonlin(pol, expData, dlg, mode, calc_method)
        % Advanced optimization using lsqnonlin
        if nargin < 4
            mode = 'all'; % Default to optimizing all parameters
        end
        if nargin < 5
            calc_method = 'TMM';
        end
        
        if isvalid(dlg)
            if strcmp(mode, 'fixed_thickness')
                dlg.Message = 'Setting up fixed-thickness optimization...';
            else
                dlg.Message = 'Setting up advanced optimization...';
            end
        end
        
        wl = get_ui_value(edWavelength, app.Defaults.wavelength);
        V = get_ui_value(edVoltage, app.Defaults.voltage);
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Get current thickness values from UI (use defaults if settings window not open)
        current_ito = get_ui_value(edITO, app.Defaults.ito);
        current_pedot = get_ui_value(edPEDOT, app.Defaults.pedot);
        
        if strcmp(mode, 'fixed_thickness')
            if strcmp(calc_method, 'Berreman')
                % Berreman: [n_te, k_te, n_tm, k_tm, angle_shift]
                lb = [1.0, 0.00, 1.0, 0.00, -3];
                ub = [1.5, 0.15, 1.5, 0.15, 3];
                x0 = [1.41, 0.05, 1.269, 0.05, 0];
                % Note: For Berreman fixed thickness, we need a special objective function
                error('Berreman fixed thickness not yet implemented for lsqnonlin');
            else
                % TMM: [n_real, kappa, angle_shift]
                lb = [1.3, 0.00, -3];    % Lower bounds
                ub = [1.5, 0.15, 3];     % Upper bounds
                x0 = [1.4, 0.05, 0];     % Initial guess
                
                % Create objective function for fixed thickness
                objective_fun = @(params) objective_function_fixed_thickness(params, current_ito, current_pedot, ...
                    pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl);
            end
            
            optimization_type = 'Fixed Thickness';
        else
            % All parameters mode
            if strcmp(calc_method, 'Berreman')
                % Berreman: [ITO_thickness, PEDOT_thickness, n_te, k_te, n_tm, k_tm, angle_shift]
                lb = [5, 40, 1.0, 0.00, 1.0, 0.00, -5];
                ub = [30, 80, 1.5, 0.25, 1.5, 0.25, 5];
                    x0 = [current_ito, current_pedot, 1.41, 0.05, 1.269, 0.05, 0];
                
                % Create objective function for Berreman
                objective_fun = @(params) objective_function_berreman(params, pol, exp_angles, exp_reflectivities, ...
                    n_bk, n_ito, n_water, degree, wl);
            else
                % TMM: [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
                lb = [5, 40, 1.0, 0.01, -5];    % Lower bounds
                ub = [30, 80, 2.0, 1.0, 5];     % Upper bounds
                x0 = [current_ito, current_pedot, 1.5, 0.5, 0]; % Use current thicknesses as starting point
                
                % Create objective function for TMM
                objective_fun = @(params) objective_function(params, pol, exp_angles, exp_reflectivities, ...
                    n_bk, n_ito, n_water, degree, wl);
            end
            
            optimization_type = 'All Parameters';
        end
        
        % Optimization options
        options = optimoptions('lsqnonlin', ...
            'Display', 'off', ...
            'MaxIterations', 1000, ...
            'MaxFunctionEvaluations', 5000, ...
            'OptimalityTolerance', 1e-8, ...
            'StepTolerance', 1e-12, ...
            'FunctionTolerance', 1e-8);
        
        if isvalid(dlg)
            dlg.Message = sprintf('Running %s optimization...', optimization_type);
        end
        
        try
            % Run optimization
            [x_opt, resnorm, residual, exitflag, output] = lsqnonlin(objective_fun, x0, lb, ub, options);
            
            % Extract results based on mode and calculation method
            if strcmp(mode, 'fixed_thickness')
                % Fixed thickness mode
                best.ito_thickness = current_ito;
                best.pedot_thickness = current_pedot;
                if strcmp(calc_method, 'Berreman')
                    % Not yet implemented
                    error('Berreman fixed thickness not implemented');
                else
                    % TMM: [n_real, kappa, angle_shift]
                    best.n_real = x_opt(1);
                    best.kappa = x_opt(2);
                    best.shift = x_opt(3);
                end
            else
                % All parameters mode
                if strcmp(calc_method, 'Berreman')
                    % Berreman: [ITO_thickness, PEDOT_thickness, n_te, k_te, n_tm, k_tm, angle_shift]
                    best.ito_thickness = x_opt(1);
                    best.pedot_thickness = x_opt(2);
                    best.n_te = x_opt(3);
                    best.k_te = x_opt(4);
                    best.n_tm = x_opt(5);
                    best.k_tm = x_opt(6);
                    best.shift = x_opt(7);
                    best.n_real = x_opt(3);  % Store n_te as n_real for compatibility
                    best.kappa = x_opt(4);   % Store k_te as kappa for compatibility
                else
                    % TMM: [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
                    best.ito_thickness = x_opt(1);
                    best.pedot_thickness = x_opt(2);
                    best.n_real = x_opt(3);
                    best.kappa = x_opt(4);
                    best.shift = x_opt(5);
                end
            end
            
            best.err = resnorm;
            best.num_valid_points = length(exp_reflectivities);
            
            % Additional optimization info
            best.optimization_method = sprintf('lsqnonlin_%s_%s', mode, lower(calc_method));
            best.calculation_method = calc_method;
            best.optimization_type = optimization_type;
            best.exitflag = exitflag;
            best.iterations = output.iterations;
            best.function_count = output.funccount;
            
            if isvalid(dlg)
                dlg.Message = sprintf('%s optimization complete (%d iterations)', optimization_type, output.iterations);
            end
            
        catch ME
            warning('%s optimization failed: %s', optimization_type, ME.message);
            if isvalid(dlg)
                dlg.Message = sprintf('%s optimization failed, falling back to grid search...', optimization_type);
            end
            best = optimize_mode_grid(pol, expData, dlg, calc_method);
        end
    end
    
    function best = optimize_mode_global(pol, expData, dlg, calc_method)
        % Global optimization using multiple starting points to avoid local minima
        if nargin < 4
            calc_method = 'TMM'; % Default to TMM if not specified
        end
        
        if isvalid(dlg)
            dlg.Message = sprintf('Setting up global optimization (%s)...', calc_method);
        end
        
        wl = get_ui_value(edWavelength, app.Defaults.wavelength);
        V = get_ui_value(edVoltage, app.Defaults.voltage);
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Get current thickness values from UI (use defaults if settings window not open)
        current_ito = get_ui_value(edITO, app.Defaults.ito);
        current_pedot = get_ui_value(edPEDOT, app.Defaults.pedot);
        
        % Parameter bounds depend on calculation method
        if strcmp(calc_method, 'Berreman')
            % For Berreman: [ITO_thickness, PEDOT_thickness, n_te, k_te, n_tm, k_tm, angle_shift]
            lb = [5, 40, 1.0, 0.00, 1.0, 0.00, -5];
            ub = [30, 80, 1.5, 0.25, 1.5, 0.25, 5];
        else
            % For TMM: [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
        lb = [5, 40, 1.0, 0.01, -5];    % Lower bounds
        ub = [30, 80, 2.0, 1.0, 5];     % Upper bounds
        end
        
        % Create objective function
        if strcmp(calc_method, 'Berreman')
            objective_fun = @(params) objective_function_berreman(params, pol, exp_angles, exp_reflectivities, ...
                n_bk, n_ito, n_water, degree, wl);
        else
        objective_fun = @(params) objective_function(params, pol, exp_angles, exp_reflectivities, ...
            n_bk, n_ito, n_water, degree, wl);
        end
        
        % Multiple starting points to explore different regions
        num_starts = 8;
        if strcmp(calc_method, 'Berreman')
            % For Berreman: [ITO, PEDOT, n_te, k_te, n_tm, k_tm, angle_shift]
            starting_points = [
                [current_ito, current_pedot, 1.41, 0.05, 1.269, 0.05, 0];
                [17.5, 60, 1.41, 0.05, 1.269, 0.05, 0];
                [10, 50, 1.35, 0.03, 1.20, 0.03, -2];
                [25, 70, 1.45, 0.10, 1.35, 0.10, 2];
                [15, 45, 1.30, 0.02, 1.15, 0.02, -1];
                [30, 75, 1.48, 0.15, 1.40, 0.15, 1];
                [20, 55, 1.38, 0.07, 1.25, 0.07, -0.5];
                [12, 65, 1.43, 0.08, 1.30, 0.08, 0.5];
            ];
        else
            % For TMM: [ITO, PEDOT, n_real, kappa, angle_shift]
        starting_points = [
            [current_ito, current_pedot, 1.5, 0.5, 0];           % Current UI values
            [17.5, 60, 1.5, 0.5, 0];                            % Center of parameter space
            [10, 50, 1.2, 0.3, -2];                             % Low thickness, low n
            [25, 70, 1.8, 0.7, 2];                              % High thickness, high n
            [15, 45, 1.1, 0.2, -1];                             % Low values
            [30, 75, 1.9, 0.8, 1];                              % High values
            [20, 55, 1.4, 0.4, -0.5];                           % Medium values
            [12, 65, 1.6, 0.6, 0.5];                            % Mixed values
        ];
        end
        
        % Optimization options for lsqnonlin
        options = optimoptions('lsqnonlin', ...
            'Display', 'off', ...
            'MaxIterations', 500, ...
            'MaxFunctionEvaluations', 2000, ...
            'OptimalityTolerance', 1e-8, ...
            'StepTolerance', 1e-12, ...
            'FunctionTolerance', 1e-8);
        
        best_error = Inf;
        best_params = [];
        best_output = [];
        
        if isvalid(dlg)
            dlg.Message = sprintf('Running global optimization (%d starting points)...', num_starts);
        end
        
        % Try optimization from each starting point
        for i = 1:num_starts
            if isvalid(dlg)
                dlg.Message = sprintf('Global optimization: Starting point %d/%d', i, num_starts);
            end
            
            try
                [x_opt, resnorm, ~, exitflag, output] = lsqnonlin(objective_fun, starting_points(i,:), lb, ub, options);
                
                if resnorm < best_error
                    best_error = resnorm;
                    best_params = x_opt;
                    best_output = output;
                end
                
            catch ME
                warning('Optimization failed for starting point %d: %s', i, ME.message);
                continue;
            end
        end
        
        if isempty(best_params)
            warning('All global optimization attempts failed, falling back to grid search...');
            best = optimize_mode_grid(pol, expData, dlg);
            return;
        end
        
        % Extract results based on calculation method
        if strcmp(calc_method, 'Berreman')
            % Berreman: [ITO_thickness, PEDOT_thickness, n_te, k_te, n_tm, k_tm, angle_shift]
            best.ito_thickness = best_params(1);
            best.pedot_thickness = best_params(2);
            best.n_real = best_params(3);  % Store n_te as n_real for compatibility
            best.kappa = best_params(4);   % Store k_te as kappa for compatibility
            best.n_te = best_params(3);
            best.k_te = best_params(4);
            best.n_tm = best_params(5);
            best.k_tm = best_params(6);
            best.shift = best_params(7);
        else
            % TMM: [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
        best.ito_thickness = best_params(1);
        best.pedot_thickness = best_params(2);
        best.n_real = best_params(3);
        best.kappa = best_params(4);
        best.shift = best_params(5);
        end
        
        best.err = best_error;
        best.num_valid_points = length(exp_reflectivities);
        
        % Additional optimization info
        best.optimization_method = sprintf('global_search_%s', lower(calc_method));
        best.calculation_method = calc_method;
        best.optimization_type = sprintf('Global Search (Multiple Starting Points, %s)', calc_method);
        best.exitflag = 1; % Success
        best.iterations = best_output.iterations;
        best.function_count = best_output.funcCount;
        best.num_starting_points = num_starts;
        
        if isvalid(dlg)
            dlg.Message = sprintf('%s Global optimization complete (%d starting points, %d iterations)', calc_method, num_starts, best_output.iterations);
        end
    end
    
    function best = optimize_mode_particleswarm(pol, expData, dlg, calc_method)
        % Professional Particle Swarm Optimization using official MATLAB toolbox
        if nargin < 4
            calc_method = 'TMM';
        end
        
        if isvalid(dlg)
            dlg.Message = sprintf('Setting up Particle Swarm optimization (%s)...', calc_method);
        end
        
        wl = str2double(get(edWavelength, 'String'));
        V = str2double(get(edVoltage, 'String'));
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Parameter bounds depend on calculation method and polarization
        if strcmp(calc_method, 'Berreman')
            if pol == 's'
                % TE (s-pol): [ITO_thickness, PEDOT_thickness, n_te, k_te, angle_shift]
                lb = [10, 40, 1.0, 0.00, -10];    % Lower bounds
                ub = [30, 80, 1.5, 0.25, 10];     % Upper bounds
                nvars = 5;
            else
                % TM (p-pol): [ITO_thickness, PEDOT_thickness, n_te, k_te, n_tm, k_tm, angle_shift]
            lb = [10, 40, 1.0, 0.00, 1.0, 0.00, -10];    % Lower bounds
            ub = [30, 80, 1.5, 0.25, 1.5, 0.25, 10];     % Upper bounds
                nvars = 7;
        end
            objective_fun = @(params) sum(objective_function_berreman(params, pol, exp_angles, exp_reflectivities, ...
                n_bk, n_ito, n_water, degree, wl).^2);
        else
            % For TMM: [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
            lb = [10, 40, 1.30, 0.01, -10];    % Lower bounds
            ub = [30, 80, 1.5, 0.2, 10];     % Upper bounds
            nvars = 5;
            objective_fun = @(params) sum(objective_function(params, pol, exp_angles, exp_reflectivities, ...
                n_bk, n_ito, n_water, degree, wl).^2);
        end
        
        % Professional Particle Swarm options
        options = optimoptions('particleswarm', ...
            'Display', 'off', ...
            'SwarmSize', 500, ...          % Number of particles (maximum exploration)
            'MaxIterations', 1000, ...      % Maximum iterations
            'MaxStallIterations', 50, ...  % Stop if no improvement
            'FunctionTolerance', 1e-8, ... % Stop if function value changes less than this
            'UseParallel', false);         % Set to true if you have parallel toolbox
        
        if isvalid(dlg)
            dlg.Message = 'Running Particle Swarm optimization...';
        end
        
        % Test if particleswarm is available
        try
            test_fun = @(x) x^2;
            particleswarm(test_fun, 1, 0, 1, optimoptions('particleswarm', 'Display', 'off', 'MaxIterations', 1));
            fprintf('Particleswarm test successful\n');
        catch test_error
            fprintf('Particleswarm test failed: %s\n', test_error.message);
            error('Particleswarm not available: %s', test_error.message);
        end
        
        try
            % Run custom Particle Swarm Optimization that tracks all local minima
            nvars = length(lb);
            [x_opt, fval, all_local_minima] = custom_particleswarm_with_tracking(objective_fun, nvars, lb, ub, options);
            
            % Extract results based on calculation method and polarization
            if strcmp(calc_method, 'Berreman')
                best.ito_thickness = x_opt(1);
                best.pedot_thickness = x_opt(2);
                
                if pol == 's'
                    % TE (s-pol): [ITO, PEDOT, n_te, k_te, angle_shift]
                    best.n_te = x_opt(3);
                    best.k_te = x_opt(4);
                    best.shift = x_opt(5);
                    best.n_tm = 1.269;  % Arbitrary placeholder (not used for TE)
                    best.k_tm = 0.05;   % Arbitrary placeholder (not used for TE)
                best.n_real = x_opt(3);  % Store n_te as n_real for compatibility
                best.kappa = x_opt(4);   % Store k_te as kappa for compatibility
                else
                    % TM (p-pol): [ITO, PEDOT, n_te, k_te, n_tm, k_tm, angle_shift]
                best.n_te = x_opt(3);
                best.k_te = x_opt(4);
                best.n_tm = x_opt(5);
                best.k_tm = x_opt(6);
                best.shift = x_opt(7);
                    best.n_real = x_opt(3);  % Store n_te as n_real for compatibility
                    best.kappa = x_opt(4);   % Store k_te as kappa for compatibility
                end
            else
                % TMM: [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
                best.ito_thickness = x_opt(1);
                best.pedot_thickness = x_opt(2);
                best.n_real = x_opt(3);
                best.kappa = x_opt(4);
                best.shift = x_opt(5);
            end
            
            best.err = fval;
            best.num_valid_points = length(exp_reflectivities);
            
            % Calculate percentage errors after optimization
            if isvalid(dlg)
                dlg.Message = sprintf('Calculating percentage errors...');
            end
            sim_reflectivities = recalculate_reflectivity_with_params(best, pol, exp_angles, calc_method, n_bk, n_ito, n_water, degree, wl);
            error_stats = calculate_percentage_errors(sim_reflectivities, exp_reflectivities);
            best.percent_error_stats = error_stats;
            
            % Additional optimization info
            best.optimization_method = sprintf('particleswarm_%s', lower(calc_method));
            best.calculation_method = calc_method;
            best.exitflag = 1; % Success
            best.iterations = options.MaxIterations;
            best.function_count = length(all_local_minima);
            best.all_local_minima = all_local_minima; % Store all local minima
            
            if isvalid(dlg)
                dlg.Message = sprintf('%s Optimization complete (%d local minima found)', calc_method, length(all_local_minima));
            end
            
            % Return the result
            return;
            
        catch ME
            warning('Particle Swarm optimization failed: %s', ME.message);
            fprintf('Error details: %s\n', ME.getReport);
            if isvalid(dlg)
                dlg.Message = sprintf('Particle Swarm failed: %s', ME.message);
            end
            % Don't automatically fall back to grid search - let user decide
            error('Particle Swarm optimization failed: %s', ME.message);
        end
    end
    
    function best = optimize_mode_particleswarm_fixed_thickness(pol, expData, dlg, calc_method)
        % Particle Swarm Optimization with fixed thickness (only optimize n, k, angle_shift)
        if nargin < 4
            calc_method = 'TMM'; % Default to TMM if not specified
        end
        
        if isvalid(dlg)
            dlg.Message = sprintf('Setting up Particle Swarm optimization (Fixed Thickness, %s)...', calc_method);
        end
        
        wl = str2double(get(edWavelength, 'String'));
        V = str2double(get(edVoltage, 'String'));
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Get current thickness values from UI (these will be FIXED)
        current_ito = get_ui_value(edITO, app.Defaults.ito);
        current_pedot = get_ui_value(edPEDOT, app.Defaults.pedot);
        
        % Parameter bounds and objective function depend on calculation method
        if strcmp(calc_method, 'Berreman')
            % Berreman fixed thickness - parameters depend on polarization
            if pol == 's'
                % TE (s-pol): Only n_te, k_te matter [n_te, k_te, angle_shift]
                lb = [1.0, 0.00, -10];    % Lower bounds
                ub = [1.5, 0.25, 10];     % Upper bounds
                nvars = 3;  % n_te, k_te, angle_shift
            else
                % TM (p-pol): Both in-plane and out-of-plane matter [n_te, k_te, n_tm, k_tm, angle_shift]
                lb = [1.0, 0.00, 1.0, 0.00, -10];    % Lower bounds
                ub = [1.5, 0.25, 1.5, 0.25, 10];     % Upper bounds
                nvars = 5;  % n_te, k_te, n_tm, k_tm, angle_shift
            end
            
            % Create objective function for Berreman fixed thickness
            objective_fun = @(params) sum(objective_function_berreman_fixed_thickness(params, current_ito, current_pedot, ...
                pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl).^2);
        else
            % TMM fixed thickness: [n_real, kappa, angle_shift]
        lb = [1.2, 0.1, -5];    % Lower bounds
        ub = [2.0, 0.5, 5];     % Upper bounds
        
            % Create objective function for TMM fixed thickness
        objective_fun = @(params) sum(objective_function_fixed_thickness(params, current_ito, current_pedot, ...
            pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl).^2);
            
            nvars = 3;  % n_real, kappa, angle_shift
        end
        
        % Professional Particle Swarm options
        options = optimoptions('particleswarm', ...
            'Display', 'off', ...
            'SwarmSize', 500, ...          % Number of particles (maximum exploration)
            'MaxIterations', 1000, ...      % Maximum iterations
            'MaxStallIterations', 50, ...  % Stop if no improvement
            'FunctionTolerance', 1e-8, ... % Stop if function value changes less than this
            'UseParallel', false);         % Set to true if you have parallel toolbox
        
        if isvalid(dlg)
            dlg.Message = 'Running Particle Swarm optimization (Fixed Thickness)...';
        end
        
        % Test if particleswarm is available
        try
            test_fun = @(x) x^2;
            particleswarm(test_fun, 1, 0, 1, optimoptions('particleswarm', 'Display', 'off', 'MaxIterations', 1));
            fprintf('Particleswarm test successful\n');
        catch test_error
            fprintf('Particleswarm test failed: %s\n', test_error.message);
            error('Particleswarm not available: %s', test_error.message);
        end
        
        try
            % Run custom Particle Swarm Optimization that tracks all local minima
            [x_opt, fval, all_local_minima] = custom_particleswarm_with_tracking(objective_fun, nvars, lb, ub, options);
            
            % Extract results (thickness values are fixed from UI)
            best.ito_thickness = current_ito;
            best.pedot_thickness = current_pedot;
            
            if strcmp(calc_method, 'Berreman')
                % Berreman - depends on polarization
                if pol == 's'
                    % TE (s-pol): [n_te, k_te, angle_shift]
                    best.n_te = x_opt(1);
                    best.k_te = x_opt(2);
                    best.shift = x_opt(3);
                    best.n_tm = 1.269;  % Arbitrary placeholder (not used for TE)
                    best.k_tm = 0.05;   % Arbitrary placeholder (not used for TE)
                    best.n_real = x_opt(1);  % Store n_te as n_real for compatibility
                    best.kappa = x_opt(2);   % Store k_te as kappa for compatibility
                else
                    % TM (p-pol): [n_te, k_te, n_tm, k_tm, angle_shift]
                    best.n_te = x_opt(1);
                    best.k_te = x_opt(2);
                    best.n_tm = x_opt(3);
                    best.k_tm = x_opt(4);
                    best.shift = x_opt(5);
                    best.n_real = x_opt(1);  % Store n_te as n_real for compatibility
                    best.kappa = x_opt(2);   % Store k_te as kappa for compatibility
                end
            else
                % TMM: [n_real, kappa, angle_shift]
            best.n_real = x_opt(1);
            best.kappa = x_opt(2);
            best.shift = x_opt(3);
            end
            
            best.err = fval;
            best.num_valid_points = length(exp_reflectivities);
            
            % Calculate percentage errors after optimization
            if isvalid(dlg)
                dlg.Message = sprintf('Calculating percentage errors...');
            end
            sim_reflectivities = recalculate_reflectivity_with_params(best, pol, exp_angles, calc_method, n_bk, n_ito, n_water, degree, wl);
            error_stats = calculate_percentage_errors(sim_reflectivities, exp_reflectivities);
            best.percent_error_stats = error_stats;
            
            % Additional optimization info
            best.optimization_method = sprintf('particleswarm_fixed_thickness_%s', lower(calc_method));
            best.calculation_method = calc_method;
            best.exitflag = 1; % Success
            best.iterations = options.MaxIterations;
            best.function_count = length(all_local_minima);
            best.all_local_minima = all_local_minima; % Store all local minima
            
            if isvalid(dlg)
                dlg.Message = sprintf('%s Fixed Thickness optimization complete (%d local minima found)', calc_method, length(all_local_minima));
            end
            
            % Return the result
            return;
            
        catch ME
            warning('Particle Swarm Fixed Thickness optimization failed: %s', ME.message);
            fprintf('Error details: %s\n', ME.getReport);
            if isvalid(dlg)
                dlg.Message = sprintf('Particle Swarm Fixed Thickness failed: %s', ME.message);
            end
            % Don't automatically fall back to grid search - let user decide
            error('Particle Swarm Fixed Thickness optimization failed: %s', ME.message);
        end
    end
    
    function error_stats = calculate_percentage_errors(sim_reflectivities, exp_reflectivities)
        % Calculate percentage error statistics
        % error_stats contains: mean_error, std_error, max_error, min_error, rmse_percent
        % Formula: (simulation - experimental) / experimental * 100%
        
        % Remove NaN values
        valid_mask = ~isnan(sim_reflectivities) & ~isnan(exp_reflectivities) & exp_reflectivities ~= 0;
        sim_valid = sim_reflectivities(valid_mask);
        exp_valid = exp_reflectivities(valid_mask);
        
        if isempty(sim_valid) || isempty(exp_valid)
            error_stats = struct();
            error_stats.mean_error_percent = NaN;
            error_stats.std_error_percent = NaN;
            error_stats.max_error_percent = NaN;
            error_stats.min_error_percent = NaN;
            error_stats.rmse_percent = NaN;
            error_stats.mae_percent = NaN;
            return;
        end
        
        % Calculate percentage errors: (sim - exp) / exp * 100
        percent_errors = ((sim_valid - exp_valid) ./ exp_valid) * 100;
        
        error_stats.mean_error_percent = mean(percent_errors);
        error_stats.std_error_percent = std(percent_errors);
        error_stats.max_error_percent = max(percent_errors);
        error_stats.min_error_percent = min(percent_errors);
        error_stats.rmse_percent = sqrt(mean(percent_errors.^2));  % Root Mean Square Error
        error_stats.mae_percent = mean(abs(percent_errors));  % Mean Absolute Error
        error_stats.num_valid_points = length(percent_errors);
    end
    
    function sim_reflectivities = recalculate_reflectivity_with_params(best, pol, exp_angles, calc_method, n_bk, n_ito, n_water, degree, wl)
        % Recalculate reflectivity using optimized parameters to get simulation values
        % This is used to calculate percentage errors
        
        if strcmp(calc_method, 'Berreman')
            % Build parameters for Berreman
            ito_thick = best.ito_thickness;
            pedot_thick = best.pedot_thickness;
            n_te = best.n_te;
            k_te = best.k_te;
            
            if pol == 's'
                % TE mode: only in-plane matters
                angle_shift = best.shift;
                n_tm = 1.269;  % Not used for TE
                k_tm = 0.05;
            else
                % TM mode: both in-plane and out-of-plane matter
                n_tm = best.n_tm;
                k_tm = best.k_tm;
                angle_shift = best.shift;
            end
            
            % Calculate permittivity tensors
            % Extract ITO refractive index components
            n_ito_real = 1.8529;
            k_ito = 0.00316;
            
            % Calculate ITO permittivity explicitly as: ε = (n² - k²) + 2nk*j
            eps_ito_real = n_ito_real^2 - k_ito^2;
            eps_ito_imag = 2 * n_ito_real * k_ito;
            eps_ito_value = eps_ito_real + eps_ito_imag * 1i;
            
            % Create ITO permittivity tensor (isotropic)
            eps_ito = eps_ito_value * eye(3);
            
            eps_te = (n_te^2 - k_te^2) + 2*n_te*k_te*1i;
            eps_tm = (n_tm^2 - k_tm^2) + 2*n_tm*k_tm*1i;
            eps_pedot = diag([eps_te, eps_te, eps_tm]);
            eps_tensors = {eps_ito, eps_pedot};
            
            n_list = [n_bk, n_ito, complex(n_te, k_te), n_water];
            d_list = [Inf, ito_thick, pedot_thick, Inf];
            
            % Apply angle transformation
            tmm_angles = -1 .* exp_angles + (273 + angle_shift);
            sim_reflectivities = nan(size(tmm_angles));
            
            for ii = 1:numel(tmm_angles)
                a = tmm_angles(ii);
                if a >= 0 && a <= 90
                    try
                        sim_reflectivities(ii) = berreman_4x4_full(pol, n_list, d_list, a*degree, wl, eps_tensors);
                    catch
                        sim_reflectivities(ii) = NaN;
                    end
                end
            end
        else
            % TMM method
            ito_thick = best.ito_thickness;
            pedot_thick = best.pedot_thickness;
            n_real = best.n_real;
            kappa = best.kappa;
            angle_shift = best.shift;
            
            n_pedot = n_real + 1i*kappa;
            n_list = [n_bk, n_ito, n_pedot, n_water];
            d_list = [Inf, ito_thick, pedot_thick, Inf];
            
            % Apply angle transformation
            tmm_angles = -1 .* exp_angles + (273 + angle_shift);
            sim_reflectivities = nan(size(tmm_angles));
            
            for ii = 1:numel(tmm_angles)
                a = tmm_angles(ii);
                if a >= 0 && a <= 90
                    try
                        data = coh_tmm(pol, n_list, d_list, a*degree, wl);
                        sim_reflectivities(ii) = data.R;
                    catch
                        sim_reflectivities(ii) = NaN;
                    end
                end
            end
        end
    end
    
    function [x_opt, fval, all_local_minima] = custom_particleswarm_with_tracking(objective_fun, nvars, lb, ub, options)
        % Custom Particle Swarm Optimization that tracks all local minima found
        % This function finds and stores all local minima, not just the global minimum
        
        % PSO parameters
        swarm_size = options.SwarmSize;
        max_iter = options.MaxIterations;
        max_stall = options.MaxStallIterations;
        func_tol = options.FunctionTolerance;
        
        % Initialize particles
        particles = zeros(swarm_size, nvars);
        velocities = zeros(swarm_size, nvars);
        personal_best = zeros(swarm_size, nvars);
        personal_best_fval = inf(swarm_size, 1);
        
        % Initialize particles randomly within bounds
        for i = 1:swarm_size
            particles(i, :) = lb + (ub - lb) .* rand(1, nvars);
            personal_best(i, :) = particles(i, :);
            personal_best_fval(i) = objective_fun(particles(i, :));
        end
        
        % Find global best
        [global_best_fval, global_best_idx] = min(personal_best_fval);
        global_best = personal_best(global_best_idx, :);
        
        % Storage for all local minima
        all_local_minima = [];
        stall_count = 0;
        
        % PSO main loop
        for iter = 1:max_iter
            % Update particles
            for i = 1:swarm_size
                % Update velocity
                w = 0.9 - 0.5 * iter / max_iter; % Inertia weight
                c1 = 2.0; % Cognitive parameter
                c2 = 2.0; % Social parameter
                
                r1 = rand(1, nvars);
                r2 = rand(1, nvars);
                
                velocities(i, :) = w * velocities(i, :) + ...
                    c1 * r1 .* (personal_best(i, :) - particles(i, :)) + ...
                    c2 * r2 .* (global_best - particles(i, :));
                
                % Update position
                particles(i, :) = particles(i, :) + velocities(i, :);
                
                % Apply bounds
                particles(i, :) = max(particles(i, :), lb);
                particles(i, :) = min(particles(i, :), ub);
                
                % Evaluate objective function
                fval = objective_fun(particles(i, :));
                
                % Update personal best
                if fval < personal_best_fval(i)
                    personal_best(i, :) = particles(i, :);
                    personal_best_fval(i) = fval;
                    
                    % Store this as a potential local minimum
                    local_min = struct();
                    local_min.parameters = particles(i, :);
                    local_min.fval = fval;
                    local_min.iteration = iter;
                    local_min.particle_id = i;
                    all_local_minima = [all_local_minima; local_min];
                end
                
                % Update global best
                if fval < global_best_fval
                    global_best = particles(i, :);
                    global_best_fval = fval;
                    stall_count = 0;
                else
                    stall_count = stall_count + 1;
                end
            end
            
            % Check convergence
            if stall_count >= max_stall
                break;
            end
        end
        
        % Return results
        x_opt = global_best;
        fval = global_best_fval;
        
        % Sort local minima by function value (best first)
        if ~isempty(all_local_minima)
            [~, sort_idx] = sort([all_local_minima.fval]);
            all_local_minima = all_local_minima(sort_idx);
        end
        
        fprintf('Custom PSO completed: %d iterations, %d local minima found\n', iter, length(all_local_minima));
        
        % Return the values
        return;
    end
    
    
    
    
    
    
    function best = optimize_mode_grid(pol, expData, dlg, calc_method)
        % Original grid search method (kept as fallback)
        if nargin < 4
            calc_method = 'TMM'; % Default to TMM if not specified
        end
        
        % Note: Grid search currently only supports TMM method
        if strcmp(calc_method, 'Berreman')
            warning('Grid search does not yet support Berreman method. Using TMM instead.');
            calc_method = 'TMM';
        end
        
        wl = str2double(get(edWavelength, 'String'));
        V = str2double(get(edVoltage, 'String'));
        
        % Thickness ranges for optimization
        % Reduced grid size to prevent memory issues (25^5 = 9.7M combinations)
        % Use 10 points instead of 25 to reduce to ~100k combinations
        ito_thickness_vals = linspace(5, 30, 10);      % ITO thickness range: 5-30 nm
        pedot_thickness_vals = linspace(40, 80, 10);    % PEDOT thickness range: 40-80 nm
        
        n_real_vals = linspace(1.0, 2, 10);
        kappa_vals   = linspace(0.01, 1, 10);
        angle_offsets = linspace(-5, 5, 10);

        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;

        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);

        best.err = Inf;
        best.n_real = NaN;
        best.kappa = NaN;
        best.shift = 0;
        best.ito_thickness = NaN;
        best.pedot_thickness = NaN;
        best.num_valid_points = 0;

        % Store all results for each shifted angle
        allResults = [];
        resultIdx = 0;

        totalIters = numel(ito_thickness_vals)*numel(pedot_thickness_vals)*numel(n_real_vals)*numel(kappa_vals)*numel(angle_offsets);
        iter = 0;
        
        % Show user the expected number of results
        if isvalid(dlg)
            dlg.Message = sprintf('Grid search: %d combinations will be generated', totalIters);
        end

        for ito_thick = ito_thickness_vals
            for pedot_thick = pedot_thickness_vals
                for nr = n_real_vals
                    for kp = kappa_vals
                        n_pedot = nr + 1i*kp;
                        n_list = [n_bk, n_ito, n_pedot, n_water];
                        d_list = [Inf, ito_thick, pedot_thick, Inf];
                        for off = angle_offsets
                            iter = iter + 1;
                            if mod(iter,1000)==0 && isvalid(dlg)
                                dlg.Message = sprintf('Running grid search... (%d/%d)', iter, totalIters);
                            end
                            tmm_angles = -1 .* exp_angles + (273 + off);
                            theo = nan(size(tmm_angles));
                            for ii = 1:numel(tmm_angles)
                                a = tmm_angles(ii);
                                if a >= 0 && a <= 90
                                    data = coh_tmm(pol, n_list, d_list, a*degree, wl);
                                    theo(ii) = data.R;
                                end
                            end
                            valid = ~isnan(theo);
                            err = sum((theo(valid) - exp_reflectivities(valid)).^2);
                            
                            % Store all results
                            resultIdx = resultIdx + 1;
                            allResults(resultIdx).ito_thickness = ito_thick;
                            allResults(resultIdx).pedot_thickness = pedot_thick;
                            allResults(resultIdx).n_real = nr;
                            allResults(resultIdx).kappa = kp;
                            allResults(resultIdx).shift = off;
                            allResults(resultIdx).error = err;
                            allResults(resultIdx).valid_points = nnz(valid);
                            
                            if err < best.err
                                best.err = err;
                                best.n_real = nr;
                                best.kappa = kp;
                                best.shift = off;
                                best.ito_thickness = ito_thick;
                                best.pedot_thickness = pedot_thick;
                                best.num_valid_points = nnz(valid);
                            end
                        end
                    end
                end
            end
        end
        
        % Store all results for export
        if pol == 's'
            app.AllOptimizationResults.TE = allResults;
        else
            app.AllOptimizationResults.TM = allResults;
        end
        
        best.optimization_method = sprintf('grid_search_%s', lower(calc_method));
        best.calculation_method = calc_method;
        best.total_combinations = totalIters;
    end
    
    function residual = objective_function(params, pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl)
        % Objective function for lsqnonlin (all parameters)
        % params = [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
        
        ito_thick = params(1);
        pedot_thick = params(2);
        n_real = params(3);
        kappa = params(4);
        angle_shift = params(5);
        
        % Calculate theoretical reflectivity
        n_pedot = n_real + 1i*kappa;
        n_list = [n_bk, n_ito, n_pedot, n_water];
        d_list = [Inf, ito_thick, pedot_thick, Inf];
        
        % Apply angle transformation
        tmm_angles = -1 .* exp_angles + (273 + angle_shift);
        theo = nan(size(tmm_angles));
        
        for ii = 1:numel(tmm_angles)
            a = tmm_angles(ii);
            if a >= 0 && a <= 90
                try
                    data = coh_tmm(pol, n_list, d_list, a*degree, wl);
                    theo(ii) = data.R;
                catch
                    theo(ii) = NaN;
                end
            end
        end
        
        % Calculate residual (difference between theory and experiment)
        valid = ~isnan(theo);
        if sum(valid) == 0
            residual = ones(size(exp_reflectivities)) * 1e6; % Large penalty for invalid calculations
        else
            residual = theo(valid) - exp_reflectivities(valid);
        end
    end
    
    function residual = objective_function_fixed_thickness(params, ito_thick, pedot_thick, pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl)
        % Objective function for lsqnonlin (fixed thickness mode)
        % params = [n_real, kappa, angle_shift]
        
        n_real = params(1);
        kappa = params(2);
        angle_shift = params(3);
        
        % Calculate theoretical reflectivity
        n_pedot = n_real + 1i*kappa;
        n_list = [n_bk, n_ito, n_pedot, n_water];
        d_list = [Inf, ito_thick, pedot_thick, Inf];
        
        % Apply angle transformation
        tmm_angles = -1 .* exp_angles + (273 + angle_shift);
        theo = nan(size(tmm_angles));
        
        for ii = 1:numel(tmm_angles)
            a = tmm_angles(ii);
            if a >= 0 && a <= 90
                try
                    data = coh_tmm(pol, n_list, d_list, a*degree, wl);
                    theo(ii) = data.R;
                catch
                    theo(ii) = NaN;
                end
            end
        end
        
        % Calculate residual (difference between theory and experiment)
        valid = ~isnan(theo);
        if sum(valid) == 0
            residual = ones(size(exp_reflectivities)) * 1e6; % Large penalty for invalid calculations
        else
            residual = theo(valid) - exp_reflectivities(valid);
        end
    end
    
    function residual = objective_function_berreman(params, pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl)
        % Objective function for Berreman 4x4 method (all parameters including thickness)
        % TE (s-pol): params = [ITO_thickness, PEDOT_thickness, n_te, k_te, angle_shift] - 5 parameters
        % TM (p-pol): params = [ITO_thickness, PEDOT_thickness, n_te, k_te, n_tm, k_tm, angle_shift] - 7 parameters
        
        ito_thick = params(1);
        pedot_thick = params(2);
        
        if pol == 's'
            % TE (s-polarization): only in-plane indices matter
            n_te = params(3);
        k_te = params(4);
            angle_shift = params(5);
            % Set arbitrary values for out-of-plane (not used for TE)
            n_tm = 1.269;
            k_tm = 0.05;
        else
            % TM (p-polarization): both in-plane and out-of-plane matter
            n_te = params(3);
            k_te = params(4);
            n_tm = params(5);
        k_tm = params(6);
        angle_shift = params(7);
        end
        
        % Calculate permittivity tensors
        % Extract ITO refractive index components
        n_ito_real = 1.8529;
        k_ito = 0.00316;
        
        % Calculate ITO permittivity explicitly as: ε = (n² - k²) + 2nk*j
        eps_ito_real = n_ito_real^2 - k_ito^2;
        eps_ito_imag = 2 * n_ito_real * k_ito;
        eps_ito_value = eps_ito_real + eps_ito_imag * 1i;
        
        % Create ITO permittivity tensor (isotropic)
        eps_ito = eps_ito_value * eye(3);
        
        % PEDOT - anisotropic (uniaxial)
        eps_te = (n_te^2 - k_te^2) + 2*n_te*k_te*1i;  % In-plane (xy)
        eps_tm = (n_tm^2 - k_tm^2) + 2*n_tm*k_tm*1i;  % Out-of-plane (z)
        eps_pedot = diag([eps_te, eps_te, eps_tm]);
        
        eps_tensors = {eps_ito, eps_pedot};
        
        % Build layer lists
        n_list = [n_bk, n_ito, complex(n_te, k_te), n_water];
        d_list = [Inf, ito_thick, pedot_thick, Inf];
        
        % Apply angle transformation
        tmm_angles = -1 .* exp_angles + (273 + angle_shift);
        theo = nan(size(tmm_angles));
        
        for ii = 1:numel(tmm_angles)
            a = tmm_angles(ii);
            if a >= 0 && a <= 90
                try
                    % Use FULL Berreman 4x4 method (proper eigenvalue decomposition)
                    theo(ii) = berreman_4x4_full(pol, n_list, d_list, a*degree, wl, eps_tensors);
                catch
                    theo(ii) = NaN;
                end
            end
        end
        
        % Calculate residual (difference between theory and experiment)
        valid = ~isnan(theo);
        if sum(valid) == 0
            residual = ones(size(exp_reflectivities)) * 1e6; % Large penalty for invalid calculations
        else
            residual = theo(valid) - exp_reflectivities(valid);
        end
    end
    
    function residual = objective_function_berreman_fixed_thickness(params, ito_thick, pedot_thick, pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl)
        % Objective function for Berreman 4x4 method (fixed thickness mode)
        % TE (s-pol): params = [n_te, k_te, angle_shift] - 3 parameters
        % TM (p-pol): params = [n_te, k_te, n_tm, k_tm, angle_shift] - 5 parameters
        
        if pol == 's'
            % TE (s-polarization): only in-plane indices matter
            n_te = params(1);
            k_te = params(2);
            angle_shift = params(3);
            % Set arbitrary values for out-of-plane (not used for TE)
            n_tm = 1.269;
            k_tm = 0.05;
        else
            % TM (p-polarization): both in-plane and out-of-plane matter
            n_te = params(1);
            k_te = params(2);
            n_tm = params(3);
            k_tm = params(4);
            angle_shift = params(5);
        end
        
        % Calculate permittivity tensors
        % Extract ITO refractive index components
        n_ito_real = 1.8529;
        k_ito = 0.00316;
        
        % Calculate ITO permittivity explicitly as: ε = (n² - k²) + 2nk*j
        eps_ito_real = n_ito_real^2 - k_ito^2;
        eps_ito_imag = 2 * n_ito_real * k_ito;
        eps_ito_value = eps_ito_real + eps_ito_imag * 1i;
        
        % Create ITO permittivity tensor (isotropic)
        eps_ito = eps_ito_value * eye(3);
        
        % PEDOT - anisotropic (uniaxial)
        eps_te = (n_te^2 - k_te^2) + 2*n_te*k_te*1i;  % In-plane (xy)
        eps_tm = (n_tm^2 - k_tm^2) + 2*n_tm*k_tm*1i;  % Out-of-plane (z)
        eps_pedot = diag([eps_te, eps_te, eps_tm]);
        
        eps_tensors = {eps_ito, eps_pedot};
        
        % Build layer lists
        n_list = [n_bk, n_ito, complex(n_te, k_te), n_water];
        d_list = [Inf, ito_thick, pedot_thick, Inf];
        
        % Apply angle transformation
        tmm_angles = -1 .* exp_angles + (273 + angle_shift);
        theo = nan(size(tmm_angles));
        
        for ii = 1:numel(tmm_angles)
            a = tmm_angles(ii);
            if a >= 0 && a <= 90
                try
                    % Use FULL Berreman 4x4 method (proper eigenvalue decomposition)
                    theo(ii) = berreman_4x4_full(pol, n_list, d_list, a*degree, wl, eps_tensors);
                catch
                    theo(ii) = NaN;
                end
            end
        end
        
        % Calculate residual (difference between theory and experiment)
        valid = ~isnan(theo);
        if sum(valid) == 0
            residual = ones(size(exp_reflectivities)) * 1e6; % Large penalty for invalid calculations
        else
            residual = theo(valid) - exp_reflectivities(valid);
        end
    end

    function [expData, ok] = load_experimental_data(pol)
        file = 'eploy_22s_Spol.xlsx'; if pol=='p', file='eploy_22s_Ppol.xlsx'; end
        if ~isfile(file)
            expData = struct(); ok=false; return;
        end
        try
            T = readtable(file);
            if width(T) < 2 || height(T) < 2
                expData = struct(); ok=false; return;
            end
            angles_orig = T{:,1};
            Rvals = T{:,2};
            mapped = -1 .* angles_orig + 273;
            [mapped_sorted, idx] = sort(mapped);
            expData.angles = mapped_sorted;
            expData.reflectivities = Rvals(idx);
            expData.original_angles = angles_orig(idx);
            ok = true;
        catch
            expData = struct(); ok=false;
        end
    end

    %====================
    % Calculation and Plotting Functions - DISABLED (graph removed)
    %====================
    % on_calculate function commented out - graph has been removed per user request
    % All plotting functions below are no longer called from the UI

    %====================
    % Computational Functions (same as Python)
    %====================
    function [angles, R] = calculate_reflectivity_vs_angle(wavelength, voltage, n_real, layer_thicknesses, polarization, delta_k, manual_k)
        degree = pi/180;
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;

        if ~isnan(manual_k)
            k_pedot = manual_k + delta_k;
        else
            k_pedot = get_k_improved(app.PEDOT, wavelength, voltage) + delta_k;
        end
        n_pedot = n_real + 1i*k_pedot;

        n_list = [n_bk, n_ito, n_pedot, n_water];
        d_list = [Inf, layer_thicknesses.ito, layer_thicknesses.pedot, Inf];

        angles = linspace(20, 89.9, 500);
        R = zeros(size(angles));
        for idx = 1:numel(angles)
            data = coh_tmm(polarization, n_list, d_list, angles(idx)*degree, wavelength);
            R(idx) = data.R;
        end
    end

    function [angles, R_base, R_mod, dR, dR_over_R] = calculate_delta_reflectivity(wavelength, voltage, n_real, layer_thicknesses, polarization, delta_k, manual_k)
        [angles, R_base] = calculate_reflectivity_vs_angle(wavelength, voltage, n_real, layer_thicknesses, polarization, 0, manual_k);
        [~, R_mod] = calculate_reflectivity_vs_angle(wavelength, voltage, n_real, layer_thicknesses, polarization, delta_k, manual_k);
        dR = R_mod - R_base;
        eps = 1e-10;
        dR_over_R = dR ./ (R_base + eps);
    end

    % All plotting functions commented out - graph has been removed per user request
    % These functions are no longer used:
    % - plot_single_graph
    % - plot_only_reflectivity  
    % - compose_title
    % - plot_experimental
    % - plot_combined_te_tm

end % <-- end of main function

