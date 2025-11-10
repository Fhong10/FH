function Epoly_35s()
% Epoly_35s - Enhanced reflectivity calculator with separate TE/TM optimization
%
% This version calculates TE and TM modes separately, allowing different angle shifts
% for each mode, using Epoly_35s_Spol.xlsx and Epoly_35s_Ppol.xlsx experimental data files.
%
% Features:
%   • Particle Swarm (All) - Optimize all parameters including thickness
%   • Particle Swarm (Fixed) - Only optimize n, k, and angle shift
%   • Global Search - Multiple lsqnonlin starting points
%   • Advanced (lsqnonlin) - Fast single optimization
%   • Separate optimization for TE and TM modes (different angle shifts allowed)
%   • Epoly_35s experimental data files (Epoly_35s_Spol.xlsx, Epoly_35s_Ppol.xlsx)
%
% Mirrors reflectivity4.py functionality with MATLAB-specific optimizations.
% -----------------------------------------------------------------------------

    %----------------------------
    % App State - Same as Python
    %----------------------------
    app = struct();
    app.Defaults.ito = 15;       % nm
    app.Defaults.pedot = 87;   % nm
    app.Defaults.wavelength = 561.0; % nm
    app.Defaults.voltage = 0.0;      % mV
    app.Defaults.n_real = 1.41;
    app.Defaults.k_manual = 0.05;
    app.Defaults.delta_k = -0.01;
    app.Defaults.n_real_te = 1.41;
    app.Defaults.n_real_tm = 1.41;
    app.Defaults.k_te = 0.05;
    app.Defaults.k_tm = 0.05;
    app.Defaults.shift_te = 0.0;
    app.Defaults.shift_tm = 0.0;
    app.OptimizedShift = 0.0; % global (single-mode) optimized shift
    app.LastOptimizationResults = struct(); % Store last optimization results for export
    app.AllOptimizationResults = struct(); % Store all optimization results for detailed export

    % Default Excel filename as requested
    app.DefaultExcel = '300rpm extinction coefficient.xlsx';
    app.PEDOT = load_pedot_data_improved(app.DefaultExcel); % improved loading

    %----------------------------
    % UI Setup with PROPER Scrolling - SCROLLABLE LEFT PANEL ONLY
    %----------------------------
    % Create a fixed-size window that fits your screen
    fig = uifigure('Name','Epoly_35s Reflectivity Calculator - Separate TE/TM Optimization',...
        'Position',[100 100 1400 800]);

    % Main grid layout - 2 columns: scrollable controls + fixed graph
    gl = uigridlayout(fig,[1 2]);
    gl.ColumnWidth = {450, '1x'};
    gl.RowHeight = {'1x'};

    %--- LEFT PANEL: SCROLLABLE CONTAINER ---
    % Create a scrollable container for the left panel
    scrollContainer = uipanel(gl,'Title','Controls');
    scrollContainer.Layout.Row = 1; 
    scrollContainer.Layout.Column = 1;
    
    % Create a scrollable panel inside the container
    scrollPanel = uipanel(scrollContainer);
    scrollPanel.Position = [0 0 430 1000]; % Make it tall to force scrolling
    scrollPanel.Scrollable = 'on';
    
    % Grid layout inside the scrollable panel
    leftgl = uigridlayout(scrollPanel,[30 1]); % 30 rows to make it tall
    leftgl.RowHeight = repmat({'fit'},1,30);
    leftgl.ColumnWidth = {'1x'};

    %--- Layer thickness panel (same as Python)
    pnlThick = uipanel(leftgl,'Title','Layer Thicknesses (nm)');
    pnlThickgl = uigridlayout(pnlThick,[4 3]); 
    pnlThickgl.ColumnWidth = {110, 90, '1x'};
    
    uilabel(pnlThickgl,'Text','ITO:','HorizontalAlignment','right');
    edITO = uieditfield(pnlThickgl,'numeric','Value',app.Defaults.ito); 
    edITO.Limits = [0 Inf];
    uilabel(pnlThickgl,'Text','');
    
    uilabel(pnlThickgl,'Text','PEDOT:PSS:','HorizontalAlignment','right');
    edPEDOT = uieditfield(pnlThickgl,'numeric','Value',app.Defaults.pedot); 
    edPEDOT.Limits = [0 Inf];
    uilabel(pnlThickgl,'Text','');
    
    lblStruct = uilabel(pnlThickgl,'Text','Glass | ITO | PEDOT:PSS | Water','FontAngle','italic');
    lblStruct.Layout.Column = [1 3];
    
    btnResetThick = uibutton(pnlThickgl,'Text','Reset to Defaults',...
        'ButtonPushedFcn',@(src,evt) reset_thickness());
    btnResetThick.Layout.Column = [1 3];

    %--- Wavelength settings (same as Python)
    pnlWl = uipanel(leftgl,'Title','Wavelength Settings');
    pnlWlgl = uigridlayout(pnlWl,[3 3]); 
    pnlWlgl.ColumnWidth = {120, 120, '1x'};
    
    uilabel(pnlWlgl,'Text','Wavelength (nm):','HorizontalAlignment','right');
    edWavelength = uieditfield(pnlWlgl,'numeric','Value',app.Defaults.wavelength); 
    edWavelength.Limits=[1 Inf];
    edWavelength.ValueChangedFcn = @(s,e) update_ri_display();
    uilabel(pnlWlgl,'Text','');
    
    btnLoadExcel = uibutton(pnlWlgl,'Text','Load PEDOT Excel...',...
        'ButtonPushedFcn',@(src,evt) on_load_excel());
    btnLoadExcel.Layout.Column = [1 2];
    
    lblVoltInfo = uilabel(pnlWlgl,'Text',''); 
    lblVoltInfo.Layout.Column = [1 3];

    %--- Voltage settings (same as Python)
    pnlV = uipanel(leftgl,'Title','Voltage Settings');
    pnlVgl = uigridlayout(pnlV,[2 3]); 
    pnlVgl.ColumnWidth = {120, 120, '1x'};
    
    uilabel(pnlVgl,'Text','Voltage (mV):','HorizontalAlignment','right');
    edVoltage = uieditfield(pnlVgl,'numeric','Value',app.Defaults.voltage);
    edVoltage.ValueChangedFcn = @(s,e) update_ri_display();
    uilabel(pnlVgl,'Text','');
    
    lblAvailV = uilabel(pnlVgl,'Text',''); 
    lblAvailV.Layout.Column = [1 3];

    %--- PEDOT refractive index panel (same structure as Python)
    pnlRI = uipanel(leftgl,'Title','PEDOT:PSS Refractive Index');
    pnlRIgl = uigridlayout(pnlRI,[15 4]);
    pnlRIgl.ColumnWidth = {90, 70, 90, '1x'};
    
    % Manual input toggle
    chkManual = uicheckbox(pnlRIgl,'Text','Manual input','Value',false,'ValueChangedFcn',@(s,e) on_manual_toggle());
    chkManual.Layout.Column = [1 4];

    % Single mode section
    uilabel(pnlRIgl,'Text','Single Mode:','FontWeight','bold'); 
    uilabel(pnlRIgl,'Text',''); 
    uilabel(pnlRIgl,'Text',''); 
    uilabel(pnlRIgl,'Text','');
    
    uilabel(pnlRIgl,'Text','Real part (n):','HorizontalAlignment','right');
    edNreal = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.n_real); 
    edNreal.ValueChangedFcn=@(s,e) update_ri_display();
    
    uilabel(pnlRIgl,'Text','Imaginary part (k):');
    edK = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.k_manual,'Editable','off');
    edK.ValueChangedFcn = @(s,e) update_ri_display();
    
    lblRIDisplay = uilabel(pnlRIgl,'Text','','FontColor',[0 0 1]); 
    lblRIDisplay.Layout.Column = [1 4];

    % Separator
    uilabel(pnlRIgl,'Text','For Combined TE/TM Plot:','FontWeight','bold'); 
    uilabel(pnlRIgl,'Text',''); 
    uilabel(pnlRIgl,'Text',''); 
    uilabel(pnlRIgl,'Text','');
    
    % TE mode section
    uilabel(pnlRIgl,'Text','TE Mode:','FontAngle','italic'); 
    uilabel(pnlRIgl,'Text',''); 
    uilabel(pnlRIgl,'Text','TM Mode:','FontAngle','italic'); 
    uilabel(pnlRIgl,'Text','');
    
    uilabel(pnlRIgl,'Text','n (TE):','HorizontalAlignment','right');
    edNte = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.n_real_te);
    uilabel(pnlRIgl,'Text','n (TM):','HorizontalAlignment','right');
    edNtm = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.n_real_tm);
    
    uilabel(pnlRIgl,'Text','k (TE):');
    edKte = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.k_te,'Editable','off');
    uilabel(pnlRIgl,'Text','k (TM):');
    edKtm = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.k_tm,'Editable','off');
    
    uilabel(pnlRIgl,'Text','Shift (TE):');
    edShiftTE = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.shift_te);
    uilabel(pnlRIgl,'Text','Shift (TM):');
    edShiftTM = uieditfield(pnlRIgl,'numeric','Value',app.Defaults.shift_tm);

    % Optimization buttons (same as Python)
    btnOptTE = uibutton(pnlRIgl,'Text','Optimize TE Only', 'ButtonPushedFcn',@(s,e) on_optimize_te());
    btnOptTM = uibutton(pnlRIgl,'Text','Optimize TM Only', 'ButtonPushedFcn',@(s,e) on_optimize_tm());
    btnOptBoth = uibutton(pnlRIgl,'Text','Optimize Both TE & TM', 'ButtonPushedFcn',@(s,e) on_optimize_both());
    
    btnOptTE.Layout.Row = 13; btnOptTE.Layout.Column = [1 2];
    btnOptTM.Layout.Row = 13; btnOptTM.Layout.Column = [3 4];
    btnOptBoth.Layout.Row = 14; btnOptBoth.Layout.Column = [1 4];
    
    % Export button
    btnExport = uibutton(pnlRIgl,'Text','Export Optimization Results', 'ButtonPushedFcn',@(s,e) on_export_results());
    btnExport.Layout.Row = 15; btnExport.Layout.Column = [1 4];

    %--- Delta k settings (same as Python)
    pnlDk = uipanel(leftgl,'Title','Delta k Settings');
    pnlDkgl = uigridlayout(pnlDk,[3 2]); 
    pnlDkgl.ColumnWidth = {180, '1x'};
    
    chkUseDk = uicheckbox(pnlDkgl,'Text','Compare with modified k','Value',true,'ValueChangedFcn',@(s,e) on_use_dk_toggle());
    
    uilabel(pnlDkgl,'Text','Δk value:','HorizontalAlignment','right');
    edDk = uieditfield(pnlDkgl,'numeric','Value',app.Defaults.delta_k);
    
    chkShowExp = uicheckbox(pnlDkgl,'Text','Show experimental data (Epoly_35s_Spol.xlsx for TE, Epoly_35s_Ppol.xlsx for TM)','Value',true);
    chkShowExp.Layout.Column = [1 2];

    %--- Polarization settings (same as Python)
    pnlPol = uipanel(leftgl,'Title','Polarization');
    bgPol = uibuttongroup(pnlPol,'SelectionChangedFcn',@(s,e) update_ri_display());
    rbTE = uiradiobutton(bgPol,'Text','TE (s-polarization)','Position',[10 10 120 22]); 
    rbTE.Value=true;
    rbTM = uiradiobutton(bgPol,'Text','TM (p-polarization)','Position',[140 10 120 22]);

    %--- Plot selection (same as Python)
    pnlPlot = uipanel(leftgl,'Title','Plot Selection');
    pnlPlotgl = uigridlayout(pnlPlot,[3 2]);
    
    ddPlot = uidropdown(pnlPlotgl,'Items',{'Reflectivity vs. Angle','ΔR vs. Angle','ΔR/R vs. Angle','ΔR/Δk vs. Angle','Combined TE/TM (all 4 curves)'},...
        'Value','Reflectivity vs. Angle');
    ddPlot.Layout.Column = [1 2];

    %--- Buttons (same as Python)
    pnlButtons = uipanel(leftgl,'Title','Actions');
    pnlButtonsgl = uigridlayout(pnlButtons,[1 2]); 
    pnlButtonsgl.ColumnWidth = {'1x','1x'};
    
    btnCalc = uibutton(pnlButtonsgl,'Text','Calculate and Plot','ButtonPushedFcn',@(s,e) on_calculate());
    btnClose = uibutton(pnlButtonsgl,'Text','Close','ButtonPushedFcn',@(s,e) close(fig));

    %--- Status (same as Python)
    pnlStatus = uipanel(leftgl,'Title','Status');
    lblStatus = uilabel(pnlStatus,'Text','Ready','FontColor',[0 0.5 0]);

    %--- RIGHT PANEL: FIXED GRAPH ---
    % The graph stays fixed and visible
    ax = uiaxes(gl); 
    ax.Layout.Row=1; 
    ax.Layout.Column=2;
    xlabel(ax,'Incident Angle (degrees)'); 
    ylabel(ax,'Reflectivity');
    grid(ax,'on');

    % Initialize display
    refresh_voltage_info();
    update_ri_display();

    %====================
    % Callbacks / helpers (same structure as Python)
    %====================
    function [oldStates] = disable_longrun_controls()
        ctrls = [btnCalc, btnOptTE, btnOptTM, btnOptBoth, btnLoadExcel];
        oldStates = arrayfun(@(b) b.Enable, ctrls, 'UniformOutput', false);
        set(ctrls, 'Enable','off');
    end
    
    function restore_controls(oldStates)
        ctrls = [btnCalc, btnOptTE, btnOptTM, btnOptBoth, btnLoadExcel];
        for i=1:numel(ctrls)
            ctrls(i).Enable = oldStates{i};
        end
    end
 
    function reset_thickness()
        edITO.Value = app.Defaults.ito;
        edPEDOT.Value = app.Defaults.pedot;
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
            lblStatus.Text = 'Loaded PEDOT Excel successfully.'; 
            lblStatus.FontColor=[0 0.5 0];
            refresh_voltage_info();
            update_ri_display();
        catch ME
            if exist('d','var') && isvalid(d), close(d); end
            if exist('oldStates','var'), restore_controls(oldStates); end
            uialert(fig,ME.message,'Error loading Excel');
            lblStatus.Text = 'Error loading Excel'; 
            lblStatus.FontColor=[0.85 0 0];
        end
    end

    function refresh_voltage_info()
        v = app.PEDOT.voltages;
        v = v(~isnan(v));
        if ~isempty(v)
            lblAvailV.Text = sprintf('Available: %.0f to %.0f mV',min(v),max(v));
        else
            lblAvailV.Text = 'No voltage data loaded';
        end
        if ~isempty(app.PEDOT.wavelengths)
            lblVoltInfo.Text = sprintf('Wavelengths: %.0f–%.0f nm  |  Columns: %d voltages',...
                min(app.PEDOT.wavelengths), max(app.PEDOT.wavelengths), numel(app.PEDOT.voltages));
        else
            lblVoltInfo.Text = '';
        end
    end

    function on_manual_toggle()
        isMan = chkManual.Value;
        edK.Editable = isMan;
        edKte.Editable = isMan;
        edKtm.Editable = isMan;
        update_ri_display();
    end

    function on_use_dk_toggle()
        use = chkUseDk.Value;
        if ~use
            if ~strcmp(ddPlot.Value,'Reflectivity vs. Angle') && ~strcmp(ddPlot.Value,'Combined TE/TM (all 4 curves)')
                ddPlot.Value = 'Reflectivity vs. Angle';
            end
        end
    end

    function update_ri_display(~,~)
        try
            wl = edWavelength.Value;
            V = edVoltage.Value;
            nreal = edNreal.Value;
            if chkManual.Value
                k = edK.Value;
                lblRIDisplay.Text = sprintf('Manual input: n = %.3f + %.6fi', nreal, k);
            else
                k = get_k_improved(app.PEDOT, wl, V);
                lblRIDisplay.Text = sprintf('At %.1f nm, %.1f mV: n = %.3f + %.6fi', wl, V, nreal, k);
            end
        catch
            lblRIDisplay.Text = '';
        end
    end

    %====================
    % Optimization Functions (same as Python)
    %====================
    function on_optimize_both()
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Optimizing','Indeterminate','on','Message','Optimizing TE...');
        try
            modeChoice = uiconfirm(fig, ['Choose how to optimize TE and TM:\n\n' ...
                '• Combined Particle Swarm (shared angle shift, minimize joint error)\n' ...
                '• Optimize Separately (choose method for TE and TM individually)'], ...
                'TE & TM Optimization', ...
                'Options', {'Combined Particle Swarm','Optimize Separately','Cancel'}, ...
                'DefaultOption', 1, 'CancelOption', 3);
            if strcmp(modeChoice,'Cancel')
                if isvalid(d), close(d); end
                restore_controls(oldStates);
                return;
            end
            
            if strcmp(modeChoice,'Combined Particle Swarm')
                lblStatus.Text='Optimizing TE & TM together...'; lblStatus.FontColor=[0 0 0.7]; drawnow;
                bestCombo = optimize_te_tm_combined_particleswarm(d);
                if isempty(bestCombo)
                    lblStatus.Text='Combined optimization cancelled or failed'; lblStatus.FontColor=[0.85 0 0];
                else
                    % Update UI with combined results
                    edITO.Value = bestCombo.ito_thickness;
                    edPEDOT.Value = bestCombo.pedot_thickness;
                    edNte.Value = bestCombo.n_real_te;
                    edKte.Value = bestCombo.kappa_te;
                    edNtm.Value = bestCombo.n_real_tm;
                    edKtm.Value = bestCombo.kappa_tm;
                    edShiftTE.Value = bestCombo.shift;
                    edShiftTM.Value = bestCombo.shift;
                    edNreal.Value = bestCombo.n_real_te;
                    edK.Value = bestCombo.kappa_te;
                    chkManual.Value = true; on_manual_toggle();
                    app.OptimizedShift = bestCombo.shift;
                    app.LastOptimizationResults.TE = bestCombo.TE;
                    app.LastOptimizationResults.TM = bestCombo.TM;
                    app.LastOptimizationResults.Combined = bestCombo;
                    
                    msg = sprintf(['Combined Particle Swarm complete:\n' ...
                        '  ITO thickness = %.2f nm\n' ...
                        '  PEDOT thickness = %.2f nm\n' ...
                        '  Angle shift = %.2f°\n\n' ...
                        '  TE:  n = %.4f + %.4fi (error = %.3e)\n' ...
                        '  TM:  n = %.4f + %.4fi (error = %.3e)\n\n' ...
                        '  Total error = %.3e'], ...
                        bestCombo.ito_thickness, bestCombo.pedot_thickness, bestCombo.shift, ...
                        bestCombo.n_real_te, bestCombo.kappa_te, bestCombo.err_te, ...
                        bestCombo.n_real_tm, bestCombo.kappa_tm, bestCombo.err_tm, ...
                        bestCombo.total_error);
                    uialert(fig, msg, 'Combined Optimization Complete');
                    lblStatus.Text = 'Combined optimization complete'; lblStatus.FontColor=[0 0.5 0];
                end
                if isvalid(d), close(d); end
                restore_controls(oldStates);
                return;
            end
            
            lblStatus.Text='Optimizing TE...'; lblStatus.FontColor=[0 0 0.7]; drawnow;
            [bestTE, hasTE] = optimize_mode('s', d, 'Optimizing TM...');
            if hasTE
                edNte.Value = bestTE.n_real;
                edKte.Value = bestTE.kappa;
                edShiftTE.Value = bestTE.shift;
                edITO.Value = bestTE.ito_thickness;
                edPEDOT.Value = bestTE.pedot_thickness;
                app.LastOptimizationResults.TE = bestTE;
            end
            d.Message = 'Optimizing TM...'; drawnow;
            [bestTM, hasTM] = optimize_mode('p', d, 'Finishing...');
            if hasTM
                edNtm.Value = bestTM.n_real;
                edKtm.Value = bestTM.kappa;
                edShiftTM.Value = bestTM.shift;
                edITO.Value = bestTM.ito_thickness;
                edPEDOT.Value = bestTM.pedot_thickness;
                app.LastOptimizationResults.TM = bestTM;
            end
            if hasTE || hasTM
                chkManual.Value = true; on_manual_toggle();
                uialert(fig, compose_results_message(bestTE,hasTE,bestTM,hasTM),'Optimization Complete');
                lblStatus.Text = 'Optimization complete'; lblStatus.FontColor=[0 0.5 0];
            else
                lblStatus.Text = 'No experimental data found'; lblStatus.FontColor=[0.85 0 0];
            end
        catch ME
            lblStatus.Text = ['Optimization error: ' ME.message]; lblStatus.FontColor=[0.85 0 0];
            uialert(fig,ME.message,'Optimization Error');
        end
        if isvalid(d), close(d); end
        restore_controls(oldStates);
    end

    function on_optimize_te()
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Optimizing','Indeterminate','on','Message','Optimizing TE...');
        try
            [best, ok] = optimize_mode('s', d, 'Finishing...');
            if ~ok
                uialert(fig,'TE experimental data (Epoly_35s_Spol.xlsx) not found or invalid.','No Data');
                lblStatus.Text='TE data not found'; lblStatus.FontColor=[0.85 0 0];
                return;
            end
            edNte.Value = best.n_real; edKte.Value = best.kappa; edShiftTE.Value = best.shift;
            edNreal.Value = best.n_real; edK.Value = best.kappa; chkManual.Value=true; on_manual_toggle();
            edITO.Value = best.ito_thickness; edPEDOT.Value = best.pedot_thickness;
            app.OptimizedShift = best.shift;
            app.LastOptimizationResults.TE = best;
            lblStatus.Text='TE optimization complete'; lblStatus.FontColor=[0 0.5 0];
            uialert(fig, format_best_message('TE',best),'TE Optimization Complete');
        catch ME
            lblStatus.Text = ['TE optimization error: ' ME.message]; lblStatus.FontColor=[0.85 0 0];
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
                uialert(fig,'TM experimental data (Epoly_35s_Ppol.xlsx) not found or invalid.','No Data');
                lblStatus.Text='TM data not found'; lblStatus.FontColor=[0.85 0 0];
                return;
            end
            edNtm.Value = best.n_real; edKtm.Value = best.kappa; edShiftTM.Value = best.shift;
            edNreal.Value = best.n_real; edK.Value = best.kappa; chkManual.Value=true; on_manual_toggle();
            edITO.Value = best.ito_thickness; edPEDOT.Value = best.pedot_thickness;
            app.OptimizedShift = best.shift;
            app.LastOptimizationResults.TM = best;
            lblStatus.Text='TM optimization complete'; lblStatus.FontColor=[0 0.5 0];
            uialert(fig, format_best_message('TM',best),'TM Optimization Complete');
        catch ME
            lblStatus.Text = ['TM optimization error: ' ME.message]; lblStatus.FontColor=[0.85 0 0];
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
            choice = uiconfirm(fig, ['What would you like to export?\n\n' ...
                '1. All Local Minima - Export all local minima found during optimization\n' ...
                '2. Best Results Only - Export only the best optimization results\n' ...
                '3. Cancel'], ...
                'Export Options', 'Options', {'All Local Minima', 'Best Results Only', 'Cancel'}, ...
                'DefaultOption', 1);
            
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
            uialert(fig,'No optimization results to export. Please run optimization first.','No Data');
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
            headers = {'Mode', 'ITO_Thickness_nm', 'PEDOT_Thickness_nm', 'n_real', 'k_imaginary', 'Angle_Shift_deg', 'Error', 'Valid_Points', 'Method'};
            allData = {};
            rowCount = 0;
            
            % Export TE result if available
            if isfield(app.LastOptimizationResults, 'TE')
                te = app.LastOptimizationResults.TE;
                rowCount = rowCount + 1;
                allData{rowCount,1} = 'TE';
                allData{rowCount,2} = te.ito_thickness;
                allData{rowCount,3} = te.pedot_thickness;
                allData{rowCount,4} = te.n_real;
                allData{rowCount,5} = te.kappa;
                allData{rowCount,6} = te.shift;
                allData{rowCount,7} = te.err;
                allData{rowCount,8} = te.num_valid_points;
                allData{rowCount,9} = te.optimization_method;
            end
            
            % Export TM result if available
            if isfield(app.LastOptimizationResults, 'TM')
                tm = app.LastOptimizationResults.TM;
                rowCount = rowCount + 1;
                allData{rowCount,1} = 'TM';
                allData{rowCount,2} = tm.ito_thickness;
                allData{rowCount,3} = tm.pedot_thickness;
                allData{rowCount,4} = tm.n_real;
                allData{rowCount,5} = tm.kappa;
                allData{rowCount,6} = tm.shift;
                allData{rowCount,7} = tm.err;
                allData{rowCount,8} = tm.num_valid_points;
                allData{rowCount,9} = tm.optimization_method;
            end
            
            % Create table and export
            T = array2table(allData, 'VariableNames', headers);
            writetable(T, fullpath);
            
            lblStatus.Text = sprintf('Results exported to %s', filename);
            lblStatus.FontColor = [0 0.5 0];
            
            uialert(fig, sprintf('Optimization results exported successfully!\nFile: %s\nTotal results: %d', fullpath, rowCount), 'Export Complete');
            
        catch ME
            lblStatus.Text = ['Export error: ' ME.message];
            lblStatus.FontColor = [0.85 0 0];
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
            
            lblStatus.Text = sprintf('All local minima exported to %s (%d solutions)', filename, rowCount);
            lblStatus.FontColor = [0 0.5 0];
            
            uialert(fig, sprintf('All local minima exported successfully!\nFile: %s\nTotal solutions: %d', fullpath, rowCount), 'Export Complete');
            
        catch ME
            lblStatus.Text = ['Export error: ' ME.message];
            lblStatus.FontColor = [0.85 0 0];
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
        if isfield(best, 'optimization_method')
            if contains(best.optimization_method, 'lsqnonlin')
                if contains(best.optimization_method, 'fixed_thickness')
                    opt_type = 'Fixed Thickness';
                    thickness_info = sprintf('(ITO=%.1f nm, PEDOT=%.1f nm - FIXED)', best.ito_thickness, best.pedot_thickness);
                else
                    opt_type = 'All Parameters';
                    thickness_info = sprintf('ITO=%.1f nm, PEDOT=%.1f nm', best.ito_thickness, best.pedot_thickness);
                end
                msg = sprintf('%s Mode (%s):\n  n = %.4f + %.4fi\n  Angle shift = %.2f°\n  %s\n  Error = %.3e\n  Iterations: %d', ...
                    mode, opt_type, best.n_real, best.kappa, best.shift, thickness_info, best.err, best.iterations);
            elseif strcmp(best.optimization_method, 'global_search')
                msg = sprintf('%s Mode (%s):\n  n = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Starting points: %d, Iterations: %d', ...
                    mode, 'Global Search', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.num_starting_points, best.iterations);
            elseif strcmp(best.optimization_method, 'particleswarm')
                msg = sprintf('%s Mode (%s):\n  n = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Iterations: %d', ...
                    mode, 'Particle Swarm', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.iterations);
            elseif strcmp(best.optimization_method, 'particleswarm_fixed_thickness')
                msg = sprintf('%s Mode (%s):\n  n = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm (FIXED)\n  PEDOT thickness = %.1f nm (FIXED)\n  Error = %.3e\n  Iterations: %d', ...
                    mode, 'Particle Swarm Fixed Thickness', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.iterations);
            else
                msg = sprintf('%s Mode (%s):\n  n = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e\n  Combinations tested: %d', ...
                    mode, 'Grid Search', best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err, best.total_combinations);
            end
        else
            % Fallback for old format
            msg = sprintf('%s Mode:\n  n = %.4f + %.4fi\n  Angle shift = %.2f°\n  ITO thickness = %.1f nm\n  PEDOT thickness = %.1f nm\n  Error = %.3e', ...
                mode, best.n_real, best.kappa, best.shift, best.ito_thickness, best.pedot_thickness, best.err);
        end
    end

    function [best, ok] = optimize_mode(pol, dlg, nextMessage)
        [expData, ok] = load_experimental_data(pol);
        if ~ok
            best = struct(); return;
        end
        
        % Choose optimization method (using uiconfirm with max 4 options)
        choice = uiconfirm(fig, sprintf(['Choose optimization method for %s mode:\n\n' ...
            '1. Particle Swarm (All) - Optimize all parameters including thickness\n' ...
            '2. Particle Swarm (Fixed) - Only optimize n, k, and angle shift\n' ...
            '3. Global Search - Multiple lsqnonlin starting points\n' ...
            '4. Advanced (lsqnonlin) - Fast single optimization'], pol), ...
            'Optimization Method', 'Options', {'Particle Swarm (All)', 'Particle Swarm (Fixed)', 'Global Search', 'Advanced (All)'}, ...
            'DefaultOption', 1);
        
        if strcmp(choice, 'Particle Swarm (All)')
            useCombinedAngle = ask_use_combined_angle();
            best = optimize_mode_particleswarm(pol, expData, dlg, useCombinedAngle);
            ok = true;
        elseif strcmp(choice, 'Particle Swarm (Fixed)')
            useCombinedAngle = ask_use_combined_angle();
            best = optimize_mode_particleswarm_fixed_thickness(pol, expData, dlg, useCombinedAngle);
            ok = true;
        elseif strcmp(choice, 'Global Search')
            best = optimize_mode_global(pol, expData, dlg);
            ok = true;
        elseif strcmp(choice, 'Advanced (All)')
            best = optimize_mode_lsqnonlin(pol, expData, dlg, 'all');
            ok = true;
        else % Grid Search (fallback)
            best = optimize_mode_grid(pol, expData, dlg);
            ok = true;
        end
        
        if isvalid(dlg), dlg.Message = nextMessage; end
    end
    
    function useCombinedAngle = ask_use_combined_angle()
        shift_preview = app.OptimizedShift;
        if isempty(shift_preview) || ~isfinite(shift_preview)
            shift_preview = 0;
        end
        prompt = sprintf(['Would you like to fix the angle shift to the current combined value?\n' ...
                          '(Current combined angle shift: %.2f°)\n\n' ...
                          'Choose "Use combined angle" to hold the shift constant,\n' ...
                          'or "Optimize angle shift" to let the optimizer fit it.'], shift_preview);
        choice = uiconfirm(fig, prompt, 'Angle Shift Option', ...
            'Options', {'Optimize angle shift','Use combined angle'}, ...
            'DefaultOption', 1, 'CancelOption', 1);
        useCombinedAngle = strcmp(choice, 'Use combined angle');
    end
    
    function best = optimize_mode_lsqnonlin(pol, expData, dlg, mode)
        % Advanced optimization using lsqnonlin
        if nargin < 4
            mode = 'all'; % Default to optimizing all parameters
        end
        
        if isvalid(dlg)
            if strcmp(mode, 'fixed_thickness')
                dlg.Message = 'Setting up fixed-thickness optimization...';
            else
                dlg.Message = 'Setting up advanced optimization...';
            end
        end
        
        wl = edWavelength.Value;
        V = edVoltage.Value;
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Get current thickness values from UI
        current_ito = edITO.Value;
        current_pedot = edPEDOT.Value;
        
        if strcmp(mode, 'fixed_thickness')
            % Fixed thickness mode: only optimize [n_real, kappa, angle_shift]
            lb = [1.3, 0.00, -3];    % Lower bounds
            ub = [1.5, 0.15, 3];     % Upper bounds
            x0 = [1.4, 0.5, 0];     % Initial guess
            
            % Create objective function for fixed thickness
            objective_fun = @(params) objective_function_fixed_thickness(params, current_ito, current_pedot, ...
                pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl);
            
            optimization_type = 'Fixed Thickness';
        else
            % All parameters mode: optimize [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
            lb = [5, 230, 1.0, 0.01, -5];    % Lower bounds
            ub = [30, 250, 2.0, 1.0, 5];     % Upper bounds
            x0 = [current_ito, current_pedot, 1.5, 0.5, 0]; % Use current thicknesses as starting point
            
            % Create objective function for all parameters
            objective_fun = @(params) objective_function(params, pol, exp_angles, exp_reflectivities, ...
                n_bk, n_ito, n_water, degree, wl);
            
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
            
            if strcmp(mode, 'fixed_thickness')
                % Extract results for fixed thickness mode
                best.ito_thickness = current_ito;
                best.pedot_thickness = current_pedot;
                best.n_real = x_opt(1);
                best.kappa = x_opt(2);
                best.shift = x_opt(3);
            else
                % Extract results for all parameters mode
                best.ito_thickness = x_opt(1);
                best.pedot_thickness = x_opt(2);
                best.n_real = x_opt(3);
                best.kappa = x_opt(4);
                best.shift = x_opt(5);
            end
            
            best.err = resnorm;
            best.num_valid_points = length(exp_reflectivities);
            
            % Additional optimization info
            best.optimization_method = sprintf('lsqnonlin_%s', mode);
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
            best = optimize_mode_grid(pol, expData, dlg);
        end
    end
    
    function best = optimize_mode_global(pol, expData, dlg)
        % Global optimization using multiple starting points to avoid local minima
        if isvalid(dlg)
            dlg.Message = 'Setting up global optimization...';
        end
        
        wl = edWavelength.Value;
        V = edVoltage.Value;
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Get current thickness values from UI
        current_ito = edITO.Value;
        current_pedot = edPEDOT.Value;
        
        % Parameter bounds [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
        lb = [5, 40, 1.0, 0.01, -5];    % Lower bounds
        ub = [30, 80, 2.0, 1.0, 5];     % Upper bounds
        
        % Create objective function
        objective_fun = @(params) objective_function(params, pol, exp_angles, exp_reflectivities, ...
            n_bk, n_ito, n_water, degree, wl);
        
        % Multiple starting points to explore different regions
        num_starts = 8;
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
        
        % Extract results
        best.ito_thickness = best_params(1);
        best.pedot_thickness = best_params(2);
        best.n_real = best_params(3);
        best.kappa = best_params(4);
        best.shift = best_params(5);
        best.err = best_error;
        best.num_valid_points = length(exp_reflectivities);
        
        % Additional optimization info
        best.optimization_method = 'global_search';
        best.optimization_type = 'Global Search (Multiple Starting Points)';
        best.exitflag = 1; % Success
        best.iterations = best_output.iterations;
        best.function_count = best_output.funcCount;
        best.num_starting_points = num_starts;
        
        if isvalid(dlg)
            dlg.Message = sprintf('Global optimization complete (%d starting points, %d iterations)', num_starts, best_output.iterations);
        end
    end
    
    function best = optimize_mode_particleswarm(pol, expData, dlg, useCombinedAngle)
        % Professional Particle Swarm Optimization using official MATLAB toolbox
        if isvalid(dlg)
            dlg.Message = 'Setting up Particle Swarm optimization...';
        end
        
        wl = edWavelength.Value;
        V = edVoltage.Value;
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Parameter bounds [ITO_thickness, PEDOT_thickness, n_real, kappa, angle_shift]
        lb = [10, 40, 1.30, 0.01, -10];    % Lower bounds
        ub = [30, 80, 1.5, 0.2, 10];       % Upper bounds
        combinedShift = app.OptimizedShift;
        if isempty(combinedShift) || ~isfinite(combinedShift)
            combinedShift = 0;
        end
        if useCombinedAngle
            lb(5) = combinedShift;
            ub(5) = combinedShift;
        end
        show_boundary_info('Particle Swarm (All parameters)', ...
            {'ITO thickness (nm)','PEDOT thickness (nm)','n_{real}','\kappa','Angle shift (deg)'}, lb, ub);
        
        % Create objective function (sum of squared errors for particleswarm)
        objective_fun = @(params) sum(objective_function(params, pol, exp_angles, exp_reflectivities, ...
            n_bk, n_ito, n_water, degree, wl).^2);
        
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
            [x_opt, fval, all_local_minima] = custom_particleswarm_with_tracking(objective_fun, 5, lb, ub, options);
            
            % Extract results
            best.ito_thickness = x_opt(1);
            best.pedot_thickness = x_opt(2);
            best.n_real = x_opt(3);
            best.kappa = x_opt(4);
            best.shift = x_opt(5);
            best.err = fval;
            best.num_valid_points = length(exp_reflectivities);
            
            % Additional optimization info
            best.optimization_method = 'particleswarm';
            best.exitflag = 1; % Success
            best.iterations = options.MaxIterations;
            best.function_count = length(all_local_minima);
            best.all_local_minima = all_local_minima; % Store all local minima
            
            if isvalid(dlg)
                dlg.Message = sprintf('Optimization complete (%d local minima found)', length(all_local_minima));
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
    
    function best = optimize_mode_particleswarm_fixed_thickness(pol, expData, dlg, useCombinedAngle)
        % Particle Swarm Optimization with fixed thickness (only optimize n, k, angle_shift)
        if isvalid(dlg)
            dlg.Message = 'Setting up Particle Swarm optimization (Fixed Thickness)...';
        end
        
        wl = edWavelength.Value;
        V = edVoltage.Value;
        
        % Physical constants
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        exp_angles = expData.original_angles;
        exp_reflectivities = expData.reflectivities(:);
        
        % Get current thickness values from UI (these will be FIXED)
        current_ito = edITO.Value;
        current_pedot = edPEDOT.Value;
        
        % Parameter bounds for fixed thickness mode [n_real, kappa, angle_shift]
        lb = [1.2, 0.1, -5];    % Lower bounds
        ub = [2.0, 0.5, 5];     % Upper bounds
        combinedShift = app.OptimizedShift;
        if isempty(combinedShift) || ~isfinite(combinedShift)
            combinedShift = 0;
        end
        if useCombinedAngle
            lb(3) = combinedShift;
            ub(3) = combinedShift;
        end
        show_boundary_info('Particle Swarm (Fixed thickness)', ...
            {'n_{real}','\kappa','Angle shift (deg)'}, lb, ub);
        
        % Create objective function for fixed thickness (sum of squared errors for particleswarm)
        objective_fun = @(params) sum(objective_function_fixed_thickness(params, current_ito, current_pedot, ...
            pol, exp_angles, exp_reflectivities, n_bk, n_ito, n_water, degree, wl).^2);
        
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
            % Run custom Particle Swarm Optimization that tracks all local minima (3 parameters only)
            [x_opt, fval, all_local_minima] = custom_particleswarm_with_tracking(objective_fun, 3, lb, ub, options);
            
            % Extract results (thickness values are fixed from UI)
            best.ito_thickness = current_ito;
            best.pedot_thickness = current_pedot;
            best.n_real = x_opt(1);
            best.kappa = x_opt(2);
            best.shift = x_opt(3);
            best.err = fval;
            best.num_valid_points = length(exp_reflectivities);
            
            % Additional optimization info
            best.optimization_method = 'particleswarm_fixed_thickness';
            best.exitflag = 1; % Success
            best.iterations = options.MaxIterations;
            best.function_count = length(all_local_minima);
            best.all_local_minima = all_local_minima; % Store all local minima
            
            if isvalid(dlg)
                dlg.Message = sprintf('Fixed Thickness optimization complete (%d local minima found)', length(all_local_minima));
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
    
    function combo = optimize_te_tm_combined_particleswarm(dlg)
        combo = [];
        [expTE, hasTE] = load_experimental_data('s');
        [expTM, hasTM] = load_experimental_data('p');
        if ~hasTE || ~hasTM
            uialert(fig, 'Combined optimization requires both TE and TM experimental datasets.', 'Missing Data', 'Icon','warning');
            return;
        end
        
        wl = edWavelength.Value;
        V = edVoltage.Value;
        n_bk = 1.518;
        n_ito = 1.8529 + 0.00316i;
        n_water = 1.333;
        degree = pi/180;
        
        lb = [10, 40, 1.30, 0.01, 1.30, 0.01, -10];
        ub = [30, 80, 1.50, 0.20, 1.50, 0.20, 10];
        
        useCombinedAngle = ask_use_combined_angle();
        combinedShift = app.OptimizedShift;
        if isempty(combinedShift) || ~isfinite(combinedShift)
            combinedShift = 0;
        end
        if useCombinedAngle
            lb(7) = combinedShift;
            ub(7) = combinedShift;
        end
        
        show_boundary_info('Particle Swarm (Combined TE & TM)', ...
            {'ITO thickness (nm)','PEDOT thickness (nm)', 'n_{real}^{TE}','\kappa^{TE}', 'n_{real}^{TM}','\kappa^{TM}','Angle shift (deg)'}, ...
            lb, ub);
        
        objective_fun = @(params) combined_error_value(params, expTE, expTM, n_bk, n_ito, n_water, degree, wl);
        
        options = optimoptions('particleswarm', ...
            'Display', 'off', ...
            'SwarmSize', 600, ...
            'MaxIterations', 1200, ...
            'MaxStallIterations', 80, ...
            'FunctionTolerance', 1e-8, ...
            'UseParallel', false);
        
        if isvalid(dlg)
            dlg.Message = 'Running Particle Swarm optimization (Combined TE & TM)...';
        end
        
        try
            [x_opt, fval, minima] = custom_particleswarm_with_tracking(objective_fun, 7, lb, ub, options);
        catch ME
            warning('Combined TE/TM Particle Swarm failed: %s', ME.message);
            uialert(fig, ME.message, 'Combined Optimization Error', 'Icon','error');
            return;
        end
        
        [~, err_te, err_tm, valid_te, valid_tm] = combined_error_value(x_opt, expTE, expTM, n_bk, n_ito, n_water, degree, wl);
        
        combo.ito_thickness = x_opt(1);
        combo.pedot_thickness = x_opt(2);
        combo.n_real_te = x_opt(3);
        combo.kappa_te = x_opt(4);
        combo.n_real_tm = x_opt(5);
        combo.kappa_tm = x_opt(6);
        combo.shift = x_opt(7);
        combo.err_te = err_te;
        combo.err_tm = err_tm;
        combo.total_error = fval;
        combo.valid_points_te = valid_te;
        combo.valid_points_tm = valid_tm;
        combo.iterations = options.MaxIterations;
        combo.all_local_minima = minima;
        combo.optimization_method = 'particleswarm_combined';
        combo.optimization_details = struct('SwarmSize', options.SwarmSize, ...
                                            'MaxIterations', options.MaxIterations, ...
                                            'MaxStallIterations', options.MaxStallIterations);
        
        combo.TE = struct('n_real', combo.n_real_te, 'kappa', combo.kappa_te, ...
            'shift', combo.shift, 'ito_thickness', combo.ito_thickness, ...
            'pedot_thickness', combo.pedot_thickness, 'err', combo.err_te, ...
            'valid_points', combo.valid_points_te, 'optimization_method', 'particleswarm_combined', ...
            'iterations', combo.iterations);
        combo.TM = struct('n_real', combo.n_real_tm, 'kappa', combo.kappa_tm, ...
            'shift', combo.shift, 'ito_thickness', combo.ito_thickness, ...
            'pedot_thickness', combo.pedot_thickness, 'err', combo.err_tm, ...
            'valid_points', combo.valid_points_tm, 'optimization_method', 'particleswarm_combined', ...
            'iterations', combo.iterations);
    end
    
    function [total_err, err_te, err_tm, valid_te, valid_tm] = combined_error_value(params, expTE, expTM, n_bk, n_ito, n_water, degree, wl)
        ito_thick = params(1);
        pedot_thick = params(2);
        n_real_te = params(3);
        kappa_te = params(4);
        n_real_tm = params(5);
        kappa_tm = params(6);
        angle_shift = params(7);
        
        [err_te, valid_te] = compute_mode_error(n_real_te, kappa_te, angle_shift, 's', expTE, ito_thick, pedot_thick, n_bk, n_ito, n_water, degree, wl);
        [err_tm, valid_tm] = compute_mode_error(n_real_tm, kappa_tm, angle_shift, 'p', expTM, ito_thick, pedot_thick, n_bk, n_ito, n_water, degree, wl);
        
        total_err = err_te + err_tm;
    end
    
    function [err, validCount] = compute_mode_error(n_real, kappa, angle_shift, pol, expDataStruct, ito_thick, pedot_thick, n_bk, n_ito, n_water, degree, wl)
        n_pedot = n_real + 1i*kappa;
        n_list = [n_bk, n_ito, n_pedot, n_water];
        d_list = [Inf, ito_thick, pedot_thick, Inf];
        
        exp_angles = expDataStruct.original_angles;
        exp_reflectivities = expDataStruct.reflectivities(:);
        
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
        valid = ~isnan(theo);
        validCount = nnz(valid);
        if validCount == 0
            err = 1e6; % large penalty
        else
            residual = theo(valid) - exp_reflectivities(valid);
            err = sum(residual.^2);
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
    
    function show_boundary_info(methodName, paramNames, lbVals, ubVals)
        lines = cell(1, numel(paramNames));
        for ii = 1:numel(paramNames)
            lines{ii} = sprintf('%s: %.3f to %.3f', paramNames{ii}, lbVals(ii), ubVals(ii));
        end
        msg = strjoin(lines, newline);
        uialert(fig, msg, [methodName ' Bounds'], 'Icon','info');
    end
    
    
    
    
    
    
    function best = optimize_mode_grid(pol, expData, dlg)
        % Original grid search method (kept as fallback)
        wl = edWavelength.Value;
        V = edVoltage.Value;
        
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
        
        best.optimization_method = 'grid_search';
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

    function [expData, ok] = load_experimental_data(pol)
        file = 'Epoly_35s_Spol.xlsx'; if pol=='p', file='Epoly_35s_Ppol.xlsx'; end
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
    % Calculation and Plotting Functions (same as Python)
    %====================
    function on_calculate()
        oldStates = disable_longrun_controls();
        d = uiprogressdlg(fig,'Title','Calculating','Indeterminate','on',...
            'Message','Computing reflectivity...');
        try
            lblStatus.Text='Calculating...'; lblStatus.FontColor=[0 0 0.7]; drawnow;
            wl = edWavelength.Value;
            V = edVoltage.Value;
            nreal = edNreal.Value;
            pol = 's'; if rbTM.Value, pol='p'; end
            layer_thicknesses = struct('ito',edITO.Value,'pedot',edPEDOT.Value);
            useDk = chkUseDk.Value;
            delta_k = edDk.Value;
            plot_choice = ddPlot.Value;
            showExp = chkShowExp.Value;
            if chkManual.Value
                kreal = edK.Value; manual_k = kreal;
            else
                kreal = get_k_improved(app.PEDOT, wl, V); manual_k = NaN;
            end

            cla(ax);

            if strcmp(plot_choice,'Combined TE/TM (all 4 curves)')
                n_te = edNte.Value; n_tm = edNtm.Value;
                if chkManual.Value
                    k_te = edKte.Value; k_tm = edKtm.Value;
                else
                    k_te = get_k_improved(app.PEDOT, wl, V);
                    k_tm = k_te;
                end
                shift_te = edShiftTE.Value;
                shift_tm = edShiftTM.Value;
                plot_combined_te_tm(ax, wl, V, n_te, k_te, n_tm, k_tm, layer_thicknesses, shift_te, shift_tm, showExp);
                lblStatus.Text='Ready'; lblStatus.FontColor=[0 0.5 0];
                if isvalid(d), close(d); end
                restore_controls(oldStates);
                return;
            end

            if useDk
                [angles, R_base, R_mod, dR, dR_over_R] = calculate_delta_reflectivity(wl, V, nreal, layer_thicknesses, pol, delta_k, manual_k);
                plot_single_graph(ax, angles, R_base, R_mod, dR, dR_over_R, wl, V, nreal, kreal, pol, plot_choice, delta_k, showExp, layer_thicknesses, app.OptimizedShift);
            else
                [angles, R_base] = calculate_reflectivity_vs_angle(wl, V, nreal, layer_thicknesses, pol, 0, manual_k);
                plot_only_reflectivity(ax, angles, R_base, wl, V, nreal, kreal, pol, showExp, layer_thicknesses, app.OptimizedShift);
            end
            lblStatus.Text='Ready'; lblStatus.FontColor=[0 0.5 0];
        catch ME
            lblStatus.Text=['Error: ' ME.message]; lblStatus.FontColor=[0.85 0 0];
            uialert(fig,ME.message,'Error during calculation');
        end
        if isvalid(d), close(d); end
        restore_controls(oldStates);
    end

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

    function plot_single_graph(axh, angles, R_base, R_mod, dR, dR_over_R, wavelength, voltage, n_real, k_real, polarization, plot_type, delta_k, show_experimental, layer_thicknesses, optimized_shift)
        cla(axh);
        hold(axh,'on');
        if strcmp(plot_type,'Reflectivity vs. Angle')
            plot(axh, angles, R_base, 'LineWidth',2,'DisplayName',sprintf('Base: n=%.3f + %.6fi',n_real,k_real));
            plot(axh, angles, R_mod, '--','LineWidth',2,'DisplayName',sprintf('Modified: n=%.3f + %.6fi',n_real,k_real+delta_k));
            if show_experimental
                plot_experimental(axh, polarization, optimized_shift);
            end
            ylabel(axh,'Reflectivity');
            title(axh, compose_title(wavelength, voltage, polarization, layer_thicknesses, optimized_shift));
            legend(axh,'Location','best');
        elseif strcmp(plot_type,'ΔR vs. Angle')
            plot(axh, angles, dR, 'LineWidth',2,'DisplayName',sprintf('ΔR (Δk = %.3f)',delta_k));
            yline(axh,0,'--','Color',[0 0 0 0.4]);
            ylabel(axh,'ΔR');
            title(axh, compose_title(wavelength, voltage, polarization, layer_thicknesses, optimized_shift, sprintf(', Δk=%.3f',delta_k)));
            legend(axh,'Location','best');
        elseif strcmp(plot_type,'ΔR/R vs. Angle')
            plot(axh, angles, dR_over_R, 'LineWidth',2,'DisplayName','ΔR/R');
            yline(axh,0,'--','Color',[0 0 0 0.4]);
            ylabel(axh,'ΔR/R');
            title(axh, compose_title(wavelength, voltage, polarization, layer_thicknesses, optimized_shift, sprintf(', Δk=%.3f',delta_k)));
            legend(axh,'Location','best');
        elseif strcmp(plot_type,'ΔR/Δk vs. Angle')
            if delta_k~=0
                plot(axh, angles, dR./delta_k, 'LineWidth',2,'DisplayName','ΔR/Δk');
                ylabel(axh,'ΔR/Δk');
                title(axh, compose_title(wavelength, voltage, polarization, layer_thicknesses, optimized_shift));
                legend(axh,'Location','best');
            else
                text(axh,0.5,0.5,'Cannot calculate ΔR/Δk when Δk = 0','Units','normalized','HorizontalAlignment','center');
            end
        end
        xlabel(axh,'Incident Angle (degrees)');
        xlim(axh,[0 90]); grid(axh,'on'); xticks(axh,0:10:90);
        hold(axh,'off');
    end

    function plot_only_reflectivity(axh, angles, R_base, wavelength, voltage, n_real, k_real, polarization, show_experimental, layer_thicknesses, optimized_shift)
        cla(axh); hold(axh,'on');
        plot(axh, angles, R_base, 'LineWidth',2,'DisplayName',sprintf('n=%.3f + %.6fi',n_real,k_real));
        if show_experimental
            plot_experimental(axh, polarization, optimized_shift);
        end
        xlabel(axh,'Incident Angle (degrees)'); ylabel(axh,'Reflectivity');
        title(axh, compose_title(wavelength, voltage, polarization, layer_thicknesses, optimized_shift));
        xlim(axh,[20 90]); grid(axh,'on'); xticks(axh,20:5:90);
        legend(axh,'Location','best');
        hold(axh,'off');
    end

    function txt = compose_title(wavelength, voltage, polarization, layer_thicknesses, optimized_shift, extra)
        if nargin<6, extra=''; end
        poltxt = 's'; if polarization=='p', poltxt='p'; end
        txt = sprintf('Reflectivity vs. Angle\nλ = %.0f nm, V = %.0f mV, %s-polarization\nITO=%.1f nm, PEDOT=%.1f nm%s',...
            wavelength, voltage, poltxt, layer_thicknesses.ito, layer_thicknesses.pedot, extra);
        if optimized_shift ~= 0
            txt = sprintf('%s\nOptimized angle shift: %.1f°',txt,optimized_shift);
        end
    end

    function plot_experimental(axh, polarization, optimized_shift)
        file = 'Epoly_35s_Spol.xlsx'; modeName='TE'; marker='o'; if polarization=='p', file='Epoly_35s_Ppol.xlsx'; modeName='TM'; marker='s'; end
        if ~isfile(file), return; end
        try
            T = readtable(file);
            if width(T) < 2, return; end
            orig = T{:,1}; Rexp = T{:,2};
            mapped = -1 .* orig + 273 + optimized_shift;
            [mapped, idx] = sort(mapped); Rexp = Rexp(idx);
            plot(axh, mapped, Rexp, '-','Marker',marker,'LineWidth',1.5,'DisplayName',sprintf('Experimental (%s)',modeName));
        catch
        end
    end

    function plot_combined_te_tm(axh, wavelength, voltage, n_te, k_te, n_tm, k_tm, layer_thicknesses, shift_te, shift_tm, showExp)
        cla(axh); hold(axh,'on');
        [angTE, Rte] = calculate_reflectivity_vs_angle(wavelength, voltage, n_te, layer_thicknesses, 's', 0, k_te);
        plot(axh, angTE, Rte, 'LineWidth',2,'DisplayName',sprintf('TE sim: n=%.3f+%.4fi',n_te,k_te));
        [angTM, Rtm] = calculate_reflectivity_vs_angle(wavelength, voltage, n_tm, layer_thicknesses, 'p', 0, k_tm);
        plot(axh, angTM, Rtm, 'LineWidth',2,'DisplayName',sprintf('TM sim: n=%.3f+%.4fi',n_tm,k_tm));
        if showExp
            if isfile('Epoly_35s_Spol.xlsx')
                T = readtable('Epoly_35s_Spol.xlsx'); if width(T)>=2
                    ang = T{:,1}; R = T{:,2};
                    mapped = -1 .* ang + 273 + shift_te;
                    [mapped, idx] = sort(mapped); R = R(idx);
                    plot(axh, mapped, R, '-o','LineWidth',1.5,'DisplayName','TE exp');
                end
            end
            if isfile('Epoly_35s_Ppol.xlsx')
                T = readtable('Epoly_35s_Ppol.xlsx'); if width(T)>=2
                    ang = T{:,1}; R = T{:,2};
                    mapped = -1 .* ang + 273 + shift_tm;
                    [mapped, idx] = sort(mapped); R = R(idx);
                    plot(axh, mapped, R, '-s','LineWidth',1.5,'DisplayName','TM exp');
                end
            end
        end
        xlabel(axh,'Incident Angle (degrees)'); ylabel(axh,'Reflectivity');
        title(axh, sprintf('Combined TE/TM Reflectivity\nλ=%.0f nm, V=%.0f mV\nITO=%.1f nm, PEDOT=%.1f nm',...
            wavelength, voltage, layer_thicknesses.ito, layer_thicknesses.pedot));
        xlim(axh,[20 90]); ylim(axh,[0 1]); grid(axh,'on'); xticks(axh,20:10:90);
        legend(axh,'Location','best');
        hold(axh,'off');
    end

end % <-- end of main function

