function varargout = Deuteros(varargin)
% DEUTEROS MATLAB code for Deuteros.fig
% Last Modified by GUIDE v2.5 10-Dec-2018 14:23:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Deuteros_OpeningFcn, ...
                   'gui_OutputFcn',  @Deuteros_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Deuteros is made visible.
function Deuteros_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for Deuteros
handles.output = hObject;
axis off

% Logo options
axes(handles.logo)

matlabImage = imread('Deuteros_logo.png');
image(matlabImage)
axis off
axis image

% Make the export Pymol_maxuptake_textbox and Pymol_maxuptake_editbox objects invisible 
set(handles.Pymol_maxuptake_editbox,'Visible','Off')
set(handles.Pymol_maxuptake_textbox,'Visible','Off')

set(handles.FilterA_dropdown,'Value',5);
set(handles.FilterB_dropdown,'Value',4);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Deuteros wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Deuteros_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% % general graphics, this will apply to any figure you open
% % (groot is the default figure object).
% set(groot, ...
% 'DefaultFigureColor', [229/255, 229/255, 229/255], ...
% 'DefaultAxesLineWidth', 0.5, ...
% 'DefaultAxesXColor', 'k', ...
% 'DefaultAxesYColor', 'k', ...
% 'DefaultAxesFontUnits', 'points', ...
% 'DefaultAxesFontSize', 8, ...
% 'DefaultAxesFontName', 'Helvetica', ...
% 'DefaultLineLineWidth', 1, ...
% 'DefaultTextFontUnits', 'Points', ...
% 'DefaultTextFontSize', 8, ...
% 'DefaultTextFontName', 'Helvetica', ...
% 'DefaultAxesBox', 'off', ...
% 'DefaultAxesTickLength', [0.02 0.025], ...
% 'DefaultGridColor', [1,1,1], ...
% 'DefaultGridAlpha', 1);
%  
% % set the tickdirs to go out - need this specific order
% set(groot, 'DefaultAxesTickDir', 'out');
% set(groot, 'DefaultAxesTickDirMode', 'manual');


% --------------------------------------------------------------------------------------------------
% 1. Import Difference and State Data
% --------------------------------------------------------------------------------------------------
 
        % ----------------------------------------------------------------------------
        % 1.1 Path to difference data file
        % ----------------------------------------------------------------------------
        function DiffDataPath_browse_Callback(~, Diff_eventdata, handles)

        check = evalin('base', 'exist(''previous_path'')' );
        
        if ~check
            disp('Previous path not found, please locate folder..')
            
            [Diff_filename, Diff_pathname] = uigetfile('*.csv','Select Difference Data ..');
            
            assignin('base','previous_path',Diff_pathname)
            disp(sprintf('Assigned %s as last accessed folder', Diff_pathname))
            disp(' ')
        else
            
            previous_path = evalin('base','previous_path')
            disp(sprintf('Previous folder %s found', previous_path))
            
            [Diff_filename, Diff_pathname] = uigetfile('*.csv','Select Difference Data ..',previous_path);

            assignin('base','previous_path',Diff_pathname)
            disp(sprintf('Assigned %s as last accessed folder', Diff_pathname))
            disp(' ')
            
        end

        Diff_fullpathname = strcat(Diff_pathname, Diff_filename);
        Diffpathdisp = set(handles.Import_differenceDataPath,'String',Diff_filename);

        set(handles.Import_textbox,'String',' ');
        set(handles.DiffDataPath_browse,'UserData',Diff_fullpathname);

        % ----------------------------------------------------------------------------
        % 1.2 Path to state data file 
        % ----------------------------------------------------------------------------
        function StateDataPath_browse_Callback(~, State_eventdata, handles)
            
        check = evalin('base', 'exist(''previous_path'')' );
        
        if ~check
            disp('Previous path not found, please locate folder..')
            
            [State_filename, State_pathname] = uigetfile('*.csv','Select Difference Data ..');
            
            assignin('base','previous_path',State_pathname)
            disp(sprintf('Assigned %s as last accessed folder', State_pathname))
            disp(' ')
        else
            
            previous_path = evalin('base','previous_path')
            disp(sprintf('Previous folder %s found', previous_path))
            
            [State_filename, State_pathname] = uigetfile('*.csv','Select Difference Data ..',previous_path);

            assignin('base','previous_path',State_pathname)
            disp(sprintf('Assigned %s as last accessed folder', State_pathname))
            disp(' ')
            
        end
            
        State_fullpathname = strcat(State_pathname, State_filename);
        Statepathdisp = set(handles.Import_stateDataPath,'String',State_filename);

        set(handles.Import_textbox,'String',' ');
        set(handles.StateDataPath_browse,'UserData',State_fullpathname);
        
        % ----------------------------------------------------------------------------
        % 1.3 Stores sequence start value
        % ----------------------------------------------------------------------------
        
        function Import_sequenceStart_Callback(hObject, eventdata, handles)
        function Import_sequenceStart_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        
        % ----------------------------------------------------------------------------
        % 1.4 Stores sequence end value
        % ----------------------------------------------------------------------------

        function Import_sequenceEnd_Callback(hObject, eventdata, handles)
        function Import_sequenceEnd_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        % ----------------------------------------------------------------------------
        % 1.5 Stores CV1 value
        % ----------------------------------------------------------------------------
        
        function Import_alpha1_Callback(hObject, eventdata, handles)
        function Import_alpha1_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        
        % --- Executes on selection change in FilterA_dropdown.
        function FilterA_dropdown_Callback(hObject, eventdata, handles)
        % hObject    handle to FilterA_dropdown (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)

        % Hints: contents = cellstr(get(hObject,'String')) returns FilterA_dropdown contents as cell array
        %        contents{get(hObject,'Value')} returns selected item from FilterA_dropdown


        % --- Executes during object creation, after setting all properties.
        function FilterA_dropdown_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to FilterA_dropdown (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called

        % Hint: popupmenu controls usually have a white background on Windows.
        %       See ISPC and COMPUTER.
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end


        
        % ----------------------------------------------------------------------------
        % 1.6 Stores CV2 value
        % ----------------------------------------------------------------------------

        function Import_alpha2_Callback(hObject, eventdata, handles)
        function Import_alpha2_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        
        
        % --- Executes on selection change in FilterB_dropdown.
        function FilterB_dropdown_Callback(hObject, eventdata, handles)
        % hObject    handle to FilterB_dropdown (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)

        % Hints: contents = cellstr(get(hObject,'String')) returns FilterB_dropdown contents as cell array
        %        contents{get(hObject,'Value')} returns selected item from FilterB_dropdown


        % --- Executes during object creation, after setting all properties.
        function FilterB_dropdown_CreateFcn(hObject, eventdata, handles)
        % hObject    handle to FilterB_dropdown (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    empty - handles not created until after all CreateFcns called

        % Hint: popupmenu controls usually have a white background on Windows.
        %       See ISPC and COMPUTER.
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        % ----------------------------------------------------------------------------
        % 1.7 Stores Nreps (Number of Repeats) value
        % ----------------------------------------------------------------------------
        
        function Import_Nreps_Callback(hObject, eventdata, handles)
        function Import_Nreps_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        % ----------------------------------------------------------------------------
        % 1.8 Button for Import Data
        % ----------------------------------------------------------------------------
        
        function [Diff_data] = Import_button_Callback(hObject, eventdata, handles)
        Diff_path = get(handles.DiffDataPath_browse,'UserData');
        State_path = get(handles.StateDataPath_browse,'UserData');

        Diff_import = csv2cell(Diff_path, 'fromfile');
        State_import = csv2cell(State_path, 'fromfile');
        
        % ------- New import method --------
        Diff_data = {};
        State_data = {};
        
        % Difference Data
        diff_num_only = Diff_import(2:end,2:end);
        diff_num_only(:,3) = {'0'};
        Diff_data.data = cell2num(diff_num_only);
        Diff_data.textdata = Diff_import(:,1);

        % State Data
        state_num_only = State_import(2:end,10:end);
        State_data.data = cell2num(state_num_only);
        
        % ----------------------------------

        % Seq Start/End Check
        seqStart = str2num(get(handles.Import_sequenceStart,'String'));
        seqEnd = str2num(get(handles.Import_sequenceEnd,'String'));

        if seqEnd < seqStart
            set(handles.seqEnd_errorbox,'String','End must be greater than start!');
        elseif isreal(seqEnd) && rem(seqEnd,1)==0
            set(handles.seqEnd_errorbox,'String',' ');
        else
            set(handles.seqEnd_errorbox,'String','Invalid last residue number!');
        end

        % Confidence Interval Check
        
        % T-distribution table values, two-tailed
        
        % Get spectrum colour
        confidence_A = get(handles.FilterA_dropdown, 'Value');
        confidence_B = get(handles.FilterB_dropdown, 'Value');

        if confidence_A == 1 % 80%, p=0.02
            CV1 = 1.886;
        elseif confidence_A == 2 % 90%, p=0.10
            CV1 = 2.920;
        elseif confidence_A == 3 % 95%, p=0.05
            CV1 = 4.303;
        elseif confidence_A == 4 % 98%, p=0.02
            CV1 = 6.965;
        elseif confidence_A == 5 % 99%, p 0.01
            CV1 = 9.925;
        elseif confidence_A == 6 % 99.5%, p=0.005
            CV1 = 14.089;
        elseif confidence_A == 7 % 99.8%, p=0.002
            CV1 = 22.327;
        elseif confidence_A == 8 % 99.9%, p=0.001
            CV1 = 31.599;
        end
        

        if confidence_B == 1 % 80%, p=0.02
            CV2 = 1.886;
        elseif confidence_B == 2 % 90%, p=0.10
            CV2 = 2.920;
        elseif confidence_B == 3 % 95%, p=0.05
            CV2 = 4.303;
        elseif confidence_B == 4 % 98%, p=0.02
            CV2 = 6.965;
        elseif confidence_B == 5 % 99%, p 0.01
            CV2 = 9.925;
        elseif confidence_B == 6 % 99.5%, p=0.005
            CV2 = 14.089;
        elseif confidence_B == 7 % 99.8%, p=0.002
            CV2 = 22.327;
        elseif confidence_B == 8 % 99.9%, p=0.001
            CV2 = 31.599;
        end

        Nrep = str2double(get(handles.Import_Nreps,'String')); % get number of repeats
        
        if Nrep == 0
            set(handles.n_errorbox,'String','Cannot have 0 repeats')
        elseif rem(Nrep,1)~=0
            set(handles.n_errorbox,'String','Must be whole number')
        elseif isnan(Nrep) == 0
            set(handles.n_errorbox,'String',' ')
        else
            set(handles.n_errorbox,'String','Enter value of n')
        end

        state_timepoints = unique(State_data.data(:,1));

        sortedStateData = sortrows(State_data.data,1);
            [~,~,uniqueIndex] = unique(sortedStateData(:,1));
            sortedStateData = mat2cell(sortedStateData,accumarray(uniqueIndex(:),1),7);


        CI_timepoints = [];
        for i=1:length(sortedStateData)
            CI_timepoints.SEM(i) = nanmean(sortedStateData{i,1}(:,5))^2;

            CI_timepoints.CV1(i) = (sqrt(CI_timepoints.SEM(i))/sqrt(Nrep))*CV1;
            CI_timepoints.CV2(i) = (sqrt(CI_timepoints.SEM(i))/sqrt(Nrep))*CV2;
        end

        % Add the sum confidence intervals including propagating error
        CI_timepoints.CV1(length(CI_timepoints.CV1)+1) = sqrt(nansum(CI_timepoints.SEM(2:end)))/sqrt(length(CI_timepoints.SEM(2:end)))*CV1;
        CI_timepoints.CV2(length(CI_timepoints.CV2)+1) = sqrt(nansum(CI_timepoints.SEM(2:end)))/sqrt(length(CI_timepoints.SEM(2:end)))*CV2;

        import.diffData = Diff_data;
        import.stateData = sortedStateData;
        import.stateDataAll = State_data;
        import.CI_timepoints = CI_timepoints;
        import.state_timepoints = state_timepoints;

        % Tell user that data has finished loading
        set(handles.Import_textbox,'String','Data loaded');

        % Change the LinearPlot section's DataType popup box values to reflect the timepoints in the data
        dataTypeTextBox = {'Coverage'; 'Redundancy'; 'Heatmap:'}; % Needs to be re-set up each time otherwise same timepoints will get re-appended
        for i=2:length(state_timepoints)+1
            if i ~= length(state_timepoints)+1
                timepoint_string = [' - Timepoint ',num2str(state_timepoints(i)), ' min'];
                dataTypeTextBox{end+1} = timepoint_string;
            elseif i == length(state_timepoints)+1
                if length(state_timepoints) > 2
                    timepoint_string = ' - Timepoint Sum';
                    dataTypeTextBox{end+1} = timepoint_string;
                end
            end
        end

        set(handles.Linear_DataType_Pop,'string',dataTypeTextBox); % Now set the new timepoints to the popup menu

        % Set enable sum on or off depending on timepoints available
        if length(state_timepoints)-1 == 1
            set(handles.Woods_EnableSum_Pop,'enable','off');
            set(handles.Woods_EnableSum_Pop,'enable','off');
        else
            set(handles.Woods_EnableSum_Pop,'enable','on');
            set(handles.Woods_EnableSum_Pop,'enable','on');
        end
        
        % Get number of residues in each peptide, -1 and -N_proline for max uptake
        difference_peptides = Diff_data.textdata(2:end,1);
        max_uptakes = [];
        for i = 1:length(difference_peptides)
            max_DU = cellfun('length',difference_peptides(i)) - count(difference_peptides(i),"P") - 1;
            max_uptakes = [max_uptakes; max_DU];
        end

        % Concatenate peptide start and end with the max uptake value and save to import button's userdata
        compact_peptide_list = horzcat(Diff_data.data(:,1:2), max_uptakes);
        import.compact_peptide_list = compact_peptide_list;
        import.filter_alphas = [CV1, CV2];
        import.filter_vals = [confidence_A, confidence_B];
        
        % Set to import variable so we can call this later when plotting
        import.linearDataType = dataTypeTextBox;
        set(hObject,'UserData',import);
        

% --------------------------------------------------------------------------------------------------
% 2. Flattened Data Maps
% --------------------------------------------------------------------------------------------------

        % ----------------------------------------------------------------------------
        % 2.1 Stores datatype options
        % ----------------------------------------------------------------------------

        function Linear_DataType_Pop_Callback(hObject, eventdata, handles)
        DataType.val = get(handles.Linear_DataType_Pop,'Value');

        if DataType.val < 3
            set(handles.Linear_ColorScale_text,'enable','off');
            set(handles.Linear_ColorScale_Pop,'enable','off');
        elseif DataType.val > 3
            set(handles.Linear_ColorScale_text,'enable','on');
            set(handles.Linear_ColorScale_Pop,'enable','on');
        end

        function Linear_DataType_Pop_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        % ----------------------------------------------------------------------------
        % 2.2 Stores ColorScale options 
        % ----------------------------------------------------------------------------
        
        function Linear_ColorScale_Pop_Callback(hObject, eventdata, handles)
        function Linear_ColorScale_Pop_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        % ----------------------------------------------------------------------------
        % 2.3 Button for generating Linear plot with selection from 2.1 and 2.2
        % ----------------------------------------------------------------------------
        
        function Linear_Plot_Callback(hObject, eventdata, handles)
        axes(handles.LinearPlot);
        cla(handles.LinearPlot);
        set(handles.LinearPlot,'Visible','on')
        cla reset

        % Get start and end sequence number from input box
        seqStart = str2double(get(handles.Import_sequenceStart,'String'));
        seqEnd = str2double(get(handles.Import_sequenceEnd,'String'));

        import = get(handles.Import_button,'UserData'); % Get import data

        timepoints = import.state_timepoints; % Get timepoints
        zero_CI = zeros(1,length(timepoints)); % Make 0 CI values to filter nothing

        % Get the selected data type from Linear Plot Options
        DataType.str = get(handles.Linear_DataType_Pop, 'String');
        DataType.val = get(handles.Linear_DataType_Pop, 'Value')

        if DataType.val == 3 % Plot Redundancy
            disp([DataType.str{DataType.val}, ' selected, select Timepoint'])
            errordlg('Select one of the timepoint options for heatmap')

        elseif DataType.val > 3 % Plot Redundancy
            disp([DataType.str{DataType.val}, ' selected'])
            tp = DataType.val-3; % Remember this does not include 0

            [peptideLengthHeight, peptideLength] = data_to_Woods(import.diffData,zero_CI(1),tp);

        else
            % Generate the matrix for the difference data
            [peptideLengthHeight, peptideLength] = data_to_Woods(import.diffData,zero_CI(1),length(zero_CI)-1);
        end

        uptakePadding = zeros(peptideLengthHeight/2,seqEnd-seqStart+1); % Make a 0 matrix for transplanting uptake values

        for i=1:2:length(peptideLength) % For each peptide, transplant into uptakePadding
            peptideNum = peptideLength(i,:);
            peptideVal = peptideLength(i+1,:);

            % Remove the NaN values
            id1 = find(isnan(peptideNum));
            id2 = find(isnan(peptideVal));
            peptideNum([id1,id2]) = [];
            peptideVal([id1,id2]) = [];

            uptakePadding((i+1)/2,peptideNum(1:end)) = peptideVal;
        end

        linearUptake = []; % Sum each column (uptake of each residue)
        redundancyVec = [];
        for i=1:length(uptakePadding)
            column = uptakePadding(:,i);
            sumUptake = mean(nonzeros(column));
            redundancySum = nnz(column);
            linearUptake = [linearUptake, sumUptake];
            redundancyVec = [redundancyVec, redundancySum];
        end

        ResidueSequence = [seqStart:seqEnd]; % Make vector of sequence to be header

        if length(ResidueSequence) ~= length(redundancyVec)
            errordlg('Data covers range outside of protein sequence given. Please check input values.')
        end

        ResidueRedundancy = vertcat(ResidueSequence, redundancyVec); % Combine to make Peptide/Redundancy vector for case 2
        ResidueUptake = vertcat(ResidueSequence, linearUptake); % Combine to make Peptide/Uptake vector for case 3

        output_struct.ResidueRedundancy = ResidueRedundancy; % Add redundancy and flattened uptake per residue and add to output structure
        output_struct.ResidueUptake = ResidueUptake;

        % Calculate the coverage for the protein
        coveredResidues = [];
        peptideResidues = peptideLength(1:2:end,:);

        for row=1:size(peptideResidues,1)
            peptiderow = peptideResidues(row,:);
            for resi=1:length(peptiderow)
                if isnan(peptiderow(resi)) == 0
                    coveredResidues = [coveredResidues, peptiderow(resi)];
                end
            end
        end

        uniqueResidues = length(unique(coveredResidues));
        expectedResidues = seqEnd-seqStart+1;
        coverage = (uniqueResidues/expectedResidues)*100;

        % ----- Plot decoration -----

        % Plot a grey box background for the theoretical length of the protein, always on
        plot(polyshape([seqStart -1 ; seqEnd -1; seqEnd 1 ; seqStart 1]),'FaceAlpha',1,'FaceColor',[0.8 0.8 0.8], 'LineStyle','none')
        hold on; box off;

        % Axis padding
        axis([seqStart-20 seqEnd+20 -3 3])

        % Axis decoration
        xlabel('Residue Number')
        set(gca,'fontsize',9)

        % Remove ticks from Y-axis
        set(gca,'YTickLabel',[ ]);
        set(gca,'Ytick',[])
        ax1 = gca;
        yruler = ax1.YRuler;
        yruler.Axle.Visible = 'off';

        % Background colour 
        set(gca,'Color',[0.95 0.95 0.95])

        % ----- Data plot -----

        % Plot depending on which case is chosen on popup menu
        if DataType.val == 1 % Plot Basic Coverage

            for i=1:2:peptideLengthHeight
                startRes = nanmin(peptideLength(i,:));
                endRes = nanmax(peptideLength(i,:));
                plot(polyshape([startRes -1 ; endRes -1; endRes 1 ; startRes 1]),'FaceAlpha',1,'FaceColor',[(65/255) (105/255) (225/255)], 'LineStyle','none')
            end

        elseif DataType.val == 2 % Plot Redundancy
            disp([DataType.str{DataType.val}, ' selected'])

            maxRedundancy = max(ResidueRedundancy(2,:));

            white = [1, 1, 1];
            purple = [0.6, 0, 0.7];
            colourStep = [linspace(white(1),purple(1),maxRedundancy)', linspace(white(2),purple(2),maxRedundancy)', linspace(white(3),purple(3),maxRedundancy)'];

            for i=1:length(ResidueRedundancy)
                if ResidueRedundancy(2,i) ~= 0
                    startRes = ResidueRedundancy(1,i)-0.5;
                    endRes = ResidueRedundancy(1,i)+0.5;
                    
                    plot(polyshape([startRes -1 ; endRes -1; endRes 1 ; startRes 1]),'FaceAlpha',1,'FaceColor',colourStep(ResidueRedundancy(2,i),:), 'LineStyle','none')
                end
            end

            cm = colormap(colourStep);
            colormap(cm)
            cb = colorbar;
            cb.Label.String = 'Redundancy';
            caxis([1 maxRedundancy]);
            cb.Ticks = [1, maxRedundancy];

        elseif DataType.val > 3 % Plot heatmap
            cla(handles.LinearPlot);

            maxUptake = max(ResidueUptake(2,:))*10000;
            minUptake = min(ResidueUptake(2,:))*10000;

            extreme = max([maxUptake, abs(minUptake)]);

            % Before colouring, need to replace the uptake value with ranking for
            % colour gradient

            % Define colours to use, and then split gradients into discrete
            % whiteRed/whiteBlue colours
            red = [1, 0.1, 0.1];
            white = [1, 1, 1];
            blue = [0.1, 0.1, 0.9];

            % Set colour scale depending on the user choice 
            ColorScale = get(handles.Linear_ColorScale_Pop,'Value');

            if ColorScale == 1 % Absolute scale, e.g. -10/+10
                whiteRed = [linspace(white(1),red(1),extreme+2)', linspace(white(2),red(2),extreme+2)', linspace(white(3),red(3),extreme+2)'];
                whiteBlue = [linspace(white(1),blue(1),extreme+2)', linspace(white(2),blue(2),extreme+2)', linspace(white(3),blue(3),extreme+2)'];

                % Make colourmap to show blue-white-red colouring
                axisVal = round(extreme/10000,1);
                colormap(redblue(100));
                cb = colorbar;
                cb.Ticks = [-axisVal, 0, axisVal];
                caxis([-axisVal, axisVal]);
                cb.Label.String = 'Uptake (Da)';


            elseif ColorScale == 2 % Normalised scale, e.g. -2/+10
                % With this scheme, the min and max uptake values are always scaled
                % to be bright red and bright blue no matter their value, 0 is
                % always white

                maxUptake = max(ResidueUptake(2,:))*10000;
                minUptake = min(ResidueUptake(2,:))*10000;

                whiteRed = [linspace(white(1),red(1),maxUptake+2)', linspace(white(2),red(2),maxUptake+2)', linspace(white(3),red(3),maxUptake+2)']; % +2 is added so that rounding max/minUptake doesn't cause bad indexing when = 0
                whiteBlue = [linspace(white(1),blue(1),abs(minUptake)+2)', linspace(white(2),blue(2),abs(minUptake)+2)', linspace(white(3),blue(3),abs(minUptake)+2)'];

                % Convert min and max to between 0-100, 50 = white
                upperVal = 50-(((maxUptake/extreme)*100)/-2);
                lowerVal = 50-(((minUptake/extreme)*100)/-2);

                cm = colormap(redblue(100)); % Make blue-white-red gradient with 100 gradations
                colormap(cm); % Apply colourmap
                caxis([round(minUptake/10000,1) round(maxUptake/10000,1)]);

                cb = colorbar; % Show colorbar
                cb.Ticks = [round(minUptake/10000,1), round(maxUptake/10000,1)];
                cb.Label.String = 'Uptake (Da)';
            end

            % Now plot neg and pos peptides separately
            for i=1:length(ResidueUptake)
                if ResidueUptake(2,i) > 0
                    colorID = int64(ResidueUptake(2,i)*10000)+1;
                    startRes = ResidueUptake(1,i)-0.5;
                    endRes = ResidueUptake(1,i)+0.5;
                    
%                     line([ResidueUptake(1,i),ResidueUptake(1,i)],[(ResidueUptake(2,i)*0)-tb_buffer,(ResidueUptake(2,i)*0)+tb_buffer],'Color',whiteRed(colorID,:),'linewidth',4);
                    plot(polyshape([startRes -1 ; endRes -1; endRes 1 ; startRes 1]),'FaceAlpha',1,'FaceColor',whiteRed(colorID,:), 'LineStyle','none')
                elseif ResidueUptake(2,i) < 0
                    colorID = int64(abs(ResidueUptake(2,i)*10000))+1;
                    startRes = ResidueUptake(1,i)-0.5;
                    endRes = ResidueUptake(1,i)+0.5;
                    
                    plot(polyshape([startRes -1 ; endRes -1; endRes 1 ; startRes 1]),'FaceAlpha',1,'FaceColor',whiteBlue(colorID,:), 'LineStyle','none')
%                     line([ResidueUptake(1,i),ResidueUptake(1,i)],[(ResidueUptake(2,i)*0)-tb_buffer,(ResidueUptake(2,i)*0)+tb_buffer],'Color',whiteBlue(colorID,:),'linewidth',4)
                else
                    startRes = ResidueUptake(1,i)-0.5;
                    endRes = ResidueUptake(1,i)+0.5;
                    
                    plot(polyshape([startRes -1 ; endRes -1; endRes 1 ; startRes 1]),'FaceAlpha',1,'FaceColor',[0.8 0.8 0.8], 'LineStyle','none')
                end
            end
        end
        
        % Plot 'Coverage = X%' text
        strmin = ['Coverage = ',num2str(round(coverage,3,'significant')),'%'];
        coverageLabel = text(0.02,0.84,strmin,'Units','normalized');
        coverageLabel.FontSize = 9;
        coverageLabel.FontWeight = 'bold';
        coverageLabel.HorizontalAlignment = 'Left';
        
        % Plot number of peptides text
        number_of_peptides = length(import.diffData.data);
        
        strmin = ['( ',num2str(number_of_peptides),' peptides )'];
        nLabel = text(0.98,0.84,strmin,'Units','normalized');
        nLabel.FontSize = 9;
        nLabel.FontWeight = 'bold';
        nLabel.HorizontalAlignment = 'right';


        set(handles.Linear_Plot,'UserData',output_struct);

        % ----------------------------------------------------------------------------
        % 2.4 Button for exporting Linear plot
        % ----------------------------------------------------------------------------
        
        function Linear_Export_Callback(hObject, eventdata, handles)
        Diff_path = get(handles.DiffDataPath_browse,'UserData');
        DataType.str = get(handles.Linear_DataType_Pop, 'String');
        DataType.val = get(handles.Linear_DataType_Pop, 'Value');
        
        if DataType.val == 1
            filename = strcat(Diff_path(1:end-4),' coverage.pdf');
        elseif DataType.val == 2
            filename = strcat(Diff_path(1:end-4),' redundancy.pdf');
        elseif DataType.val > 2
            timepoint_string = DataType.str{DataType.val};
            filename = strcat(Diff_path(1:end-4),' heatmap ',timepoint_string(13:end),'.pdf');
        end
        
        export_fig(handles.LinearPlot,filename)
        
        
% --------------------------------------------------------------------------------------------------
% 3. Woods Plot
% --------------------------------------------------------------------------------------------------

        % ----------------------------------------------------------------------------
        % 3.1 Data type: Binary or Gradation
        % ----------------------------------------------------------------------------
        
        function Woods_DataType_Pop_Callback(hObject, eventdata, handles)
        function Woods_DataType_Pop_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end

        % ----------------------------------------------------------------------------
        % 3.2 Filter using CV 1 or 2
        % ----------------------------------------------------------------------------
        
        function Woods_FilterUsing_Pop_Callback(hObject, eventdata, handles)
        function Woods_FilterUsing_Pop_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        
        % ----------------------------------------------------------------------------
        % 3.3 Show SUM data or not
        % ----------------------------------------------------------------------------
        
        function Woods_EnableSum_Pop_Callback(hObject, eventdata, handles)
        function Woods_EnableSum_Pop_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        
        % ----------------------------------------------------------------------------
        % 3.4 Button for plotting Woods plot
        % ----------------------------------------------------------------------------
        
        function Woods_Plot_Callback(hObject, eventdata, handles)
%         WoodsUIPanel = uipanel('parent',gcf,...
%                                'Position',[0.02742073693230506 0.010945273631840797 0.7437874892887748 0.6388059701492538],...
%                                'BackgroundColor',[1,1,1],...
%                                'HighlightColor',[1,1,1],...
%                                'BorderType','None',...
%                                'Clipping','on',...
%                                'tag','WoodsPlotUIPanel');
                           
        import = get(handles.Import_button,'UserData'); % Import the data

        seqStart = str2num(get(handles.Import_sequenceStart,'String')); % Get start and end numbers
        seqEnd = str2num(get(handles.Import_sequenceEnd,'String'));
        
        % ----- Plotting -----
        
        arrayfun(@delete,findobj(gcf,'type','axes','tag','WoodsPlot')); % clear axes each time you press plot

        % colours to use
        grey = [150/255 150/255 150/255];
        white = [255/255 255/255 255/255];
        blueWhite = [220/255 220/255 255/255];
        blue = [126/255 167/255 207/255];
        redWhite = [255/255 220/255 220/255];
        red = [244/255 98/255 100/255];

        % Get things for plotting
        timepoints = import.state_timepoints; % Get timepoints

        confidenceIntervals = import.CI_timepoints; % Get confidence intervals
        confidenceInterval_CV1 = confidenceIntervals.CV1;
        confidenceInterval_CV2 = confidenceIntervals.CV2;

        data = import.diffData; % Get difference data
        
        % Generate the sum column 
        lengthData = length(data.data(1,:));
        for k=1:length(data.data)
            peptideSum = sum(data.data(k,4:end));
            data.data(k,lengthData+1) = peptideSum;
        end
        
        % Get absolute max differential uptake value for yaxis scaling
        ymax = max(max(abs(data.data(:,4:end-1))));
        ymax = max([ymax, confidenceInterval_CV1, confidenceInterval_CV2]); % Also take into account CI values for axis determination
        
        ymax_sum = max(max(abs(data.data(:,4:end))));
        ymax_sum = max([ymax_sum, confidenceInterval_CV1, confidenceInterval_CV2]);
        
        ymax_round = round(ymax/0.5)*0.5; % get nearest 0.5 value
        buffer = 0.5; % increase margin spacing by factor, default 0.2
        ymax_sum_round = round(ymax_sum/0.5)*0.5;
%         buffer_sum = ymax_sum_round*0.2;

        % Get values for other plot options
        enable_sum = get(handles.Woods_EnableSum_Pop,'Value'); % Check if enable sum is enabled
        filter_using = get(handles.Woods_FilterUsing_Pop,'Value'); % Check which filter to use 
        DataType = get(handles.Woods_DataType_Pop,'Value'); % Check data type to plot

        % PyMOL_struct = {} % Set up empty cell to put data in for exporting

        peptide_data = []; % Set up empty matrix to put data in for exporting
        max_uptake_over_all_timepoints = []; % Collect all the max uptakes to PyMOL output scaling later

        % -------------------------------------------------
        % Pack some subpanels for plotting each Woods plot
        % -------------------------------------------------

        rows = 2; % default is always two rows
        cols = round((length(timepoints)/2),0); % round number of timepoints up to the nearest even number, inclusive of sum

        p = panel(handles.WoodsUIPanel); % create a new panel
        p.pack() % first create a global panel inside the uipanel
        p(1).pack(rows,cols) % pack in 2 rows with a suitable number of columns
        
        curr_timepoint = 1; % start counting at first timepoint
        
        A_val = get(handles.FilterA_dropdown, 'Value');
        B_val = get(handles.FilterB_dropdown, 'Value');
        choices = get(handles.FilterA_dropdown, 'String'); % choice of confidence interval filters
        
        for r=1:rows
            for c=1:cols
                i = curr_timepoint; 

                if curr_timepoint < length(timepoints)
                    p(1,r,c).select() % Select the correct subplot 
                    hold on; grid on; % turn the plot formatting on
                    
                    set(gca,'color',[229/255, 229/255, 229/255]);
                    
                    hAxes = gca;     %Axis handle
                    hAxes.XRuler.Axle.LineStyle = 'none';  
                    hAxes.YRuler.Axle.LineStyle = 'none';
                    hAxes.GridColor = [1,1,1];
                    hAxes.GridAlpha = 1;
                    set(gca,'TickDir','out');
                    
                    title({'';'';['Exposure Time ' num2str(timepoints(curr_timepoint+1)) ' min']},'fontweight','bold')
                    
                    
                    % disp(sprintf('row number %i, col number %i, title %s', r, c, num2str(timepoints(curr_timepoint+1))))
                    
                    line([seqStart seqEnd],[0 0],'Color','white','LineStyle','-','LineWidth',1) % baseline
                    CV1_pos = dashline([seqStart seqEnd],[confidenceInterval_CV1(i+1) confidenceInterval_CV1(i+1)],1.5,1,1.5,1,'LineWidth',1.5,'Color','black');
                    CV1_neg = dashline([seqStart seqEnd],[-confidenceInterval_CV1(i+1) -confidenceInterval_CV1(i+1)],1.5,1,1.5,1,'LineWidth',1.5,'Color','black'); % CV1 neg
                    CV2_pos = dashline([seqStart seqEnd],[confidenceInterval_CV2(i+1) confidenceInterval_CV2(i+1)],0.5,0.5,0.5,0.5,'LineWidth',1.5,'Color','black'); % CV2 pos
                    CV2_neg = dashline([seqStart seqEnd],[-confidenceInterval_CV2(i+1) -confidenceInterval_CV2(i+1)],0.5,0.5,0.5,0.5,'LineWidth',1.5,'Color','black'); % CV2 neg

                    % Get confidenceFilter depending on 
                    if filter_using == 1
                        confidence_A = choices(A_val);
                        
                        confidenceFilter = confidenceInterval_CV1(i+1);
                        CV_text = strcat(confidence_A,'%: 0±',num2str(round(confidenceInterval_CV1(i+1),2)),' Da');
                    elseif filter_using == 2
                        confidence_B = choices(B_val);
                        
                        confidenceFilter = confidenceInterval_CV2(i+1);
                        CV_text = strcat(confidence_B,'%: 0±',num2str(round(confidenceInterval_CV2(i+1),2)),' Da');
                    end
                    
                    [peptideLengthHeight, peptideLength, confidenceMatHighHeight, confidenceMatHigh, confidenceMatLowHeight, confidenceMatLow] = data_to_Woods(data,confidenceFilter,i);

                    output_struct.peptide_data{i,1} = [peptideLength]; % Add the peptide data to the peptide_data list
                    output_struct.peptide_data{i,2} = [confidenceMatHigh];
                    output_struct.peptide_data{i,3} = [confidenceMatLow];

                    peptideSequence = peptideLength(1:2:end,:); % Get all residues
                    peptideUptake = peptideLength(2:2:end,:); % Get all uptakes

                    maxUptake = max(max(peptideUptake)); % Calculate the max and min uptake 
                    minUptake = min(min(peptideUptake));

                    % Calculate the most extreme value and append it to the max_uptake_over_all_timepoints
                    extreme = max([maxUptake, abs(minUptake)]);
                    max_uptake_over_all_timepoints = [max_uptake_over_all_timepoints; extreme];

                    peptide_linewidths = ymax/40;

                    if DataType == 1 % binary colouration for deprot/prot peptides

                        % -- For plotting number of each type of peptide in legend
                        % Count number of peptides in each category, put inside
                        % each if switch to avoid accumulating peptides
                        nonsig_peptides = 0;
                        sigpos_peptides = 0;
                        signeg_peptides = 0;

                        % Plot Woods peptides for each timepoint
                        for j=1:2:peptideLengthHeight % nonsignificant peptides under CI, grey
                            %plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',peptide_linewidths,'Color',grey)
                            
                            startRes = nanmin(peptideLength(j,:)); endRes = nanmax(peptideLength(j,:));
                            startUptake = peptideLength(j+1,1)-peptide_linewidths; endUptake = peptideLength(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',grey, 'LineStyle','none')
                            nonsig_peptides = nonsig_peptides+1;
                        end
                        for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                            startRes = nanmin(confidenceMatHigh(j,:)); endRes = nanmax(confidenceMatHigh(j,:));
                            startUptake = confidenceMatHigh(j+1,1)-peptide_linewidths; endUptake = confidenceMatHigh(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',red, 'LineStyle','none')
                            
                            %plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',peptide_linewidths,'Color',red)
                            sigpos_peptides = sigpos_peptides+1;
                        end
                        for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                            startRes = nanmin(confidenceMatLow(j,:)); endRes = nanmax(confidenceMatLow(j,:));
                            startUptake = confidenceMatLow(j+1,1)-peptide_linewidths; endUptake = confidenceMatLow(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',blue, 'LineStyle','none')
                            
                            %plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',peptide_linewidths,'Color',blue);
                            signeg_peptides = signeg_peptides+1;
                        end


                    elseif DataType == 2 % graded coloration for deprot/prot peptides

                        % -- For plotting number of each type of peptide in legend
                        % Count number of peptides in each category, put inside
                        % each if switch to avoid accumulating peptides
                        nonsig_peptides = 0;
                        sigpos_peptides = 0;
                        signeg_peptides = 0;

                        extreme = max([maxUptake, abs(minUptake)])*10000; % x10000 scaling for colour values

                        whiteRed = [linspace(redWhite(1),red(1),extreme+1)', linspace(redWhite(2),red(2),extreme+1)', linspace(redWhite(3),red(3),extreme+1)']; % +1 is added so that rounding max/minUptake doesn't cause bad indexing
                        whiteBlue = [linspace(blueWhite(1),blue(1),extreme+1)', linspace(blueWhite(2),blue(2),extreme+1)', linspace(blueWhite(3),blue(3),extreme+1)'];

                        for j=1:2:peptideLengthHeight % nonsignificant peptides under CI, grey
                            startRes = nanmin(peptideLength(j,:)); endRes = nanmax(peptideLength(j,:));
                            startUptake = peptideLength(j+1,1)-peptide_linewidths; endUptake = peptideLength(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',grey, 'LineStyle','none')
                            
                            %plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',peptide_linewidths,'Color',grey)
                            nonsig_peptides = nonsig_peptides+1;
                        end
                        for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                            scaleColor = int16(confidenceMatHigh(j+1,:)*10000);
                            
                            startRes = nanmin(confidenceMatHigh(j,:)); endRes = nanmax(confidenceMatHigh(j,:));
                            startUptake = confidenceMatHigh(j+1,1)-peptide_linewidths; endUptake = confidenceMatHigh(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',whiteRed(scaleColor(1,1),:), 'LineStyle','none')
                            
                            %plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',peptide_linewidths,'Color',whiteRed(scaleColor(1,1),:))
                            sigpos_peptides = sigpos_peptides+1;
                        end
                        for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                            scaleColor = abs(int16(confidenceMatLow(j+1,:)*10000));
                            
                            startRes = nanmin(confidenceMatLow(j,:)); endRes = nanmax(confidenceMatLow(j,:));
                            startUptake = confidenceMatLow(j+1,1)-peptide_linewidths; endUptake = confidenceMatLow(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',whiteBlue(scaleColor(1,1),:), 'LineStyle','none')
                            %plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',peptide_linewidths,'Color',whiteBlue(scaleColor(1,1),:));
                            signeg_peptides = signeg_peptides+1;
                        end

                        extreme = extreme/10000;
                    end

                    % Number of peptides
                    number_significant_positive = sigpos_peptides;
                    number_significant_negative = signeg_peptides;
                    number_nonsignificant = nonsig_peptides-sigpos_peptides-signeg_peptides;

                    % Plotting invisible markers for the peptide categories
                    h(1) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','white');
                    h(2) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','red');
                    h(3) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','blue');
                    h(4) = scatter(NaN,NaN,'filled','s','MarkerFaceColor',grey);

                    CI_legend = legend([h(2), h(3), h(4)],...
                                        {['Deprotected (',num2str(number_significant_positive),')'],...
                                        ['Protected (',num2str(number_significant_negative),')'],...
                                        ['Not significant (',num2str(number_nonsignificant),')']},...
                                        'location', 'northeast',...
                                        'fontSize',8);
                    legend('boxoff')
                    text(0.02,0.99,['\it',CV_text],'Units','Normalized','FontSize',8); % display the confidence interval

                    % Y-axis options for normal plots
                    ylow = round(-ymax_round-(ymax_round*buffer),0);
                    yhigh = round(ymax_round+(ymax_round*buffer),0);

                    % Ticks
                    tick_increments = round(((yhigh-ylow)/5)*2)/2;
                    ylim([ylow-tick_increments,yhigh+tick_increments])
                    ax = gca; ax.YAxis.TickLabelFormat = '%,.1f';
                    yticks([ylow-tick_increments:tick_increments:yhigh+tick_increments]);
                    
                    curr_timepoint = curr_timepoint+1; % move onto next timepoint
                    
                elseif (curr_timepoint == length(timepoints)) && (enable_sum == 2)

                    p(1,r,c).select() % Select the correct subplot 
                    hold on; grid on; % turn the plot formatting on
                    
                    set(gca,'color',[229/255, 229/255, 229/255]);
                    
                    hAxes = gca;     %Axis handle
                    hAxes.XRuler.Axle.LineStyle = 'none';  
                    hAxes.YRuler.Axle.LineStyle = 'none';
                    hAxes.GridColor = [1,1,1];
                    hAxes.GridAlpha = 1;
                    set(gca,'TickDir','out');
                    
                    title({'';'';['Exposure Time Sum']},'fontweight','bold')
                    
                    % Plot lines for CI and 0 baseline
                    line([seqStart seqEnd],[0 0],'Color','black','LineStyle','-','LineWidth',1) % baseline
                    CV1_pos = dashline([seqStart seqEnd],[confidenceInterval_CV1(i+1) confidenceInterval_CV1(i+1)],1.5,1,1.5,1,'LineWidth',1.5,'Color','black');
                    CV1_neg = dashline([seqStart seqEnd],[-confidenceInterval_CV1(i+1) -confidenceInterval_CV1(i+1)],1.5,1,1.5,1,'LineWidth',1.5,'Color','black'); % CV1 neg
                    CV2_pos = dashline([seqStart seqEnd],[confidenceInterval_CV2(i+1) confidenceInterval_CV2(i+1)],0.5,0.5,0.5,0.5,'LineWidth',1.5,'Color','black'); % CV2 pos
                    CV2_neg = dashline([seqStart seqEnd],[-confidenceInterval_CV2(i+1) -confidenceInterval_CV2(i+1)],0.5,0.5,0.5,0.5,'LineWidth',1.5,'Color','black'); % CV2 neg

                    
                    % Get confidenceFilter depending on 
                    if filter_using == 1
                        confidence_A = choices(A_val);
                        
                        confidenceFilter = confidenceInterval_CV1(i+1);
                        CV_text = strcat(confidence_A,'%: 0±',num2str(round(confidenceInterval_CV1(i+1),2)),' Da');
                    elseif filter_using == 2
                        confidence_B = choices(B_val);
                        
                        confidenceFilter = confidenceInterval_CV2(i+1);
                        CV_text = strcat(confidence_B,'%: 0±',num2str(round(confidenceInterval_CV2(i+1),2)),' Da');
                    end

                    [peptideLengthHeight, peptideLength, confidenceMatHighHeight, confidenceMatHigh, confidenceMatLowHeight, confidenceMatLow] = data_to_Woods(data,confidenceFilter,i);

                    output_struct.peptide_data{i,1} = [peptideLength]; % Add the peptide data to the peptide_data list
                    output_struct.peptide_data{i,2} = [confidenceMatHigh];
                    output_struct.peptide_data{i,3} = [confidenceMatLow];

                    peptideSequence = peptideLength(1:2:end,:); % Get all residues
                    peptideUptake = peptideLength(2:2:end,:); % Get all uptakes

                    maxUptake = max(max(peptideUptake)); % Calculate the max and min uptake 
                    minUptake = min(min(peptideUptake));

                    % Calculate the most extreme value and append it to the max_uptake_over_all_timepoints
                    extreme = max([maxUptake, abs(minUptake)]);
                    max_uptake_over_all_timepoints = [max_uptake_over_all_timepoints; extreme];
                    
                    peptide_linewidths = ymax_sum/50;

                    if DataType == 1

                        % -- For plotting number of each type of peptide in legend
                        % Count number of peptides in each category, put inside
                        % each if switch to avoid accumulating peptides
                        nonsig_peptides = 0;
                        sigpos_peptides = 0;
                        signeg_peptides = 0;

                        % Plot Woods peptides for each timepoint
                        for j=1:2:peptideLengthHeight % nonsignificant peptides under CI, grey
                            %plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',peptide_linewidths,'Color',grey)
                            
                            startRes = nanmin(peptideLength(j,:)); endRes = nanmax(peptideLength(j,:));
                            startUptake = peptideLength(j+1,1)-peptide_linewidths; endUptake = peptideLength(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',grey, 'LineStyle','none')
                            nonsig_peptides = nonsig_peptides+1;
                        end
                        for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                            startRes = nanmin(confidenceMatHigh(j,:)); endRes = nanmax(confidenceMatHigh(j,:));
                            startUptake = confidenceMatHigh(j+1,1)-peptide_linewidths; endUptake = confidenceMatHigh(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',red, 'LineStyle','none')
                            
                            %plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',peptide_linewidths,'Color',red)
                            sigpos_peptides = sigpos_peptides+1;
                        end
                        for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                            startRes = nanmin(confidenceMatLow(j,:)); endRes = nanmax(confidenceMatLow(j,:));
                            startUptake = confidenceMatLow(j+1,1)-peptide_linewidths; endUptake = confidenceMatLow(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',blue, 'LineStyle','none')
                            
                            %plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',peptide_linewidths,'Color',blue);
                            signeg_peptides = signeg_peptides+1;
                        end

                    elseif DataType == 2

                        % -- For plotting number of each type of peptide in legend
                        % Count number of peptides in each category, put inside
                        % each if switch to avoid accumulating peptides
                        nonsig_peptides = 0;
                        sigpos_peptides = 0;
                        signeg_peptides = 0;

                        extreme = max([maxUptake, abs(minUptake)])*10000;

                        % split colour gradients into 10,000x the max uptake value to allow referencing using the uptake value 
                        whiteRed = [linspace(redWhite(1),red(1),extreme+1)', linspace(redWhite(2),red(2),extreme+1)', linspace(redWhite(3),red(3),extreme+1)']; % +1 is added so that rounding max/minUptake doesn't cause bad indexing
                        whiteBlue = [linspace(blueWhite(1),blue(1),extreme+1)', linspace(blueWhite(2),blue(2),extreme+1)', linspace(blueWhite(3),blue(3),extreme+1)'];

                        % Plot sum data for each type of peptide, with colour according to the uptake value
                        for j=1:2:peptideLengthHeight % nonsignificant peptides under CI, grey
                            startRes = nanmin(peptideLength(j,:)); endRes = nanmax(peptideLength(j,:));
                            startUptake = peptideLength(j+1,1)-peptide_linewidths; endUptake = peptideLength(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',grey, 'LineStyle','none')
                            
                            %plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',peptide_linewidths,'Color',grey)
                            nonsig_peptides = nonsig_peptides+1;
                        end
                        for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                            scaleColor = int16(confidenceMatHigh(j+1,:)*10000);
                            
                            startRes = nanmin(confidenceMatHigh(j,:)); endRes = nanmax(confidenceMatHigh(j,:));
                            startUptake = confidenceMatHigh(j+1,1)-peptide_linewidths; endUptake = confidenceMatHigh(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',whiteRed(scaleColor(1,1),:), 'LineStyle','none')
                            
                            %plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',peptide_linewidths,'Color',whiteRed(scaleColor(1,1),:))
                            sigpos_peptides = sigpos_peptides+1;
                        end
                        for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                            scaleColor = abs(int16(confidenceMatLow(j+1,:)*10000));
                            
                            startRes = nanmin(confidenceMatLow(j,:)); endRes = nanmax(confidenceMatLow(j,:));
                            startUptake = confidenceMatLow(j+1,1)-peptide_linewidths; endUptake = confidenceMatLow(j+1,1)+peptide_linewidths;

                            plot(polyshape([startRes startUptake ; endRes startUptake; endRes endUptake ; startRes endUptake]),'FaceAlpha',1,'FaceColor',whiteBlue(scaleColor(1,1),:), 'LineStyle','none')
                            %plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',peptide_linewidths,'Color',whiteBlue(scaleColor(1,1),:));
                            signeg_peptides = signeg_peptides+1;
                        end

                        extreme = extreme/10000;
                    end 

                    % Y-axis options for Sum Plot
                    ylow = round(-ymax_sum_round-(ymax_sum_round*buffer),0);
                    yhigh = round(ymax_sum_round+(ymax_sum_round*buffer),0);

                    % Ticks
                    tick_increments = round(((yhigh-ylow)/7)*2)/2;
                    ax = gca; ax.YAxis.TickLabelFormat = '%,.1f';
                    yticks([ylow:tick_increments:yhigh]);
                    ylim([ylow,yhigh])

                    % Number of peptides
                    number_significant_positive = sigpos_peptides;
                    number_significant_negative = signeg_peptides;
                    number_nonsignificant = nonsig_peptides-sigpos_peptides-signeg_peptides;

                    % Plotting invisible markers for the peptide categories
                    h(1) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','white');
                    h(2) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','red');
                    h(3) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','blue');
                    h(4) = scatter(NaN,NaN,'filled','s','MarkerFaceColor',grey);
                    
                    CI_legend = legend([h(2), h(3), h(4)],...
                                        {['Deprotected (',num2str(number_significant_positive),')'],...
                                        ['Protected (',num2str(number_significant_negative),')'],...
                                        ['Not significant (',num2str(number_nonsignificant),')']},...
                                        'location', 'northeast',...
                                        'fontSize',8);
                    legend('boxoff')
                    
                    text(0.02,0.99,['\it',CV_text],'Units','Normalized','FontSize',8); % display the confidence interval

                    curr_timepoint = curr_timepoint+1; % move onto next timepoint
                    
                end
                
                % Set x-axis, if start residue is less than 50, set to 0
                if seqStart < 50
                    xlim([0 seqEnd])
                else
                    xlim([seqStart seqEnd])
                end
                
                % Axis labelling
                p.xlabel('Residue number');
                p.ylabel('Deuterium uptake (Da)');
                
                xlabh = get(gca,'xlabel');
                set(xlabh,'Units','normalized');
                set(xlabh,'Position',get(xlabh,'Position') - [0 0.005 0]);
                
                ylabh = get(gca,'ylabel');
                set(ylabh,'Units','normalized');
                set(ylabh,'Position',get(ylabh,'Position') - [0.007 0 0]);
                
                xh = get(gca,'xlabel'); % handle to the label object
                ph = get(xh,'position'); % get the current position property
                ph(2) = 4*ph(2) ;        % double the distance, 
                                       % negative values put the label below the axis
                set(xh,'position',ph,'fontweight','bold');   % set the new position
                
                xh = get(gca,'ylabel'); % handle to the label object
                set(xh,'fontweight','bold');   % set the new position

            end
        end


        output_struct.max_uptake_over_all_timepoints = max(max_uptake_over_all_timepoints);

        % Data cursor options
        datacursormode on

        dcm_obj = datacursormode(gcf);
        set(dcm_obj,'UpdateFcn',@linearplotdatacursor,'DisplayStyle','datatip');

        set(handles.Woods_Export,'UserData',output_struct);
        
        % -----------------------------------------
        % Generate peptide list for exporting later
        % -----------------------------------------
        
        filtered_peptides.headers = {'Start','End','Max Uptake','Type'}; % Headings for the csv file
        
        for i=1:length(timepoints) % Extend headers by number of timepoints + sum
            if i ~= length(timepoints)
            filtered_peptides.headers{end+1} = ['Abs DU t=',num2str(timepoints(i+1)),' (Da)'];
            else
                filtered_peptides.headers{end+1} = ['Sum DU (Da)'];
            end
        end
        
        for i=1:length(timepoints) % Extend headers by RFU headers
            if i ~= length(timepoints)
            filtered_peptides.headers{end+1} = ['RFU t=',num2str(timepoints(i+1)),' (%)'];
            else
                filtered_peptides.headers{end+1} = ['Sum RFU (%)'];
            end
        end

        deprot_peptidelists = []; % Contains the start and end residues of each peptide across all timepoints for prot and deprot
        prot_peptidelists = [];

        per_timepoint_peptidelists = {};
        [times, ~] = size(output_struct.peptide_data);
        
        for timepoint = 1:times

            % Deprotected peptides
            deprotected_data = output_struct.peptide_data{timepoint,2};
            [rows, ~] = size(deprotected_data);

            peptidelist = [];
            for row = 1:2:rows
                peptide_row = deprotected_data(row,:);
                peptide_start = peptide_row(1,1);
                peptide_end = max(peptide_row);

                if ~isnan(peptide_start) && ~isnan(peptide_end)
                    startEnd = [peptide_start, peptide_end];
                    deprot_peptidelists = [deprot_peptidelists; startEnd]; % populate the peptide lists
                    peptidelist = [peptidelist; startEnd];
                end
            end
            per_timepoint_peptidelists{timepoint,1} = peptidelist;

            % Protected peptides
            protected_data = output_struct.peptide_data{timepoint,3};
            [rows, ~] = size(protected_data);

            peptidelist = [];
            for row = 1:2:rows
                peptide_row = protected_data(row,:);
                peptide_start = peptide_row(1,1);
                peptide_end = max(peptide_row);

                if ~isnan(peptide_start) && ~isnan(peptide_end)
                    startEnd = [peptide_start, peptide_end];
                    prot_peptidelists = [prot_peptidelists; startEnd];
                    peptidelist = [peptidelist; startEnd];
                end
            end
            per_timepoint_peptidelists{timepoint,2} = peptidelist;

        end
        
        % if peptide in difference_peptide_list is in deprot_peptidelists, extract and append label
        difference_peptide_list = import.diffData.data;

        filtered_peptides.deprot_peptides = [];
        for peptide_row = 1: length(difference_peptide_list)
            peptide = difference_peptide_list(peptide_row,:);
            difference_start_end = peptide(1:2);

            check_match = intersect(difference_start_end, deprot_peptidelists);
            if ~isempty(check_match)

                uptake_sum = sum(peptide(4:end));
                peptide(1,end+1) = uptake_sum;

                for timepoint = 1:length(per_timepoint_peptidelists)
                    if exist('per_timepoint_peptidelists(timepoint,1)','var')
                        peptide_list = cell2mat(per_timepoint_peptidelists(timepoint,1));

                        check_if_match = length(intersect(peptide(1:2),peptide_list));

                        if check_if_match ~= 2
                            peptide(timepoint+3) = 0;
                        end
                    end
                end

                if sum(peptide(4:end)) ~= 0
                    peptide = num2cell(peptide);
                    peptide{3} = 'Deprotected';
                    filtered_peptides.deprot_peptides = [filtered_peptides.deprot_peptides; peptide];
                end
            end
        end
        
        filtered_peptides.prot_peptides = [];
        for peptide_row = 1: length(difference_peptide_list)
            peptide = difference_peptide_list(peptide_row,:);
            difference_start_end = peptide(1:2);

            check_match = intersect(difference_start_end, prot_peptidelists);
            if ~isempty(check_match)
                
                uptake_sum = sum(peptide(4:end)); % calculate sum uptake 
                peptide(1,end+1) = uptake_sum;

                for timepoint = 1:length(per_timepoint_peptidelists)
                    if exist('per_timepoint_peptidelists(timepoint,2)','var')
                        peptide_list = cell2mat(per_timepoint_peptidelists(timepoint,2));

                        check_if_match = length(intersect(peptide(1:2),peptide_list));

                        if check_if_match ~= 2
                            peptide(timepoint+3) = 0;
                        end
                    end
                end

                if sum(peptide(4:end)) ~= 0
                    peptide = num2cell(peptide);
                    peptide{3} = 'Protected';
                    filtered_peptides.prot_peptides = [filtered_peptides.prot_peptides; peptide];
                end
                
            end
        end
        
        % if any uptake is below the confidence filter, change to 0 
        if filter_using == 1
            for i = 1:length(confidenceInterval_CV1)-1
                
                if ~isempty(filtered_peptides.prot_peptides)
                    % Loop over all protected peptides, change all outside of CI filter to 0
                    for j = 1:length(filtered_peptides.prot_peptides(:,3+i))
                        if filtered_peptides.prot_peptides{j,3+i} > -confidenceInterval_CV1(i+1)
                            filtered_peptides.prot_peptides{j,3+i} = 0;
                        end
                    end
                end
                
                % Loop over all deprotected peptides, change all outside of CI filter to 0
                if ~isempty(filtered_peptides.deprot_peptides)
                    for j = 1:length(filtered_peptides.deprot_peptides(:,3+i))
                        if filtered_peptides.deprot_peptides{j,3+i} < confidenceInterval_CV1(i+1)
                            filtered_peptides.deprot_peptides{j,3+i} = 0;
                        end
                    end
                end
            end
        elseif filter_using == 2
            for i = 1:length(confidenceInterval_CV2)-1
                
                % Loop over all protected peptides, change all outside of CI filter to 0
                if ~isempty(filtered_peptides.prot_peptides)
                    for j = 1:length(filtered_peptides.prot_peptides(:,3+i))
                        if filtered_peptides.prot_peptides{j,3+i} > -confidenceInterval_CV2(i+1) || isnan(filtered_peptides.prot_peptides{j,3+i})
                            filtered_peptides.prot_peptides{j,3+i} = 0;
                        end
                    end
                end
                
                % Loop over all deprotected peptides, change all outside of CI filter 2 to 0
                if ~isempty(filtered_peptides.deprot_peptides)
                    for j = 1:length(filtered_peptides.deprot_peptides(:,3+i))
                        if filtered_peptides.deprot_peptides{j,3+i} < confidenceInterval_CV2(i+1) || isnan(filtered_peptides.deprot_peptides{j,3+i})
                            filtered_peptides.deprot_peptides{j,3+i} = 0;
                        end
                    end
                end
            end
        end
        
        % Remove any empty uptake rows from filtered_peptides
        if ~isempty(filtered_peptides.deprot_peptides)
            [rows,~] = size(filtered_peptides.deprot_peptides);
            temp = [];
            for i = 1:rows
                uptake_row = filtered_peptides.deprot_peptides(:,4:end);
                if sum([uptake_row{i,:}]) ~= 0
                    temp = [temp; filtered_peptides.deprot_peptides(i,:)];
                end
            end
            
            filtered_peptides.deprot_peptides = temp; % overwrite original with temp cell array
        end
        
        % repeat for protected peptides
        if ~isempty(filtered_peptides.prot_peptides)
            [rows,~] = size(filtered_peptides.prot_peptides);
            temp = [];
            for i = 1:rows
                uptake_row = filtered_peptides.prot_peptides(:,4:end);
                if sum([uptake_row{i,:}]) ~= 0
                    temp = [temp; filtered_peptides.prot_peptides(i,:)];
                end
            end

            filtered_peptides.prot_peptides = temp; % overwrite original with temp cell array
        end

        % For each peptide in the export data file, get the max uptake and calculate RFU
        compact_peptide_list = import.compact_peptide_list;
        if ~isempty(filtered_peptides.prot_peptides)
            filtered_peptides.prot_peptides_RFU = absolute_to_RFU(filtered_peptides.prot_peptides,compact_peptide_list);
        end
        if ~isempty(filtered_peptides.deprot_peptides)
            filtered_peptides.deprot_peptides_RFU = absolute_to_RFU(filtered_peptides.deprot_peptides,compact_peptide_list);
        end
        
        set(hObject,'UserData',filtered_peptides);
        
        % ----------------------------------------------------------------------------
        % 3.5 Woods plot datatips customisation
        % ----------------------------------------------------------------------------
        
        function txt = linearplotdatacursor(empt,event_obj)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        txt = {['Residue: ',num2str(pos(1))],...
                  ['Uptake: ',num2str(pos(2)), ' Da']};
        
        % ----------------------------------------------------------------------------
        % 3.6 Button for exporting Woods plot
        % ----------------------------------------------------------------------------
        
        function Woods_Export_Callback(hObject, eventdata, handles)
        Diff_path = get(handles.DiffDataPath_browse,'UserData');
        woods_filename = strcat(Diff_path(1:end-4), '_WoodsPlot.pdf');

        p = findobj(gcf,'tag','WoodsUIPanel'); % Find the axes object in the GUI
        export_fig(p,woods_filename)
        
        % ----------------------------------------------------------------------------
        % 3.7 Button for significant peptides
        % ----------------------------------------------------------------------------
        
        function Generate_filtered_peptides_Callback(hObject, eventdata, handles)
        filename = get(handles.DiffDataPath_browse,'UserData');
        filtered_peptides = get(handles.Woods_Plot,'UserData'); % Import filtered peptide lists
        
        if ~isempty(filtered_peptides.deprot_peptides) && ~isempty(filtered_peptides.prot_peptides)
            peptide_output_block = [filtered_peptides.headers; filtered_peptides.deprot_peptides_RFU];
            peptide_output_block = [peptide_output_block; filtered_peptides.prot_peptides_RFU];
            
        elseif isempty(filtered_peptides.deprot_peptides)
            peptide_output_block = [filtered_peptides.headers; filtered_peptides.prot_peptides_RFU];
            
        elseif isempty(filtered_peptides.prot_peptides)
            peptide_output_block = [filtered_peptides.headers; filtered_peptides.deprot_peptides_RFU];
        end
        
        % Export into csv file
        peptideList_filename = strcat(filename(1:end-4), '_peptide_list.csv');
        
        dlmcell(peptideList_filename,peptide_output_block,',');

% --------------------------------------------------------------------------------------------------
% 4. Export Options
% --------------------------------------------------------------------------------------------------

        % ----------------------------------------------------------------------------
        % 4.1 Store chain string for PyMOL script
        % ----------------------------------------------------------------------------
        
        function Pymol_pdbchain_editbox_Callback(hObject, eventdata, handles)
        function Pymol_pdbchain_editbox_CreateFcn(hObject, eventdata, handles)
        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
            set(hObject,'BackgroundColor','white');
        end
        
                % ----------------------------------------------------------------------------
                % 4.2.1 Coverage
                % ----------------------------------------------------------------------------
                
                function Pymol_coverage_checkbox_Callback(hObject, eventdata, handles)
        
                        % ----------------------------------------------------------------------------
                        % 4.2.1.1 Selection of colour for residues with coverage
                        % ----------------------------------------------------------------------------

                        function Pymol_coverage_dropdown_Callback(hObject, eventdata, handles)
                        function Pymol_coverage_dropdown_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end

                        % ----------------------------------------------------------------------------
                        % 4.2.1.2 Selection of colour for residues without coverage
                        % ----------------------------------------------------------------------------

                        function Pymol_nocoverage_dropdown_Callback(hObject, eventdata, handles)
                        function Pymol_nocoverage_dropdown_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end
                
                % ----------------------------------------------------------------------------
                % 4.2.2 Redundancy
                % ----------------------------------------------------------------------------
                
                function Pymol_redundancy_checkbox_Callback(hObject, eventdata, handles)

                        % ----------------------------------------------------------------------------
                        % 4.2.2.1 Selection of which spectrum to use for redundancy
                        % ----------------------------------------------------------------------------

                        function Pymol_spectrum_dropdown_Callback(hObject, eventdata, handles)
                        function Pymol_spectrum_dropdown_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end

                % ----------------------------------------------------------------------------
                % 4.2.3 Absolute uptake
                % ----------------------------------------------------------------------------
                
                function Pymol_uptake_checkbox_Callback(hObject, eventdata, handles)

                        % ----------------------------------------------------------------------------
                        % 4.2.3.1 Type of data to plot: absolute uptake in Da or relative uptake as %
                        % ----------------------------------------------------------------------------
                        
                        function Pymol_uptakedata_dropdown_Callback(hObject, eventdata, handles)
                        function Pymol_uptakedata_dropdown_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end
                        
                        % ----------------------------------------------------------------------------
                        % 4.2.3.2 Scale: select Two, Four or Six colour gradations
                        % ----------------------------------------------------------------------------
                        
                        % Two colours will plot all protected uptake as blue, all deprotected as red
                        
                        % Four colours will find most extreme uptake value (prot or deprot), and apply minor/major
                        % colouration as -100:-50 = dark blue, -50:0 = light blue, 0:50 = light red, 50:100 as dark red
                        % where 100 is most extreme uptake
                        
                        % Six colours will use 1/3 gradations instead of 1/2 as four colours.
                        % -100:-66 = dark blue, -66:-33 = medium blue, -33:0 = light blue
                        % 0:33 = light red, 33:66 = medium red, 66:100 = dark red
                    
                        function Pymol_scale_dropdown_Callback(hObject, eventdata, handles)
                        function Pymol_scale_dropdown_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end

                        % ----------------------------------------------------------------------------
                        % 4.2.3.3 Scale range: select automatic colour scaling or user-defined
                        % ----------------------------------------------------------------------------

                        function Pymol_scalerange_dropdown_Callback(hObject, eventdata, handles)
                            
                        % If scale_range is set to automatic, do nothing, if user-defined, show max_uptake box
                        scale_range = get(handles.Pymol_scalerange_dropdown,'Value');
                        if scale_range == 2
                            set(handles.Pymol_maxuptake_editbox,'Visible','On')
                            set(handles.Pymol_maxuptake_textbox,'Visible','On')
                        else
                            set(handles.Pymol_maxuptake_editbox,'Visible','Off')
                            set(handles.Pymol_maxuptake_textbox,'Visible','Off')
                        end
                            
                        function Pymol_scalerange_dropdown_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end
                        
                        % ----------------------------------------------------------------------------
                        % 4.2.3.4 Max uptake: if 4.2.3.3 is user-defined, stores value for max uptake
                        % ----------------------------------------------------------------------------

                        function Pymol_maxuptake_editbox_Callback(hObject, eventdata, handles)
                        function Pymol_maxuptake_editbox_CreateFcn(hObject, eventdata, handles)
                        if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                            set(hObject,'BackgroundColor','white');
                        end

        % ----------------------------------------------------------------------------
        % 4.3 Button for generating PyMOL script
        % ----------------------------------------------------------------------------
        
        function Button_pymol_gen_Callback(hObject, eventdata, handles)
        filename = get(handles.DiffDataPath_browse,'UserData'); % get filename
        import = get(handles.Import_button,'UserData');
        PDBchain = string(get(handles.Pymol_pdbchain_editbox,'String')); % get chain of PDB from user
        filtered_peptides = get(handles.Woods_Plot,'UserData'); % Import filtered peptide lists
        enable_sum = get(handles.Woods_EnableSum_Pop,'Value'); % Is sum data turned on?
        linearPlot_data = get(handles.Linear_Plot,'UserData'); % get data from linear plot
        
        % Format chain string
        if strlength(PDBchain) > 1
            % Format string as "(chain A | chain B)"
            PDBchain_split = strsplit(PDBchain,',');
            
            list_of_chains = [];
            for i=1:strlength(PDBchain)-floor((strlength(PDBchain)/2))
                if i ~= strlength(PDBchain)-floor((strlength(PDBchain)/2))
                    chain_str = strcat({'chain '},PDBchain_split(i),{' |'});
                    list_of_chains = [list_of_chains chain_str];
                else
                    chain_str = strcat({'chain '},PDBchain_split(i));
                    list_of_chains = [list_of_chains chain_str];
                end
                
                
            end
            
            list_of_chains = strjoin(list_of_chains);
            PDBchain = strcat('(',list_of_chains,')');
        else
            PDBchain = strcat({'chain '},PDBchain);
        end

        
        % What parameters were used for CI calculation, and which CV1/CV2 was applied?
        filter_using = get(handles.Woods_FilterUsing_Pop,'Value');
        filter_choices = get(handles.Woods_FilterUsing_Pop,'String');
        
        filter_alphas = import.filter_alphas; % returns alpha values
        filter_vals = import.filter_vals; % returns choices
        choices = get(handles.FilterA_dropdown, 'String'); % choice of confidence interval filters

        % Get confidenceFilter depending on 
        CV_val_Str = filter_choices(filter_using);
        CV_val = strcat(choices(filter_vals(filter_using)),{'%, '},num2str(filter_alphas(filter_using)));
        
        nReps = string(get(handles.Import_Nreps,'String'));
        
        % Text to go on top of all pml files
        pml_header = {'# Make prettier';
                     'util.perfomance(0);';
                     'as cartoon;';
                     'set cartoon_oval_length, 0.9;';
                     'set opaque_background, off;';
                     'set cartoon_loop_radius, 0.25;';
                     'space pymol;';
                     'bg_color black;';
                     'set spec_reflect, 0;';
                     'set ray_shadow, 0;'};
                 
        % Append log text to each file for user
        Diff_path = string(get(handles.DiffDataPath_browse,'UserData'));
        State_path = string(get(handles.StateDataPath_browse,'UserData'));
        diff_filename_delim = strsplit(Diff_path,'.');
        
        log_datetime = string(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm'));
        log_text = [strjoin(['# File was generated on ', log_datetime],'');
                    strjoin(['#    Difference file: ', Diff_path],'');
                    strjoin(['#    State file: ', State_path],''); ' ';
                    '# Confidence interval parameters ';
                    strjoin(['#    Data filtered using: ', CV_val_Str, ' (', CV_val, ')'],'');
                    strjoin(['#    Number of repeats: ', nReps],'');
                    ' '];
                    
        % Get list of residues with and without coverage
        residues_coverage = []; residues_nocoverage = [];
        for i = 1:length(linearPlot_data.ResidueRedundancy)
            if linearPlot_data.ResidueRedundancy(2,i) > 0
                residues_coverage = [residues_coverage; linearPlot_data.ResidueRedundancy(1,i)];
            else 
                residues_nocoverage = [residues_nocoverage; linearPlot_data.ResidueRedundancy(1,i)];
            end
        end
        
        % Determine which Data type checkboxes are checked
        coverage_checkbox = get(handles.Pymol_coverage_checkbox, 'Value');
        redundancy_checkbox = get(handles.Pymol_redundancy_checkbox, 'Value');
        uptake_checkbox = get(handles.Pymol_uptake_checkbox, 'Value');
        
            % Options for coverage
            % --------------------
            
            % Get colour to use for coverage/no coverage peptides
            solid_colours = {'tv_red','tv_blue','tv_green','brightorange','magenta','teal',...
                             'black','white','grey'};
            
            % Export coverage to pml if checkbox is enabled
            if coverage_checkbox == 1
                coverage_colour_dropdown = solid_colours(get(handles.Pymol_coverage_dropdown, 'Value'));
                nocoverage_colour_dropdown = solid_colours(get(handles.Pymol_nocoverage_dropdown, 'Value'));

                % Base cell array for adding residue colours to
                coverage_textblock = {' '; strjoin(['# Colour ', coverage_colour_dropdown,', residues with coverage'],'')};

                % Append all residues with coverage
                for r = 1:length(residues_coverage)
                    coverage_textblock = [coverage_textblock; strjoin(['color ', coverage_colour_dropdown, ', ', PDBchain, ' & i. ', mat2str(residues_coverage(r)),';'],'')];
                end

                % Append all residues without coverage
                coverage_textblock = [coverage_textblock; ' '; strjoin(['# Colour ', nocoverage_colour_dropdown,', residues with no coverage'],'')];
                for r = 1:length(residues_nocoverage)
                    coverage_textblock = [coverage_textblock; strjoin(['color ', nocoverage_colour_dropdown, ', ', PDBchain, ' & i. ', mat2str(residues_nocoverage(r)),';'],'')];
                end

                % Finally add scene button for coverage
                coverage_textblock = [log_text; pml_header; coverage_textblock; ' '; 'scene Coverage, store;'];

                % Get difference file name and print to pml file
                coverage_filename = strcat(diff_filename_delim(1:end-1), '_coverage.pml');

                coverage_print = fopen(coverage_filename,'w');
                fprintf(coverage_print, '%s\n',coverage_textblock{:});
                fclose(coverage_print);
            end
            
            % Overwrite for grey/white colouration to everything after customisable part
            % Base cell array for adding residue colours to
            coverage_textblock = {' '; '# Colour white, residues with coverage'};

            % Append all residues with coverage

            for r = 1:length(residues_coverage)
                coverage_textblock = [coverage_textblock; strjoin(['color white, ', PDBchain, ' & i. ', mat2str(residues_coverage(r)),';'],'')];
            end

            % Append all residues without coverage
            coverage_textblock = [coverage_textblock; ' '; '# Colour grey, residues with no coverage'];
            for r = 1:length(residues_nocoverage)
                coverage_textblock = [coverage_textblock; strjoin(['color grey, ', PDBchain, ' & i. ', mat2str(residues_nocoverage(r)),';'],'')];
            end
            
            % Options for redundancy
            % ----------------------
            
            % Get spectrum colour
            spectrum_colours = get(handles.Pymol_spectrum_dropdown, 'Value');
            if redundancy_checkbox == 1
                if spectrum_colours == 1
                    spectrum = 'yellow_white_red';
                    color_length = 3;
                elseif spectrum_colours == 2
                    spectrum = 'yellow_white_blue';
                    color_length = 3;
                elseif spectrum_colours == 3
                    spectrum = 'magenta_white_cyan';
                    color_length = 3;
                elseif spectrum_colours == 4
                    spectrum = 'magenta_white_green';
                    color_length = 3;
                elseif spectrum_colours == 5
                    spectrum = 'magenta_white_yellow';
                    color_length = 3;
                elseif spectrum_colours == 6
                    spectrum = 'green_white_magenta';
                    color_length = 3;
                elseif spectrum_colours == 7
                    spectrum = 'red_blue';
                    color_length = 2;
                elseif spectrum_colours == 8
                    spectrum = 'red_cyan';
                    color_length = 2;
                elseif spectrum_colours == 9
                    spectrum = 'red_green';
                    color_length = 2;
                elseif spectrum_colours == 10
                    spectrum = 'red_yellow';
                    color_length = 2;
                elseif spectrum_colours == 11
                    spectrum = 'blue_green';
                    color_length = 2;
                elseif spectrum_colours == 12
                    spectrum = 'blue_magenta';
                    color_length = 2;
                elseif spectrum_colours == 13
                    spectrum = 'blue_red';
                    color_length = 2;
                elseif spectrum_colours == 14
                    spectrum = 'blue_yellow';
                    color_length = 2;
                elseif spectrum_colours == 15
                    spectrum = 'cyan_red';
                    color_length = 2;
                elseif spectrum_colours == 16
                    spectrum = 'cyan_magenta';
                    color_length = 2;
                elseif spectrum_colours == 17
                    spectrum = 'yellow_cyan_white';
                    color_length = 3;
                elseif spectrum_colours == 18
                    spectrum = 'green_blue';
                    color_length = 2;
                elseif spectrum_colours == 19
                    spectrum = 'green_magenta';
                    color_length = 2;
                elseif spectrum_colours == 20
                    spectrum = 'green_red';
                    color_length = 2;
                elseif spectrum_colours == 21
                    spectrum = 'yellow_red';
                    color_length = 2;
                elseif spectrum_colours == 22
                    spectrum = 'yellow_blue';
                    color_length = 2;
                elseif spectrum_colours == 23
                    spectrum = 'yellow_cyan';
                    color_length = 2;
                elseif spectrum_colours == 24
                    spectrum = 'yellow_green';
                    color_length = 2;
                elseif spectrum_colours == 25
                    spectrum = 'yellow_magenta';
                    color_length = 2;
                elseif spectrum_colours == 26
                    spectrum = 'magenta_blue';
                    color_length = 2;
                elseif spectrum_colours == 27
                    spectrum = 'magenta_cyan';
                    color_length = 2;
                elseif spectrum_colours == 28
                    spectrum = 'magenta_green';
                    color_length = 2;
                elseif spectrum_colours == 29
                    spectrum = 'magenta_yellow';
                    color_length = 2;
                elseif spectrum_colours == 30
                    spectrum = 'rainbow';
                    color_length = 1;
                elseif spectrum_colours == 31
                    spectrum = 'rainbow2';
                    color_length = 1;
                elseif spectrum_colours == 32
                    spectrum = 'rainbow2_rev';
                    color_length = 1;
                elseif spectrum_colours == 33
                    spectrum = 'rainbow_cycle';
                    color_length = 1;
                elseif spectrum_colours == 34
                    spectrum = 'rainbow_cycle';
                    color_length = 1;
                elseif spectrum_colours == 35
                    spectrum = 'rainbow_rev';
                    color_length = 1;
                end
                    
                % Create textblock for redundancy colouring
                redundancy_textblock = {' '; strjoin(['# Set b-factor of residue and apply ', string(spectrum),' spectrum'],'')};
                for i = 1: length(linearPlot_data.ResidueRedundancy)
                    redundancy_textblock = [redundancy_textblock; strjoin(['alter i. ', string(linearPlot_data.ResidueRedundancy(1,i)),' & ', PDBchain,', b=', string(linearPlot_data.ResidueRedundancy(2,i)), ';'],'')];
                end
                
                % Add spectrum b line, check colour length first
                if color_length == 1 || color_length == 2
                    redundancy_textblock = [redundancy_textblock; ' '; strjoin(['spectrum b, ', spectrum,', minimum=0, maximum=',string(max(linearPlot_data.ResidueRedundancy(2,:))), ';'],'')];
                else
                    redundancy_textblock = [redundancy_textblock; ' '; strjoin(['spectrum b, ', spectrum,', minimum=,', string(-(max(linearPlot_data.ResidueRedundancy(2,:)))),', maximum=',string(max(linearPlot_data.ResidueRedundancy(2,:))), ';'],'')];
                end
                
                % Finally add scene button for coverage
                redundancy_textblock = {[log_text; pml_header; redundancy_textblock; ' '; 'scene Redundancy, store;']};
                
                % Get difference filename and print to pml file 
                redundancy_filename = strcat(diff_filename_delim(1:end-1), '_redundancy.pml');
                
                redundancy_print = fopen(redundancy_filename,'w');
                fprintf(redundancy_print, '%s\n',redundancy_textblock{:});
                fclose(redundancy_print);
            end
            
            
            
            % Options for uptake
            % ---------------------------
            
            % Get options for dropdowns 
            if uptake_checkbox == 1
                
                % get number of timepoints, default state_timepoints includes 0, so -1 if no SUM, leave if SUM
                if enable_sum == 2
                    number_of_timepoints = length(import.state_timepoints);
                else
                    number_of_timepoints = length(import.state_timepoints)-1;
                end
                
                % Get options from the export uptake dropdown menus
                datatype = get(handles.Pymol_uptakedata_dropdown,'Value'); % 1 = absolute, 2 = relative %
                scale_size = get(handles.Pymol_scale_dropdown, 'Value'); % 1 = two, 2 = four, 3 = six
                scale_range = get(handles.Pymol_scalerange_dropdown,'Value'); % 1 = auto, 2 = user_defined
                
                % Absolute uptake data
                % --------------------
                % filtered_peptides.prot_peptides and filtered_peptides.deprot_peptides
                
                % Relative uptake data
                % --------------------
                % filtered_peptides.prot_peptides_RFU and filtered_peptides.deprot_peptides_RFU
                
                % Controls whether Absolute or RFU data is used for next parts
                if datatype == 1
                    if ~isempty(filtered_peptides.prot_peptides) % check if there are protected peptides
                        protected_peptides = [cell2mat(filtered_peptides.prot_peptides(:,1:2)), cell2mat(filtered_peptides.prot_peptides(:,3+1:end))];
                        
                        if enable_sum == 1 % Delete sum data if user does not need Sum  
                            residue_ranges_prot = protected_peptides(:,1:end-1);
                        else
                            residue_ranges_prot = protected_peptides(:,1:end);
                        end
                    end
                    
                    if ~isempty(filtered_peptides.deprot_peptides) % check if there are deprotected peptides
                        deprotected_peptides = [cell2mat(filtered_peptides.deprot_peptides(:,1:2)), cell2mat(filtered_peptides.deprot_peptides(:,3+1:end))];
                        
                        if enable_sum == 1 % Delete sum data if user does not need Sum  
                            residue_ranges_deprot = deprotected_peptides(:,1:end-1);
                        else
                            residue_ranges_deprot = deprotected_peptides(:,1:end);
                        end
                    end
                    
                    datatype_str = "absolute";
                    
                elseif datatype == 2
                    if enable_sum == 1 % Delete sum data if user does not need Sum 
                        if ~isempty(filtered_peptides.prot_peptides) % check if there are protected peptides
                            protected_peptides = [cell2mat(filtered_peptides.prot_peptides_RFU(:,1:2)), cell2mat(filtered_peptides.prot_peptides_RFU(:,4+number_of_timepoints+1+1:end))];
                            residue_ranges_prot = protected_peptides(:,1:end-1);
                        end
                        
                        if ~isempty(filtered_peptides.deprot_peptides) % check if there are protected peptides
                            deprotected_peptides = [cell2mat(filtered_peptides.deprot_peptides_RFU(:,1:2)), cell2mat(filtered_peptides.deprot_peptides_RFU(:,4+number_of_timepoints+1+1:end))];
                            residue_ranges_deprot = deprotected_peptides(:,1:end-1);
                        end
                        
                        datatype_str = "relative";
                    else
                        if ~isempty(filtered_peptides.prot_peptides) % check if there are protected peptides
                            protected_peptides = [cell2mat(filtered_peptides.prot_peptides_RFU(:,1:2)), cell2mat(filtered_peptides.prot_peptides_RFU(:,4+number_of_timepoints+1:end))];
                            residue_ranges_prot = protected_peptides(:,1:end);
                        end
                        
                        if ~isempty(filtered_peptides.deprot_peptides) % check if there are protected peptides
                            deprotected_peptides = [cell2mat(filtered_peptides.deprot_peptides_RFU(:,1:2)), cell2mat(filtered_peptides.deprot_peptides_RFU(:,4+number_of_timepoints+1:end))];
                            residue_ranges_deprot = deprotected_peptides(:,1:end);
                        end
                        
                        datatype_str = "relative";
                    
                    end
                end
                
                if scale_range == 1 % automatic
                    % Determine most extreme values to use for automatic cutoff
                    if ~isempty(filtered_peptides.deprot_peptides) && ~isempty(filtered_peptides.prot_peptides)
                        highest_val = max([max(max(abs(residue_ranges_prot(:,3:end)))), max(max(abs(residue_ranges_deprot(:,3:end))))]);
                    elseif ~isempty(filtered_peptides.deprot_peptides)
                        highest_val = max(max(abs(residue_ranges_deprot(:,3:end))));
                    elseif ~isempty(filtered_peptides.prot_peptides)
                        highest_val = max(max(abs(residue_ranges_prot(:,3:end)))); % fixed bug here 19102018
                    end
                    
                    scale_str = "automatic_scale";
                    
                elseif scale_range == 2 % manual
                    user_val = get(handles.Pymol_maxuptake_editbox,'String');
                    
                    % determine if user_val is a number, and if +/-
                    num_check = length(str2double(user_val));
                    if num_check == 0
                        errordlg('Max uptake must be float or integer','Invalid input');
                    end
                    
                    % make user-defined value absolute to standardize input
                    highest_val = abs(str2double(user_val));
                    scale_str = "manual_scale";
                end
                
                % pre-calc the half and third of max thresholds
                highest_val_half = highest_val/2;
                highest_val_third = highest_val/3;

                if scale_size == 1 % Two colours

                    two_color_textblock = {' '; '# Colour all significantly protected peptides blue and deprotected peptides red'};
                    
                    if ~isempty(filtered_peptides.deprot_peptides)
                        [rows_d, ~] = size(residue_ranges_deprot);
                    end
                    if ~isempty(filtered_peptides.prot_peptides)
                        [rows_p, ~] = size(residue_ranges_prot);
                    end
                    
                    for t=1:number_of_timepoints
                        
                        if ~isempty(filtered_peptides.deprot_peptides)
                            residue_ranges_deprot = sortrows(residue_ranges_deprot, 2+t,'ascend');
                        end
                        if ~isempty(filtered_peptides.prot_peptides)
                            residue_ranges_prot = sortrows(residue_ranges_prot, 2+t,'descend');
                        end
                    
                        % check if enable sum is on, and if the loop is going over SUM data
                        if (enable_sum == 2) && (t == number_of_timepoints)
                            timepoint = "SUM";
                        else
                            timepoint = string(import.state_timepoints(t+1));
                        end

                        % Add text for which timepoint
                        two_color_textblock = [two_color_textblock; ' '; strjoin(['# Peptides for timepoint ',timepoint],''); coverage_textblock; ' '];

                        %prot text first
                        if ~isempty(filtered_peptides.prot_peptides)
                            for r=1:rows_p
                                if residue_ranges_prot(r,t+2) ~= 0
                                    two_color_textblock = [two_color_textblock; strjoin(['color tv_blue, ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                end
                            end
                        end

                        %deprot text
                        if ~isempty(filtered_peptides.deprot_peptides)
                            for r=1:rows_d
                                if residue_ranges_deprot(r,t+2) ~= 0
                                    two_color_textblock = [two_color_textblock; strjoin(['color tv_red, ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
                                end
                            end
                        end

                        % Add scene button
                        two_color_textblock = [two_color_textblock; strjoin(['scene Timepoint_', timepoint, ', store;'],'')];

                    end

                    % Add log and header
                    two_color_textblock = [log_text; pml_header; two_color_textblock];

                    % Get difference filename and print to pml file 
                    suffix = strjoin(['_uptake_', datatype_str, '_two_colour_', scale_str,'.pml'],'');
                    uptake_filename = strcat(diff_filename_delim(1:end-1), suffix);

                    uptake_print = fopen(uptake_filename,'w');
                    fprintf(uptake_print, '%s\n',two_color_textblock{:});
                    fclose(uptake_print);


                elseif scale_size == 2 % Four colours

                    four_color_textblock = {' '; '# Colour all significantly protected peptides blue and deprotected peptides red';
                                            strjoin(['# Threshold of +/-', string(highest_val_half), ' applied for minor/major differences (Four colours)'],'')};

                    if ~isempty(filtered_peptides.deprot_peptides)
                        [rows_d, ~] = size(residue_ranges_deprot);
                    end
                    if ~isempty(filtered_peptides.prot_peptides)
                        [rows_p, ~] = size(residue_ranges_prot);
                    end

                    for t=1:number_of_timepoints
                        
                        if ~isempty(filtered_peptides.deprot_peptides)
                            residue_ranges_deprot = sortrows(residue_ranges_deprot, 2+t,'ascend')
                        end
                        if ~isempty(filtered_peptides.prot_peptides)
                            residue_ranges_prot = sortrows(residue_ranges_prot, 2+t,'descend');
                        end

                        if (enable_sum == 2) && (t == number_of_timepoints)
                            timepoint = "SUM";
                        else
                            timepoint = string(import.state_timepoints(t+1));
                        end

                        % Add text for which timepoint
                        four_color_textblock = [four_color_textblock; ' '; strjoin(['# Peptides for timepoint ',timepoint],''); coverage_textblock; ' '];

                        %prot text first
                        if ~isempty(filtered_peptides.prot_peptides)
                            for r=1:rows_p
                                uptake = residue_ranges_prot(r,t+2);
                                if (uptake < 0) && (uptake > -highest_val_half) % minor protection
                                    four_color_textblock = [four_color_textblock; strjoin(['color lightblue, ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                elseif (uptake < -highest_val_half) % major protection
                                    four_color_textblock = [four_color_textblock; strjoin(['color tv_blue, ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                end
                            end
                        end

                        % deprotected
                        if ~isempty(filtered_peptides.deprot_peptides)
                            for r=1:rows_d
                                uptake = residue_ranges_deprot(r,t+2);
                                if (uptake > 0) && (uptake < highest_val_half) % minor deprotection
                                    four_color_textblock = [four_color_textblock; strjoin(['color salmon, ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
                                elseif (uptake > highest_val_half) % major protection
                                    four_color_textblock = [four_color_textblock; strjoin(['color tv_red, ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
                                end
                            end
                        end

                        % Add scene button
                        four_color_textblock = [four_color_textblock; strjoin(['scene Timepoint_', timepoint, ', store;'],'')];
                    end

                    % Add log and header
                    four_color_textblock = [log_text; pml_header; four_color_textblock];

                    % Get difference filename and print to pml file 
                    suffix = strjoin(['_uptake_', datatype_str, '_four_colour_', scale_str,'.pml'],'');
                    uptake_filename = strcat(diff_filename_delim(1:end-1), suffix);

                    uptake_print = fopen(uptake_filename,'w');
                    fprintf(uptake_print, '%s\n',four_color_textblock{:});
                    fclose(uptake_print);

                elseif scale_size == 3 % Six colours

                    six_color_textblock = {' '; '# Colour all significantly protected peptides blue and deprotected peptides red';
                                            strjoin(['# Threshold of +/-', string(highest_val_third), ' and ', string(2*highest_val_third), ' applied for minor/medium/major differences (Six colours)'],'')};

                    if ~isempty(filtered_peptides.deprot_peptides)
                        [rows_d, ~] = size(residue_ranges_deprot);
                    end
                    if ~isempty(filtered_peptides.prot_peptides)
                        [rows_p, ~] = size(residue_ranges_prot);
                    end

                    for t=1:number_of_timepoints
                        
                        if ~isempty(filtered_peptides.deprot_peptides)
                            residue_ranges_deprot = sortrows(residue_ranges_deprot, 2+t,'ascend');
                        end
                        if ~isempty(filtered_peptides.prot_peptides)
                            residue_ranges_prot = sortrows(residue_ranges_prot, 2+t,'descend');
                        end

                        if (enable_sum == 2) && (t == number_of_timepoints)
                            timepoint = "SUM";
                        else
                            timepoint = string(import.state_timepoints(t+1));
                        end

                        % Add text for which timepoint
                        six_color_textblock = [six_color_textblock; ' '; strjoin(['# Peptides for timepoint ',timepoint],''); coverage_textblock; ' '];

                        %prot text first
                        if ~isempty(filtered_peptides.prot_peptides)
                            for r=1:rows_p
                                uptake = residue_ranges_prot(r,t+2);
                                if (uptake < 0) && (uptake > -highest_val_third) % minor protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color lightblue, ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];

                                elseif (uptake < -highest_val_third) && (uptake > -2*highest_val_third) % medium protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color teal, ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];

                                elseif (uptake < -2*highest_val_third) % major protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color tv_blue, ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                end
                            end
                        end

                        % deprotected
                        if ~isempty(filtered_peptides.deprot_peptides)
                            for r=1:rows_d
                                uptake = residue_ranges_deprot(r,t+2);
                                if (uptake > 0) && (uptake < highest_val_third) % minor deprotection
                                    six_color_textblock = [six_color_textblock; strjoin(['color salmon, ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];

                                elseif (uptake > highest_val_third) && (uptake < 2*highest_val_third) % medium deprotection
                                    six_color_textblock = [six_color_textblock; strjoin(['color brightorange, ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];

                                elseif (uptake > 2*highest_val_third) % major protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color tv_red, ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
                                end
                            end
                        end

                        % Add scene button
                        six_color_textblock = [six_color_textblock; strjoin(['scene Timepoint_', timepoint, ', store;'],'')];
                    end

                    % Add log and header
                    six_color_textblock = [log_text; pml_header; six_color_textblock];

                    % Get difference filename and print to pml file 
                    suffix = strjoin(['_uptake_', datatype_str, '_six_colour_', scale_str,'.pml'],'');
                    uptake_filename = strcat(diff_filename_delim(1:end-1), suffix);

                    uptake_print = fopen(uptake_filename,'w');
                    fprintf(uptake_print, '%s\n',six_color_textblock{:});
                    fclose(uptake_print);
                end
            end % if uptake_checkbox == 1
                
