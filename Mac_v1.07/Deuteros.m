function varargout = Deuteros(varargin)
% DEUTEROS MATLAB code for Deuteros.fig
% Last Modified by GUIDE v2.5 13-Jun-2018 13:49:23

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

% Logo options
axes(handles.logo)
matlabImage = imread('Deuteros_logo.png');
image(matlabImage)
axis off
axis image

% Make the export Pymol_maxuptake_textbox and Pymol_maxuptake_editbox objects invisible 
set(handles.Pymol_maxuptake_editbox,'Visible','Off')
set(handles.Pymol_maxuptake_textbox,'Visible','Off')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Deuteros wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Deuteros_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------------------------------------
% 1. Import Difference and State Data
% --------------------------------------------------------------------------------------------------

        % ----------------------------------------------------------------------------
        % 1.1 Path to difference data file
        % ----------------------------------------------------------------------------
        
        function DiffDataPath_browse_Callback(~, Diff_eventdata, handles)
            
        persistent previous_path
        
        % Check if a previous file was already opened, if yes, start uigetfile there
        if exist('previous_path','var') == 1
            defaultFileName = fullfile(previous_path, '*.csv'); % if a previous directory was opened, start there with browse button
            [Diff_filename, Diff_pathname] = uigetfile(defaultFileName,'Select Difference Data ..'); % get full path to the .csv file 
            previous_path = Diff_pathname; % store this as the most recently accessed path
        
        else
            [Diff_filename, Diff_pathname] = uigetfile('*.csv','Select Difference Data ..');
            previous_path = Diff_pathname;
        end

        Diff_fullpathname = strcat(Diff_pathname, Diff_filename);
        Diffpathdisp = set(handles.Import_differenceDataPath,'String',Diff_filename);

        set(handles.Import_textbox,'String',' ');
        set(handles.DiffDataPath_browse,'UserData',Diff_fullpathname);

        % ----------------------------------------------------------------------------
        % 1.2 Path to state data file 
        % ----------------------------------------------------------------------------
        function StateDataPath_browse_Callback(~, State_eventdata, handles)
            
        persistent previous_path
        
        % Check if a previous file was already opened, if yes, start uigetfile there
        if exist('previous_path','var')
            defaultFileName = fullfile(previous_path, '*.csv'); % if a previous directory was opened, start there with browse button
            [State_filename, State_pathname] = uigetfile({'*.csv'},'Select State Data ..');  % get full path to the .csv file 
            previous_path = State_pathname; % store this as the most recently accessed path
        else

            [State_filename, State_pathname] = uigetfile({'*.csv'},'Select State Data ..'); 
            previous_path = State_pathname;
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
        
        % ----------------------------------------------------------------------------
        % 1.6 Stores CV2 value
        % ----------------------------------------------------------------------------

        function Import_alpha2_Callback(hObject, eventdata, handles)
        function Import_alpha2_CreateFcn(hObject, eventdata, handles)
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
        CV1 = str2double(get(handles.Import_alpha1,'String'));
        CV2 = str2double(get(handles.Import_alpha2,'String'));
        Nrep = str2double(get(handles.Import_Nreps,'String'));

        if Nrep == 0
            set(handles.n_errorbox,'String','Cannot have 0 repeats')
        elseif rem(Nrep,1)~=0
            set(handles.n_errorbox,'String','Must be whole number')
        elseif isnan(Nrep) == 0
            set(handles.n_errorbox,'String',' ')
        else
            set(handles.n_errorbox,'String','Enter value of n')
        end

        if isnan(CV1) == 0
            set(handles.CV1_errorbox,'String',' ')
        elseif isnan(CV1) == 1
            set(handles.Import_alpha1,'String','0')
            set(handles.CV1_errorbox,'String',' ')
            CV1 = str2double(get(handles.Import_alpha1,'String'));
        else
            set(handles.CV1_errorbox,'String','Enter value of a1')
        end

        if isnan(CV2) == 0
            set(handles.CV2_errorbox,'String',' ')
        elseif isnan(CV2) == 1
            set(handles.Import_alpha2,'String','0')
            set(handles.CV1_errorbox,'String',' ')
            CV2 = str2double(get(handles.Import_alpha2,'String'));
        else
            set(handles.CV2_errorbox,'String','Enter value of a1')
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
        DataType.val = get(handles.Linear_DataType_Pop, 'Value');

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
        plot([seqStart,seqEnd],[0,0],'LineWidth',25,'Color',[0.85 0.85 0.85]);
        hold on; box off;

        % Axis padding
        axis([seqStart-20 seqEnd+20 -0.01 0.01])

        % Axis decoration
        xlabel('Residue Number')
        set(gca,'fontsize',11)

        % Remove ticks from Y-axis
        set(gca,'YTickLabel',[ ]);
        set(gca,'Ytick',[])
        ax1 = gca;
        yruler = ax1.YRuler;
        yruler.Axle.Visible = 'off';

        % Background colour 
        set(gca,'Color',[0.95 0.95 0.95])
        
        % Set top/bottom buffer for data/grey background
        tb_buffer = 0.0035;

        % ----- Data plot -----

        % Plot depending on which case is chosen on popup menu
        if DataType.val == 1 % Plot Basic Coverage

            for i=1:2:peptideLengthHeight
                plot(peptideLength(i,:),peptideLength(i+1,:)*0, 'linewidth',25,'Color',[(65/255) (105/255) (225/255)])
            end

        elseif DataType.val == 2 % Plot Redundancy
            disp([DataType.str{DataType.val}, ' selected'])

            maxRedundancy = max(ResidueRedundancy(2,:));

            white = [1, 1, 1];
            purple = [0.6, 0, 0.7];
            colourStep = [linspace(white(1),purple(1),maxRedundancy)', linspace(white(2),purple(2),maxRedundancy)', linspace(white(3),purple(3),maxRedundancy)'];

            for i=1:length(ResidueRedundancy)
                if ResidueRedundancy(2,i) ~= 0
                    line([ResidueRedundancy(1,i),ResidueRedundancy(1,i)],[(ResidueRedundancy(2,i)*0)-tb_buffer,(ResidueRedundancy(2,i)*0)+tb_buffer],'Color',colourStep(ResidueRedundancy(2,i),:),'linewidth',4)
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
                    line([ResidueUptake(1,i),ResidueUptake(1,i)],[(ResidueUptake(2,i)*0)-tb_buffer,(ResidueUptake(2,i)*0)+tb_buffer],'Color',whiteRed(colorID,:),'linewidth',4);
                elseif ResidueUptake(2,i) < 0
                    colorID = int64(abs(ResidueUptake(2,i)*10000))+1;
                    line([ResidueUptake(1,i),ResidueUptake(1,i)],[(ResidueUptake(2,i)*0)-tb_buffer,(ResidueUptake(2,i)*0)+tb_buffer],'Color',whiteBlue(colorID,:),'linewidth',4)
                end
            end
        end
        
        % Plot 'Coverage = X%' text
        strmin = ['Coverage = ',num2str(round(coverage,3,'significant')),'%'];
        coverageLabel = text(0.02,0.84,strmin,'Units','normalized');
        coverageLabel.FontSize = 11;
        coverageLabel.FontWeight = 'bold';
        coverageLabel.HorizontalAlignment = 'Left';
        
        % Plot number of peptides text
        number_of_peptides = length(import.diffData.data);
        
        strmin = ['( ',num2str(number_of_peptides),' peptides )'];
        nLabel = text(0.98,0.84,strmin,'Units','normalized');
        nLabel.FontSize = 11;
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
        WoodsUIPanel = uipanel('parent',gcf,'Position',[-8.163265306122449E-4 -8.802816901408451E-4 0.8146938775510204 0.6575704225352113],'BackgroundColor',[1,1,1],'HighlightColor',[1,1,1],'BorderType','None','Clipping','On','tag','WoodsPlotUIPanel');
        import = get(handles.Import_button,'UserData'); % Import the data

        seqStart = str2num(get(handles.Import_sequenceStart,'String')); % Get start and end numbers
        seqEnd = str2num(get(handles.Import_sequenceEnd,'String'));
        
        % ----- Plotting -----
        
        arrayfun(@delete,findobj(gcf,'type','axes','tag','WoodsPlot')); % clear axes each time you press plot

        % colours to use
        grey = [150/255 150/255 150/255];
        white = [255/255 255/255 255/255];
        blueWhite = [220/255 220/255 255/255];
        blue = [20/255 20/255 255/255];
        redWhite = [255/255 220/255 220/255];
        red = [255/255 20/255 20/255];

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
        ymax_sum = max(max(abs(data.data(:,4:end))));
        
        ymax_round = round(ymax/0.5)*0.5;
        buffer = ymax_round*0.1;
        ymax_sum_round = round(ymax_sum/0.5)*0.5;
        buffer_sum = ymax_sum_round*0.1;

        % Get values for other plot options
        enable_sum = get(handles.Woods_EnableSum_Pop,'Value'); % Check if enable sum is enabled
        filter_using = get(handles.Woods_FilterUsing_Pop,'Value'); % Check which filter to use 
        DataType = get(handles.Woods_DataType_Pop,'Value'); % Check data type to plot

        % PyMOL_struct = {} % Set up empty cell to put data in for exporting

        peptide_data = []; % Set up empty matrix to put data in for exporting
        max_uptake_over_all_timepoints = []; % Collect all the max uptakes to PyMOL output scaling later

        for i=1:length(timepoints)

            if i ~= length(timepoints)
                axes('Parent',WoodsUIPanel,'Units', 'normalized', 'Position', [0 0 1 1]);
                subplot(2,round((length(timepoints)/2),0),i,'Parent',WoodsUIPanel,'Tag','WoodsPlot')
                
                % --- Subplot controls --- 
                pos = get(gca, 'Position'); % [x y width height]
                pos(3) = 0.38; % manually change width
                pos(4) = 0.37; % manually change height
                set(gca, 'Position', pos)
                % ------------------------ 
    
                title({'';'';['Exposure Time ' num2str(timepoints(i+1)) ' min']})
                xlabel('Residue Number')
                ylabel('Deuterium Uptake (Da)')
                set(gca,'FontSize',11)
                box on; hold on;

                xlim([seqStart seqEnd]) % Set the x-limit here

                % Plot lines for CI and 0 baseline
                line([seqStart seqEnd],[0 0],'Color','black','LineStyle','-','LineWidth',1) % baseline
                CV1_pos = line([seqStart seqEnd],[confidenceInterval_CV1(i+1) confidenceInterval_CV1(i+1)],'Color','black','LineStyle',':','LineWidth',2.5); % CV1 pos
                CV1_neg = line([seqStart seqEnd],[-confidenceInterval_CV1(i+1) -confidenceInterval_CV1(i+1)],'Color','black','LineStyle',':','LineWidth',2.5); % CV1 neg
                CV2_pos = line([seqStart seqEnd],[confidenceInterval_CV2(i+1) confidenceInterval_CV2(i+1)],'Color','black','LineStyle','--','LineWidth',1.5); % CV2 pos
                CV2_neg = line([seqStart seqEnd],[-confidenceInterval_CV2(i+1) -confidenceInterval_CV2(i+1)],'Color','black','LineStyle','--','LineWidth',1.5); % CV2 neg

                CV1_text = ['CI1: 0±',num2str(round(confidenceInterval_CV1(i+1),2)),' Da'];
                CV2_text = ['CI2: 0±',num2str(round(confidenceInterval_CV2(i+1),2)),' Da'];

                % Get confidenceFilter depending on 
                if filter_using == 1
                    confidenceFilter = confidenceInterval_CV1(i+1);
                elseif filter_using == 2
                    confidenceFilter = confidenceInterval_CV2(i+1);
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

                if DataType == 1 % binary colouration for deprot/prot peptides
                    
                    % -- For plotting number of each type of peptide in legend
                    % Count number of peptides in each category, put inside
                    % each if switch to avoid accumulating peptides
                    nonsig_peptides = 0;
                    sigpos_peptides = 0;
                    signeg_peptides = 0;

                    % Plot Woods peptides for each timepoint
                    for j=1:2:peptideLengthHeight % nonsignificant peptides under CI, grey
                        hold on
                        grid on
                        box on
                        plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',3,'Color',grey)
                        nonsig_peptides = nonsig_peptides+1;
                    end
                    for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                        hold on
                        grid on
                        plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',3,'Color',red)
                        sigpos_peptides = sigpos_peptides+1;
                    end
                    for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                        hold on
                        grid on
                        plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',3,'Color',blue);
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
                        hold on
                        grid on
                        box on
                        plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',3,'Color',grey)
                        nonsig_peptides = nonsig_peptides+1;
                    end
                    for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                        hold on
                        grid on

                        scaleColor = int16(confidenceMatHigh(j+1,:)*10000);
                        plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',3,'Color',whiteRed(scaleColor(1,1),:))
                        sigpos_peptides = sigpos_peptides+1;
                    end
                    for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                        hold on
                        grid on

                        scaleColor = abs(int16(confidenceMatLow(j+1,:)*10000));
                        plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',3,'Color',whiteBlue(scaleColor(1,1),:));
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

                % Axis options
                ylim([-ymax_round-buffer_sum,ymax_round+buffer_sum])
                ax = gca; ax.YAxis.TickLabelFormat = '%,.1f';

                % Plot legends
                CI_leg = legend([CV1_pos, CV2_pos, h(1), h(2), h(3), h(4)],{CV1_text,CV2_text,'',['Deprotected (',num2str(number_significant_positive),')'],['Protected (',num2str(number_significant_negative),')'],['Not significant (',num2str(number_nonsignificant),')']},'Orientation','Vertical','location','northeastoutside');
                title(CI_leg,'Legend')
                legend('boxoff')
                CI_leg.FontSize = 11;

            elseif (i == length(timepoints)) && (enable_sum == 2)

                Woods_ax = axes('Parent',WoodsUIPanel,'Units', 'normalized', 'Position', [0 0 1 1]);
                subplot(2,round((length(timepoints)/2),0),i,'Parent',WoodsUIPanel,'Tag','WoodsPlot')
                
                % --- Subplot controls --- 
                pos = get(gca, 'Position'); % [x y width height]
                pos(3) = 0.38; % manually change width; 0.38
                pos(4) = 0.37; % manually change height
                set(gca, 'Position', pos)
                % ------------------------ 
                
                title({'';'';['Exposure Time Sum']})
                xlabel('Residue Number')
                ylabel('Deuterium Uptake (Da)')
                set(gca,'FontSize',11)
                box on; hold on;

                xlim([seqStart seqEnd]) % Set the x-limit here

                % Plot lines for CI and 0 baseline
                line([seqStart seqEnd],[0 0],'Color','black','LineStyle','-','LineWidth',1) % baseline
                CV1_pos = line([seqStart seqEnd],[confidenceInterval_CV1(i+1) confidenceInterval_CV1(i+1)],'Color','black','LineStyle',':','LineWidth',2); % CV1 pos
                CV1_neg = line([seqStart seqEnd],[-confidenceInterval_CV1(i+1) -confidenceInterval_CV1(i+1)],'Color','black','LineStyle',':','LineWidth',2); % CV1 neg
                CV2_pos = line([seqStart seqEnd],[confidenceInterval_CV2(i+1) confidenceInterval_CV2(i+1)],'Color','black','LineStyle','--','LineWidth',1); % CV2 pos
                CV2_neg = line([seqStart seqEnd],[-confidenceInterval_CV2(i+1) -confidenceInterval_CV2(i+1)],'Color','black','LineStyle','--','LineWidth',1); % CV2 neg

                CV1_text = ['CI1: 0±',num2str(round(confidenceInterval_CV1(i+1),2)),' Da'];
                CV2_text = ['CI2: 0±',num2str(round(confidenceInterval_CV2(i+1),2)),' Da'];

                % Get confidenceFilter depending on 
                if filter_using == 1
                    confidenceFilter = confidenceInterval_CV1(i+1);
                elseif filter_using == 2
                    confidenceFilter = confidenceInterval_CV2(i+1);
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

                if DataType == 1
                    
                    % -- For plotting number of each type of peptide in legend
                    % Count number of peptides in each category, put inside
                    % each if switch to avoid accumulating peptides
                    nonsig_peptides = 0;
                    sigpos_peptides = 0;
                    signeg_peptides = 0;

                    % Plot Woods peptides for each timepoint
                    for j=1:2:peptideLengthHeight % nonsignificant peptides under CI, grey
                        hold on
                        grid on
                        box on
                        plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',3,'Color',grey)
                        nonsig_peptides = nonsig_peptides+1;
                    end
                    for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                        hold on
                        grid on
                        plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',3,'Color',red)
                        sigpos_peptides = sigpos_peptides+1;
                    end
                    for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                        hold on
                        grid on
                        plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',3,'Color',blue);
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
                        hold on
                        grid on
                        box on
                        plot(peptideLength(j,:),peptideLength(j+1,:),'LineWidth',3,'Color',grey)
                        nonsig_peptides = nonsig_peptides+1;
                    end
                    for j=1:2:confidenceMatHighHeight % significant peptides over CI, positive, red
                        hold on
                        grid on

                        scaleColor = int64(confidenceMatHigh(j+1,:)*10000);
                        plot(confidenceMatHigh(j,:),confidenceMatHigh(j+1,:),'LineWidth',3,'Color',whiteRed(scaleColor(1,1),:))
                        sigpos_peptides = sigpos_peptides+1;
                    end
                    for j=1:2:confidenceMatLowHeight % significant peptides over CI, negative, blue
                        hold on
                        grid on

                        scaleColor = abs(int64(confidenceMatLow(j+1,:)*10000));
                        plot(confidenceMatLow(j,:),confidenceMatLow(j+1,:),'LineWidth',3,'Color',whiteBlue(scaleColor(1,1),:));
                        signeg_peptides = signeg_peptides+1;
                    end

                    extreme = extreme/10000;
                end

                % Y-axis options for Sum Plot
                ylim([-ymax_sum_round-buffer_sum,ymax_sum_round+buffer_sum])
                ax = gca; ax.YAxis.TickLabelFormat = '%,.1f';
                
                % Number of peptides
                number_significant_positive = sigpos_peptides;
                number_significant_negative = signeg_peptides;
                number_nonsignificant = nonsig_peptides-sigpos_peptides-signeg_peptides;
                
                % Plotting invisible markers for the peptide categories
                h(1) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','white');
                h(2) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','red');
                h(3) = scatter(NaN,NaN,'filled','s','MarkerFaceColor','blue');
                h(4) = scatter(NaN,NaN,'filled','s','MarkerFaceColor',grey);

                % Plot legends
                CI_leg = legend([CV1_pos, CV2_pos, h(1), h(2), h(3), h(4)],{CV1_text,CV2_text,'',['Deprotected (',num2str(number_significant_positive),')'],['Protected (',num2str(number_significant_negative),')'],['Not significant (',num2str(number_nonsignificant),')']},'Orientation','Vertical','location','northeastoutside','Parent',WoodsUIPanel);
                title(CI_leg,'Legend')
                legend('boxoff')
                CI_leg.FontSize = 11;
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

        p = findobj(gcf,'tag','WoodsPlotUIPanel'); % Find the axes object in the GUI
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
        
        % What parameters were used for CI calculation, and which CV1/CV2 was applied?
        CV_filter = get(handles.Woods_FilterUsing_Pop,'Value');
        
        if CV_filter == 1
            CV_val = string(get(handles.Import_alpha1,'String'));
            CV_val_Str = "CV1";
        elseif CV_filter == 2
            CV_val = string(get(handles.Import_alpha2,'String'));
            CV_val_Str = "CV2";
        end
        
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
                    coverage_textblock = [coverage_textblock; strjoin(['color ', coverage_colour_dropdown, ', chain ', PDBchain, ' & i. ', mat2str(residues_coverage(r)),';'],'')];
                end

                % Append all residues without coverage
                coverage_textblock = [coverage_textblock; ' '; strjoin(['# Colour ', nocoverage_colour_dropdown,', residues with no coverage'],'')];
                for r = 1:length(residues_nocoverage)
                    coverage_textblock = [coverage_textblock; strjoin(['color ', nocoverage_colour_dropdown, ', chain ', PDBchain, ' & i. ', mat2str(residues_nocoverage(r)),';'],'')];
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
                coverage_textblock = [coverage_textblock; strjoin(['color white, chain ', PDBchain, ' & i. ', mat2str(residues_coverage(r)),';'],'')];
            end

            % Append all residues without coverage
            coverage_textblock = [coverage_textblock; ' '; '# Colour grey, residues with no coverage'];
            for r = 1:length(residues_nocoverage)
                coverage_textblock = [coverage_textblock; strjoin(['color grey, chain ', PDBchain, ' & i. ', mat2str(residues_nocoverage(r)),';'],'')];
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
                    redundancy_textblock = [redundancy_textblock; strjoin(['alter i. ', string(linearPlot_data.ResidueRedundancy(1,i)),' & chain ', PDBchain,', b=', string(linearPlot_data.ResidueRedundancy(2,i)), ';'],'')];
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
                        highest_val = max(max(abs(residue_ranges_deprot(:,3:end))));
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
                                    two_color_textblock = [two_color_textblock; strjoin(['color tv_blue, chain ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                end
                            end
                        end

                        %deprot text
                        if ~isempty(filtered_peptides.deprot_peptides)
                            for r=1:rows_d
                                if residue_ranges_deprot(r,t+2) ~= 0
                                    two_color_textblock = [two_color_textblock; strjoin(['color tv_red, chain ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
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
                                    four_color_textblock = [four_color_textblock; strjoin(['color lightblue, chain ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                elseif (uptake < -highest_val_half) % major protection
                                    four_color_textblock = [four_color_textblock; strjoin(['color tv_blue, chain ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                end
                            end
                        end

                        % deprotected
                        if ~isempty(filtered_peptides.deprot_peptides)
                            for r=1:rows_d
                                uptake = residue_ranges_deprot(r,t+2);
                                if (uptake > 0) && (uptake < highest_val_half) % minor deprotection
                                    four_color_textblock = [four_color_textblock; strjoin(['color salmon, chain ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
                                elseif (uptake > highest_val_half) % major protection
                                    four_color_textblock = [four_color_textblock; strjoin(['color tv_red, chain ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
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
                                    six_color_textblock = [six_color_textblock; strjoin(['color lightblue, chain ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];

                                elseif (uptake < -highest_val_third) && (uptake > -2*highest_val_third) % medium protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color teal, chain ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];

                                elseif (uptake < -2*highest_val_third) % major protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color tv_blue, chain ', PDBchain ,' & resi ',string(residue_ranges_prot(r,1)),':', string(residue_ranges_prot(r,2)),';'],'')];
                                end
                            end
                        end

                        % deprotected
                        if ~isempty(filtered_peptides.deprot_peptides)
                            for r=1:rows_d
                                uptake = residue_ranges_deprot(r,t+2);
                                if (uptake > 0) && (uptake < highest_val_third) % minor deprotection
                                    six_color_textblock = [six_color_textblock; strjoin(['color salmon, chain ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];

                                elseif (uptake > highest_val_third) && (uptake < 2*highest_val_third) % medium deprotection
                                    six_color_textblock = [six_color_textblock; strjoin(['color brightorange, chain ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];

                                elseif (uptake > 2*highest_val_third) % major protection
                                    six_color_textblock = [six_color_textblock; strjoin(['color tv_red, chain ', PDBchain ,' & resi ',string(residue_ranges_deprot(r,1)),':', string(residue_ranges_deprot(r,2)),';'],'')];
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
                
