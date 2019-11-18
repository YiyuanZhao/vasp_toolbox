clear variables;
%% Read NOSCAR File
dirname = uigetdir('C:\');
if dirname == 0
    f = warndlg('Selection cenceled.','Warning');
    clear variables;
    return
end
filename = [dirname,'\DOSCAR'];
startRow = 6;
formatSpec = '%11s%12s%16s%18s%[^\n\r]';
fileID = fopen(filename,'r');
if fileID == -1
    f = errordlg('NOSCAR file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
% 将非数值文本替换为 NaN。
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4]
    % 将输入元胞数组中的文本转换为数值。已将非数值文本替换为 NaN。
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % 创建正则表达式以检测并删除非数值前缀和后缀。
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % 在非千位位置中检测到逗号。
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % 将数值文本转换为数值。
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 查找非数值元胞
raw(R) = {NaN}; % 替换非数值元胞
DOSCAR = cell2mat(raw);
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

%% Read EIGENVALUE File
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [6, Inf];
opts.Delimiter = " ";
opts.VariableNames = ["order", "energy", "eigenvalue", "Var4"];
opts.SelectedVariableNames = ["order", "energy", "eigenvalue"];
opts.VariableTypes = ["double", "double", "double", "string"];
opts = setvaropts(opts, [1, 2], "TrimNonNumeric", true);
opts = setvaropts(opts, [1, 2], "ThousandsSeparator", ",");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
filename = [dirname,'\EIGENVAL'];
if ~exist(filename,'file')
    f = errordlg('EIGENVAL file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end
EIGENVAL = readtable(filename, opts);
clear opts
%% Read KPOINTS File
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [5, Inf];
opts.Delimiter = [" ", "!"];
opts.VariableNames = ["x", "y", "z", "symmetrydot"];
opts.VariableTypes = ["double", "double", "double", "categorical"];
opts = setvaropts(opts, 3, "TrimNonNumeric", true);
opts = setvaropts(opts, 3, "ThousandsSeparator", ",");
opts = setvaropts(opts, 4, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

filename = [dirname,'\KPOINTS'];
if ~exist(filename,'file')
    f = errordlg('KPOINTS file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end
kpoints = readtable(filename, opts);
kpoints(isnan(kpoints.x),:) = [];
clear opts

opts = delimitedTextImportOptions("NumVariables",1);
opts.DataLines = [2, 2];
opts.Delimiter = " ";
opts.VariableNames = "pointNumber";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

kpointsExtra = readtable(filename, opts);
kpoints.pointNumber = zeros(length(kpoints.x),1);
kpoints.pointNumber(1) = kpointsExtra.pointNumber;
clear opts dirname filename kpointsExtra
%% DOS Plot
fig1 = figure();
plot(DOSCAR(2:end,1)-DOSCAR(1,4),DOSCAR(2:end,2));
set(gca,'XLim',[-8 8]);

%% Band Plot
tmp.kDistance.length = length(kpoints.x)./2;
tmp.plotUpperBound = 10;
tmp.plotLowerBound = -10;
% Generate x-axis scaling parameter
for i = 1: tmp.kDistance.length
    tmp.highSymmetryDotArray(i) = kpoints.symmetrydot(2*i - 1);
    tmp.kDistance.x(i) = abs( kpoints.x(2*i) - kpoints.x(2*i-1));
    tmp.kDistance.y(i) = abs( kpoints.y(2*i) - kpoints.y(2*i-1));
    tmp.kDistance.z(i) = abs( kpoints.z(2*i) - kpoints.z(2*i-1));
    tmp.kDistance.r(i) = sqrt(tmp.kDistance.x(i).^2 + tmp.kDistance.y(i).^2 + tmp.kDistance.z(i)^2);
end

tmp.netot = EIGENVAL.order(1);
tmp.nktot = EIGENVAL.energy(1);
tmp.nbands = EIGENVAL.eigenvalue(1);

tmp.highSymmetryDotArray(tmp.kDistance.length + 1) = kpoints.symmetrydot(end);
tmp.scalingParameter = 2 * tmp.kDistance.r;
tmp.removeLength = length(kpoints.symmetrydot)/2 - 1;
tmp.criticalPoint = calc_criticalPoint(kpoints.pointNumber(1),tmp.nktot,tmp.scalingParameter,tmp.removeLength);

fig2 = figure();
hold on

for i = 1: tmp.nbands
    for j = 1: tmp.nktot
        bandplot.order(j) = calc_order(kpoints.pointNumber(1), tmp.nktot,j, tmp.scalingParameter, tmp.criticalPoint);
        bandplot.energy{i}(j) = EIGENVAL.energy( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) ) - DOSCAR(1, 4);
        bandplot.eigenvalue{i}(j) = EIGENVAL.eigenvalue( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) );
        if mod(j, kpoints.pointNumber(1)) == 0 && j ~=  tmp.nktot
            bandplot.energy{i}(j) = NaN;
            bandplot.eigenvalue{i}(j) = NaN;
        end
    end
    tmp.removePosition = find(isnan(bandplot.energy{i}));
    tmp.removeLength = length(tmp.removePosition);
    bandplot.energy{i}(tmp.removePosition) = [];
    bandplot.eigenvalue{i}(tmp.removePosition) = [];
    bandplot.order(tmp.removePosition) = [];
    plot(bandplot.order(1: end),bandplot.energy{i});
end
% Generate high symmetric point line
for i = 1: length(tmp.scalingParameter)
    bandplot.xSymmetryGrid = [tmp.criticalPoint(i), tmp.criticalPoint(i)];
    bandplot.ySymmetryGrid = [tmp.plotLowerBound, tmp.plotUpperBound];
    plot(bandplot.xSymmetryGrid,bandplot.ySymmetryGrid,'k');
end

hold off
xticks([0, tmp.criticalPoint]);
xticklabels(tmp.highSymmetryDotArray);
axis([0,tmp.criticalPoint(end),tmp.plotLowerBound,tmp.plotUpperBound]);
clear i j tmp
%% Inline Function
function criticalPoint = calc_criticalPoint(pointNumber, nktot, scalingParameter, removeLength)
sum = 0;
criticalPoint = zeros(1,nktot/pointNumber);
for i = 1: length(criticalPoint)
    if i == 1
        criticalPoint(i) = sum + pointNumber * scalingParameter(i);
        sum = criticalPoint(i);
    else
        criticalPoint(i) = sum + (pointNumber - 1) * scalingParameter(i);
        sum = criticalPoint(i);
    end    
end

end

function order = calc_order(pointNumber, nktot, j, scalingParameter, criticalPoint)
levelIndex = floor(j/pointNumber);
if floor(levelIndex) == 0
    order = j * scalingParameter(1);
else
    if mod(j,pointNumber) == 0
        order = criticalPoint(levelIndex);
    else
        order = criticalPoint(levelIndex) + (j - 1 - pointNumber * levelIndex) * scalingParameter(levelIndex + 1);
    end
end

end
