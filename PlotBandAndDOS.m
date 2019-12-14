clear variables;
%% Read INCAR File
%Open dialog GUI
dirname = uigetdir('C:\');
%dirname = 'D:\tmp';
if dirname == 0
    f = warndlg('Selection cenceled.','Warning');
    clear variables;
    return
end

filename = [dirname,'\INCAR'];
if ~exist(filename,'file')
    f = errordlg('INCAR file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end

opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [1, Inf];
opts.Delimiter = "=";
opts.VariableNames = ["ALGO", "Fast"];
opts.VariableTypes = ["char", "char"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
INCAR = readtable(filename, opts);
INCAR = table2cell(INCAR);
% Delete empty elements & annotations
numIdx{1} = cellfun('isempty',INCAR(:, 2)) == 1;
INCAR(numIdx{1}, :) = [];
numIdx{2} = ~cellfun ( 'isempty', regexp(INCAR(:,1),"^(!)$") );
INCAR(numIdx{2}, :) = [];
numIdx{3} = ~cellfun ( 'isempty', regexp(INCAR(:,1),"note") );
INCAR(numIdx{3}, :) = [];
for numIdx = 1: length(INCAR)
    INCAR{numIdx,2} = regexprep(INCAR{numIdx,2},'(\s)+!.*','');
end
INCAR(:, 3) = INCAR(:, 2);
numIdx = cellfun(@(x) ~isnan(str2double(x)), INCAR(:, 2));
% Structure of INCAR cell: INCAR(:, 1) = Control tag names, INCAR(:,
% 2) = Value(double format), INCAR(:, 3) = Value (char format)
INCAR(numIdx, 2) = cellfun(@(x) {str2double(x)}, INCAR(numIdx, 2));
% Spin polarized control
if  ~prod( cellfun('isempty', regexp(INCAR(:,1), "(?!\!)^(\s)+MAGMOM") ) ) == 1
    tmp.spinPolarizedFlag = 1;
elseif  ~prod( cellfun('isempty', regexp(INCAR(:,1), "(?!\!)^(\s)+ISPIN") ) ) == 1
    numIdx = ~cellfun ( 'isempty', regexp(INCAR(:,1),"ISPIN") );
    if INCAR{numIdx, 2} == 2
        tmp.spinPolarizedFlag = 1;
    else
        tmp.spinPolarizedFlag = 0;
    end
else
    tmp.spinPolarizedFlag = 0;
end
% Read LORBIT flag (without annotated)
if ~prod( cellfun('isempty', regexp(INCAR(:,1), "(?!\!)^(\s)+LORBIT") ) ) == 1
    numIdx = ~cellfun ( 'isempty', regexp(INCAR(:,1),"LORBIT") );
    tmp.lorbit = INCAR{numIdx, 2};
else
    tmp.lorbit = 0;
end
% Read NEDOS flag (without annotated)
if ~prod( cellfun('isempty', regexp(INCAR(:,1), "(?!\!)^(\s)+NEDOS") ) ) == 1
    numIdx = ~cellfun ( 'isempty', regexp(INCAR(:,1),"NEDOS") );
    tmp.nedos = INCAR{numIdx, 2};
else
    tmp.nedos = 301;
end
clear opts numIdx

%% Read POSCAR File
filename = [dirname,'\POSCAR'];
if ~exist(filename,'file')
    f = errordlg('POSCAR file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end

opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [1, Inf];
opts.Delimiter = ["\t", " "];
opts.VariableNames = ["col1", "col2", "col3", "col4"];
opts.VariableTypes = ["string", "string", "string", "string"];
opts = setvaropts(opts, [1, 2, 3, 4], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
POSCARtemp = readtable(filename, opts);
POSCARtemp = table2array(POSCARtemp);
POSCARtemp(ismissing(POSCARtemp)) = "";
for i = 1:length(POSCARtemp)
    tmp.numIdx = find( strlength(POSCARtemp(i, :)) );
    POSCAR.stringMat(i, 1: length(tmp.numIdx)) = POSCARtemp(i, tmp.numIdx);
end
POSCAR.stringMat(ismissing(POSCAR.stringMat)) = "";
POSCAR.element.name = POSCAR.stringMat(6, :);
POSCAR.element.name( POSCAR.element.name=='' )= [];
POSCAR.element.num = str2double (POSCAR.stringMat(7, :) );
POSCAR.element.num(isnan(POSCAR.element.num)) = [];
clear opts i POSCARtemp

%% Read DOSCAR File
filename = [dirname,'\DOSCAR'];
fileID = fopen(filename,'r');
if fileID == -1
    f = errordlg('NOSCAR file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end

switch tmp.spinPolarizedFlag
    case 0
        opts = delimitedTextImportOptions("NumVariables", 4);
        try
            opts.DataLines = [6, 6+tmp.nedos];
        catch
            opts.DataLines = [6, 307];
        end
        opts.Delimiter = " ";
        opts.VariableNames = ["energy", "DOS", "intDOS", "fermiLevel"];
        opts.VariableTypes = ["double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        DOSCAR.dos = readtable(filename, opts);
        clear opts
    case 1
        opts = delimitedTextImportOptions("NumVariables", 5);
        try
            opts.DataLines = [6, 6+tmp.nedos];
        catch
            opts.DataLines = [6, 307];
        end
        opts.Delimiter = " ";
        opts.VariableNames = ["energy", "DOS_Up", "DOS_Down", "intDOS_Up", "intDOS_Down"];
        opts.VariableTypes = ["double", "double", "double", "double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
        DOSCAR.dos = readtable(filename, opts);
        clear opts
end
DOSCAR.fermiLevel = DOSCAR.dos(1, 4);
DOSCAR.fermiLevel = table2array(DOSCAR.fermiLevel);
DOSCAR.dos(1, :) = [];

switch tmp.lorbit
    case 10
        switch tmp.spinPolarizedFlag
            case 0
                opts = delimitedTextImportOptions("NumVariables", 4);
                opts.DataLines = [tmp.nedos+7, Inf];
                opts.Delimiter = " ";
                opts.VariableNames = ["energy", "sDOS", "pDOS", "dDOS"];
                opts.VariableTypes = ["double", "double", "double", "double"];
                opts.ExtraColumnsRule = "ignore";
                opts.EmptyLineRule = "read";
                opts.ConsecutiveDelimitersRule = "join";
                opts.LeadingDelimitersRule = "ignore";
                DOSCAR.projectedDOStemp = readtable(filename, opts);
                clear opts
            case 1
                opts = delimitedTextImportOptions("NumVariables", 7);
                opts.DataLines = [tmp.nedos+7, Inf];
                opts.Delimiter = " ";
                opts.VariableNames = ["energy", "sDOS_Up", "sDOS_Down", "pDOS_Up", "pDOS_Down", "dDOS_Up", "dDOS_Down"];
                opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];
                opts = setvaropts(opts, [6, 7], "TrimNonNumeric", true);
                opts = setvaropts(opts, [6, 7], "ThousandsSeparator", ",");
                opts.ExtraColumnsRule = "ignore";
                opts.EmptyLineRule = "read";
                opts.ConsecutiveDelimitersRule = "join";
                opts.LeadingDelimitersRule = "ignore";
                DOSCAR.projectedDOStemp = readtable(filename, opts);
                clear opts
        end
    case 11
        switch tmp.spinPolarizedFlag
            case 0
                opts = delimitedTextImportOptions("NumVariables", 10);
                opts.DataLines = [tmp.nedos+7, Inf];
                opts.Delimiter = " ";
                opts.VariableNames = ["energy", "s_DOS", "py_DOS", "pz_DOS", "px_DOS", "dxy_DOS", "dyz_DOS", "dz2_r2__DOS", "dxz_DOS", "dx2_y2_DOS"];
                opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                opts = setvaropts(opts, [4, 5, 6, 7, 8, 9, 10], "TrimNonNumeric", true);
                opts = setvaropts(opts, [4, 5, 6, 7, 8, 9, 10], "ThousandsSeparator", ",");
                opts.ExtraColumnsRule = "ignore";
                opts.EmptyLineRule = "read";
                opts.ConsecutiveDelimitersRule = "join";
                opts.LeadingDelimitersRule = "ignore";
                DOSCAR.projectedDOStemp = readtable(filename, opts);
                clear opts
            case 1
                opts = delimitedTextImportOptions("NumVariables", 19);
                opts.DataLines = [tmp.nedos+7, Inf];
                opts.Delimiter = " ";
                opts.VariableNames = ["energy", "s_DOS_Up", "s_DOS_Down", "py_DOS_Up", "py_DOS_Down", "pz_DOS_Up", "pz_DOS_Down", "px_DOS_Up", "px_DOS_Down", "dxy_DOS_Up", "dxy_DOS_Down", "dyz_DOS_Up", "dyz_DOS_Down", "dz2_r2__DOS_Up", "dz2_r2__DOS_Down", "dxz_DOS_Up", "dxz_DOS_Down", "x2_y2_DOS_Up", "x2_y2_DOS_Down"];
                opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
                opts = setvaropts(opts, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "TrimNonNumeric", true);
                opts = setvaropts(opts, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], "ThousandsSeparator", ",");
                opts.ExtraColumnsRule = "ignore";
                opts.EmptyLineRule = "read";
                opts.ConsecutiveDelimitersRule = "join";
                opts.LeadingDelimitersRule = "ignore";
                DOSCAR.projectedDOStemp = readtable(filename, opts);
                clear opts
        end               
end

if tmp.lorbit ~= 0
    for index = 1: height(DOSCAR.projectedDOStemp) / (tmp.nedos + 1)
    tmp.lowerbound =  (tmp.nedos + 1)*index - 300;
    tmp.upperbound =  (tmp.nedos + 1)*index;
    DOSCAR.projectedDOSforEachAtom{index} = DOSCAR.projectedDOStemp(tmp.lowerbound:tmp.upperbound, :);
    end
end

clear index
%% Read EIGENVALUE File
filename = [dirname,'\EIGENVAL'];
if ~exist(filename,'file')
    f = errordlg('EIGENVAL file not found','File Error');
    set(f,'fontsize',40);
    clear variables;
    return
end

switch tmp.spinPolarizedFlag
    case 0
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
    case 1
        opts = delimitedTextImportOptions("NumVariables", 25);
        opts.DataLines = [6, Inf];
        opts.Delimiter = " ";
        opts.VariableNames = ["order", "energyUp", "energyDown", "eigenvalueUp", "eigenvalueDown", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25"];
        opts.SelectedVariableNames = ["order", "energyUp", "energyDown", "eigenvalueUp", "eigenvalueDown"];
        opts.VariableTypes = ["double", "double", "double", "double", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char"];
        opts = setvaropts(opts, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, [1, 2], "TrimNonNumeric", true);
        opts = setvaropts(opts, [1, 2], "ThousandsSeparator", ",");
        opts = setvaropts(opts, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25], "EmptyFieldRule", "auto");
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts.ConsecutiveDelimitersRule = "join";
        opts.LeadingDelimitersRule = "ignore";
end
EIGENVAL = readtable(filename, opts);
clear opts

%% Read KPOINTS File
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [5, Inf];
opts.Delimiter = [" ", "!", "\t"];
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
% fig1 = figure();
% plot(DOSCAR.dos.energy-DOSCAR.fermiLevel,DOSCAR.dos.DOS_Up,DOSCAR.dos.energy-DOSCAR.fermiLevel,-DOSCAR.dos.DOS_Down);
% set(gca,'XLim',[-8 8]);

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

tmp = getNproperties(tmp,EIGENVAL);
tmp.highSymmetryDotArray(tmp.kDistance.length + 1) = kpoints.symmetrydot(end);
tmp.scalingParameter = 2 * tmp.kDistance.r;
tmp.removeLength = length(kpoints.symmetrydot)/2 - 1;
tmp.criticalPoint = calc_criticalPoint(kpoints.pointNumber(1),tmp.nktot,tmp.scalingParameter,tmp.removeLength);

fig2 = figure();
hold on
switch tmp.spinPolarizedFlag
    case 0
        [bandplot, tmp] = plotBandLineWithoutSpin(tmp,EIGENVAL,kpoints,DOSCAR.fermiLevel);   
    case 1
         [bandplot, tmp] = plotBandLineWithSpin(tmp,EIGENVAL,kpoints,DOSCAR.fermiLevel);
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
% clear i j tmp
%% Inline Function
function tmp = getNproperties(tmp,EIGENVAL)
switch tmp.spinPolarizedFlag
    case 0
        tmp.netot = EIGENVAL.order(1);
        tmp.nktot = EIGENVAL.energy(1);
        tmp.nbands = EIGENVAL.eigenvalue(1);  
    case 1
        tmp.netot = EIGENVAL.order(1);
        tmp.nktot = EIGENVAL.energyUp(1);
        tmp.nbands = EIGENVAL.energyDown(1);  
end

end

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

function [bandplot, tmp] = plotBandLineWithoutSpin(tmp,EIGENVAL,kpoints,fermiLevel)
for i = 1: tmp.nbands
    for j = 1: tmp.nktot
        bandplot.order(j) = calc_order(kpoints.pointNumber(1), tmp.nktot,j, tmp.scalingParameter, tmp.criticalPoint);
        bandplot.energy{i}(j) = EIGENVAL.energy( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) ) - fermiLevel;
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

end

function [bandplot, tmp] = plotBandLineWithSpin(tmp,EIGENVAL,kpoints,fermiLevel)
for i = 1: tmp.nbands
    for j = 1: tmp.nktot
        bandplot.order(j) = calc_order(kpoints.pointNumber(1), tmp.nktot,j, tmp.scalingParameter, tmp.criticalPoint);
        bandplot.energyUp{i}(j) = EIGENVAL.energyUp( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) ) - fermiLevel;
        bandplot.energyDown{i}(j) = EIGENVAL.energyDown( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) ) - fermiLevel;
        bandplot.eigenvalueUp{i}(j) = EIGENVAL.eigenvalueUp( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) );
        bandplot.eigenvalueDown{i}(j) = EIGENVAL.eigenvalueDown( (tmp.nbands + 2) * j + i - (tmp.nbands - 1) );
        if mod(j, kpoints.pointNumber(1)) == 0 && j ~=  tmp.nktot
            bandplot.energyUp{i}(j) = NaN;
            bandplot.energyDown{i}(j) = NaN;
            bandplot.eigenvalueUp{i}(j) = NaN;
            bandplot.eigenvalueDown{i}(j) = NaN;
        end
    end
    tmp.removePosition = find(isnan(bandplot.energyUp{i}));
    tmp.removeLength = length(tmp.removePosition);
    bandplot.energyUp{i}(tmp.removePosition) = [];
    bandplot.energyDown{i}(tmp.removePosition) = [];
    bandplot.eigenvalueUp{i}(tmp.removePosition) = [];
    bandplot.eigenvalueDown{i}(tmp.removePosition) = [];
    bandplot.order(tmp.removePosition) = [];
    plot(bandplot.order(1: end),bandplot.energyUp{i},'r--');
    plot(bandplot.order(1: end),bandplot.energyDown{i},'b');
end
legend({'SpinUp','SpinDown'},'AutoUpdate','off','Location','southeast');
end