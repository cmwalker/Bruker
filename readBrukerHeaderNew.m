function struct = readBrukerHeader(filename)

filename;
% Open the File
fid = fopen(filename, 'r', 'l');
% Read in data to a cell delimited by ## and $$
C = textscan(fid,'%s','Delimiter',{'##','$$'});
tmp = C{1,1}; %Strip the wrapper cell
% Remove empy cell entries
for i = numel(tmp):-1:1
if isempty(tmp{i})
tmp(i) = [];
end
end
for i = 1:numel(tmp)
    [names{i} values{i}] = strtok(tmp{i}, '=');
end
fclose(fid);



filename;
fid = fopen(filename, 'r', 'l');
if (fid==-1)
    error('File not found in readBrukerHeader');
end

rawdata = fread(fid, Inf, 'char=>char')';
temp = rawdata;
index = 1;

while ~isempty(temp)
    [rawitems{index} temp] = strtok(temp, ['##'; '$$']);
     index = index+1;
end

index2=1;

for index = 1:length(rawitems)
    currstring = rawitems{index};
    if(length(currstring)>=5)
            if(strcmp(currstring(1:5), ' @vis'))
                continue;
            end
    end
    if(length(currstring)>=4)
        if(strcmp(currstring(1:4), 'END='))
            break;
        end
    end
    % KAM add for error with PV6 method headers
    if(length(currstring)>=20)
        if(strcmp(currstring(1:20), 'PVM_StartupShimList='))
            break;
        end
    end
    % end add
    [name value] = strtok(currstring, '=');
    if(isempty(value))
        continue;
    end
    %Now, we parse value...
    %There are two possibilities.  Either the value is an array, or it
    %isn't.
    if(strcmp(value, '='))
        eval(['struct.' name '=[];']);
        output.(name) = [];
        continue;
    end
    %Remove all newlines
    value(find(value==sprintf('\n')))='';
    %Excise the newline and =
    value = value(2:end);
    if((length(value)>=2) && strcmp(value(1:2),'( '))
        %Handles array stuff
        [size values] = strtok(value(2:end), ')');
        %Again, excise the paren
        values = values(2:end);
        size = str2num(size);
        if(numel(size)==1)
            size = [size 1];
        end
        num = str2num(values);
        if(~isempty(num))
            cmdstr = ['struct.' name '= reshape(num, size);'];
        else
            cmdstr = '';
            num_parens = length(find(values=='('));
            if isempty(values)
                continue
            end
            if((values(1)=='<' || values(2)=='<') && values(end)=='>')
                if(values(1)=='<')
                    cmdstr = ['struct.' name '=values(2:end-1);'];
                else
                    cmdstr = ['struct.' name '=values(3:end-1);'];
                end
            elseif(num_parens==prod(size))
                for index3=1:size
                    [thisval values] = strtok(values, ')');
                    comp_cmdstr{index3} = ['struct.' name '{' num2str(index3) '} = ''' thisval(3:end) ''';'];
                end
                cmdstr = 'for index3=1:size, eval(comp_cmdstr{index3});, end';
            else
                for index3=1:size
                    [thisval values] = strtok(values);
                    comp_cmdstr{index3} = ['struct.' name '{' num2str(index3) '} = ''' thisval ''';'];
                end
                cmdstr = 'for index3=1:size, eval(comp_cmdstr{index3});, end';
            end
        end
        eval(cmdstr);
    else
        %Single-value stuff
        num = str2num(value);
        if(~isempty(num))
            output.(name) = num;
            cmdstr = ['struct.' name '=' value ';'];
        else
            output.(name) = value;
            cmdstr = ['struct.' name '= ''' value ''';'];
        end
        eval(cmdstr);
    end
end
fclose(fid);
