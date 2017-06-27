function rawlist = getrawdata(pth,type)
%GETRAWDATA Digests a full folder structure worth of data
%   Detailed explanation goes here

rawlist = struct('name',{{}},'group',{{{}}},'data',{{}},'program',{{}});

files = dir(pth);
elim = [];
foldr = [];
for i = 1:length(files)
	if files(i).name(1) == '.'
		elim = [elim i];
	elseif isdir(strcat(pth,files(i).name)) == 1
		elim = [elim i];
		foldr = [foldr i];
	elseif files(i).bytes == 0
		elim = [elim i];
	elseif strcmp(files(i).name,'outfile.BbB')
		elim = [elim i];
	end
end
foldr = files(foldr);
files(elim) = [];

ndex = 0;
for i = 1:length(files)
	rawlist.name(ndex+1) = {files(i).name};
    rawlist.group(ndex+1) = {strcat('Data',num2str(i,4))};
    rawlist.program(ndex+1) = {local_MedPCprog(strcat(pth,files(i).name))};
    if strcmp(type,'randycode')
		rawlist.data(ndex+1) = {local_randycode(strcat(pth,files(i).name))};
	elseif strcmp(type,'tempblock')
		rawlist.data(ndex+1) = {local_tempblock(strcat(pth,files(i).name))};
    elseif strcmp(type,'balsammatch')
        rawlist.data(ndex+1) = {local_balsammatch(strcat(pth,files(i).name))};
    end
	ndex = ndex+1;
end

for i = 1:length(foldr)
    out = getrawdata(strcat(pth,foldr(i).name,'/'),type);
    odex = length(out.name);
    for j = 1:odex
        gdex = length(out.group(:,j));
        rawlist.group(1:gdex,ndex+j) = out.group(:,j);
        rawlist.group(gdex+1,ndex+j) = {foldr(i).name};
    end
    rawlist.name(ndex+1:ndex+odex) = out.name;
	rawlist.data(ndex+1:ndex+odex) = out.data;
	rawlist.program(ndex+1:ndex+odex) = out.program;
	ndex = ndex+odex;
end


end

function prog = local_MedPCprog(pth)
%LOCAL_MEDPCPROG Designed to handle generic randycode data
%   Detailed explanation goes here
	fid = fopen(pth);
	q = textscan(fid, '%s');
    fclose(fid);
    ddex = find(ismember(q{1,1}, 'MSN:')==1) + 1;
    prog = q{1,1}(ddex);
end

function data = local_randycode(pth)
%LOCAL_RANDYCODE Designed to handle generic randycode data
%   Detailed explanation goes here
	fid = fopen(pth);
    q = textscan(fid, '%s');
    fclose(fid);
    ddex = find(ismember(q{1,1}, 'W:')==1) + 1;
    d = q{1,1}(ddex:length(q{1,1}));
    for i = 1:length(d)
        if d{i}(length(d{i})) == ':'
            ph(i) = (1==0);
        else
            ph(i) = (1==1);
        end
    end
    ph = find(ph);
    for i = 1:length(ph)
        data(i,1) = str2num(d{ph(i)});
        data(i,2) = mod(data(i,1),10000);
        data(i,1) = (data(i,1)-data(i,2))./10000000;
    end
end

function data = local_tempblock(pth)
%LOCAL_TEMPBLOCK Designed to handle and index the temporal blocking data
%   Detailed explanation goes here
	fid = fopen(pth);
	q = textscan(fid, '%s');
    fclose(fid);
    edex = find(ismember(q{1,1}, 'U:')==1) + 1;
    ddex = find(ismember(q{1,1}, 'V:')==1) + 1;
    e = q{1,1}(edex:ddex-2);
    d = q{1,1}(ddex:length(q{1,1}));
    ph = find((mod(1:length(d),6))~=1);
    for i = 1:length(ph)
        data(i,1) = str2num(d{ph(i)});
        data(i,2) = str2num(e{ph(i)});
    end
end

function data = local_balsammatch(pth)
%LOCAL_BALSAMMATCH Designed to handle and index the matching rat data collected
%in the Balsam Lab 2009-Present
%   Detailed explanation goes here
	fid = fopen(pth);
	q = textscan(fid, '%s');
    fclose(fid);
    ddex = find(ismember(q{1,1}, 'D:')==1) + 1;
    edex = find(ismember(q{1,1}, 'F:')==1) - 1;
    if isempty(edex)
        edex = length(q{1,1});
    end
    d = q{1,1}(ddex:edex);
    resp = find((mod(1:length(d),4))==2);
    rein = find((mod(1:length(d),4))==3);
    rtim = find((mod(1:length(d),4))==0);
    for i = 1:length(resp)
        data(i,1) = str2num(d{resp(i)});
        data(i,2) = str2num(d{rein(i)});
        data(i,3) = str2num(d{rtim(i)});
    end
    try
        if length(data) == 0
            pth
        end
    catch
        pth
    end
end


