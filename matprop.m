function const = matprop(varargin)
% {pname,pref}, {sname,sref}, {mname,mref}, cref, dbfile
pname = varargin{1}{1};
pref = varargin{1}{2};
sname = varargin{2}{1};
sref = varargin{2}{2};
if (nargin < 4)
    cref = pref;
else
    if isempty(varargin{4})
        cref = pref;
    else
        cref = varargin{4};
    end %if
end %if
if (nargin < 5)
    dbfile = 'material_database.mat';
else
    dbfile = varargin{5};
end %if
if (nargin < 3)
    mname = 'Air';
    mref = 'std';
else
    if isempty(varargin{3})
        mname = 'Air';
        mref = 'std';
    else
        mname = varargin{3}{1};
        mref = varargin{3}{2};
    end %if
end %if
load(dbfile,'-mat','material_prop','contact_prop','medium_prop');
mp = material_properties();
pprop = mp.particle_prop;
sprop = mp.substrate_prop;
cprop = mp.contact_prop;
mprop = mp.medium_prop;
const = [];
%% particle and substrate properties assaign
for i=1:size(material_prop,2)
    iname = strcmpi(material_prop{1,i}(:,1),'Name');
    iref = strcmpi(material_prop{1,i}(:,1),'Reference');
    if (strcmpi(material_prop{1,i}{iname,2},pname)) &&...
            (strcmpi(material_prop{1,i}{iref,2},pref))
        ip = i;
    end %if
    if (strcmpi(material_prop{1,i}{iname,2},sname)) &&...
            (strcmpi(material_prop{1,i}{iref,2},sref))
        is = i;
    end %if
end %i
% particle
for i=1:size(material_prop{1,ip},1)
    vari = find(strcmpi(material_prop{1,ip}{i,1},pprop(:,1)));
    if ~isempty(vari)
        var = pprop{vari,4};
        const.(var) = material_prop{1,ip}{i,2};
    end %if
end %i
% substrate
for i=1:size(material_prop{1,is},1)
    vari = find(strcmpi(material_prop{1,is}{i,1},sprop(:,1)));
    if ~isempty(vari)
        var = sprop{vari,4};
        const.(var) = material_prop{1,is}{i,2};
    end %if
end %i
%% contact
for i=1:size(contact_prop,2)
    ipname = strcmpi(contact_prop{1,i}(:,1),'Particle Name');
    isname = strcmpi(contact_prop{1,i}(:,1),'Substrate Name');
    imname = strcmpi(contact_prop{1,i}(:,1),'Medium Name');
    iref = strcmpi(contact_prop{1,i}(:,1),'Reference');
    if ((strcmpi(contact_prop{1,i}{ipname,2},pname)) &&...
            (strcmpi(contact_prop{1,i}{isname,2},sname)) &&...
            (strcmpi(contact_prop{1,i}{imname,2},mname))) &&...
            (strcmpi(contact_prop{1,i}{iref,2},cref))
        ic = i;
    end %if
end %i
for i=1:size(contact_prop{1,ic},1)
    vari = find(strcmpi(contact_prop{1,ic}{i,1},cprop(:,1)));
    if ~isempty(vari)
        var = cprop{vari,4};
        const.(var) = contact_prop{1,ic}{i,2};
    end %if
end %i
%% medium
for i=1:size(medium_prop,2)
    iname = strcmpi(medium_prop{1,i}(:,1),'Name');
    iref = strcmpi(medium_prop{1,i}(:,1),'Reference');
    if (strcmpi(medium_prop{1,i}{iname,2},mname)) &&...
            (strcmpi(medium_prop{1,i}{iref,2},mref))
        im = i;
    end %if
end %i
for i=1:size(medium_prop{1,im},1)
    vari = find(strcmpi(medium_prop{1,im}{i,1},mprop(:,1)));
    if ~isempty(vari)
        var = mprop{vari,4};
        const.(var) = medium_prop{1,im}{i,2};
    end %if
end %i
%% Check missing properties
icheck = find(isfield(const,pprop(:,4)) == 0);
if ~isempty(icheck)
    disp('The following Particle property is missing:');
    disp(pprop(icheck,1));
end %if
icheck = find(isfield(const,sprop(:,4)) == 0);
if ~isempty(icheck)
    disp('The following Substrate property is missing:');
    disp(sprop(icheck,1));
end %if
icheck = find(isfield(const,mprop(:,4)) == 0);
if ~isempty(icheck)
    disp('The following Medium property is missing:');
    disp(mprop(icheck,1));
end %if
icheck = find(isfield(const,cprop(:,4)) == 0);
if ~isempty(icheck)
    disp('The following Contact property is missing:');
    disp(cprop(icheck,1));
end %if
% %% Calculate Composite Young's Modulus
const.K = comp_modulus(const.E1,const.nu1,const.E2,const.nu2);

