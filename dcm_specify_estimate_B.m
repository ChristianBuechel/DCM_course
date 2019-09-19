function dcm_specify_estimate

% Interface for stepping the user through creating a DCM
% FORMAT DCM = spm_dcm_specify_ui(SPM,xY)
%
% SPM      - SPM structure from SPM.mat
% xY       - (optional) VOI structures to be inserted into the DCM
%
% DCM      - DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify_ui.m 6399 2015-04-07 16:05:10Z guillaume $


%==========================================================================
% Outputs
%==========================================================================


hostname = char(getHostName(java.net.InetAddress.getLocalHost));
switch hostname
    case 'mahapralaya'
        base               = '';
    case 'REM'
        base               = 'c:\Users\buechel\Documents\SPM_Kurse\SPM_2019\';
    otherwise
        error('Only hosts REM and mahapralaya accepted');
end


all_sub   = [1];

ana_name     = 'analysis_2019_DCM';

clear subj;
label        = 'Kurs_full';

%We have those regions
%1 VOI_V1_1
%2 VOI_V5_1
%3 VOI_SPC_1

%and those inputs
%1 photic
%2 motion
%3 attention
region_fnames = {'VOI_V1_1.mat','VOI_V5_1.mat','VOI_SPC_1.mat'};

for g=1:size(region_fnames,2)
    clear a;
    a = load([base filesep ana_name filesep region_fnames{g}]);
    xY(g) = a.xY;
end

use_regions = [1 2 3];

%use_inputs  = {[1], [2], [3],  [2 3]  , [1 2 3]}; %inputs can be grouped and combined
%use_weights = {[1]',[1]',[1]', [1 -1]', [1 1 1]'}; %finally we will do a use_inputs * use_weights so can weight everything flexibly
%use_names   = {'photic','motion','attention','motXatt','phot_mot_att'};

use_inputs  = {[1], [2], [3]}; %inputs can be grouped and combined
use_weights = {[1]',[1]',[1]'}; %finally we will do a use_inputs * use_weights so can weight everything flexibly
use_names   = {'photic','motion','attention'};


m           = size(use_regions,2);

%clear xY;
%load([base filesep volunteer filesep ana_name filesep 'VOIs.mat']);
clear SPM;
load([base filesep ana_name filesep 'SPM.mat']);

%-Get cell array of region structures
%--------------------------------------------------------------------------
% add missing fields in xY
%         for i=1:size(xY,2)
%             xY(i).Sess = 1;
%             xY(i).str  = [num2str(xY(i).spec) 'mm ' xY(i).def];
%             xY(i).X0   = [SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]) SPM.xX.K(xY(i).Sess).X0];%X0 contains hpass filter
%         end

%==========================================================================
% Inputs
%==========================================================================

%-Get (nc) 'causes' or inputs U
%--------------------------------------------------------------------------
nsess  = size(SPM.Sess,2);
U.dt   = SPM.Sess(1).U(1).dt; %assume dt is the same for all sessions
U.name = [];
for i=1:size(use_inputs,2)
    U.name{end + 1} = use_names{i};
end
U.u    = [];
Uu     = [];
end_u  = size(SPM.Sess(1).U(1).u,1);
for i=1:size(use_inputs,2)
    Uu             = [Uu cat(2,SPM.Sess(1).U(use_inputs{i}).u) * use_weights{i}];
end

Uu       = Uu(33:end_u,:);
U.u      = [U.u; Uu];

nc       = size(U.u,2);
%==========================================================================
% Timings
%==========================================================================
%-VOI timings
%--------------------------------------------------------------------------
RT     = SPM.xY.RT;
t0     = SPM.xBF.T0;
t      = SPM.xBF.T;
T0     = RT * t0 / t;
delays = repmat(T0,1,m);

%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------
TE    = 0.04;

%==========================================================================
% Model options
%==========================================================================

options.nonlinear  = 0;
options.two_state  = 0;
options.stochastic = 0;
options.centre     = 0;




%==========================================================================
% Graph connections
%==========================================================================
a      = ones(m,m);
a(3,1) = 0;
a(1,3) = 0;

%V1 <--> V5
%V5 <--> SPC


b        = zeros(m,m,nc);
b(2,1,2) = 1; %motion modulating V1 --> V5

b_non_diag = a - diag(diag(a));
non_zero = find(b_non_diag);

all_b = [];
n_b = [1:size(non_zero,1)];
for bb = n_b
    my_b = nchoosek(non_zero,bb);
    for xx = 1:size(my_b,1)
        b_new = zeros(size(a));
        b_new(my_b(xx,:)) = 1;
        all_b{end+1} = b_new;
    end
end


v1v5  = [0 1 0;1 0 0;0 0 0];
v5spc = [0 0 0;0 0 1;0 1 0];

family.names = {'V1V5','V5SPC','both'};

for models = 1:size(all_b,2)
    
    b(:,:,3) = all_b{models}; % here we only change which connection is modualted by attention
    
    
    family.partition(models) = 0;
    
    if any(any(all_b{models}&v1v5)) 
        family.partition(models) = family.partition(models)+1;
    end
    
    if any(any(all_b{models}&v5spc)) 
        family.partition(models) = family.partition(models)+2;
    end
    
    c      = zeros(m,nc);
    
    c(1,1) = 1;  %photic into V1
    
    d     = zeros(m,m,0);
    
    %==========================================================================
    % Response
    %==========================================================================
    
    %-Response variables & confounds (NB: the data have been whitened)
    %--------------------------------------------------------------------------
    n     = length(use_regions);             % number of regions
    v     = length(xY(1).u);                 % number of time points
    clear Y;
    Y.dt  = SPM.xY.RT;
    Y.X0  = xY(1).X0;
    %Y.X0  = [];
    for i = 1:n
        Y.y(:,i)  = xY(use_regions(i)).u;
        Y.name{i} = xY(use_regions(i)).name;
    end
    
    %-Error precision components (one for each region) - i.i.d. (because of W)
    %--------------------------------------------------------------------------
    Y.Q        = spm_Ce(ones(1,n)*v);
    
    %==========================================================================
    % DCM structure
    %==========================================================================
    clear DCM;
    %-Store all variables in DCM structure
    %--------------------------------------------------------------------------
    DCM.a       = a;
    DCM.b       = b;
    DCM.c       = c;
    DCM.d       = d;
    DCM.U       = U;
    DCM.Y       = Y;
    DCM.xY      = xY(use_regions);
    DCM.v       = v;
    DCM.n       = n;
    DCM.TE      = TE;
    DCM.delays  = delays;
    DCM.options = options;
    
    %[DCM] = spm_dcm_estimate(DCM);
    %save([base filesep ana_name filesep  'DCM_' label '_' sprintf('%3.3d',models) '.mat'],'DCM');
end
save(['family.mat'],'family');
