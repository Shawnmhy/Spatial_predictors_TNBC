%% Nuclei Annotation
%% --------------------------------------------------------------
%% Created by: Haoyang Mi - Johns Hopkins - 04/26/2022
%% --------------------------------------------------------------
%% Description:
%%% Convert mask from label matrix to boundaries
%%% Images and masks courtesy of Keren et al at Stanford University
%% --------------------------------------------------------------
%% Input:
%%% wd: working directory should contain the following subfolders
%%%     + Mask
%% --------------------------------------------------------------

%% # list all sample names

% get the sample file names
mydir = pwd;
idcs   = strfind(pwd,'\'); % change '/' to '\' for windows OS
sampledir = mydir(1:idcs(end)-1);
%%
% list all samples names
allSample = dir(fullfile(sampledir, 'Images', '*_VesselMask.tiff'));

% get cell mask specifically

% create directory to save all overlays (sample level)
if (~exist(fullfile(sampledir, 'Vessel_Boundary'), 'dir'))
    mkdir(fullfile(sampledir, 'Vessel_Boundary'))
end


% for loop to go through every sample and list all ROIs

for i = 1:length(allSample)
    
    % make a folder to save overlay images
    % i = 1
    
    % get the real sample name at the corresponding index
    Sample = allSample(i).name;
    
    
            
    % remove the .tiff suffix
    prefix = '.tiff';
    realname = strcat(strrep(Sample,prefix,''), '_boundary');
        
        
    % current TIFF mask
    mask = imread(fullfile(sampledir,'Images', Sample), 1);
    
    % get the boundary mask
    [B,L,N,A] = bwboundaries(mask); 
    %figure; imshow(mask); hold on; 

 % Loop through object boundaries
    
    
    % boundary list
    bdryList_all = [];
    
    if ~isempty(B)
        
        % initialize
        id = 1;
        ring = 1;

        comb = {};

        for m = 1:N
            % Boundary m is the parent of a hole if the k-th column 
            % of the adjacency matrix A contains a non-zero element 
            
                         
            % parent object
            [patr, patc] = size(B{m});
                
            % draw id    

            ring_vec_pat = ring * ones(patr, 1);
                
            % hole indicator
            hole_ind_pat = zeros(patr, 1);
            
            % there could be no hole
            pat_hole_merge = horzcat(B{m}, ring_vec_pat, hole_ind_pat);
            
            % if this is a polygon with hole
            if (nnz(A(:,m)) > 0) 
                % Loop through the children of boundary m and combine them
                % to parent
   
                
                
                for l = find(A(:,m))' 
                    
                    % draw id, parent is 1, holes start from 2
                    ring = ring + 1;

                    [subr, subc] = size(B{l});
              
                    %Add a new field to contain the original position as I need this for  
                    %drawing the orginal borders (which are correct)
                
                    ring_vec = ring * ones(subr, 1);
                    %B{l} = {};
                                   
                    % add converted child boundary to parent in a for loop
                    %fliped = flip(B{l}, 1);
                    %parent = [parent;fliped];
                    
                           
                    % hole indicator
                    hole_ind_child = ones(subr, 1);
                    
                    % flip the child boundary
                    fliped = flip(B{l}, 1);
                    
                    hole_bdry = horzcat(fliped, ring_vec, hole_ind_child);
                    
                    pat_hole_merge = vertcat(pat_hole_merge, hole_bdry);
                    
   
                end
                

                comb{id} = pat_hole_merge;
                id = id + 1;
                ring = ring + 1;
            
            % if this is a polygon without hole
            else
                comb{id} = pat_hole_merge;
                ring = ring + 1;
                id = id + 1;

            end
        end
    
        % transpose
        comb = comb';
        
        for idx = 1:length(comb) % loop through each piece
            % boundary list of each piece
        
            bdlist = comb{idx};
        
            % bdlist size    
            [row, col] = size(bdlist);
    
            % bdlist index vector   
            idx_vec = idx * ones(row, 1);
    
            % append    
            bdryList = horzcat(idx_vec, bdlist(:,3), bdlist(:,2), bdlist(:,1), bdlist(:,4));
    
            % vertical bind
            bdryList_all = vertcat(bdryList_all, bdryList);
        end
        
        namesplit = split(realname, '_');
        
        % new name for morphometric table
        newName_boundary = [namesplit{1}, '_', namesplit{2}, '_', namesplit{3}, '.csv'];
    
        % add name
        bdryList_all = array2table(bdryList_all);

        bdryList_all.Properties.VariableNames = {'id', 'ring', 'x', 'y', 'hole'};

        writetable(bdryList_all,fullfile(sampledir,'Vessel_Boundary', newName_boundary))
    end
end
