function peptide_list_max_uptake = absolute_to_RFU(peptide_list, compact_peptide_list)
% Takes as input a formatted peptide list and converts absolute uptake values to RFU based on the peptide sequence

% Determine number of peptides
[num_peptides, peptide_data] = size(peptide_list);

% Loop over peptides and cross-reference with the max uptake list
peptide_uptakes = [];
for i = 1:num_peptides
    ref = cell2mat(peptide_list(i,1:2));
    
    for j = 1:length(compact_peptide_list)
        comp = compact_peptide_list(j,:);
        
        if ref(1,1) == comp(1,1) && ref(1,2) == comp(1,2)
            peptide_uptakes = [peptide_uptakes; comp];
        end
    end
end
    
% peptide_uptakes contains [seqStart seqEnd maxUptake] for each peptide

% Now make a new cell array and divide each absolute DU with max Uptake for RFU
peptide_RFU = [];
absolutes = peptide_list(:,4:end);

for i = 1:num_peptides
    abs = cell2mat(absolutes(i,:));
    maxUptake = peptide_uptakes(i,3);
    
    peptide_RFU = [peptide_RFU; (abs/maxUptake)*100]; %Store the RFU value
end

% Concatenate the RFU data onto the same peptide_list and return
peptide_list_max_uptake = [peptide_list(:,1:2), num2cell(peptide_uptakes(:,3)), peptide_list(:,3:end), num2cell(peptide_RFU)];