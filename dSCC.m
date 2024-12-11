function dSCC = dSCC(chromosome_file, alpha, method_type)
    % dSCC: Calculate the dSCC value for a given chromosome
    % INPUTS:
    % chromosome_file: Path to the chromosome contact map (square matrix) file
    % alpha: Scaling factor for freq2dist
    % method_type: ShNeigh method type (1 or 2)
    % OUTPUT:
    % dSCC: Spearman correlation coefficient between original and reconstructed distances

    % Load the chromosome contact map
    heatmapMatrix = load(chromosome_file);

    % Run ShNeigh to get reconstructed coordinates
    [posMatrix, ~] = ShNeigh(heatmapMatrix, method_type);
    reconstructed_coords = posMatrix(:, 2:4);

    % Compute reconstructed distances
    reconstructed_distances = pdist(reconstructed_coords);
    reconstructed_flattened = squareform(reconstructed_distances);
    reconstructed_flattened = reconstructed_flattened(:);

    % Filter the original heatmap to match valid nodes
    selbin = sum(heatmapMatrix > 0) > 2;
    filtered_heatmap = heatmapMatrix(selbin, selbin);
    original_distances = freq2dist(filtered_heatmap, alpha, 1.0);
    original_flattened = original_distances(:);

    % Compute Spearman correlation
    dSCC = corr(original_flattened, reconstructed_flattened, 'Type', 'Spearman');

    % Display the result
    disp(['dSCC value for ', chromosome_file, ': ', num2str(dSCC)]);
end
