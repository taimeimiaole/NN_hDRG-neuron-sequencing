# title: Source code for figure 5b
# author: Hanying Yan
# date: 2024/08/01


################ 1. Setup R env & load cell boundary annotation data ################
{
  rm(list=ls())
  
  library(rgeos)
  library(ggplot2)
  
  load("Figure5B_raw_data.Rdata")
  # json_drg1-drg4 are list of annotated cells in 4 sections and meta is the metadata for all annotated cells
}

################ 2. Define functions for hypothesis test ################ 
{
  find_centroid_area = function(cell_id, json_sample){
    #'@description The function calculates the centroid, area and radius of a cell
    #'@param cell_id The cell id of a cell
    #'@param json_sample A list of annotated cell boundaries of a sample
    #'@return A vector of extracted cell info 
    
    cell_boundaries =  list(x = json_sample[[cell_id]]$x,y = json_sample[[cell_id]]$y)
    coords = matrix(unlist(cell_boundaries), ncol = 2, byrow = FALSE)
    polygon = SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))
    area = gArea(polygon)
    centroid = coordinates(gCentroid(polygon))
    radius = sqrt(area/pi)
    polygon_cell_id = c(cell_id, centroid, area, radius)
    
    return(polygon_cell_id)
  }
  
  find_centroid_area_all = function(json_sample, meta, cell_type_col){
    
    #'@description The function calculates all cell info from a sample
    #'@param json_sample A list of cell boundaries of a sample
    #'@param meta The meta data for all annotated cells
    #'@param cell_type_col The column in metadata defines the cell type
    #'@return A data frame of extracted cell info
    
    cell_ids = names(json_sample)
    centroid_area_df = sapply(cell_ids, function(x) find_centroid_area(x, json_sample))
    
    #clean the data frame with colnames
    rownames(centroid_area_df) = c("cell_id", 'x', 'y', 'area', 'radius')
    centroid_area_df = data.frame(t(centroid_area_df))
    centroid_area_df$x = as.numeric(centroid_area_df$x)
    centroid_area_df$y = as.numeric(centroid_area_df$y)
    centroid_area_df$point_type = 'centroid'
    
    #add additional cell info from meta
    meta_subset = meta[meta$cell_id %in% cell_ids, c(cell_type_col, 'cell_id')]
    centroid_area_df = merge(centroid_area_df, meta_subset)
    rownames(centroid_area_df) = centroid_area_df$cell_id
    return(centroid_area_df)
  }
  
  cal_cell_type_stats = function(cell_type_sample, cell_type_col, meta_sample, total_area_sample){
    #'@description The function calculates cell statistics of a cell type
    #'@param cell_type_sample The specific cell type
    #'@param cell_type_col The column in metadata defines the cell type
    #'@param meta_sample The meta data for all annotated cells in a sample
    #'@param total_area_sample The total area of a sample
    #'@return A vector of cell stats under the null hypothesis
    
    # Find the cell ids of the specific cell type 
    cell_type_meta_sample = meta_sample[grep(cell_type_sample, meta_sample[, cell_type_col]), ]
    cell_type_cell_ids_sample= cell_type_meta_sample$cell_id
    
    # Calculate unit area and radius of a cell under the null hypothesis 
    n_cell_type_sample = length(cell_type_cell_ids_sample)
    unit_cell_type_area_sample = total_area_sample/n_cell_type_sample
    unit_cell_type_radius_sample = sqrt(unit_cell_type_area_sample/pi)
    
    cell_type_stat_sample = c(n_cell_type_sample, unit_cell_type_area_sample, unit_cell_type_radius_sample)
    return(cell_type_stat_sample)
  }
  
  cal_cell_type_stats_all = function(cell_types_sample, cell_type_col, meta_sample, total_area_sample){
    #'@description The function calculates cell statistics of all cell types in a sample
    #'@param cell_types_sample The unique cell types in a sample
    #'@param cell_type_col The column in metadata defines the cell type
    #'@param meta_sample The meta data for all annotated cells in a sample
    #'@param total_area_sample The total area of a sample
    #'@return A data frame of cell stats for each sample under the null hypothesis
    
    cell_type_stat_sample_df = sapply(cell_types_sample, 
                                      function(x) cal_cell_type_stats(x, cell_type_col, meta_sample, total_area_sample))
    rownames(cell_type_stat_sample_df) = c("n", 'unit_area', 'unit_radius')
    cell_type_stat_sample_df = data.frame(t(cell_type_stat_sample_df))
    return(cell_type_stat_sample_df)
  }
  
  cal_centroid_distance = function(centroid_area_df){
    #'@description The function calculates Euclidian distance between every two cells
    #'@param centroid_area_df A data frame of cell info in a sample
    #'@return A matrix of cell distance for each pair of cells of the same cell type
    
    n=nrow(centroid_area_df)
    dist_matrix = matrix(NA, nrow=n, ncol=n)
    for (i in 1:n){
      for (j in 1:n){
        dist_matrix[i,j] = sqrt((centroid_area_df[i,'x']-centroid_area_df[j,'x'])^2 +
                                  (centroid_area_df[i,'y']-centroid_area_df[j,'y'])^2)
      }
    }
    
    dist_df = data.frame(dist_matrix)
    colnames(dist_df) = rownames(centroid_area_df)
    return(dist_df)
  }
  
  calc_chisq = function(cell_type_sample, centroid_sample, cell_type_stat_df_sample, diameter, meta_sample, cell_type_col){
    #'@description The function calculates cell statistics of a cell type
    #'@param cell_types_sample The unique cell types in a sample
    #'@param centroid_sample The centroid info all cells in the sample
    #'@param cell_type_stat_df_sample The statistics of cells for each cell type
    #'@param diameter The fixed diameter in micron used to search for neighbor cells
    #'@param meta_sample The meta data for all annotated cells in a sample
    #'@param cell_type_col The column in metadata defines the cell type
    #'@return A data frame of chisq test results
    
    # Calculate centroid distance matrix
    cell_type_centroid_sample = centroid_sample[centroid_sample[,cell_type_col]==cell_type_sample, ]
    cell_type_dist_df_sample = cal_centroid_distance(cell_type_centroid_sample)
    
    # Get the number of observed cells and expected cells in the neighborhood under the null hypothesis
    n_observed_radius = diameter/2
    n_expected_vector = pi * (n_observed_radius^2) / cell_type_stat_df_sample[cell_type_sample, 'unit_area']
    n_observed_vector = colSums(cell_type_dist_df_sample <= n_observed_radius)
    
    # Calculate chisq and p value
    X = (n_observed_vector-n_expected_vector)^2/n_expected_vector
    dof = 1
    pvalue = 1 - pchisq(X, df = dof)
    
    cell_type_chisq_df_sample = data.frame(cell_type = rep(cell_type_sample, length(n_observed_vector)),
                                           cell_id = names(n_observed_vector),
                                           n_observed = n_observed_vector, n_expected = n_expected_vector,
                                           chisq = X, pvalue = pvalue)
    return(cell_type_chisq_df_sample)
  }
  
  calc_chisq_all = function(centroid_sample, cell_type_stat_df_sample, diameter, meta_sample, cell_type_col){
    #'@description The function calculates cell statistics of all cell types in a sample
    #'@param centroid_sample The centroid info all cells in the sample
    #'@param cell_type_stat_df_sample The statistics of cells for each cell type
    #'@param diameter The fixed diameter in micron used to search for neighbor cells
    #'@param meta_sample The meta data for all annotated cells in a sample
    #'@param cell_type_col The column in metadata defines the cell type
    #'@return A data frame of chisq test results of all cell types
    
    cell_types_sample = rownames(cell_type_stat_df_sample)
    cell_type_chisq_list_sample = sapply(cell_types_sample, 
                                         function(x) calc_chisq(x, centroid_sample, cell_type_stat_df_sample, 
                                                                diameter, meta_sample, cell_type_col), simplify = FALSE)
    
    cell_type_chisq_df_sample = do.call(rbind, cell_type_chisq_list_sample)
    return(cell_type_chisq_df_sample)
  }
  
  uniform_dist_chisq_test= function(meta, name, total_area_sample, cell_type_col, json_sample, diameter){
    #'@description The function perform uniform distribution hypothesis test
    #'@param meta The meta data for all annotated cells
    #'@param name The name of the sample
    #'@param total_area_sample The total area of a sample
    #'@param cell_type_col The column in metadata defines the cell type
    #'@param json_sample A list of annotated cell boundaries of a sample
    #'@param diameter The fixed diameter in micron used to search for neighbor cells
    #'@return A list of data frame which contains: centroid info of all cells and
    #'        cells stats, chisq test results for each cell type
    
    # Extract all unique cell types in the sample
    meta_sample = subset(meta, meta$cell_id %in% names(json_sample))
    cell_types_sample = unique(meta_sample[,cell_type_col])
    
    # Calculate info of all cells in the sample
    centroid_sample = find_centroid_area_all(json_sample, meta, cell_type_col)
    
    # Find statistics of cells for each cell type
    cell_type_stat_df_sample = cal_cell_type_stats_all(cell_types_sample, cell_type_col, meta_sample, total_area_sample)
    
    # Perform chisq test and return all results
    cell_type_chisq_df_sample = calc_chisq_all(centroid_sample, cell_type_stat_df_sample, diameter, meta_sample, cell_type_col)
    return(list(centroid=centroid_sample, stats = cell_type_stat_df_sample, chisq = cell_type_chisq_df_sample))
  }
}


################ 3. Define functions to generate figure 5b  ################
{
  plot_group_percentage_cell_type = function(chisq_sample, bw, type){
    #'@description The function generate the p value distribution plot 
    #'@chisq_sample The data frame of p values calculated
    #'@param bw The bin width of hist plot
    #'@param type The specific cell type
    #'@return p-value distribution plot of a cell type

    chisq_sample_cell_type = subset(chisq_sample, cell_type==type)
    total = dim(chisq_sample_cell_type)[1]
    
    # Refine cell type name 
    if (type=="hAd.LTMR"){
      type = expression("hA"*delta*".LTMR")
    }
    if (type=="hAb.LTMR"){
      type = expression("hA"*beta*".LTMR")
    }
    if (type=="hPEP.PIEZOh"){
      type = expression("hPEP.PIEZO"^"h")
    }
    
    # Generate and return p-value distribution plot
    p = ggplot(chisq_sample_cell_type, aes(x=pvalue)) +   
      geom_histogram(aes(y=stat(count)/sum(stat(count)), fill=sample), breaks=seq(0,1,bw)) +
      scale_x_continuous(breaks=seq(0,1,bw)) + 
      labs(title = type, x = "pvalue", y = "Percentage") +
      theme(axis.text.x = element_text(angle = 45), legend.position = "none") + 
      annotate(geom="text", label = paste0("Total: ", total), x = Inf, y = Inf, hjust = 1, vjust = 1) +
      scale_fill_manual(values = c('Section 1' = '#F8766D', 'Section 2' = '#7CAE00',
                                   'Section 3' = '#00BFC4', 'Section 4' = '#C77CFF'))
  }    
  
  plot_group_percentage_cell_type_all = function(chisq_sample, bw, cell_types_order) {
    #'@description The function generate the p value distribution plot for all cell types
    #'@chisq_sample The data frame of p values calculated
    #'@param bw The bin width of hist plot
    #'@param cell_types_order The cell types order used to arrange plots
    #'@return p-value distribution plot of all cell type
    
    test_list = lapply(cell_types_order, function(x) plot_group_percentage_cell_type(chisq_sample, bw, x))
    ps = ggpubr::ggarrange(plotlist=test_list, nrow = 4, ncol = 4,common.legend = TRUE, legend="right")
    return(ps)
  }
}


################ 4. Perform uniform distribution hypothesis test & export plot ################ 
{
  ### total section area calculated by ImageJ
  total_area1 = 7424858
  total_area2 = 4484872
  total_area3 = 3924645
  total_area4 = 3906858
  
  ### apply functions to perform chisq test
  chisq_d500_drg1 = uniform_dist_chisq_test(meta, 'drg1', total_area1, 'cell_type_manual', json_drg1, 500)
  chisq_d500_drg2 = uniform_dist_chisq_test(meta, 'drg2', total_area2, 'cell_type_manual', json_drg2, 500)
  chisq_d500_drg3 = uniform_dist_chisq_test(meta, 'drg3', total_area3, 'cell_type_manual', json_drg3, 500)
  chisq_d500_drg4 = uniform_dist_chisq_test(meta, 'drg4', total_area4, 'cell_type_manual', json_drg4, 500)
  
  ### Combine chisq test results from all 4 samples and exclude cell type 'unidentified'
  chisq_d500_drgs = rbind(chisq_d500_drg1$chisq, chisq_d500_drg2$chisq, chisq_d500_drg3$chisq, chisq_d500_drg4$chisq)
  chisq_d500_drgs$sample = c(rep('Section 1', 436), rep('Section 2', 375), rep('Section 3', 298), rep('Section 4', 306))
  chisq_d500_drgs_no_unidentified = subset(chisq_d500_drgs, cell_type != 'unidentified')
  
  ### Generate plot
  cell_types_order = names(table(chisq_d500_drgs_no_unidentified$cell_type))[c(16,4,5,6,15,12,13,11,9,8,10,7,2,1,14,3)]
  # provide cell type order to align figure 5a
  ps = plot_group_percentage_cell_type_all(chisq_d500_drgs_no_unidentified, 0.1, cell_types_order)
  ggsave("Figure5B.svg", ps, bg = "white", dpi = 400, width = 300, height = 300, units='mm') 
}

