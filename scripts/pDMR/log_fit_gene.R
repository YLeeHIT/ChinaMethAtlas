fit_genes <- function(gene_list, final_HA, output_path,
                      fit_model = "logistic",
                      log_transform = FALSE,
                      outFlag = TRUE,
                      width = 8,
                      height = 6) 
{

# ==========================================
# Gene-wise Methylation Curve Fitting Function
#
# Arguments:
#   gene_list     : vector of gene names to fit
#   final_HA      : data frame containing methylation data
#   output_path   : output file path for summary results
#   fit_model     : fitting model type ("logistic", "exponential", or "linear")
#   log_transform : whether to apply log-transformation to methylation values
#   outFlag       : whether to remove outliers during fitting
#   width, height : figure size for output plots
#
# Output:
#   A summary table of fitting parameters and model performance.
# ==========================================

    final_result <- data.frame()  # initialize storage

    # Validate model type
    if (!(fit_model %in% c("logistic", "exponential", "linear"))) {
        stop("Invalid fit_model. Choose from 'logistic', 'exponential', or 'linear'.")
    }

    # Iterate over each gene
    for (gene in gene_list) {
        current_flag <- outFlag
        current_fit_coe <- 1
        success <- FALSE

        # Try up to 3 attempts with adaptive parameters
        for (attempt in 1:3) {
            tryCatch({
            # ----- Select fitting model -----
                if (fit_model == "logistic") {
                    result <- fun_fit_methylation_logistic(
                        final_HA = final_HA,
                        gene_name = gene,
                        output_dir = dirname(output_path),
                        Flag = current_flag,
                        fit_coe = current_fit_coe,
                        log_transform = log_transform,
                        width = width,
                        height = height                                                          
                        )
                            
                } else if (fit_model == "exponential") {
                    result <- fun_fit_methylation_exponential(
                        final_HA = final_HA,
                        gene_name = gene,
                        output_dir = dirname(output_path),
                        Flag = current_flag,
                        log_transform = log_transform,
                        start_a = 1,
                        start_b = 0.1,
                        start_c = 0,
                        width = width,
                        height = height                                                                                   
                        )
                            
                } else if (fit_model == "linear") {
                    result <- fun_fit_methylation_linear(
                        final_HA = final_HA,
                        gene_name = gene,
                        output_dir = dirname(output_path),
                        Flag = current_flag,
                        log_transform = log_transform,
                        width = width,
                        height = height                                                                             
                        )    
                }

            # ----- Store success record -----
            result$Status <- "Success"
            final_result <- rbind(final_result, result)
            success <- TRUE
            cat("Gene", gene, "processed successfully on attempt", attempt, "\n")
            break  # exit retry loop on success
                                
        }, error = function(e) {
            # ----- Handle fitting failure -----
            cat("Error in processing gene:", gene, "on attempt", attempt, "\n")

            # Adjust parameters dynamically
            if (attempt == 1) {
                current_flag <- FALSE        # try without outlier removal    
            } else if (attempt == 2) {
                current_fit_coe <- 0.5       # reduce fitting coefficient    
            } else if (attempt == 3) {
                current_fit_coe <- 2         # increase fitting coefficient    
            }     
        })        
    }

    # ----- If still failed after 3 attempts -----
    if (!success) {
        cat("Gene", gene, "failed after 3 attempts.\n")

    failed_record <- switch(
        fit_model,
        "logistic" = data.frame(
            Gene = gene,
            ymin = NA, ymax = NA, k = NA, x0 = NA,
            p_ymin = NA, p_ymax = NA, p_k = NA, p_x0 = NA,
            R2 = NA, F_value = NA, p_value = NA,
            Status = "Failed"
            ),
        "exponential" = data.frame(
            Gene = gene,
            a = NA, b = NA, c = NA,
            p_a = NA, p_b = NA, p_c = NA,
            R2 = NA, F_value = NA, p_value = NA,
            Status = "Failed"                                                             
            ),
        "linear" = data.frame(
            Gene = gene,
            Intercept = NA, p_Intercept = NA,
            Slope = NA, p_Slope = NA,
            R2 = NA, p_Value = NA,
            Status = "Failed"
            )                         
        )
        final_result <- rbind(final_result, failed_record)
        
    }
}

