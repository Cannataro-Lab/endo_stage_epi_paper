
# 
# compound_obj <- compound
# cesa <- cesa_samples_by_groups
# stages_name = "Stage"

variant_prevalence_function <- function(compound_obj, cesa, stages_name = "Stage"){
  # IMPORTANT:
  # assumes every sample has the same coverage
  
  prev_list <- vector(mode = "list", length = length(compound_obj))
  names(prev_list) <- compound_obj$compounds$compound_name
  stages <- unique(cesa$samples[,..stages_name])
  stages <- as.matrix(stages)[,1] # pull out just unique stage names 
  
  samples_df <- as.data.frame(cesa$samples)
  
  for(compound_name in names(prev_list)){
    
    this_compound_df <- data.frame(compound = compound_name, stages, tumors_with = rep(NA,length(stages)), 
                                                                   tumors_without= rep(NA,length(stages)),
                                   total_samples = rep(NA,length(stages))) 
    
    samples_with_compound <- compound_obj$samples_with[[compound_name]] 
    
    for(stage_name in this_compound_df$stages){
      
      tumors_in_stage <- samples_df[samples_df[,stages_name] == stage_name,"Unique_Patient_Identifier"]
      
      this_compound_df[this_compound_df$stages == stage_name,"total_samples"] <- length(tumors_in_stage)
      
      this_compound_df[this_compound_df$stages == stage_name,"tumors_with"] <- length(intersect(tumors_in_stage,samples_with_compound))
      
      this_compound_df[this_compound_df$stages == stage_name,"tumors_without"] <- 
        this_compound_df[this_compound_df$stages == stage_name,"total_samples"] - this_compound_df[this_compound_df$stages == stage_name,"tumors_with"]
      
    }
    
    prev_list[[compound_name]] <- this_compound_df
    
  }
  
  return(do.call(what = rbind, args = prev_list))
  
}






