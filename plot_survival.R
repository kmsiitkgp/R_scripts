#!/usr/bin/env Rscript

# BiocManager::install("SigCheck")

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")

# Graph plotting packages
library("ggplot2")
library("grid")

# Specialized Graph plotting packages
library("survival")
library("survminer")
library("SigCheck")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

## NOTE: Some predefined columns name MUST be present in the meta/clinical data
## "Sample_id" : this column MUST contain sample names
## "Time:      : this column MUST contain survival duration in days
## "Status     : this column MUST contain Dead/Alive status (Alive = 0, Dead = 1)
## "Sex        : this column MUST contain sex/Sex "Male" or "Female"

## Declare column in metadata to subset samples and its values
## If plotting the whole data, set subset_group <- NA and subset_value <- NA
## If metadata has column named "Race" with "Asian", "African", "American" 
## values and you want to plot only Asian and African origin patients, define
## subset_group <- "Race" and subset_value <- c("Asian", "African").
subset_group <- "Sex"   #"Received.platinum" # "Sample.collected.pre.platinum" 
subset_value <- c("Male", "Female")  #c("Y", "N")  #c("Y", "N", "NE")
#subset_group <- "Stage"   #"Received.platinum" # "Sample.collected.pre.platinum" 
#subset_value <- c("T1", "T2", "T3", "T4")  #c("Y", "N")  #c("Y", "N", "NE")
subset_group <- NA
subset_value <- NA

## "plot_by": define parameter to plot survival.
## To plot by gene expression, set as "Expression" (default). 
## To plot by "Race" etc, a column named "Race" MUST exist in metadata.
plot_by <- "Expression"

## "split_by": define parameter to split data. 
## If you want to plot curves of gene X for male and female patients, then a 
# column named "Sex" MUST exist in metadata.
split_by <- subset_group  #NA

## If you want to plot curves of gene X for male and female patients SEPARATELY,
## then set combine_plot <- FALSE
## If "split_by" <- NA, then combine_plot variable is ignored and its 
## value doesnt matter
combine_plot <- FALSE 

## If plotting by expression, expression cutoffs will be calculated to stratify
## patients into HIGH vs LOW groups. Define if you want to calculate a single 
## cutoff for all patients in subset_groups or individual cutoffs for each
## subset_group
multiple_cutoff <- TRUE

## Choose a stratify_criteria to set cutoff for stratifying samples
## "m" : stratify using median cutoff (top 50% vs bottom 50%)
## "t" : stratify using tertile expression (top 33% vs bottom 33%)
## "q" : stratify using quartile expression (top 25% vs bottom 25%)
## "o" : stratify using optimal cutoff calculated by survminer
## "th": stratify using thirds cutoff
stratify_criteria <- "o"

# Define reference level for calculating Hazard Ratio (HR)
# If you are plotting ONLY by expression, reference is set to "LOW" (default)
# If you are plotting by expression independent favtor like Race or Sex, set
# appropriate reference
reference <- "NA"
reference <- dplyr::if_else(plot_by == "Expression", "LOW", reference)

## Indicate if you want to see confidence intervals in the plot
confidence_interval <- FALSE

## Indicate the title of legend
legend_title <- "Expression"  #Sex #Race

## indicate if you want to plot the risk table
plot_risk_table <- FALSE

## indicate if you want to plot the survival curve
plot_curve <- TRUE  #TRUE

## Indicate if you want to plot all quartiles or only HIGH vs LOW
all_quartiles <- FALSE

## Indicate if you are plotting a gene signature or single gene
gene_signature <- TRUE

# # VERY VERY IMPORTANT: ggsurvplot() labels the groups in alphabetical order. So, 
# # when we want to use custom labels, initialize them in alphabetical order. 
# # Eg: c("High", "Low") instead of  c("Low, "High")
legend_label <- c("Female", "Male")           # Reference is female
color_palette <- c("#F6D2E0", "#C8E7F5")      # 1st color~female, 2nd color~male
color_palette <- c("orchid2", "dodgerblue2")  # 1st color~female, 2nd color~male
color_palette <- c("#EE7AE9", "#1C86EE")      # 1st color~female, 2nd color~male

legend_label <- c("High", "Low")              # Reference is Low
color_palette <- c("#DB6D00", "#490092")      # 1st color~high, 2nd color~low

## Store path of parent directory i.e. root directory for the project
parent_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"

## Declare global variables for survival curves. 
variable_x <- "Time"
variable_y <- "Status" 

#******************************************************************************#
#                STEP 1: CREATE FUNCTIONS TO PLOT SURVIVAL CURVE               #
#******************************************************************************#

wrangle_data <- function(stratify_criteria){
  
  stats <- list("gene" = c(), 
                "group" = c(),
                "lower_cutoff" = c(),
                "middle_cutoff" = c(),
                "upper_cutoff" = c(),
                "HR" = c(), 
                "CI_lower" = c(), 
                "CI_upper" = c(), 
                "pvalue" = c())
  
  # STEP 1: Calculate cutoffs
  # If cutoffs need to be calculated for each group, subset the expr_df and pass
  # it to calculate_cutoffs(). Else, pass entire expr_df to calculate_cutoffs()
  if (multiple_cutoff == TRUE & !is.na(split_by)) {
    
    # Create empty dataframe to store results of calculate_cutoffs() for each group
    survival_data <- data.frame(model = " ")
    
    # Calculate cutoffs for each group in split_by
    for (x in (expr_df %>% dplyr::distinct(get(split_by)))[[1]]) {
      mat <- expr_df %>% 
        dplyr::filter(get(split_by) == x) %>% 
        dplyr::select("Sample_id", "GEO_ID", "Sex", "Time", "Status", all_of(gene))
      mat <- calculate_cutoffs(mat, x)
      
      # Save the data from output of calculate_cutoffs()
      survival_data <- dplyr::bind_rows(survival_data, mat[[1]])
      stats$gene <- c(stats$gene, mat[[2]]$gene)
      stats$group <- c(stats$group, mat[[2]]$group)
      stats$lower_cutoff <- c(stats$lower_cutoff, mat[[2]]$lower)
      stats$middle_cutoff <- c(stats$middle_cutoff, mat[[2]]$middle)
      stats$upper_cutoff <- c(stats$upper_cutoff, mat[[2]]$upper)
    }
    
    # Populate the model variable by concatenating "Expression" and "split_by"
    survival_data <- survival_data %>%
      dplyr::mutate(model = paste0(Expression, "_", get(split_by))) %>%
      dplyr::filter(!is.na(get(split_by)))
  } else{
    mat <- expr_df
    x <- "NA"
    mat <- calculate_cutoffs(mat, x)
    
    # Save the data from output of calculate_cutoffs()
    survival_data <- mat[[1]]
    stats$gene <- c(stats$gene, mat[[2]]$gene)
    stats$group <- c(stats$group, mat[[2]]$group)
    stats$lower_cutoff <- c(stats$lower_cutoff, mat[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, mat[[2]]$middle)
    stats$upper_cutoff <- c(stats$upper_cutoff, mat[[2]]$upper)
    
    # Rename the column "Expression" to "model"
    survival_data <- survival_data %>%
      dplyr::rename(model = Expression)
  }
  
  # STEP 2: Calculate survival stats
  # If each group has to be plotted in separate plots, subset the survival_data
  # and pass it to plot_survival(). Else, pass entire survival_data to 
  # plot_survival().
  if (combine_plot == "FALSE") {
    if (!is.na(split_by)){
      for (x in (expr_df %>% dplyr::distinct(get(split_by)))[[1]]) {
        s_data <- survival_data %>% dplyr::filter(get(split_by) == x)
        cox_stats <- plot_survival(s_data, x)
        
        # Save the data from output of plot_survival()
        stats$HR <- c(stats$HR, cox_stats$HR )
        stats$CI_lower <- c(stats$CI_lower, cox_stats$CI_lower)
        stats$CI_upper <- c(stats$CI_upper, cox_stats$CI_upper)
        stats$pvalue <- c(stats$pvalue, cox_stats$pvalue)
      }
    } else {
      s_data <- survival_data
      x <- "NA"
      cox_stats <- plot_survival(s_data, x)
      
      # Save the data from output of plot_survival()
      stats$HR <- c(stats$HR, cox_stats$HR )
      stats$CI_lower <- c(stats$CI_lower, cox_stats$CI_lower)
      stats$CI_upper <- c(stats$CI_upper, cox_stats$CI_upper)
      stats$pvalue <- c(stats$pvalue, cox_stats$pvalue)
    }
  }  else{
    s_data <- survival_data
    x <- "NA"
    cox_stats <- plot_survival(s_data, x)
    
    # Save the data from output of plot_survival()
    stats$HR <- c(stats$HR, cox_stats$HR )
    stats$CI_lower <- c(stats$CI_lower, cox_stats$CI_lower)
    stats$CI_upper <- c(stats$CI_upper, cox_stats$CI_upper)
    stats$pvalue <- c(stats$pvalue, cox_stats$pvalue)
  }
  
  return(list(survival_data, stats))
}

calculate_cutoffs <- function(df, group){
  
  # Identify upper & lower cutoffs based on stratify_criteria
  #*************************Split samples by median**************************#
  if(stratify_criteria == "m"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[3]]
    cutoff_upper_end <- quartiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #****************Split samples into top and bottom tertiles****************#
  if(stratify_criteria == "t"){
    tertiles <- stats::quantile(x = df[[gene]],
                                probs = c(0, 0.33, 0.66, 1),
                                na.rm=TRUE)
    
    cutoff_lower_end <- tertiles[[2]]
    cutoff_upper_end <- tertiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #***************Split samples into top and bottom quartiles****************#
  if(stratify_criteria == "q"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[2]]
    cutoff_upper_end <- quartiles[[4]]
    cutoff_middle <- quartiles[[3]]
  }
  
  #*********************Split expression range by thirds*********************#
  if(stratify_criteria == "th"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    iqr <- stats::IQR(x = df[[gene]],
                      na.rm=TRUE)
    
    # Normal range of expression values lie between cutoff_lower & cutoff_upper
    cutoff_upper <- quartiles[[4]]+1.5*iqr
    cutoff_lower <- dplyr::if_else(quartiles[[1]]-1.5*iqr > 0, quartiles[[1]]-1.5*iqr, 0)
    
    # Based on normal range of expression, identify onethird & twothird cutoff
    cutoff_lower_end <- cutoff_lower + (cutoff_upper-cutoff_lower)/3
    cutoff_upper_end <- cutoff_lower + (cutoff_upper-cutoff_lower)*2/3
    cutoff_middle <- "NA"
  }
  
  #***************Split expression range using optimal cutoff****************#
  if(stratify_criteria == "o"){
    
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    # Sometimes quartiles will look like: 
    # 0%       25%      50%      75%     100% 
    # 0.000000 0.000000 0.000000 0.000000 3.495493 
    # In such cases, surv_cutpoint() will fail. So, we add extra if() here.
    if (quartiles[[4]] > quartiles[[2]]){
      res.cut <- survminer::surv_cutpoint(data = df,
                                          time = "Time",
                                          event = "Status",
                                          variables = gene)
      
      cutoff_lower_end <- res.cut$cutpoint$cutpoint
      cutoff_upper_end <- res.cut$cutpoint$cutpoint
      cutoff_middle <- "NA"
    } else{
      #cat("Surv cutpoint unable to detect optimum cutoff")
      cutoff_lower_end <- "NA"
      cutoff_upper_end <- "NA"
      cutoff_middle <- "NA"
    }
  }
  
  # Categorize the data based on above cutoffs
  if (all_quartiles == "TRUE"){
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::if_else(get(gene) > cutoff_upper_end, "HIGH", 
                                                dplyr::if_else(get(gene) < cutoff_lower_end, "LOW",
                                                               dplyr::if_else(get(gene) < cutoff_middle, "MED_LOW",
                                                                              "MED_HIGH"))))
  } else{
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::if_else(get(gene) > cutoff_upper_end, "HIGH", 
                                                dplyr::if_else(get(gene) < cutoff_lower_end, "LOW", "MID"))) %>%
      dplyr::filter(Expression != "MID")
  }
  
  # # Print the cutoffs
  # cat("\nGene         :", gene)
  # cat("\nGroup        :", group)
  # cat("\nLower cutoff :", cutoff_lower_end)
  # cat("\nUpper cutoff :", cutoff_upper_end)
  # cat("\nMiddle cutoff:", cutoff_middle)
  
  # Create a list to store cutoff values
  ls <- list("group" = c(), 
             "gene" = c(), 
             "lower" = c(), 
             "upper" = c(), 
             "middle" = c())
  
  ls$group <- c(group)
  ls$gene <- c(gene)
  ls$lower <- c(cutoff_lower_end)
  ls$upper <- c(cutoff_upper_end)
  ls$middle <- c(cutoff_middle)
  
  # Return the df
  return(list(df, ls))
  
}

plot_survival <- function(survival_data, group){
  
  # If all samples belong to one group (like LOW or HIGH or males or female),
  # then quit the function as comparison cannot be done
  if (nrow(survival_data %>% dplyr::count(model)) > 1){
    
    # # Remove rows with NA or negative values in Time
    # survival_data <- survival_data %>% dplyr::filter(Time > 0)
    
    # Create a survival object where Alive = 0, Dead = 1
    survival_object <- survival::Surv(time = survival_data$Time,
                                      event = survival_data$Status,
                                      type = "right",
                                      origin = 0)
    
    # Create a formula for plotting survival curve
    survival_formula <- survival_object ~ model
    
    # Create a fit for survival curve. survfit() gives error in ggsurvplot(). Use surv_fit()
    survival_curve <- survminer::surv_fit(formula = survival_formula,
                                          data = survival_data,
                                          type = "kaplan-meier",
                                          group.by = NULL,
                                          match.fd = FALSE)
    
    # Check summary of the survival curve with time duration of our interest
    #cat("\nRange of survival (days):", range(survival_data$Time, na.rm=TRUE), "\n")
    base::summary(survival_curve, times = base::seq(from = floor(range(survival_data$Time, na.rm=TRUE)[[1]]/500),
                                                    to = ceiling(range(survival_data$Time, na.rm=TRUE)[[2]]/500)*500,
                                                    by = 500))
    
    # Create a Cox model for the survival curve and calculate stats
    cox_model <- survival::coxph(formula = survival_formula,
                                 data = survival_data)
    #print(summary(cox_model))
    cat("\n")
    
    # Calculate HR, 95% CI for HR, p-val
    # NOTE: Variable mentioned here is the numerator in h1(t)/h0(t).
    # The reference variable h0(t) will not be mentioned in co-efficients.
    # Make sure this is not the reference level i.e. low expression. If this is
    # the reference, then we need to reverse the HR ratio, legend labels
    #print(names(cox_model$coefficients))  
    
    # If all samples belong to more than 2 groups (like LOW, MID, HIGH), then
    # we cannot have survival stats. SO, we set them to 0.
    if (nrow(survival_data %>% dplyr::count(model)) == 2){
      if (stringr::str_detect(names(cox_model$coefficients), reference)){
        HR <- round(exp(-cox_model$coefficients[[1]]), 2)
        CI <- round(exp(-confint(cox_model)), 2)
        CI_1 <- CI[1]
        CI[1] <- CI[2]
        CI[2] <- CI_1
        p_val <- survminer::surv_pvalue(fit = survival_curve,
                                        stratify_criteria = "survdiff",
                                        test.for.trend = FALSE,
                                        combine = FALSE)
      } else {
        HR <- round(exp(cox_model$coefficients[[1]]),2)
        CI <- round(exp(confint(cox_model)),2)
        p_val <- survminer::surv_pvalue(fit = survival_curve,
                                        stratify_criteria = "survdiff",
                                        test.for.trend = FALSE,
                                        combine = FALSE)
      }
    } else {
      HR <- 0
      CI <- c(0, 0)
      p_val <- c(0, 0)
    }
    
    # Plot the survival curve using survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: a survival curve and a risk table
    # Saving it using cowplot() first and then using ggsave() works nicely as
    # compared to saving directly using ggsave()
    if(plot_curve =="TRUE"){
      # Plot the p-value
      grob1 <- grobTree(textGrob(label = paste0("p = ", round(p_val[[2]],2)),
                                 x = 0.50,
                                 y = 0.95,
                                 hjust = 0,
                                 gp = grid::gpar(fontfamily="Times", fontface="bold", col="black", fontsize=20)))
      
      # Plot the HR-value
      grob2 <- grobTree(textGrob(label = paste0("HR = ", round(HR,1), " [", round(CI[1],1), ", ", round(CI[2],1), "]"),
                                 x = 0.50,
                                 y = 0.87,
                                 hjust = 0,
                                 gp = grid::gpar(fontfamily="Times", fontface="bold", col="black", fontsize=20)))
      
      # Plot the survival curve
      survival_plot <- survminer::ggsurvplot(fit = survival_curve,
                                             conf.int = confidence_interval,
                                             legend.title = legend_title,
                                             legend.labs=legend_label,
                                             palette=color_palette,
                                             size = 2,
                                             censor.size = 9) %++%
        ggplot2::annotation_custom(grob1) %++%
        ggplot2::annotation_custom(grob2) %++%
        ggplot2::labs(title = toupper(dplyr::if_else(gene=="combined.exp", "", gene)),
                      #fill = legend_title,
                      x = "Time (Days)",
                      y = "Survival Probability")  %++%
        ggplot2::theme(plot.title =   element_text(family="sans", face="bold", colour="black", size=25, hjust=0.5),
                       plot.caption = element_text(family="sans", face="bold", colour="black", size=25, hjust=0),
                       axis.title.x = element_text(family="sans", face="bold", colour="black", size=25, hjust=0.5),
                       axis.title.y = element_text(family="sans", face="bold", colour="black", size=25, hjust=0.5),
                       axis.text.x =  element_text(family="sans", face="bold", colour="black", size=20, hjust=0.5),
                       axis.text.y =  element_text(family="sans", face="bold", colour="black", size=20, hjust=1),
                       legend.title = element_text(family="sans", face="bold", colour="black", size=25, hjust=0.5),
                       legend.text =  element_text(family="sans", face="bold", colour="black", size=25, hjust=0.5),
                       #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                       legend.position = "bottom",   #c(0.78, 0.9)
                       legend.justification = "center",
                       legend.direction = "horizontal")
      
      # Save the plot
      ggplot2::ggsave(filename = paste0(gse, "_", group, "_", stratify_criteria, "_", gene, ".pdf"),
                      plot = last_plot(),
                      device = "pdf",
                      path = paste0(parent_path, "Surv_plots/"),
                      scale = 1,
                      width = 6,
                      height = 7,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = NULL)
    }
    
    # Plot the risk table  using survminer::ggrisktable()
    if (plot_risk_table == "TRUE"){
      risk_table <- survminer::ggrisktable(fit = survival_curve,
                                           risk.table.type = "absolute",
                                           legend.labs=legend_label,
                                           fontsize = 7) %++%
        ggplot2::labs(title = paste0("Number at Risk"),
                      fill = legend_title,
                      x = "Time (Days)",
                      y = legend_title)  %++%
        ggplot2::theme(plot.title =   element_text(family="sans", face="bold", colour="black", size=25, hjust = 0.5),
                       plot.caption = element_text(family="sans", face="bold", colour="black", size=25, hjust = 0),
                       axis.title.x = element_text(family="sans", face="bold", colour="black", size=25, hjust = 0.5),
                       axis.title.y = element_text(family="sans", face="bold", colour="black", size=25, hjust = 0.5),
                       axis.text.x =  element_text(family="sans", face="plain", colour="black", size=20, hjust = 0.5),
                       #axis.text.y =  element_text(family="sans", face="bold", colour="black", size=20, hjust = 0.5),
                       #legend.title = element_text(family="sans", face="bold", colour="black", size=20, hjust = 0.5),
                       #legend.text =  element_text(family="sans", face="bold", colour="black", size=20, hjust = 0.5),
                       #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                       #legend.position = c(0.78, 0.9),
                       #legend.justification = "center",
                       legend.direction = "vertical")
      
      # Save the plot
      ggplot2::ggsave(filename = paste0("Risk Table ", gse, "_", group, "_", stratify_criteria, "_", gene, ".pdf"),
                      plot = last_plot(),
                      device = "pdf",
                      path = paste0(parent_path, "Surv_plots/"),
                      scale = 1,
                      width = 6,
                      height = 2.2,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = NULL)
    }
  } else {
    HR <- 0
    CI <- c(0, 0)
    p_val <- c(0, 0)
  }
  
  # Create a list to store survival stats
  ls <- list("group" = c(), "HR" = c(), "CI_lower" = c(), "CI_upper" = c(), "pvalue" =c())
  
  ls$group <- c(group)
  ls$HR <- c(HR)
  ls$CI_lower <- c(CI[1])
  ls$CI_upper <- c(CI[2])
  ls$pvalue <- c(p_val[[2]])
  
  return(ls)
}

# We need to calculate a combined expression value of all genes in the gene
# signature for each sample. Next, we use surv_cutpoint() to classify samples 
# into high and low groups and plot survival curves.

# There are several approaches to calculate the combined expression value. While
# normal z-scaling is logical, https://doi.org/10.1186/gb-2006-7-10-r93 
# recommends a more advanced z-scaling which is implemented below. Refer
# "Measuring gene set expression" in "Materials & stratify_criteria" section.

advanced_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  ix <- is.element(toupper(rownames(eset)), toupper(gset))
  cat(sum(ix))
  
  if (sum(ix)>0){
    avg_gset <- base::apply(X=eset[ix,], MARGIN=2, FUN=mean, na.rm=TRUE)
    avg_all <- base::apply(X=eset, MARGIN=2, FUN=mean, na.rm=TRUE)
    sd_all <- base::apply(X=eset, MARGIN=2, FUN=sd, na.rm=TRUE)
    z <- (avg_gset - avg_all)*sqrt(sum(ix))/sd_all
  } else{
    z <- NA
  }
  
  return(z)
}

normal_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  eset <- eset[gset,]
  a <- t(scale(t(eset)))
  
  z <- colSums(a, na.rm=TRUE) 
  
  return(z)
}

#******************************************************************************#
#                     STEP 2: PLOT SURVIVAL CURVES BY SEX                      #
#******************************************************************************#

# Read survival data
expr_df <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "BBN original figs/1A_data.xlsx"))

gse <- "Mice Expt"
sex <- "Male vs Female"
i <-  ""

# Plot the gene signature
summary <- plot_survival(stratify_criteria)     

#******************************************************************************#
#                    STEP 3: PLOT SURVIVAL CURVES FOR GENES                    #
#******************************************************************************#

# Define the gse of the project
gse <- "TCGA_BLCA"         # DESeq2 normalized counts
gse <- "Blaveri"           # already log2 transformed & median centered
gse <- "mskcc"             # already log2 transformed
gse <- "gse13507"          # already log2 transformed
gse <- "gse31684"          # already log2 transformed
gse <- "gse32894"          # already log2 transformed
#gse <- "Imvigor210"

wb <- openxlsx::createWorkbook()
for (gse in c("TCGA_BLCA", "Blaveri", "gse13507", "mskcc", "gse32894", "gse31684")){
  
  # Import read_data
  read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_Normalized.xlsx"))
  colnames(read_data)[1] <- "SYMBOL"
  
  # Import meta_data and subset meta_data if needed
  meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_Metadata.xlsx"))  
  #%>% dplyr::filter(grepl(pattern = "T2|T3|T4", x = Stage))
  
  # Reformat metadata 
  meta_data <- meta_data %>% 
    dplyr::mutate(GEO_ID = make.names(names = GEO_ID)) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("GEO_ID", .keep_all = TRUE)
  
  if (!is.na(subset_group)){
    meta_data <- meta_data %>% 
      dplyr::filter(get(subset_group) %in% subset_value)
  }
  
  # Reformat read data
  normalized_counts <- read_data %>%
    dplyr::mutate(SYMBOL = make.names(names = SYMBOL)) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
  normalized_counts <- normalized_counts[,intersect(make.names(meta_data$GEO_ID), colnames(normalized_counts))]
  normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]
  
  # Perform log1p transformation on the subsetted data if not done before
  # Perform median centering for each gene across samples as median expression is 
  # more robust to outliers as compared to mean expression.
  if (gse == "TCGA_BLCA"){
    normalized_counts <- log(1+normalized_counts, base=2)
  }
  if (gse != "Blaveri"){
    t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
    normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
  }
  
  #*****************USE THIS SECTION FOR PLOTTING INDIVIDUAL GENES***************#
  
  #wb <- openxlsx::createWorkbook()
  if (gene_signature == FALSE){
    
    # Create a list of genes for which survival curves need to be plotted
    plot_genes <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "j.xlsx"), sheet = "Non-epithelial DEGs") %>%
      dplyr::filter(padj < 0.05 & log2FoldChange <= -0.58) %>%
      dplyr::distinct_at("Human.homolog") %>%
      unlist(use.names = FALSE)
    plot_genes <- setdiff(plot_genes, "None")
    
    # Merge expression data with survival data
    expr_df <- normalized_counts %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("GEO_ID") %>%
      dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID"))
    
    # Create a list to store survminer cutoffs, coxph stats, etc..
    stats <- list("gene" = c(),
                  "group" = c(),
                  "lower_cutoff" = c(),
                  "middle_cutoff" = c(),
                  "upper_cutoff" = c(),
                  "HR" = c(),
                  "CI_lower" = c(),
                  "CI_upper" = c(),
                  "pvalue" = c())
    
    # Plot survival curves
    for (gene in intersect(unique(plot_genes), rownames(normalized_counts))) {
      plot_curve <- FALSE
      summary <- wrangle_data(stratify_criteria)
      stats$gene          <- c(stats$gene, summary[[2]]$gene)
      stats$group         <- c(stats$group, summary[[2]]$group)
      stats$lower_cutoff  <- c(stats$lower_cutoff, summary[[2]]$lower)
      stats$middle_cutoff <- c(stats$middle_cutoff, summary[[2]]$middle)
      stats$upper_cutoff  <- c(stats$upper_cutoff, summary[[2]]$upper)
      stats$HR            <- c(stats$HR, summary[[2]]$HR )
      stats$CI_lower      <- c(stats$CI_lower, summary[[2]]$CI_lower)
      stats$CI_upper      <- c(stats$CI_upper, summary[[2]]$CI_upper)
      stats$pvalue        <- c(stats$pvalue, summary[[2]]$pvalue)
    }
    
    stats_df <- data.frame(stats)
    
    # Calculate padj for male and female patients separately
    stats_male <- stats_df %>% dplyr::filter(group == "Male")
    stats_male$padj <- stats::p.adjust(p = stats_male$pvalue, method = "fdr", n = length(stats_male$pvalue))
    stats_female <- stats_df %>% dplyr::filter(group == "Female")
    stats_female$padj <- stats::p.adjust(p = stats_female$pvalue, method = "fdr", n = length(stats_female$pvalue))
    stats_df <- dplyr::bind_rows(stats_male, stats_female)
    
    # Save the results
    openxlsx::addWorksheet(wb, sheetName = "7B")
    openxlsx::writeData(wb, sheet = "7B", x = stats_df, rowNames = FALSE)
  }
  #openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Individual_gene_stats.xlsx"), overwrite = TRUE)
  
  #*****************USE THIS SECTION FOR PLOTTING GENE SIGNATURES**************#
  
  if (gene_signature == TRUE){
    plot_curve <- TRUE

    m <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Individual_gene_stats_final.xlsx"), sheet = "7B") %>%
    #m <- stats_df %>%
      dplyr::filter(group == "Male" & padj < 0.05 & HR >= 1) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE)

    f <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Individual_gene_stats_final.xlsx"), sheet = "7B") %>%
    #f <- stats_df %>%
      dplyr::filter(group == "Female" & padj < 0.05 & HR >= 1) %>%
      dplyr::select(gene) %>%
      unlist(use.names = FALSE)

    common <- base::intersect(m,f)

    plot_genes <- base::setdiff(f, common)
    # plot_genes <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "j.xlsx"), sheet = "Signature") %>%
    #   dplyr::select("Male.Epithelial") %>%
    #   unlist(use.names = FALSE)
    
    plot_genes <- c("DDR2", "MDH2", "PGK1", "PGAM1", "ALDOA", "GPI")

    # Calculate z-score using the function described in above paper
    expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
    
    # Merge expression data with survival data
    expr_df <- expr_df %>%
      data.frame() %>%
      dplyr::rename(combined.exp = identity(1)) %>%
      tibble::rownames_to_column("GEO_ID") %>%
      dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID"))
    
    gene <- "combined.exp"
    if (nrow(expr_df > 0)){
      summary <- wrangle_data(stratify_criteria)
    }
    
    surv_df <- summary[[1]]
    stats_df <- as.data.frame(summary[[2]])
    
    # Save the results
    openxlsx::addWorksheet(wb, sheetName = paste0("Summary_",gse))
    openxlsx::writeData(wb, sheet = paste0("Summary_",gse), x = stats_df, rowNames = FALSE)
    openxlsx::addWorksheet(wb, sheetName = gse)
    openxlsx::writeData(wb, sheet = gse, x = surv_df, rowNames = FALSE)
  }
}
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "7B_summary.xlsx"), overwrite = TRUE)

#******************************************************************************#
#                      STEP 4: CHECK IF SIGNATURE IS GOOD                      #
#******************************************************************************#

# Create an AnnotatedDataFrame object for meta data
p_data <- meta_data
rownames(p_data) <- make.names(p_data$Sample_id)
p_data <- Biobase::AnnotatedDataFrame(data = p_data,
                                      varMetadata = data.frame("labelDescription" = colnames(p_data))) 

# Create an AnnotatedDataFrame object for features
f_data <- data.frame("SYMBOL" = rownames(normalized_counts))
rownames(f_data) <- make.names(f_data$SYMBOL)
f_data <- Biobase::AnnotatedDataFrame(data = f_data,
                                      varMetadata = data.frame("labelDescription" = colnames(f_data))) 

# Create an ExpressionSet object for read data
e_data <- as.matrix(normalized_counts)
eset <- Biobase::ExpressionSet(assayData = e_data,
                               phenoData = p_data,
                               featureData = f_data,
                               annotation = "custom")

# NOTE: classes parameter in sigCheck() = Status column of varLabels(eset)
# NOTE: survival parameter in sigCheck() = Time column of varLabels(eset)
varLabels(eset)
eset$Status
eset$Time
# NOTE: annotation parameter in sigCheck() = SYMBOL column of fvarLabels(eset)
# NOTE: plot_genes and genes in SYMBOL column must have some overlap
fvarLabels(eset)   

# Perform anlysis for male and female
# sigCheck is not good at classifying samples into optimal groups. So, we 
# manually classify the samples using survminer and import the classification
# into sigCheck() using scoreMethod and threshold parameters


#*****************USE THIS SECTION FOR PLOTTING GENE SIGNATURES****************#
p <- c()

for (i in 1:100){
  
  plot_genes <- base::sample(x=rownames(normalized_counts), size=64, replace = FALSE)
  
  # Calculate z-score using the function described in above paper
  expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
  
  # Merge expression data with survival data
  expr_df <- expr_df %>%
    data.frame() %>%
    dplyr::rename(combined.exp = identity(1)) %>%
    tibble::rownames_to_column("GEO_ID") %>%
    dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID"))
  
  gene <- "combined.exp"
  if (nrow(expr_df > 0)){
    summary <- wrangle_data(stratify_criteria)
  }
  
  surv_df <- summary[[1]] %>%
    dplyr::filter(Status == 0 | Status == 1) %>%
    dplyr::mutate(Expression = stringi::stri_replace_all_regex(str = Expression,
                                                               pattern = c("HIGH", "LOW"),
                                                               replacement = c("1", "0"),
                                                               vectorise_all = FALSE)) %>%
    dplyr::filter(Sex == "Male")
  
  sigCheck_score <- function(eset){
    e <- surv_df  %>% 
      dplyr::select(Expression) %>%
      unlist(., use.names = FALSE) %>%
      as.numeric()
    
    return(e)
  }
  
  # Format the object to remove sample not present in surv_df
  # Also, remove samples which have status other than 0 or 1.
  eset_subset <- eset[, eset$GEO_ID %in% surv_df$GEO_ID]
  
  # Create a SigCheck object
  check <- sigCheck(expressionSet = eset_subset, 
                    classes = "Status", 
                    survival = "Time",
                    signature = plot_genes,
                    annotation = "SYMBOL",
                    scoreMethod = sigCheck_score(eset),
                    threshold = sum(sigCheck_score(eset))/length(sigCheck_score(eset)))
  
  p <- c(p, check@survivalPval)
  
  
  # sigCheckRandom(check = check,
  #                iterations=100)
}

#******************************************************************************#