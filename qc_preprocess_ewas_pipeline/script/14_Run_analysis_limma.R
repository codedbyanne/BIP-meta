
# LOAD LIBRARY
suppressMessages({
  library(readxl)
})


### DECIDE WHICH SCRIPT TO RUN DEPENDING ON THE NR OF MALES AND FEMALES AFTER QC
samplesheet <- as.data.frame(read_xlsx("output/Tables/samplesheet_after_QC.xlsx"))

if (multiple_samples == "yes"){
  samplesheet <- samplesheet[which(samplesheet$Baseline_sample == "Y"),]
}

nr_males <- nrow(samplesheet[which(samplesheet$Sex == "M"), ])
nr_females <- nrow(samplesheet[which(samplesheet$Sex == "F"), ])

if (nr_males >= 35 & nr_females >= 35){
  
  
  if (tissue_type == "blood"){
    
    
    print("Run sex-specific case-control analyses.")
    print("1) Run females")
    source("script/14_CaseControl_Analysis_Females_limma_blood.R")
    print("2) Run males")
    source("script/14_CaseControl_Analysis_Males_limma_blood.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4_CaseControl_sex-specific.Rmd", 
                      output_file = paste0(dir_gen,"Report_Part5_CaseControl_Analysis_sex-specific.html"),
                      params = list(dir_gen = dir_gen))
    
  }
    
  
  if (tissue_type == "saliva"){
    
    print("Run sex-specific case-control analyses.")
    print("1) Run females")
    source("script/14_CaseControl_Analysis_Females_limma_saliva.R")
    print("2) Run males")
    source("script/14_CaseControl_Analysis_Males_limma_saliva.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4_CaseControl_sex-specific.Rmd", 
                      output_file = paste0(dir_gen,"Report_Part5_CaseControl_Analysis_sex-specific.html"),
                      params = list(dir_gen = dir_gen))
  }
  
} else {
  
  if (tissue_type == "blood"){
    
    print("Running sex-adjusted case-control analysis.")
    
    source("script/14_CaseControl_Analysis_Females_Males_limma_blood.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4_CaseControl_sex-adjusted.Rmd", 
                      output_file = paste0(dir_gen,"Report_Part5_CaseControl_Analysis_sex_adjusted.html"),
                      params = list(dir_gen = dir_gen)
    )
  }
  
  if (tissue_type == "saliva"){
    
    print("Running sex-adjusted case-control analysis.")
    
    source("script/14_CaseControl_Analysis_Females_Males_limma_saliva.R")
    
    # CREATE A REPORT
    print("Create report")
    source("2_Paths.R")
    rmarkdown::render("script/Report_Part4_CaseControl_sex-adjusted.Rmd", 
                      output_file = paste0(dir_gen,"Report_Part5_CaseControl_Analysis_sex_adjusted.html"),
                      params = list(dir_gen = dir_gen)
    )
  }
}