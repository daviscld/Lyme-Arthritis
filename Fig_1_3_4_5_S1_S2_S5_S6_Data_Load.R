#Data cleaning and processing
data <- as.data.frame(read_excel("../SerologyData_LymeArthritis.xlsx"),
                      header = TRUE)

data_corrected <- data[1:84,] #exclude irrelevant or pbs samples for now
jf_ind <- which(data_corrected$Fluid == "Joint Fluid")
s_ind <- which(data_corrected$Fluid == "Serum")
data_corrected$Sample[jf_ind] <- paste0(data_corrected$Sample[jf_ind], "jf")
data_corrected$Sample[s_ind] <- paste0(data_corrected$Sample[s_ind], "s")

data_new <- as.data.frame(read_excel("../Lyme_SteereBowman_Demographic.xlsx"),
                          header = FALSE)[,3:34]

#Create interpretable column names
tmp <- c("Serum_ID", "JF_ID", "SampleDate","Age","Sex","HLADRB1_type","HLADRB1_hirisk",
         "filler","Arth_onset_date","PO_abx_start_date","duration_abx_sample","mo_on_abx", "resolve_po_abx",
         "filler2","time_sample_preIV","IV_abx_start_date","sample_during_IV","time_sample_postIV",
         "resolve_IV_abx", "filler3", "DMARD_start_date","duration_DMARD_to_arth_resolve","total_arth_dur",
         "arth_resolution_date","synovec_date","IgG_Bb_ELISA","Bb_JF_PCR","JF_WBC","JF_perc_polys",
         "JF_perc_lymphs","ESR","CRP")

colnames(data_new) <- tmp
drop <- c("filler","filler2","filler3")
data_new <- data_new[3:nrow(data_new),-which(colnames(data_new) %in% drop)]

#Initial df of interest: sex, age, HLADR status (high-risk and also actual)
keep <- c("Serum_ID","JF_ID", "Sex","Age","HLADRB1_type","HLADRB1_hirisk","total_arth_dur")
dem_data <- data_new[,which(colnames(data_new) %in% keep)]
tmp <- c(dem_data$Serum_ID, dem_data$JF_ID)
dem_data <- rbind(dem_data[,3:7], dem_data[,3:7])
dem_data$ID <- tmp
dem_data <- dem_data[-which(is.na(dem_data$ID)),]

# there are 3 individuals who were apparently included in the refractory group but who never
#received IV abx and went straight from oral abx to prolonged NSAID course (notably NOT a DMARD).
#These were very early patients so the clinical algorithm / treatment pathway was not yet
#established and that's why they didn't get IV abx and DMARD. They PROBABLY should be in the
#refractory group but not totally clear they wouldnt have responded to IV abx.
#These patients were thus excluded.
drop_pt <- c("GANC01","GANC02","GANC03","GALA24s",
             "GANC04", "GANC05", "GANC06", "GANC07", "GANC08",
             "GALA28s","GALA29s","GALA24jf","GALA28jf","GALA29jf")
dem_data <- dem_data[-which(dem_data$ID %in% drop_pt),]
data_corrected <- data_corrected[-which(data_corrected$Sample %in% drop_pt),]
data_corrected <- merge(data_corrected, dem_data, by.x = "Sample", by.y = "ID",
                        all.x = TRUE)

#Create vector of fluid types
fluid <- data_corrected$Fluid
fluid <-factor(fluid)

#Create vector of clinical phenotype
phenotype <- as.character(data_corrected$Disease)
phenotype <- map(strsplit(phenotype, " "), 2)
phenotype <- as.factor(unlist(phenotype))
data_corrected$Disease <- phenotype

#Create vector of patient ID
patientID <- data_corrected$Sample

#Create vector of age
age <- as.numeric(data_corrected$Age)

#Create vector of sex
sex <- as.factor(data_corrected$Sex)

#Create vector of HLADRB1 type
HLADRB1_type <- as.factor(data_corrected$HLADRB1_type)

#Create vector of HLADRB1 risk
HLADRB1_hirisk <- as.factor(data_corrected$HLADRB1_hirisk)

#Create vector of total arthritis duration
total_arth_dur <- as.numeric(data_corrected$total_arth_dur)

#use 2x PBS control to correct
pbs <- 2*as.numeric(colMeans(data[which(data[,1] == "PBS"),4:288]))
pbs_adcp <- 2*as.numeric(data[which(data[,1] == "PBS")[1],289:296])
pbs <- append(pbs, pbs_adcp)
exclude <- c("Sample","Disease", "Fluid", "Sex","Age","HLADRB1_type",
             "HLADRB1_hirisk","total_arth_dur")
`%!in%` <- Negate(`%in%`)
num_ind <- which(colnames(data_corrected) %!in% exclude)
data_corrected <- sweep(data_corrected[,num_ind],2,pbs,"-")
data_corrected[data_corrected < 0] <- 0
data_corrected[is.na(data_corrected)] <- 0

#log10 transform
logdata_corrected <- data_corrected
logdata_corrected[logdata_corrected == 0] <- 1
logdata_corrected <- log10(logdata_corrected)
logdata_corrected[,which(grepl("ADCP",colnames(logdata_corrected)))] <- data_corrected[,which(grepl("ADCP",colnames(data_corrected)))]
data_corrected <- logdata_corrected

#Remove certain antigens
remove_ag <- c("EBOV", "HA1", "Tetanus", "MMP10")

for (ag in remove_ag){
  data_corrected <- data_corrected[,-which(grepl(ag, colnames(data_corrected)))]
}

#Make all apoB100 names the same
colnames(data_corrected) <- gsub("ApoB 100|apoB 100|ApoB100|apoB100", "apoB", colnames(data_corrected))

# Remove features with too much missing data (>50% are 0)
# Make all polar plots and antigen-specific activity plots BEFORE doing this
remove_ind <- list()
n = 0
for (column in colnames(data_corrected)){
  n = n+1
  values <- data_corrected[,n]
  if (length(values[which(values < 0.2)])/length(values) > 0.499)
    remove_ind <- append(remove_ind, n)
}
data_corrected <- data_corrected[,-unlist(remove_ind)]

#Pull out identifying number from Excel sheet subject name
identity <- rownames(data_corrected)

patientID_multilevel <- gsub("jf|s","", patientID)

#Create data frame with corresponding identity, strain, and resistor status
df_id <- data.frame(identity = identity, fluid = fluid, phenotype = phenotype,
                    age = age, sex = sex, HLADRB1_type = HLADRB1_type,
                    HLADRB1_hirisk = HLADRB1_hirisk, total_arth_dur = total_arth_dur,
                    sample = patientID)

df_features <- data.frame(name = colnames(data_corrected))
df_features$label <- factor(df_features$name)

my_colors <- list(fluid = c("Joint Fluid" = '#6666ae', "Serum" = '#cb9ac6'),
                  phenotype = c("Refractory" = '#363795', "Responsive" = '#41bb93'))


