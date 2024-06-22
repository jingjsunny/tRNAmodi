library(readxl)
setwd("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/Input MS data/")
file_list <- list.files()
file_length <- length (file_list)
#Load modification list
fc_mod_list <- read.csv("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/mod_list_fc.csv")
fc_mod_list <- fc_mod_list [,2]
Plate_condition <- read.csv("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/Plate condition.csv")
PA14_library <- read_excel("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/PA14 library plate designations.xlsx")
w <- 1
repeat {
  #Read batch table csv file
  MS <- read.csv(file_list[w], stringsAsFactors = FALSE)
  basename <- gsub(" MS_all.csv","",file_list[w])
  #Get row names
  rownames(MS) <- gsub ("\\d+\\_", "", MS[,1])
  #Cleanup data
  MS <- MS[-c(1), -c(1:4)]
  #Prepare column name
  n <- as.integer(seq (2, ncol(MS), by = 5))
  mod_name <- colnames(MS)
  mod_name <- gsub(".Results","",mod_name[n-1])
  #Extract MS and UV signal
  MS_UV_signal <- sapply(MS[,n],as.numeric)
  MS_UV_signal[is.na(MS_UV_signal)]<-0
  colnames(MS_UV_signal) <- mod_name
  rownames(MS_UV_signal) <- rownames(MS)
  
  ###Get samples with abnormal input amount, UV mean <0.2 of the mean of the whole plate
  m <- ncol(MS_UV_signal)
  rA_mean <- mean(MS_UV_signal[,(m-3)],na.rm=TRUE)
  rG_mean <- mean(MS_UV_signal[,(m-1)],na.rm=TRUE)
  rU_mean <- mean(MS_UV_signal[,m],na.rm=TRUE)
  UV_ab_rows <- (MS_UV_signal[, (m-3)]< rA_mean*0.2) | (MS_UV_signal[, (m-1)]< rG_mean*0.2) |(MS_UV_signal[, m]< rU_mean*0.2)
  UV_ab <- MS_UV_signal[UV_ab_rows, ,drop = FALSE]
  #Get a full list of samples with abnormal inputs
  temp1 <- data.frame(Plate = rep (basename,nrow(UV_ab)), Well = rownames(UV_ab), UV_ab)
  if (w == 1) {
    UV_ab_all <- temp1
  }
  if (w > 1) {
    UV_ab_all <- rbind(UV_ab_all, temp1)
  }
  #remove abnormal data
  MS_UV_signal <- MS_UV_signal[!UV_ab_rows,]

  ###Calculate normalized peak area
  N_peakarea_all <- MS_UV_signal[,1:(m-4)]/rowSums(MS_UV_signal[,c((m-3),(m-1),m)])
  #Get a full table of normalized peak areas
  temp2 <- data.frame(Plate = rep (basename,nrow(N_peakarea_all)), Well = rownames(N_peakarea_all), N_peakarea_all)
  if (w == 1) {
    Peakarea_all <- temp2
  }
  if (w > 1) {
    Peakarea_all <- rbind(Peakarea_all, temp2)
  }
  
  ###Calculate fold change and combine into a full table
  N_peakarea <- N_peakarea_all[,fc_mod_list]
  FC <- N_peakarea
  j <- 0
  repeat {
    # Define the row indices for the current plate
    start_row <- 12 * j + 1
    if (start_row > nrow(N_peakarea)) { 
      break
    }
    end_row <- min(12 * j + 12, nrow(N_peakarea))
    # Extract the data for the current plate
    filtered_data <- N_peakarea[start_row:end_row,]
    # Calculate initial averages, excluding NAs
    initial_averages <- colMeans(filtered_data, na.rm = TRUE)
    # Filter out values more than twice and less than half the initial average
    for (t in 1:ncol(filtered_data)) {
      large_values_indices <- which(filtered_data[, t] > 2 * initial_averages[t])
      filtered_data[large_values_indices, t] <- NA
    }
    processed_averages <- colMeans(filtered_data, na.rm = TRUE)
    for (t in 1:ncol(filtered_data)) {
      small_values_indices <- which(filtered_data[, t] < 0.5 * processed_averages[t])
      filtered_data[small_values_indices, t] <- NA
    }
    # Recompute the averages, excluding the large values
    final_averages <- colMeans(filtered_data, na.rm = TRUE)
    # Update the FC values
    for (i in 1:nrow(filtered_data)) {
      FC[start_row + i - 1, ] <- N_peakarea[start_row + i - 1, ] / final_averages
    }
    
    # Increment j and check if all plates have been processed
    j <- j + 1
    if (j > 7 ) {
      break
    }
  }
  #Get a full list of fold changes
  temp3 <- data.frame(Plate = rep (basename,nrow(FC)), Well = rownames(FC), FC)
  if (w == 1) {
    FC_all <- temp3
  }
  if (w > 1) {
    FC_all <- rbind(FC_all, temp3)
  }
  
  ### report fold change >2 or <0.5 data
  FC_F2_list <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(FC_F2_list) <- c('Plate','Well','Modification','FC')
  i <- 1
  n <- 1
  repeat {
    j <- 1
    repeat{
      if(!is.na(FC[i,j])){
        if (FC[i,j] < 0.5 | FC[i,j] > 2){
          FC_F2_list[n,1] <- basename
          FC_F2_list[n,2] <- rownames(FC)[i]
          FC_F2_list[n,3] <- colnames(FC)[j]
          FC_F2_list[n,4] <- FC[i,j]
          n <- n+1
        }
      }
      j <- j+1
      if (j > ncol(FC)) {break}
    }
    i <- i+1
    if (i > nrow(FC)){break}
  }
  #Get a full list of fold changes
  if (w == 1) {
    FC_F2_full_list <- FC_F2_list
  }
  if (w > 1) {
    FC_F2_full_list <- rbind(FC_F2_full_list, FC_F2_list)
  }
  #Get a full list of fold changes with all data
  if(length(unique(FC_F2_list[,2])) == 1) {
    temp4 <- data.frame(Plate = basename, Well = unique(FC_F2_list[,2]),t(FC[c(unique(FC_F2_list[,2])),]))
  }
  if(length(unique(FC_F2_list[,2])) > 1) {
  temp4 <- data.frame(Plate = rep (basename,length(unique(FC_F2_list[,2]))), Well = unique(FC_F2_list[,2]), FC[c(unique(FC_F2_list[,2])),])
  }
  if (w == 1) {
    FC_F2_data <- temp4
  }
  if (w > 1) {
    FC_F2_data <- rbind(FC_F2_data, temp4)
  }
  w <- w+1
  if (w > file_length) {break}
}

###remove empty well in the library that filled with WT strain
ID1 <- paste(Plate_condition$Plate,Plate_condition$Note)
Peakarea_all <- Peakarea_all[!paste(Peakarea_all$Plate,Peakarea_all$Well)%in%ID1,]
FC_all <- FC_all[!paste(FC_all$Plate,FC_all$Well)%in%ID1,]
FC_F2_data <- FC_F2_data[!paste(FC_F2_data$Plate,FC_F2_data$Well)%in%ID1,]
FC_F2_full_list <- FC_F2_full_list[!paste(FC_F2_full_list$Plate,FC_F2_full_list$Well)%in%ID1,]
UV_ab_all <- UV_ab_all[!paste(UV_ab_all$Plate,UV_ab_all$Well)%in%ID1,]

###Annotate each mutant with gene information
PA14_library$new_label <- gsub ("PAMr_nr_mas_0", "Plate ", PA14_library$`PA14 NR Set Plate (PO_PlateLabel96)`)
PA14_library$new_label <- gsub ("PAMr_nr_mas_", "Plate ", PA14_library$new_label)
PA14_library$new_label <- gsub ("ExMr_nr_mas_01_", "Plate ExMr_", PA14_library$new_label)
PA14_library$new_label <- gsub ("PAPho_nr_mas_01_1", "Plate PAPho", PA14_library$new_label)
#Annotate Peakarea_all
Peakarea_all <- cbind(Peakarea_all, 
                      Gene_ID = rep(NA, nrow(Peakarea_all)), 
                      Gene_name = rep(NA, nrow(Peakarea_all)), 
                      PAO1_ID = rep(NA, nrow(Peakarea_all)),
                      Blast = rep(NA, nrow(Peakarea_all)))
i <- 1
repeat {      
  location <- which(Peakarea_all$Plate == c(PA14_library[i,25]) & Peakarea_all$Well == c(PA14_library[i,2]))
  Peakarea_all[location, 45:48] <- PA14_library[i, c(9, 14, 13,23)]
  i <- i+1
  if (i > nrow(PA14_library)){break}
}
#Annotate FC_all
FC_all <- cbind(FC_all, 
                Gene_ID = rep(NA, nrow(FC_all)), 
                Gene_name = rep(NA, nrow(FC_all)), 
                PAO1_ID = rep(NA, nrow(FC_all)),
                Blast = rep(NA, nrow(FC_all)))
i <- 1
repeat {      
  location <- which(FC_all$Plate == c(PA14_library[i,25]) & FC_all$Well == c(PA14_library[i,2]))
  FC_all[location, 33:36] <- PA14_library[i, c(9, 14, 13,23)]
  i <- i+1
  if (i > nrow(PA14_library)){break}
}
#Annotate FC_F2_data
FC_F2_data <- cbind(FC_F2_data, 
                    Gene_ID = rep(NA, nrow(FC_F2_data)), 
                    Gene_name = rep(NA, nrow(FC_F2_data)), 
                    PAO1_ID = rep(NA, nrow(FC_F2_data)),
                    Blast = rep(NA, nrow(FC_F2_data)))
i <- 1
repeat {      
  location <- which(FC_F2_data$Plate == c(PA14_library[i,25]) & FC_F2_data$Well == c(PA14_library[i,2]))
  FC_F2_data[location, 33:36] <- PA14_library[i, c(9, 14, 13,23)]
  i <- i+1
  if (i > nrow(PA14_library)){break}
}

#Annotate UV_ab_all
UV_ab_all <- cbind(UV_ab_all, 
                    Gene_ID = rep(NA, nrow(UV_ab_all)), 
                    Gene_name = rep(NA, nrow(UV_ab_all)), 
                    PAO1_ID = rep(NA, nrow(UV_ab_all)),
                   Blast = rep(NA, nrow(UV_ab_all)))
i <- 1
repeat {      
  location <- which(UV_ab_all$Plate == c(PA14_library[i,25]) & UV_ab_all$Well == c(PA14_library[i,2]))
  UV_ab_all[location, 49:52] <- PA14_library[i, c(9, 14, 13,25)]
  i <- i+1
  if (i > nrow(PA14_library)){break}
}

###report mutants that has more than 10 modifications with FC >1.5 or <0.7, those mutants need to be revisited manually.
filter_condition <- function(row) {
  sum(row > 1.5 | row < 0.7)
}
rows_to_keep <- apply(FC_all[, 3:32], 1, filter_condition) > 10
abnormal_FC <- FC_all[rows_to_keep, ]

setwd("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/Output")
write.csv(Peakarea_all,"Peakarea_all.csv",row.names = FALSE)
write.csv(FC_all,"FC_all.csv",row.names = FALSE)
write.csv(FC_F2_data,"FC_F2_data.csv",row.names = FALSE)
write.csv(FC_F2_full_list,"FC_F2_full_list.csv",row.names = FALSE)
write.csv(UV_ab_all,"UV_ab_all.csv",row.names = FALSE)
