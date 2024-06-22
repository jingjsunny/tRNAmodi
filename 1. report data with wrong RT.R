setwd("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/MS_all_0317/")
#Read batch table csv file
csv_file_name <- list.files(pattern = "\\.csv$", full.names = TRUE)
MS <- read.csv(csv_file_name, stringsAsFactors = FALSE)
#Get modification inclusion list
mod_list <- read.csv("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/mod_list_fc.csv")
mod_list <- mod_list [,2]
#Create a data to store abnormal RT list
RT_ab_list <- data.frame(matrix(ncol=4, nrow=0))
colnames(RT_ab_list) <- c('Plate','Well','Modification','RT')
#Get a table of retention time
Modi_columns <- grep(".Results", names(MS))
RT_list <- MS[,Modi_columns]
rownames(RT_list) <- MS[,1]
colnames(RT_list) <- gsub(".Results","",colnames(RT_list))
RT_list <- RT_list[-1,mod_list]
basename <- gsub(" MS_all.csv","",file_list[w])
#Calculate mean RT for each modification
  for (i in 1:ncol(RT_list))
  {
    RT_list[,i] <- as.numeric(RT_list[,i])
  }
  RT_mean <- colMeans(RT_list)
  #Get abnormal retention time
  delta_RT <- 0.2 #setup delta RT here
  for (i in 1:ncol(RT_list)){
    abnormal_indices <- which(RT_list[,i] > (RT_mean[i]+delta_RT) | RT_list[,i] < (RT_mean[i]-delta_RT))
    if (length(abnormal_indices) > 0) {
      n <- length(abnormal_indices)
      df <- data.frame(Plate = rep(basename,n),
                         Well = rownames(RT_list)[abnormal_indices],
                         Modification = rep (colnames(RT_list)[i],n), 
                         RT = RT_list[abnormal_indices,i])
      RT_ab_list <- rbind (RT_ab_list, df)
    }
  }

setwd("~/Zoho WorkDrive (AMR IRG)/My Folders/1. Running project/4. tRNA screening Pseudomonas/Data analysis for publication 20240303/")
write.csv(RT_ab_list,"RT_ab_list.csv")
