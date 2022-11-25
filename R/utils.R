
extract_data <- function(folder_name){
  #
  #' @description This function aims at extracting all the datas from
  #' a specific slim data repository.
  #'
  #
  #' @param folder_name The name of the slim data folder
  #
  #' @return data object containing all the relevant datas
  
  # At the end we will return the object return data
  # Depending on the presence of specific file in the foldern, return_data
  # will be updated with the relevant object
  return_data = list()
  
  # In this part we will create all possible file name
  # Step 1 and step 2 corresponds to simulation on which we decided
  # to extract data at 2 different time point
  
  #genome, m2 and m3 has been created if we have already converted vcf file into genome
  #as this is a time consuming process, we prefer to obtain the matrix from the read.table function
  #rather than going through all the vcf conversion process
  filename_genome <- paste(folder_name, "genome", sep="")
  filename_m2 <- paste(folder_name, "m2", sep="")
  filename_m3 <- paste(folder_name, "m3", sep="")
  
  
  filename_genome_step1 <- paste(folder_name, "genome_step1.vcf", sep="")

  filename_m2_step1 <- paste(folder_name, "mutationm2_step1.txt", sep="")

  
  filename_m3_step1 <- paste(folder_name, "mutationm3_step1.txt", sep="")

  
  filename_position_step1 <- paste(folder_name, "position_ind_step1", sep="")

  
  filename_var1_step1 <- paste(folder_name, "var1_step1", sep="")
  filename_var1_pred <- paste(folder_name, "var1_pred", sep="")
  
  filename_var2_step1 <- paste(folder_name, "var2_step1", sep="")
  filename_var2_pred <- paste(folder_name, "var2_pred", sep="")
  
  filename_fitness_step1 <- paste(folder_name, "fitness_step1", sep="")
  # Fitness pred is the predicted fitness if there is no adaptation at all
  # i.e what would be the fitness of the individual of step 1 if it faced
  # the environment of step 2
  filename_fitness_pred <- paste(folder_name, "fitness_pred", sep="")
  
  # As all the files are not always present in all the simulation (for instance
  # obtaining data of 2 time point is not always relevant depending on what
  # you want to do), for each filename, we check if the file exist, if it exists
  # we add the content of the file to the return_data object
  
  if (file.exists(filename_genome)){
    print("Obtaining genome from previous conversion")
    genome_maf <- read.table(filename_genome)
    
    print("End of genome Extraction @genome_step1 to get it")
    return_data <- c(return_data, genome_step1 = list(genome_maf))
    
    nb_individuals <- dim(genome_maf)[1]
    return_data <- c(return_data, nb_individuals = nb_individuals)
  }
  else {
    print("Extracting genome step 1 from VCF file ...")
    # For the moment we directly convert VCF files into matrix of SNPs data
    # At some point, maybe we need to collect all the information of the VCF file
    # in order to perform alignment between step 1 and step 2
    genome_vcf_step1 <- read.table(file = filename_genome_step1, header=FALSE)
    genome_simu_step1 = genome_vcf_step1[, 10:dim(genome_vcf_step1)[2]]
    genome_simu_step1 = data.frame(t(genome_simu_step1), stringsAsFactors=FALSE)
    print("Converting genome to 0, 1 and 2 values ...")
    print("Step 1 ...")
    genome_simu_step1[genome_simu_step1=="0|0"] = 0
    print("Step 1 OK")
    print("Step 2 ...")
    genome_simu_step1[genome_simu_step1=="0|1"] = 1
    print("Step 2 OK")
    print("Step 3 ...")
    genome_simu_step1[genome_simu_step1=="1|0"] = 1
    print("Step 3 OK")
    print("Step 4 ...")
    genome_simu_step1[genome_simu_step1=="1|1"] = 2
    print("Step 4 OK")
    genome_simu_step1[] <- lapply(genome_simu_step1, as.numeric)
    
    maf <- maf_filtering(genome_simu_step1, 0.02)
    genome_maf <- maf$genome_maf
    above_maf <- maf$above_maf
    
    write.table(genome_maf, filename_genome)
    
    print("End of genome Extraction @genome_step1 to get it")
    return_data <- c(return_data, genome_step1 = list(genome_maf))
    
    nb_individuals <- dim(genome_simu_step1)[1]
    return_data <- c(return_data, nb_individuals = nb_individuals)
    
  }
  
  if (file.exists(filename_m2)){
    print("Obtaining m2 from previous conversion")
    index_m2_step1 <- read.table(filename_m2)[,1]
    
    print("$m2_step1 to access it")
    return_data <- c(return_data, m2_step1 = list(index_m2_step1))
  }else{
    if (file.exists(filename_m2_step1)){
      print("Extracting mutationm2_step1 ...")
      mutationm2_step1 <- read.delim(file=filename_m2_step1, sep=" ", header=FALSE)
      index_m2_step1 = seq(1,dim(genome_vcf_step1)[1])[ genome_vcf_step1[,2] %in% (mutationm2_step1[,7] + 1)]
      print("$m2_step1 to access it")
      # Reindexation
      index_m2_step1 <- snp_reindexation(above_maf, index_m2_step1)
      write.table(index_m2_step1, filename_m2)
      return_data <- c(return_data, m2_step1 = list(index_m2_step1))
    }
    
  }
  
  if (file.exists(filename_m3)){
    print("Obtaining m3 from previous conversion")
    index_m3_step1 <- read.table(filename_m3)[,1]
    
    print("$m3_step1 to access it")
    return_data <- c(return_data, m3_step1 = list(index_m3_step1))
  }else{
    if (file.exists(filename_m3_step1)){
      print("Extracting mutationm3_step1 ...")
      mutationm3_step1 <- read.delim(file=filename_m3_step1, sep=" ", header=FALSE)
      index_m3_step1 = seq(1,dim(genome_vcf_step1)[1])[ genome_vcf_step1[,2] %in% (mutationm3_step1[,7] + 1)]
      print("$m3_step1 to access it")
      # Reindexation
      index_m3_step1 <- snp_reindexation(above_maf, index_m3_step1)
      write.table(index_m3_step1, filename_m3)
      return_data <- c(return_data, m3_step1 = list(index_m3_step1))
    }
  }

  
  
  if (file.exists(filename_position_step1)){
    print("Extracting position step 1 of individuals ...")
    pos_indiv_step1 = read.table(file = filename_position_step1, header=FALSE)
    pos_step1 <- convert_pos_indiv(pos_indiv_step1)
    print("$pos_step1 to access it")
    return_data <- c(return_data, pos_step1 = list(pos_step1))
  }
  

  
  if (file.exists(filename_var1_step1)){
    print("Extracting var1_step1 .. ")
    var1_step1 = read.table(file = filename_var1_step1, header=FALSE)
    print("$var1_step1 to acces it")
    return_data <- c(return_data, var1_step1 = list(var1_step1))
  }
  

  
  if (file.exists(filename_var1_pred)){
    print("Extracting var1_pred .. ")
    var1_pred = read.table(file = filename_var1_pred, header=FALSE)
    print("$var1_pred to acces it")
    return_data <- c(return_data, var1_pred = list(var1_pred))
  }
  
  
  if (file.exists(filename_var2_step1)){
    print("Extracting var2_step1 .. ")
    var2_step1 = read.table(file = filename_var2_step1, header=FALSE)
    print("$var1_step2 to acces it")
    return_data <- c(return_data, var2_step1 = list(var2_step1))
  }
  

  
  if (file.exists(filename_var2_pred)){
    print("Extracting var2_pred .. ")
    var2_pred = read.table(file = filename_var2_pred, header=FALSE)
    print("$var2_pred to acces it")
    return_data <- c(return_data, var2_pred = list(var2_pred))
  }
  
  
  if (file.exists(filename_fitness_step1)){
    print("Extracting fitness step 1..")
    fitness_step1 = read.table(file = filename_fitness_step1, header=FALSE)
    print("$fitness_step1 to access it")
    return_data <- c(return_data, fitness_step1 = list(fitness_step1))
  }
  
  
  if (file.exists(filename_fitness_pred)){
    print("Extracting fitness pred")
    fitness_pred = read.table(file = filename_fitness_pred, header=FALSE)
    print("$fitness_pred to access it")
    return_data <- c(return_data, fitness_pred = list(fitness_pred))
  }
  
  return(return_data)
}


maf_filtering <- function(genome, maf_threshold){
  upper_maf_threshold <- 1 - maf_threshold
  colmean_genome <- colMeans(genome) / 2
  above_maf <- (colmean_genome > maf_threshold) & (colmean_genome < upper_maf_threshold)
  genome_maf <- genome[,above_maf]
  
  return(list(genome_maf = genome_maf, above_maf=above_maf))
}

convert_pos_indiv <- function(pos_indiv){
  
  nb_individuals <- dim(pos_indiv)[1] / 2
  index_x <- 2*seq(0,(nb_individuals-1)) + 1
  index_y <- 2*(seq(0,(nb_individuals-1)) + 1)
  
  x_pos <- pos_indiv[index_x,1]
  y_pos <- pos_indiv[index_y,1]
  return(list(x=x_pos, y=y_pos))
}

align_genome <- function(genome_vcf_1, genome_vcf_2){
  mut_id_1 <- sapply(seq(1,dim(genome_vcf_1)[1]), function(x) strsplit(as.character(genome_vcf_1[,8]),";")[[c(x,1)]])
  mut_id_2 <- sapply(seq(1,dim(genome_vcf_2)[1]), function(x) strsplit(as.character(genome_vcf_2[,8]),";")[[c(x,1)]])
  
  return(list(mut1 = mut_id_1 %in% mut_id_2, mut2 = mut_id_2 %in% mut_id_1))
}

snp_reindexation <- function(truefalselist, snps_set){
  new_causal_set <- c()
  for (snps in snps_set){
    if (truefalselist[snps]){
      new_causal_set <- c(new_causal_set, sum(truefalselist[1:snps]))
    }
  }
  return(new_causal_set)
}



pos_to_pop <- function(x_pos, y_pos, vertical_length, vertical_step, horizontal_length, horizontal_step){
  
  
  length_vertical_step <- vertical_length / (vertical_step - 1)
  y_pos <- y_pos - (length_vertical_step / 2)
  
  
  # We handle edge case where y_pos is out of its bound
  y_pos <- max(0, y_pos)
  y_pos <- min(y_pos, vertical_length - (length_vertical_step + 1e-3))
  y_comp <- y_pos %/% (length_vertical_step)
  
  length_horizontal_step <- horizontal_length / (horizontal_step - 1)
  x_pos <- x_pos - (length_horizontal_step / 2)
  
  x_pos <- max(0, x_pos)
  x_pos <- min(x_pos, horizontal_length - (length_horizontal_step + 1e-3))
  
  x_comp = x_pos %/% length_horizontal_step + 1
  
  
  return((horizontal_step-2) * y_comp + x_comp)
}

get_pop <- function(pos, vertical_length, vertical_step, horizontal_length, horizontal_step){
  nb_ind <- length(pos$x)
  pop <- c()
  for (i in seq(1,nb_ind)){
    x_pos <- pos$x[i]
    y_pos <- pos$y[i]
    pop <- c(pop, pos_to_pop(x_pos, y_pos, vertical_length, vertical_step, horizontal_length, horizontal_step))
  }
  
  return(pop)
}

get_mean_by_pop <- function(var, indiv_pop){
  
  nb_pop <- length(table(indiv_pop))
  if (is.null(dim(var))){
    mean_var <- sapply(seq(1,nb_pop), function(x) mean(var[indiv_pop==x]))
    return(mean_var)
  }else{
    if (dim(var)[2] == 1){
      mean_var <- sapply(seq(1,nb_pop), function(x) mean(var[indiv_pop==x,]))
      return(mean_var)
    }
    mean_var <- sapply(seq(1,nb_pop), function(x) colMeans(var[indiv_pop==x,]))
    return(t(mean_var))
  }
  
}

sample_individuals <- function(data, nb_ind_sampling, pop){
  min_nb_ind <- min(table(pop))
  if (nb_ind_sampling > min_nb_ind){
    stop("nb_ind is greater than nb of individuals in some population, it is thus impossible to subsample")
  }
  
  nb_pop <- length(unique(pop))
  
  # sample genome, var1
  genome <- c()
  var1 <- c()
  var1_pred <- c()
  var2_exists <- ('var2_step1' %in% names(data))
  
  if (var2_exists){
    var2 <- c()
    var2_pred <- c()
  }
  
  fitness <- c()
  fitness_pred <- c()
  pos <- list()
  pos$x <- c()
  pos$y <- c()
  nb_ind_par_pop <- table(pop)
  
  new_pop <- c()
  
  for (i in seq(1, nb_pop)){
    current_nb_ind <- nb_ind_par_pop[i]
    current_genome <- data$genome_step1[pop==i,]
    current_var1 <- as.matrix(data$var1_step1[pop==i,])
    current_var1_pred <- as.matrix(data$var1_pred[pop==i,])
    current_fitness <- as.matrix(data$fitness_step1[pop==i,])
    current_fitness_pred <- as.matrix(data$fitness_pred[pop==i,])
    current_pos_x <- data$pos_step1$x[pop==i]
    current_pos_y <- data$pos_step1$y[pop==i]
    current_sampling <- sort(sample(seq(1,current_nb_ind), nb_ind_sampling))
    
    genome <- rbind(genome, current_genome[current_sampling,])
    var1 <- rbind(var1, as.matrix(current_var1[current_sampling,]))
    var1_pred <- rbind(var1_pred, as.matrix(current_var1_pred[current_sampling,]))
    fitness <- rbind(fitness, as.matrix(current_fitness[current_sampling,]))
    fitness_pred <- rbind(fitness_pred, as.matrix(current_fitness_pred[current_sampling,]))
    pos$x <- rbind(pos$x, current_pos_x[current_sampling])
    pos$y <- rbind(pos$y, current_pos_y[current_sampling])
    
    if (var2_exists){
      current_var2 <- as.matrix(data$var2_step1[pop==i,])
      current_var2_pred <- as.matrix(data$var2_pred[pop==i,])
      
      var2 <- rbind(var2, as.matrix(current_var2[current_sampling,]))
      var2_pred <- rbind(var2_pred, as.matrix(current_var2_pred[current_sampling,]))
    }
    
    new_pop <- c(new_pop, rep(i, nb_ind_sampling))
    
  }
  
  data$genome_step1 <- genome
  data$var1_step1 <- var1
  data$var1_pred <- var1_pred
  data$pos_step1 <- pos
  data$fitness_step1 <- fitness
  data$fitness_pred <- fitness_pred
  
  if (var2_exists){
    data$var2_step1 <- var2
    data$var2_pred <- var2_pred
  }
  
  return(list(data=data, pop=new_pop))
  
}

freq_by_pop <- function(genome, indiv_pop){
  nb_pop <- length(table(indiv_pop))
  genome_freq <- c()
  for (i in seq(1,nb_pop)){
    current_genome <- genome[indiv_pop == i,]
    current_genome_freq <- colMeans(current_genome)/2
    genome_freq <- rbind(genome_freq, current_genome_freq)
  }
  return(genome_freq)
  
}

get_fitness_variation <- function(data, pop){
  fitness_variation <- data$fitness_pred - data$fitness_step1
  mean_fit_var <- get_mean_by_pop(fitness_variation, pop)
  return(mean_fit_var)
}


