paquetes <- c("terra", "sp", "sf", "dplyr", "ranger", "tuneRanger", 
              "caret", "pROC", "pdp", "tidyr", "dismo", "usdm", "blockCV")

# Función para instalar y cargar automáticamente
instalar_y_cargar <- function(paq) {
  if (!require(paq, character.only = TRUE)) {
    message("Instalando paquete: ", paq, "...")
    install.packages(paq, dependencies = TRUE, repos = "https://cran.rstudio.com/")
    library(paq, character.only = TRUE)
  } else {
    message("✓ ", paq, " ya está instalado")
  }
}

cat("=== Verificando e instalando paquetes necesarios ===\n\n")
sapply(paquetes, instalar_y_cargar)

cat("\n¡Todos los paquetes están listos!\n")


rm(paquetes)
rm(instalar_y_cargar)
#### Función para cargar y combinar archivos TIF ######
combine_tif_files <- function(folder_path) {
  # Obtener la lista de archivos TIF en la carpeta
  tif_files <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  # Verificar si hay archivos TIF
  if (length(tif_files) == 0) {
    stop("No se encontraron archivos TIF en la carpeta especificada.")
  }
  
  # Cargar todos los archivos TIF como una lista de SpatRaster
  raster_list <- lapply(tif_files, terra::rast)
  
  # Obtener los nombres de los archivos (sin la ruta ni la extensión)
  layer_names <- tools::file_path_sans_ext(basename(tif_files))
  
  # Combinar los rasters en un solo SpatRaster
  combined_raster <- terra::rast(raster_list)
  
  # Asignar los nombres originales a las capas
  names(combined_raster) <- layer_names
  
  # Devolver el raster combinado
  return(combined_raster)
}


###### preparar dataframe para SDM 1 especie #########
SDM_prep<- function(occurrences,env,dist_bkg = 500000, seed = 666){
  
  data = occurrences
  ## occs to sf
  occs.sf <- sf::st_as_sf(data, coords = c("lon","lat"), crs = terra::crs(env))
  
  ##change projection for calculus
  eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)
  
  #### making a buffer ####
  occs.buf <- sf::st_buffer(occs.sf, dist = dist_bkg) |> 
    sf::st_union() |> 
    sf::st_sf() |>
    sf::st_transform(crs = terra::crs(env))
  
  ### cutting the raster with the buffer ###
  
  env.buf <- terra::crop(env, occs.buf)
  
  env.buf <- terra::mask(env.buf, occs.buf)
  
  ### to dataframe
  env.buf.df<-terra::as.data.frame(env.buf, xy = T)
  names(env.buf.df)[1:2]<-c("lon","lat")
  
  ####### making a mask ########
  
  mask <- env.buf[[1]]  * 0 + 1
  
  clls <- terra::extract(mask, data[,c("lon","lat")], cells = TRUE)
  
  
  vect <- duplicated(clls$cell)
  
  data$cell<-clls$cell
  
  dupv <- data[vect,] ## duplicados
  nodp <- data[!vect,] ##no duplicados
  
  
  dupv <- dupv[,c('lat', 'lon','cell')] %>% mutate(class = 'Duplicado')
  nodp <- nodp[,c('lat', 'lon','cell')] %>% mutate(class = 'Limpio')
  
  clls <- as_tibble(clls)
  
  ##juntamos todo de nuevo
  occr.2 <- rbind(nodp,dupv)
  
  ##numero de datos duplicados y datos limpios
  table(occr.2$class)
  
  occr.2 <- mutate(occr.2, class = factor(class, levels = c('Limpio', 'Duplicado')))
  
  occr.2 = subset(occr.2,occr.2$class == "Limpio")
  
  
  clls <- terra::extract(mask, occr.2[,c("lon","lat")], cells = TRUE)
  
  mask[clls$cell] <- NA 
  
  ##sampling pseudo-aubsences (backgrounds)
  set.seed(seed)
  bckg <- as_tibble(terra::spatSample(x = mask, 
                                      size = nrow(clls), 
                                      na.rm = TRUE, 
                                      replace = FALSE, 
                                      as.df = TRUE, 
                                      xy = TRUE, 
                                      warn = TRUE,
                                      values = TRUE))
  
  colnames(bckg)[1:2]<-c("lon","lat")
  
  ### combine pseudo-absences, real presences and simulated presences 
  occr.2<-mutate(occr.2, status = "presence", p = 1)
  
  occr.2 <- occr.2[,-c(3,4)] 
  
  
  bckg<-bckg[,-3]
  
  bckg<-mutate(bckg,status="absence" ,p = 0)
  
  occr.2<-rbind(occr.2,bckg)
  
  model_data <- terra::extract(env.buf, occr.2[,c("lon","lat")], df = TRUE, ID = FALSE)
  
  model_data
  
  ### VIF
  vif=vifstep(model_data,th=10)
  goodvars=vif@results$Variables
  
  env.buf <- subset(env.buf, goodvars) ## only VIF<10
  
  ####### extracting data again #######
  
  model_data <- terra::extract(env.buf, occr.2[,c("lon","lat")], df = TRUE, ID = FALSE)
  
  model_data<-cbind(occr.2,model_data)
  
  model_data$p<-as.factor(model_data$p)
  model_data<- na.omit(model_data)
  
  print(vif)
  return(list(SDM_df = model_data, env = env.buf))
}

###### calcular autocorrelacion espacial ######

spat_ran<- function(env_buff, samples = 5000, seed = 666){
  #### SDM_buff es un raster resultado de la funcion SDM_prep
  
  #### spatial autocorrelation 
  set.seed(666)
  au_corr <- cv_spatial_autocor(r = env_buff, 
                                num_sample = 5000)
  
  spat_ran<- median(au_corr$plots$barchart$data$range)  ### spatial block size for cross validation
  
  if (spat_ran < 100000) {
    print("We recommend a spatial range of 100000")
  } else {
    print(glue("We recommend a spatial range of {spat_ran}"))
  }
}

###### validacion cruzada ########

SDM_cross<-function(SDM_df,env_buff,spat_range = 100000,folds = 4,seed = 666){
  
  points<-st_as_sf(SDM_df,coords = c("lon", "lat"), crs = terra::crs(env_buff))
  
  ######## spatial blocks #########
  set.seed(seed)
  cross_valid <- cv_spatial(
    x = points,
    column = "p", # columna binaria (binary or multi-class)
    r = env_buff,
    k = folds, # numero de folds
    size = spat_range, # tamaño de los bloques
    selection = "random", # blockes aleatorios
    iteration = 500, # find evenly dispersed folds
    progress = T, # turn off progress bar
    biomod2 = FALSE, 
    raster_colors = terrain.colors(10, rev = TRUE), # options from cv_plot for a better colour contrast
    plot = F
    ) 
  
  return(cross_valid)
}

##### preparacion de modelos #######



RF_SDM<-function(data,folds,number_trees = 1000,seed = 666){
  
  ### lista de modelos
  model_list<-vector("list", length(folds)) ##folders
  
  
  ### testing
  test_table<-as.data.frame(data[,4])
  colnames(test_table)[1]<-"p"
  test_table$predict<-NA
  
  ##### random forest tuning #####
  #### we are going to tuning with all data
  
  #model_data <- na.omit(model_data)
  
  task=makeClassifTask(data = as.data.frame(data[,-c(1:3)]),target = "p")
  task$type
  set.seed(123)
  tun<-tuneRanger(task = task,measure = list(acc),num.trees = number_trees,
                  iters = 100)
  
  ###hyperparameters
  tun$recommended.pars
  tun$results
  tun<-tun$recommended.pars
  tun
  mtry<-tun$mtry
  sample<-tun$sample.fraction
  node_size<-tun$min.node.size
  
  tune_params<-list(mtry = mtry,
                    sample = sample,
                    node_size = node_size)
  
  model_data=data[,-c(1,2,3)]
  
  ####the SDM (cross validation) ######
  set.seed(seed)
  for(k in seq_len(length(folds))){
    # extracting the training and testing indices
    # this way works with folds_list list (but not folds_ids)
    trainSet <- unlist(folds[[k]][1]) # training set indices; first element
    testSet <- unlist(folds[[k]][2]) # testing set indices; second element
    rf <- ranger(p ~ ., na.omit(model_data[trainSet, ]), 
                 importance = "permutation",
                 num.trees = number_trees, ##number of trees
                 mtry = mtry, ##tree size
                 min.node.size = node_size,
                 sample.fraction = sample, ###sample fraction
                 classification = T,
                 probability = T) # model fitting on training set
    test_table$predict[testSet] <- predict(rf, model_data[testSet, ])$predictions[, 2] # predict the test set
    model_list[[k]]<-rf
  }
  
  test_table$fold<-NA
  
  ####adding folds to the test dataframe
  
  for (i in seq_along(folds)) {
    # extract rows of the fods
    fold_indices <- unlist(folds[[i]][2])
    
    # vrify that fold_indices is a numeric vector
    if (!is.numeric(fold_indices)) {
      stop("fold_indices are not a numeric vector")
    }
    
    # set the fold number to the rows
    test_table$fold[fold_indices] <- i
  }
  
  test_table$status<-model_data$status
  
  auc_all_data<-NULL
  
  for(i in 1:4){
    auc_all_data[i]<-auc(roc(subset(test_table,test_table$fold==i)$p,as.numeric(round(subset(test_table,test_table$fold==i)$predict))))
  }
  
  models<-list(models = model_list,
               tune = tune_params,
               auc = auc_all_data,
               test_table = test_table)
}




####### confussion matrix 4 plots ########

SDM_conf_4plots<-function(lista_SDM,name = "SDM_confussion_plot.png"){
  
  auc_all_data<-lista_SDM$auc
  
  png(glue::glue("{name}.png"), width = 10, height = 7, res = 300, units = "in")
  
  # Set up a 2x2 plotting grid with adjusted inner margins
  par(mfrow = c(2, 2), mar = c(4, 6, 4, 2) + 0.1, oma = c(0, 0, 0, 0))  # Increased left margin (mar[2])
  
  # Loop through the 4 folds
  for (k in 1:4) {
    # Subset data for the current fold
    actual <- lista_SDM$test_table$p[lista_SDM$test_table$fold == k]
    predicted <- lista_SDM$test_table$predict[lista_SDM$test_table$fold == k]
    
    # Create confusion matrix
    all_classes <- levels(factor(lista_SDM$test_table$p))  # All possible classes
    conf_matrix <- table(
      Actual = factor(actual, levels = all_classes),
      Predicted = factor(round(predicted,0), levels = all_classes)
    )
    
    # Plot the confusion matrix as a heatmap
    image(
      x = 1:ncol(conf_matrix),
      y = 1:nrow(conf_matrix),
      z = t(conf_matrix),  # Transpose for correct orientation
      col = gray.colors(10, start = 1, end = 0),  # Grayscale for heatmap
      xlab = "Predicted",
      ylab = "Label",
      main = paste("Fold", k),
      xaxt = "n", yaxt = "n"  # Suppress default axes
    )
    
    # Add axis labels (class names)
    axis(1, at = 1:ncol(conf_matrix), labels = colnames(conf_matrix), las = 2)
    axis(2, at = 1:nrow(conf_matrix), labels = rownames(conf_matrix), las = 2)
    
    # Add text to show the counts in the confusion matrix
    for (i in 1:nrow(conf_matrix)) {
      for (j in 1:ncol(conf_matrix)) {
        text(j, i, conf_matrix[i, j], col = "red", font = 2)  # Red text for counts
      }
    }
    
    test_acc <- round(auc_all_data[[k]],2)
    
    mtext(
      text = sprintf("Test AUC: %.2f", test_acc),
      side = 3, line = 1, outer = FALSE, cex = 0.8, adj = 0
    )
  }
  
  dev.off()
  
}

##### importancia relativa promediada y normalizada ######


SDM_imp <- function(RF_list, name = "variable_importance_plot"){
  ##### RF_list es una lista producto de RF_SDM
  
  importance_list <- lapply(RF_list, function(model) {
    as.data.frame(model$variable.importance) %>%
      tibble::rownames_to_column(var = "variable") %>%
      rename(importance = 2)
  })
  
  importance_list <- lapply(RF_list, function(model) {
    imp <- model$variable.importance
    imp_normalized <- imp / sum(imp)  # Normalizar
    as.data.frame(imp_normalized) %>%
      tibble::rownames_to_column(var = "variable") %>%
      rename(importance = 2)
  })
  
  # Paso 2: Combinar las importancias en un solo dataframe
  importance_df <- bind_rows(importance_list, .id = "model") %>%
    group_by(variable) %>%
    summarise(
      mean_importance = mean(importance),
      sd_importance = sd(importance),
      # Calcular intervalo de confianza (95%)
      lower_ci = mean_importance - 1.96 * sd_importance / sqrt(10),  # 10 modelos
      upper_ci = mean_importance + 1.96 * sd_importance / sqrt(10)
    ) %>%
    # Ordenar de mayor a menor por importancia promedio
    arrange((mean_importance)) %>%
    # Convertir la variable a factor para mantener el orden en el gráfico
    mutate(variable = factor(variable, levels = variable))
  
  # Ordenar de mayor a menor importancia
  importance_df <- importance_df[order(importance_df$mean_importance, decreasing = TRUE), ]
  importance_df$variable <- factor(importance_df$variable, levels = importance_df$variable)
  
  # Posiciones en el eje Y (1 a n, de arriba hacia abajo)
  y_pos <- length(importance_df$variable):1
  
  png(glue::glue("{name}.png"), units = "in", width = 10, height = 7, res = 300)
  # Gráfico horizontal
  plot(importance_df$mean_importance, y_pos,
       xlim = c(0, max(importance_df$upper_ci) * 1.1),  # espacio extra a la derecha
       ylim = c(0.5, length(y_pos) + 0.5),
       yaxt = "n", xlab = "Relative Importance", ylab = "",
       pch = 19, col = "steelblue", cex = 1.3)
  
  # Bigotes (intervalos de confianza)
  segments(importance_df$lower_ci, y_pos, importance_df$upper_ci, y_pos,
           col = "darkgray", lwd = 2)
  
  # Puntos centrales
  points(importance_df$mean_importance, y_pos, pch = 19, col = "steelblue", cex = 1.3)
  
  # Etiquetas de variables en el eje Y
  axis(2, at = y_pos, labels = importance_df$variable, las = 1, cex.axis = 0.9)
  
  # Línea vertical en 0
  abline(v = 0, col = "gray", lty = 2)
  # Opcional: valores numéricos al lado de los puntos
  text(importance_df$mean_importance, y_pos,
       labels = round(importance_df$mean_importance, 3),
       pos = 4, cex = 0.85, col = "black", font = 1)
  dev.off()
  
  return(importance_df)
}


SDM_pdp<- function(RF_list,RF_df, pred_var = NULL, dividir = 1, real_var_name = NULL){
  
  pdp_list <- vector("list", length(RF_list))
  
  for (i in 1:length(RF_list) ){
    pd = partial(RF_list[[i]],
                 pred.var = pred_var,
                 train = RF_df[,-c(1:3)],
                 type = "classification",
                 prob = TRUE,
                 which.class=1)
    
    pdp_list[[i]]<-pd
  }
  
  pdp<- do.call(rbind, pdp_list)
  
  summary_pdp<- pdp %>%
    group_by(pdp[,1]) %>% # Agrupa por valores de la variable independiente
    summarise(
      y_promedio = mean(yhat),
      y_sd = sd(yhat), # Desviación estándar
      ci_lower = y_promedio - qt(0.975, n()-1) * y_sd / sqrt(n()),
      ci_upper = y_promedio + qt(0.975, n()-1) * y_sd / sqrt(n())
    )
  
  colnames(summary_pdp)[1] <- real_var_name
  
  summary_pdp[,1] <-  summary_pdp[,1]/dividir
  return(summary_pdp)
}

##### prediciendo nicho en el espacio #########

SDM_spatial_pred<-function(RF_list, occurrences, spat_rast, pol){
  ## RF_list is a list product of RF_SDM
  ## occurrences, there are a product of SDM_prep
  ## spat_rast is a raster in a terra format (needs the same variables 
  ## of your models)
  ## pol is a polygon to cut where the models should predict
  
  ## the result is a list with two dataframes: 
  ## a DF predicting in the all raster
  ## a DF predicting only in the raster with a positive MESS (MESS >= 0)
  ## a polygon of the MESS with only negative values
  
  ## in this manner, you can see where your model is extrapolating and 
  ## be cautious with your interpretations
  ## this do not make a map, because each person can custom it as they want
  
  ##### first cut with the polygon #######
  poly_rast <- terra::crop(spat_rast, pol)
  poly_rast <- terra::mask(spat_rast, pol)
  
  #### mess and mess mask ######
  mess_map <- dismo::mess(x = raster::stack(poly_rast), v = occurrences[,-c(1:4)], full = FALSE)
  mess_map<-terra::rast(mess_map)
  mess_mask <- mess_map >= 0
  
  
  #### second cut rast #####
  poly_rast_2 <- mask(poly_rast, mess_mask, maskvalues = FALSE)
  
  #### poly mess
  poly_mess <- as.polygons(mess_map < 0,
                           na.rm = TRUE,
                           aggregate = TRUE,
                           dissolve = TRUE)
  
  poly_mess<- st_as_sf(poly_mess)
  
  poly_mess <- poly_mess %>% 
    filter(mess == 1)
  
  #### df 1 with all the polygon ######
  df_1<-terra::as.data.frame(poly_rast, xy = T)
  colnames(df_1)[1:2] <- c('lon', 'lat')
  
  
  #### df 2 with mess polygon
  df_2<-terra::as.data.frame(poly_rast_2, xy = T)
  colnames(df_2)[1:2] <- c('lon', 'lat')
  
  ####predicting niche 1 #####
  
  preds_1<-matrix(NA,ncol=length(RF_list$models),nrow=nrow(df_1))
  
  preds_1<-as.data.frame(preds_1)
  
  set.seed(666)
  
  for(i in 1:4){
    model_predictions<-predict(RF_list$models[[i]],df_1[,-c(1,2)])$predictions[, 2]
    preds_1[,i]<-model_predictions
  }
  
  mean_preds_1 <- apply(preds_1, 1, mean)
  
  sd_model_1 <- apply(preds_1, 1, sd)
  
  coords_1<- df_1[,c(1,2)]
  
  predictions_df<-data.frame(lon = coords_1$lon,
                             lat = coords_1$lat, 
                             probs = mean_preds_1,
                             sd = sd_model_1)
  #### prediciendo nicho poligono 2 ######
  
  preds_2<-matrix(NA,ncol=length(RF_list$models),nrow=nrow(df_2))
  
  preds_2<-as.data.frame(preds_2)
  
  set.seed(666)
  
  for(i in 1:4){
    model_predictions<-predict(RF_list$models[[i]],df_2[,-c(1,2)])$predictions[, 2]
    preds_2[,i]<-model_predictions
  }
  
  
  mean_preds_2 <- apply(preds_2, 1, mean)
  
  sd_model_2 <- apply(preds_2, 1, sd)
  
  coords_2<- df_2[,c(1,2)]
  
  predictions_df_2<-data.frame(lon = coords_2$lon,
                               lat = coords_2$lat, 
                               probs = mean_preds_2,
                               sd = sd_model_2)
  
  list_predictions<-list(predictions = predictions_df,
                         mess_predictions = predictions_df_2,
                         mess_polygon = poly_mess,
                         df_all = df_1,
                         df_mess = df_2)
  
  return(list_predictions)
}