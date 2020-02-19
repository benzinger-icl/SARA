# SARA Classification Model by Lauren Koenig
# 8/15/2019
# Takes volumes and cortical thicknesses as calculated by FreeSurfer
# Corrects volumes based on Intracranial volumes regression approach
# Calculates SARA predicted probability of AD, and groups based upon diagnoses given to previous clinical patients
# Optionally creates a csv of the ICV-corrected volumes and thicknesses, along with SARA probability and grouping
# Optionally creates a graph for each subject describing the SARA probability and grouping (for clinical use)
# The SARA probability indicates the liklihood 0-100% probability) that the observed pattern of atrophy indicates clinical Alzheimer disease
# The SARA grouping is intended to aid in interpretation by indicating if the SARA probability falls within 80% or 90% sensitivity or specificity 
# Intended primarily for patients age 55-85, younger patients more likely to have false negatives
# make graphs only if single patient (can have multiple MR's)
#only uses MRs within age range to calculate slope

SARA <- function(patient_data_location, normals_data_location, patient_data_source = "Patient", output_directory=NA, make_csv=TRUE, make_graphs=FALSE, 
                              ID_name = "CLINICALDATA.ID", age_label = "Age", 
                              IntraCranialVol = "IntraCranialVol",
                              Left.Amygdala_volume = "Left.Amygdala_volume", Right.Amygdala_volume="Right.Amygdala_volume",
                              Left.Hippocampus_volume="Left.Hippocampus_volume", Right.Hippocampus_volume="Right.Hippocampus_volume",
                              Left.Inf.Lat.Vent_volume="Left.Inf.Lat.Vent_volume", Right.Inf.Lat.Vent_volume="Right.Inf.Lat.Vent_volume",
                              lh_entorhinal_thickness="lh_entorhinal_thickness", rh_entorhinal_thickness="rh_entorhinal_thickness",
                              lh_inferiorparietal_thickness="lh_inferiorparietal_thickness", rh_inferiorparietal_thickness="rh_inferiorparietal_thickness") 
  {

  #########################
  # SET UP
  #########################
  #read in R environment for normals data
  load(normals_data_location)
  
  #Read in the data
  patient <- read.csv(as.character(patient_data_location)) 
  
  #rename subject ID label if necessary
  if ( !("CLINICALDATA.ID" %in% colnames(patient) )) {
    colnames(patient)[colnames(patient)==ID_name] <- "CLINICALDATA.ID"
  }  
  #rename Age label if necessary
  if (!("Age" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==age_label] <- "Age"
  }
  #rename region names if necessary
  if (!("IntraCranialVol" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==IntraCranialVol] <- "IntraCranialVol"
  }
  if (!("Left.Amygdala_volume" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==Left.Amygdala_volume] <- "Left.Amygdala_volume"
  }
  if (!("Right.Amygdala_volume" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==Right.Amygdala_volume] <- "Right.Amygdala_volume"
  }
  if (!("Left.Hippocampus_volume" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==Left.Hippocampus_volume] <- "Left.Hippocampus_volume"
  }
  if (!("Right.Hippocampus_volume" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==Right.Hippocampus_volume] <- "Right.Hippocampus_volume"
  }
  if (!("Left.Inf.Lat.Vent_volume" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==Left.Inf.Lat.Vent_volume] <- "Left.Inf.Lat.Vent_volume"
  }
  if (!("Right.Inf.Lat.Vent_volume" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==Right.Inf.Lat.Vent_volume] <- "Right.Inf.Lat.Vent_volume"
  }
  if (!("lh_entorhinal_thickness" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==lh_entorhinal_thickness] <- "lh_entorhinal_thickness"
  }
  if (!("rh_entorhinal_thickness" %in% colnames(patient) )){
    colnames(patient)[colnames(patient)==rh_entorhinal_thickness] <- "rh_entorhinal_thickness"
  }
  if (!("lh_inferiorparietal_thickness" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==lh_inferiorparietal_thickness] <- "lh_inferiorparietal_thickness"
  }
  if (!("rh_inferiorparietal_thickness" %in% colnames(patient) ) ){
    colnames(patient)[colnames(patient)==rh_inferiorparietal_thickness] <- "rh_inferiorparietal_thickness"
  }
  
  
  
  
  
  if ( !("CLINICALDATA.ID" %in% colnames(normal_aging_means) )) {
    colnames(normal_aging_means)[colnames(normal_aging_means)==ID_name] <- "CLINICALDATA.ID"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==ID_name] <- "CLINICALDATA.ID"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==ID_name] <- "CLINICALDATA.ID"
  }
  if ( !("Age" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==age_label] <- "Age"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==age_label] <- "Age"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==age_label] <- "Age"
  }
  if (!("IntraCranialVol" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==IntraCranialVol] <- "IntraCranialVol"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==IntraCranialVol] <- "IntraCranialVol"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==IntraCranialVol] <- "IntraCranialVol"
  }
  if ( !("Left.Amygdala_volume" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==Left.Amygdala_volume] <- "Left.Amygdala_volume"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==Left.Amygdala_volume] <- "Left.Amygdala_volume"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==Left.Amygdala_volume] <- "Left.Amygdala_volume"
  }
  if (  !("Right.Amygdala_volume" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==Right.Amygdala_volume] <- "Right.Amygdala_volume"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==Right.Amygdala_volume] <- "Right.Amygdala_volume"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==Right.Amygdala_volume] <- "Right.Amygdala_volume"
  }
  if (!("Left.Hippocampus_volume" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==Left.Hippocampus_volume] <- "Left.Hippocampus_volume"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==Left.Hippocampus_volume] <- "Left.Hippocampus_volume"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==Left.Hippocampus_volume] <- "Left.Hippocampus_volume"
  }
  if (  !("Right.Hippocampus_volume" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==Right.Hippocampus_volume] <- "Right.Hippocampus_volume"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==Right.Hippocampus_volume] <- "Right.Hippocampus_volume"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==Right.Hippocampus_volume] <- "Right.Hippocampus_volume"
  }
  if (  !("Left.Inf.Lat.Vent_volume" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==Left.Inf.Lat.Vent_volume] <- "Left.Inf.Lat.Vent_volume"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==Left.Inf.Lat.Vent_volume] <- "Left.Inf.Lat.Vent_volume"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==Left.Inf.Lat.Vent_volume] <- "Left.Inf.Lat.Vent_volume"
  }
  if (!("Right.Inf.Lat.Vent_volume" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==Right.Inf.Lat.Vent_volume] <- "Right.Inf.Lat.Vent_volume"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==Right.Inf.Lat.Vent_volume] <- "Right.Inf.Lat.Vent_volume"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==Right.Inf.Lat.Vent_volume] <- "Right.Inf.Lat.Vent_volume"
  }
  if (!("lh_entorhinal_thickness" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==lh_entorhinal_thickness] <- "lh_entorhinal_thickness"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==lh_entorhinal_thickness] <- "lh_entorhinal_thickness"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==lh_entorhinal_thickness] <- "lh_entorhinal_thickness"
  }
  if ( !("rh_entorhinal_thickness" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==rh_entorhinal_thickness] <- "rh_entorhinal_thickness"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==rh_entorhinal_thickness] <- "rh_entorhinal_thickness"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==rh_entorhinal_thickness] <- "rh_entorhinal_thickness"
  }
  if (  !("lh_inferiorparietal_thickness" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==lh_inferiorparietal_thickness] <- "lh_inferiorparietal_thickness"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==lh_inferiorparietal_thickness] <- "lh_inferiorparietal_thickness"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==lh_inferiorparietal_thickness] <- "lh_inferiorparietal_thickness"
  }
  if (  !("rh_inferiorparietal_thickness" %in% colnames(normal_aging_means) )){
    colnames(normal_aging_means)[colnames(normal_aging_means)==rh_inferiorparietal_thickness] <- "rh_inferiorparietal_thickness"
    colnames(normal_aging_sd)[colnames(normal_aging_sd)==rh_inferiorparietal_thickness] <- "rh_inferiorparietal_thickness"
    colnames(ICV_correction_coeff)[colnames(ICV_correction_coeff)==rh_inferiorparietal_thickness] <- "rh_inferiorparietal_thickness"
  }
  
  
  
  
  
  
  
  #prevents scientific notation
  options(scipen=999)
  
  
  
  
  
  
  ########################################
  # Intracranial Volume Correction 
  ########################################
  
  
  #separate volume and thickness measurements
  patient_thickness <- patient[grep("thickness", names(patient), value=FALSE, ignore.case=T)]
  patient_volume <- patient[grep("thickness", names(patient), value=FALSE, invert=TRUE, ignore.case=T)]
  
  #remove columns that aren't in both patient and normals (easy way to clean up / removes most demographic data)
  patient_thickness <- cbind(patient_thickness[, which(colnames(patient_thickness)%in% colnames(normal_aging_means))])
  patient_volume <- cbind(patient_volume[, which(colnames(patient_volume)%in% colnames(normal_aging_means))])
  
  #remove extra columns from volumes that aren't volumes (will need to adjust if particular csv file has additional columns added to it)
  patient_volume <- patient_volume[, !grepl(("CLINICALDATA.ID|Age|IntraCranialVol|WM_hypointensities|non_WM_hypointensities"), colnames(patient_volume), ignore.case = T)]

  #Correct for head volume size
  for (i in colnames(patient_volume)) { #goes through each ROI in patient_volume (which contains only volumes)
   patient_volume[i] <- data.frame( as.matrix(patient_volume[i]) - as.numeric(ICV_correction_coeff[i])*(patient$IntraCranialVol - Normals_ICV_mean )) #normalizes all the patient data based on ICV
    }
  

  patient_sub <- patient_volume
  
  #add non corrected columns back in to the dataframe
  for (i in colnames(patient[!(colnames(patient)%in% colnames(patient_volume))])) {
    patient_sub[i] <- patient[i]
  }
  

  
  ########################################
  #  Calculate Z-scores (if no graphs)
  ######################################## 
  if (make_csv == T && make_graphs == F) {  #if make graphs is false haven't yet calculated the z-scores, if true, don't want to repeat the loop
    
    #Set normals data's min and max ages
    n_min <- ceiling(min(normal_aging_means$Age, na.rm = TRUE))
    n_max <- floor(max(normal_aging_means$Age, na.rm = TRUE))
    
    #Set up if patient is within the age range covered by super normal cohort
    #creates true/false factor for each session if age at session is within normals (for z score and graph axis)
    patient$ages_within_normals <- patient$Age
    for (i in 1:length(patient$Age)){
      if ((patient$Age[i] >= min(normal_aging_means$Age, na.rm = TRUE)) && (patient$Age[i] <= max(normal_aging_means$Age, na.rm = TRUE))){
        patient$ages_within_normals[i] <- TRUE
      }
      else{
        patient$ages_within_normals[i] <- FALSE
      }
    }
    
    #Loop that goes through each region in subset matrix defined above
    for (Q in 1:(length(patient_sub))) {      
      
      ROI_name <- as.character(colnames(patient_sub[Q])) #still grabs correct column from both dataframes even if the column order gets messed up
      if (all(!grepl(ROI_name, c("CLINICALDATA.ID", "Age", "SARA_prediction", "SARA_category", "IntraCranialVol") ) )) {
        
        #Smooth population SD data with loess
        sd_fit <- loess(as.matrix(normal_aging_sd[ROI_name]) ~ c(n_min:n_max), degree=1, span=0.7)

        #Calculate z scores and associates percentiles for each of patient's MR data points
        region_z_score = matrix(nrow=1, ncol=length(patient$Age))
        for (i in 1:length(patient$Age)) {
          
          if (patient$ages_within_normals[i] == TRUE){
            #sometimes it doesn't want to grab the exact age for some reason? so I'm using >=, so will grab the next 0.01 age up and use that mean if that happens (usually <0.001 difference in z)
            region_z_score[i] = ((patient_sub[i,Q] - normal_aging_means[normal_aging_means$Age >= round(patient$Age[i], digits=2),ROI_name][1]) / #rounds to nearest 0.01 age for mean (because the saved means are for every 0.01 age)
              predict(sd_fit, newdata=patient$Age[i]))
             }else { 
            region_z_score[i] = NA  #leaves z-score blank if patient age is outside age range of normal group
          }
          patient_sub[paste(ROI_name,"_zscore", sep = "")] <- data.frame(t(region_z_score))  
        } #end of zscore calculation for loop going through each age
      }  #end of if ROI is actually a region loop
    } # end of Q for loop going through each region
  } #end of if statement selecting for csv T graphs F
  
  
  
  
  
  
  
  
  
  
  
  ########################################
  #  GRAPH Patient Volumes on Normal Aging Line (and calculate z-scores)
  ########################################
  
  if (make_graphs == T) { #making graphs is optional, intended for use by clinician
    packages <- c("ggplot2", "ggrepel")
    for( i in packages ){
      #  require returns TRUE invisibly if it was able to load package
      if( ! require( i , character.only = TRUE ) ){
        #  If package was not able to be loaded then re-install
        install.packages( i , dependencies = TRUE )
        #  Load package after installing
        require( i , character.only = TRUE )
      }
    }
    
    
    #Set working directory to where you want the png (graphs) and csv (table) files to go
    if (!is.na(output_directory)){
      setwd(output_directory) ###adjust to your file location where you want the output to go
    }
    
    
    #Set normals data's min and max ages
    n_min <- ceiling(min(normal_aging_means$Age, na.rm = TRUE))
    n_max <- floor(max(normal_aging_means$Age, na.rm = TRUE))
    
    #Set up if patient is within the age range covered by super normal cohort
    #creates true/false factor for each session if age at session is within normals (for z score and graph axis)
    patient$ages_within_normals <- patient$Age
    for (i in 1:length(patient$Age)){
      if ((patient$Age[i] >= min(normal_aging_means$Age, na.rm = TRUE)) && (patient$Age[i] <= max(normal_aging_means$Age, na.rm = TRUE))){
        patient$ages_within_normals[i] <- TRUE
      }
      else{
        patient$ages_within_normals[i] <- FALSE
      }
    }
    
    
    #Determines if using longitudinal data, specifically if 2+ data points within the normals age 
    if ( sum(patient$ages_within_normals) > 1 ){
      is_longitudinal <- TRUE
    } else {
      is_longitudinal <- FALSE
    }
    
    #Set up function to annotate graphs cleanly
    annotate_text2 <- function(label, x, y, facets=NULL, hjust=0, vjust=0, color='black', alpha=NA,
                               family=thm$text$family, size=thm$text$size, fontface=1, lineheight=1.0,
                               box_just=ifelse(c(x,y)<0.5,0,1), margin=unit(size/2, 'pt'), thm=theme_get()) {
      x <- scales::squish_infinite(x)
      y <- scales::squish_infinite(y)
      data <- if (is.null(facets)) data.frame(x=NA) else data.frame(x=NA, facets)
      
      tg <- grid::textGrob(
        label, x=0, y=0, hjust=hjust, vjust=vjust,
        gp=grid::gpar(col=alpha(color, alpha), fontsize=size, fontfamily=family, fontface=fontface, lineheight=lineheight)
      )
      ts <- grid::unit.c(grid::grobWidth(tg), grid::grobHeight(tg))
      vp <- grid::viewport(x=x, y=y, width=ts[1], height=ts[2], just=box_just)
      tg <- grid::editGrob(tg, x=ts[1]*hjust, y=ts[2]*vjust, vp=vp)
      inner <- grid::grobTree(tg, vp=grid::viewport(width=unit(1, 'npc')-margin*2, height=unit(1, 'npc')-margin*2))
      
      layer(
        data = NULL,
        stat = StatIdentity,
        position = PositionIdentity,
        geom = GeomCustomAnn,
        inherit.aes = TRUE,
        params = list(
          grob=grid::grobTree(inner), 
          xmin=-Inf, 
          xmax=Inf, 
          ymin=-Inf, 
          ymax=Inf
        )
      )
    }
    
    
    for (Q in 1:(length(patient_sub))) {      #Loop that goes through each region in subset matrix defined above
      
      ROI_name <- as.character(colnames(patient_sub[Q])) #still grabs correct column from both dataframes even if the column order gets messed up
      if (all(!grepl(ROI_name, c("CLINICALDATA.ID", "Age", "SARA_prediction", "SARA_category", "IntraCranialVol") ) )){
        
        #Smooth population data with loess
        smoothFit <- loess( as.matrix(normal_aging_means[ROI_name]) ~ normal_aging_means$Age, degree=1, span=0.05 ) #locally fits a polynomial surface
        smoothPred <- predict( smoothFit, newdata=n_min:n_max, se=TRUE ) #uses loess fit to predict data for each age, used for graph
        
        sd_fit <- loess(as.matrix(normal_aging_sd[ROI_name]) ~ c(seq(n_min,n_max, by=1)), degree=1, span=0.7)
        sd_pred <- predict(sd_fit, newdata=seq(n_min,n_max, by=1))

        plotFrame <- data.frame( Age=c(n_min:n_max), fit=smoothPred$fit, sd=sd_pred ) #combines all the model's data into single dataframe
        
        #Calculate z scores and associates percentiles for each of patient's MR data points
        region_z_score = matrix(nrow=1, ncol=length(patient$Age))
        region_percentile = matrix(ncol=length(patient$Age))
        
        #Calculate slopes for the patient and the normal group within the same age range as patient's data
        if (is_longitudinal == 1) {
          patient_longitudinal <- patient_sub[patient$ages_within_normals==1,] #makes dataset excluding any data points outside the normals age range
          patient_lm <- round( lm( as.matrix(patient_longitudinal[ROI_name]) ~ patient_longitudinal$Age )$coefficients, digits=4) #linear model of patient's data by age
          patient_slope <- paste("Patient Change per year:", patient_lm["patient_longitudinal$Age"] ) #taking slope of model and adding descriptive text
          normals_at_pat_age <- subset(plotFrame, Age <= max(patient_longitudinal$Age, na.rm = TRUE) & Age >= min(patient_longitudinal$Age, na.rm = TRUE)) #creating dataframe with the modeled normal data from patient's age at first MR to age at last MR
          normals_lm <- round( lm( normals_at_pat_age$fit ~ normals_at_pat_age$Age )$coefficients, digits=4) #making a linear model of the normals data just within the age range of the patient's longitudinal data
          normals_slope <- paste("Normal Group Change per year:", normals_lm["normals_at_pat_age$Age"])
        }

        for (i in 1:length(patient$Age)){
          
          if (patient$ages_within_normals[i] == TRUE){
            #sometimes it doesn't want to grab the exact age for some reason? so I'm using >=, so will grab the next 0.01 age up and use that mean if that happens (usually <0.001 difference in z)
            region_z_score[i] = ((patient_sub[i,Q] - normal_aging_means[normal_aging_means$Age >= round(patient$Age[i], digits=2),ROI_name][1]) / #rounds to nearest 0.01 age for mean (because the saved means are for every 0.01 age)
                                   predict(sd_fit, newdata=patient$Age[i]))
            region_percentile[i] = round( pnorm(region_z_score[i])*100, digits = 0 ) #Converts z score into a percentile
            
            if (region_percentile[i] != "NaN" && as.numeric(region_percentile[i]) < 1) { region_percentile[i]="1" } #defines percentile if below limit
            
            if (region_percentile[i] != "NaN" && as.numeric(region_percentile[i])>99) { region_percentile[i]="99" } #defines percentile if above limit
            
            if (region_percentile[i] != "NaN" && (as.numeric(region_percentile[i])>=1) && (as.numeric(region_percentile[i])<=99)){
              #print(paste("percentile is ", region_percentile[i])) #prints out the percentile to command line
              region_percentile[i] = paste(", ", region_percentile[i]," Percentile", sep="") #add a comma to the version being used in the display
            }
            
          } else { 
            region_percentile[i] = ""  #leaves percentile blank if patient is outside age range of normal group
          }
        }
        
        patient_sub[paste(ROI_name,"_zscore", sep = "")] <- data.frame(t(region_z_score))
        
        ########################################
        # Graph Making
        #########################################
        
        #ggplot2 options
        options <- theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18),
                         axis.title.x=element_text(size=18,face="bold"),
                         axis.title.y=element_text(size=18,angle=90,face="bold"),
                         legend.text=element_text(size=18),
                         legend.title=element_blank(),legend.position="right",
                         plot.title=element_text(size=18,face="bold",hjust = 0.5),legend.title.align=0.5,
                         plot.subtitle=element_text(size=16))
        
        #Make Plot elements from normals data
        region_plot <- ggplot(plotFrame,aes(x=Age,y=fit))
        region_smoothLine <- geom_line(aes(x=Age,y=fit,colour="Normals"),size=1.25)
        region_smoothRibbon2 <- geom_ribbon(aes(x=Age,ymax=fit+2*sd,ymin=fit-2*sd, fill="2 Std Dev"),alpha=0.1)
        region_smoothRibbon <- geom_ribbon(aes(x=Age,ymax=fit+sd,ymin=fit-sd, fill="1 Std Dev"),alpha=0.5)
        
        
        #Make Plot elements from patient data
        region_scatterNorm <- geom_point(data=patient_sub, aes(x=Age,y=eval(parse(text = ROI_name)), group="Patient", colour="Patient"), size=4.0) #graphs patient's datapoints
        region_lineNorm <- geom_line( data=patient_sub, aes(x=Age, y=eval(parse(text = ROI_name)), group="Patient", colour="Patient"),size=1.2) #connects patient's datapoints with lines
        
        
        #Add text on patient's data points and calculated percentiles
        #Set up the label
        
        display_label=paste(as.matrix(format(round(patient_sub[ROI_name],digits=2), nsmall=2)), t(region_percentile), sep="")
        
        #Use ggrepel to display the label so it doesn't get cut off
        region_volume_text <- geom_text_repel(data=patient_sub, mapping=aes(Age, eval(parse(text = ROI_name)), label=display_label),
                                              size=6, point.padding=unit(.5, "lines"), box.padding = unit(0.5, "lines"), force=5 )
        
        
        
        
        #Make Aesthetic Plot elements
        region_xlabel <- xlab("Age") #x label, doesn't change
        if (grepl("thickness", ROI_name, ignore.case=T)) {
          region_ylabel <- ylab(expression(bold(paste('Thickness (mm)'))))
        } else {
          region_ylabel <- ylab(expression(bold(paste('Normalized Volume  ',(mm^{3})))))
        }
        
        
        
        #Set y max
        a <- (plotFrame[order(plotFrame$fit + 2.5*plotFrame$sd, decreasing = T),][1,]$fit + 2.5*plotFrame[order(plotFrame$fit + 2.5*plotFrame$sd, decreasing = T),][1,]$sd) #normals max + 2.5SD (at the max)
        b <- (max(patient_sub[ROI_name], na.rm = TRUE) + 0.5*mean(plotFrame$sd, na.rm = TRUE)) # patient's y max + 0.5 average SD
        
        if (a > b) {  y_max <- a } else{ y_max <- b }
        
        #Set y min
        c <- (plotFrame[order(plotFrame$fit - 2.5*plotFrame$sd),][1,]$fit - 2.5*plotFrame[order(plotFrame$fit - 2.5*plotFrame$sd),][1,]$sd) #lowest normals data - 2.5 SD (at the min) 
        d <- (min( patient_sub[ROI_name], na.rm = TRUE) - 0.5*mean(plotFrame$sd, na.rm = TRUE)) #patient's data - 0.5 average SD
        
        if ( c < d  ) { y_min <- c } else {  y_min <- d }
        
        
        #set title
        c <- strsplit(gsub("\\.|_", " ", ROI_name, ignore.case = FALSE, perl = FALSE,
                           fixed = FALSE, useBytes = FALSE), " ")[[1]] 
        region_title <- ggtitle( paste(toupper(substring(c, 1,1)), substring(c, 2),
                                       sep="", collapse=" "), subtitle = NULL)
        
        
        
        #set title: if longitudinal include subtitle with slopes
        if (is_longitudinal == 1) { 
          c <- strsplit(gsub("\\.|_", " ", ROI_name, ignore.case = FALSE, perl = FALSE,
                             fixed = FALSE, useBytes = FALSE), " ")[[1]] 
          region_title <- ggtitle( paste(toupper(substring(c, 1,1)), substring(c, 2),
                                         sep="", collapse=" "), subtitle = paste(patient_slope,"\n",normals_slope))  
        } else {
          c <- strsplit(gsub("\\.|_", " ", ROI_name, ignore.case = FALSE, perl = FALSE,
                             fixed = FALSE, useBytes = FALSE), " ")[[1]] 
          region_title <- ggtitle( paste(toupper(substring(c, 1,1)), substring(c, 2),
                                         sep="", collapse=" "), subtitle = NULL)
          }
        
        region_customY <- scale_y_continuous(limits=c(y_min,y_max))
        
        #set x-axis values
        x_min <- n_min - 1
        x_max <- n_max
        
        if (max(patient$Age) >= x_max) {
          # if the max patient age is over the normals max age, use the patient age as the x axis max
          x_max <- max(patient$Age, na.rm = TRUE) 
        }
        if (min(patient$Age, na.rm = TRUE) <= n_min){
          # if the min patient age is below the normals min age, then use the patient age as the x axis min
          x_min <- min(patient$Age, na.rm = TRUE) 
        }
        
        region_customX <- scale_x_continuous(limits=c(5*floor(x_min/5), 5*ceiling(x_max/5)), breaks=(seq(from = (10*floor(x_min/10)), to = 10*ceiling(x_max/10), by = 10)))
        
        #set legend labels - do not mess with, reorders labels without clear cause  
        legend  <- scale_colour_manual(name="", labels=c("Average", "Patient"),  
                                       values = c("1 Std Dev" = "#F8766D", "Normals"="red3", "Patient"="#673dbd", "2 Std Dev"="lightblue"))  
        
        #Combine the different elements of the plot
        Fig <- (region_plot + theme_bw() + options  + region_smoothRibbon2 + region_smoothRibbon  + region_smoothLine
                + region_lineNorm + region_scatterNorm + region_volume_text 
                + region_xlabel + region_ylabel + region_title + region_customY + region_customX + 
                  legend)
        
        #Fig
        
        #make full filename
        file_name_full <- paste(patient_data_source, "_", ROI_name, ".png", sep="") #combines subject number and region name to make file name
        png(filename= file_name_full ,width=950,height=480) #makes png file
        print(Fig) #adds Fig to png file
        dev.off() #Closes the graphs generated by print(Fig)
        
      } #end of loop that starts in Model Fitting section
      
      
    }
  }
  
  
  ########################################
  #  Apply SARA Classification Model
  ########################################
  #combine volumes and average thicknesses across hemispheres
  Amygdala <- patient_sub$Left.Amygdala_volume + patient_sub$Right.Amygdala_volume
  Hippocampus <- patient_sub$Left.Hippocampus_volume + patient_sub$Right.Hippocampus_volume
  InfLatVent <- patient_sub$Left.Inf.Lat.Vent_volume + patient_sub$Right.Inf.Lat.Vent_volume
  Entorhinal_thickness <- (patient_sub$lh_entorhinal_thickness + patient_sub$rh_entorhinal_thickness)/2
  InferiorParietal_thickness <- (patient_sub$lh_inferiorparietal_thickness + patient_sub$rh_inferiorparietal_thickness)/2
  
  # Calculate SARA classification score
  patient_sub$SARA_prediction <- (exp(1)^( 14.4694972676 + -0.0005619816*Amygdala + -0.0006719734*Hippocampus +  0.0002468182*InfLatVent + -0.6800497755*Entorhinal_thickness + -3.2822467954*InferiorParietal_thickness) )  / 
    (1 + exp(1)^(14.4694972676 + -0.0005619816*Amygdala + -0.0006719734*Hippocampus +  0.0002468182*InfLatVent + -0.6800497755*Entorhinal_thickness + -3.2822467954*InferiorParietal_thickness) )
  
  #determine category which is based on a cohort of previous patients
  patient_sub$SARA_category <- NA
  for (i in 1:length(patient_sub$SARA_prediction)){
    if (patient_sub$SARA_prediction[i] < 0.29) {patient_sub$SARA_category[i] <- "Includes <10% of those with AD" }
    if (patient_sub$SARA_prediction[i] >= 0.29 & patient_sub$SARA_prediction[i] < 0.45) {patient_sub$SARA_category[i] <- "Includes <20% of those with AD" }
    if (patient_sub$SARA_prediction[i] >=0.45 & patient_sub$SARA_prediction[i] <=0.61) {patient_sub$SARA_category[i] <- "Patient Results Are Equivocal" }
    if (patient_sub$SARA_prediction[i] <= 0.82 & patient_sub$SARA_prediction[i] > 0.61) {patient_sub$SARA_category[i] <- "Includes <20% of those with Non-AD" }
    if (patient_sub$SARA_prediction[i] > 0.82) {patient_sub$SARA_category[i] <- "Includes <10% of those with Non-AD" }
  }
  
  
  
  
  ########################################
  #  GRAPH SARA Classification RESULTS
  ########################################
  if (make_graphs == T) { #making graphs is optional, intended for use by clinician
    
    
    for (i in 1:length(patient_sub$CLINICALDATA.ID)) {
      #make graph
      boxplot_data <- c(0,0.29, 0.45, 0.61, 0.82, 1)
      
      Fig <- ggplot(data.frame(boxplot_data), aes(x=boxplot_data, y=0) ) +
        #description and results above the color bar
        annotate("text", x = 0.5, y = 1.85, size=5, label ="The weighted combination of hippocampal volume, amygdala volume, the inferior portion of the lateral ventricles,\nentorhinal cortex thickness, and inferior parietal cortex thickness determines to what extent a patient's atrophy is 'AD-like'") +
        
        #color for each category section of colorbar
        geom_rect(  aes(xmin = 0, xmax = 0.29),  ymin = 0.5, ymax = 1, alpha = 0.2, fill="#018571", size=1) + 
        geom_rect(  aes(xmin = 0.29, xmax = 0.45),  ymin = 0.5, ymax = 1, alpha = 0.2, fill="#80cdc1", size=1) + 
        geom_rect(  aes(xmin = 0.45, xmax = 0.61),  ymin = 0.5, ymax = 1, alpha = 0.2, fill="grey", size=1) + 
        geom_rect(  aes(xmin = 0.61, xmax = 0.82),  ymin = 0.5, ymax = 1, alpha = 0.2, fill="#dfc27d", size=1) + 
        geom_rect(  aes(xmin = 0.82, xmax = 1),  ymin = 0.5, ymax = 1, alpha = 0.1, fill="#a6611a", size=1) + 
        
        #patient's line and category
        annotate("text", x = 0.5, y = 1.45, size=7, label=patient_sub$SARA_category[i], col="#1a0f2f" ) +
        annotate("segment",x=patient_sub$SARA_prediction[i],xend=patient_sub$SARA_prediction[i], y=0.4,yend=1.1, size=2, col="#673dbd")  + #patient line
        annotate("text", size=8, x = patient_sub$SARA_prediction[i], y = 1.2, col="#673dbd", label = ((round(patient_sub$SARA_prediction[i],2)))) + 
        
        #category labeling of color bar
        annotate("text", x = 0.91, y = 0.75, label ="Includes <10%\nof those\nwith Non-AD", size=6  ) + 
        annotate("text", x = 0.145, y = 0.75, label ="Includes <10%\nof those\nwith AD", size=6 ) + 
        annotate("text", x = 0.37, y = 0.75, label ="Includes <20%\nof those\nwith AD", size=6  ) + 
        annotate("text", x = 0.53, y = 0.75, label ="Equivocal", size=6  ) + 
        annotate("text", x = 0.715, y = 0.75, label ="Includes <20%\nof those\nwith Non-AD", size=6  ) + 
        
        #lines that extend below color bar
        annotate("segment",x=0,xend=0, y=0.1,yend=1, size=1, col="black")  +
        annotate("segment",x=0.29,xend=0.29, y=0.1,yend=1, size=1, col="black")  +
        annotate("segment",x=0.45,xend=0.45, y=0.1,yend=1, size=1, col="black")  +
        annotate("segment",x=0.61,xend=0.61, y=0.1,yend=1, size=1, col="black")  +
        annotate("segment",x=0.82,xend=0.82, y=0.1,yend=1, size=1, col="black")  +
        annotate("segment",x=1,xend=1, y=0.1,yend=1, size=1, col="black")  +
        
        #labelling of the lines that extend below color bar
        geom_text(aes(label = boxplot_data), col="black", size=6)   +
        
        #graph details
        scale_x_continuous(limits = c(0,1)) +
        scale_y_continuous(limits = c(0,2)) +
        theme(panel.background = element_blank(),
              axis.text = element_blank(), 
              axis.title.x=element_text(size=18),
              axis.ticks = element_blank(),
              plot.title=element_text(size=22,face="bold",hjust = 0.5),
              legend.title.align=0.5 )  +
        labs(title="SARA Classification Model", x = "Probability of Having Alzheimer Pattern of Atrophy", y="")
      
      #make file
      file_name_full <- paste(patient$CLINICALDATA.ID[i], "_SARA_classification", ".png", sep="") #combines subject ID and region name to make file name
      png(filename= file_name_full ,width=950,height=480) #makes png file
      print(Fig) #adds Fig to png file
      dev.off() #Closes the graphs generated by print(Fig)
    }
  }
  
  
  
  ########################################
  #  WRITE CSV FILE WITH ALL DATA
  ########################################
  if (make_csv == T) {
    if (!is.na(output_directory)){
      setwd(output_directory) ###adjust to your file location where you want the output to go
    }
    #rename ID if necessary
    if (ID_name != "CLINICALDATA.ID"){
      colnames(patient_sub)[colnames(patient_sub)=="CLINICALDATA.ID"] <- ID_name
    }
    #rename Age if necessary
    if (age_label != "Age"){
      colnames(patient_sub)[colnames(patient_sub)=="Age"] <- age_label
    }
    
    
    #rename region names if necessary (z-scores too)
    if (IntraCranialVol != "IntraCranialVol"){
      colnames(patient_sub)[colnames(patient_sub)=="IntraCranialVol"] <- IntraCranialVol
    }
    if (Left.Amygdala_volume != "Left.Amygdala_volume"){
      colnames(patient_sub)[colnames(patient_sub)=="Left.Amygdala_volume"] <- Left.Amygdala_volume
      colnames(patient_sub)[colnames(patient_sub)=="Left.Amygdala_volume_zscore"] <- paste(Left.Amygdala_volume,"_zscore", sep = "")
    }
    if (Right.Amygdala_volume != "Right.Amygdala_volume"){
      colnames(patient_sub)[colnames(patient_sub)=="Right.Amygdala_volume"] <- Right.Amygdala_volume 
      colnames(patient_sub)[colnames(patient_sub)=="Right.Amygdala_volume_zscore"] <-  paste(Right.Amygdala_volume,"_zscore", sep = "")
    }
    
    if (Left.Hippocampus_volume != "Left.Hippocampus_volume"){
      colnames(patient_sub)[colnames(patient_sub)=="Left.Hippocampus_volume"] <- Left.Hippocampus_volume
      colnames(patient_sub)[colnames(patient_sub)=="Left.Hippocampus_volume_zscore"] <-  paste(Left.Hippocampus_volume,"_zscore", sep = "")
    }
    if (Right.Hippocampus_volume != "Right.Hippocampus_volume"){
      colnames(patient_sub)[colnames(patient_sub)=="Right.Hippocampus_volume"] <- Right.Hippocampus_volume
      colnames(patient_sub)[colnames(patient_sub)=="Right.Hippocampus_volume_zscore"] <-  paste(Right.Hippocampus_volume,"_zscore", sep = "")
    }
    
    if (Left.Inf.Lat.Vent_volume != "Left.Inf.Lat.Vent_volume"){
      colnames(patient_sub)[colnames(patient_sub)=="Left.Inf.Lat.Vent_volume"] <- Left.Inf.Lat.Vent_volume
      colnames(patient_sub)[colnames(patient_sub)=="Left.Inf.Lat.Vent_volume_zscore"] <-  paste(Left.Inf.Lat.Vent_volume,"_zscore", sep = "")
    }
    if (Right.Inf.Lat.Vent_volume != "Right.Inf.Lat.Vent_volume"){
      colnames(patient_sub)[colnames(patient_sub)=="Right.Inf.Lat.Vent_volume"] <- Right.Inf.Lat.Vent_volume
      colnames(patient_sub)[colnames(patient_sub)=="Right.Inf.Lat.Vent_volume_zscore"] <-  paste(Right.Inf.Lat.Vent_volume,"_zscore", sep = "")
    }
    
    if (lh_entorhinal_thickness != "lh_entorhinal_thickness"){
      colnames(patient_sub)[colnames(patient_sub)=="lh_entorhinal_thickness"] <- lh_entorhinal_thickness
      colnames(patient_sub)[colnames(patient_sub)=="lh_entorhinal_thickness_zscore"] <-  paste(lh_entorhinal_thickness,"_zscore", sep = "")
    }
    if (rh_entorhinal_thickness != "rh_entorhinal_thickness"){
      colnames(patient_sub)[colnames(patient_sub)=="rh_entorhinal_thickness"] <- rh_entorhinal_thickness
      colnames(patient_sub)[colnames(patient_sub)=="rh_entorhinal_thickness_zscore"] <-  paste(rh_entorhinal_thickness,"_zscore", sep = "")
    }
    
    if (lh_inferiorparietal_thickness != "lh_inferiorparietal_thickness"){
      colnames(patient_sub)[colnames(patient_sub)=="lh_inferiorparietal_thickness"] <- lh_inferiorparietal_thickness
      colnames(patient_sub)[colnames(patient_sub)=="lh_inferiorparietal_thickness_zscore"] <-  paste(lh_inferiorparietal_thickness,"_zscore", sep = "")
    }
    if (rh_inferiorparietal_thickness != "rh_inferiorparietal_thickness"){
      colnames(patient_sub)[colnames(patient_sub)=="rh_inferiorparietal_thickness"] <- rh_inferiorparietal_thickness
      colnames(patient_sub)[colnames(patient_sub)=="rh_inferiorparietal_thickness_zscore"] <-  paste(rh_inferiorparietal_thickness,"_zscore", sep = "")
    }
    filename=paste(patient_data_source, "_ICVcorrected_regions_with_SARA_", Sys.Date(), ".csv", sep="") #making csv to write data into
    write.table(file=filename,as.matrix(patient_sub),col.names=TRUE, row.names=FALSE,sep=",",dec=".") #writing data into csv file
    
  }
  
} #end of function
