library(Hmisc)

get_data<- function(){
    # Read in google sheet
    gs4_deauth()
    df<- read_sheet("https://docs.google.com/spreadsheets/d/1raDMYU-MpSikywNyi1FEHR0wJPOAOaXbtwwMaVPOSQ/edit#gid=0")
    
    # Remove empty rows
    df <- df %>% drop_na(author)
    df <- df %>% drop_na(n)
  
    # Get correct formats
    df <- df %>% mutate_at(vars(age, age_onset, duration_illness_years, male_percent, duration_treatment,white_matter_sd, thalamus_mn), ~as.numeric(as.character(.))) 
    df$notes <- as.character(df$notes)
    df$link_to_paper <- as.character(df$link_to_paper)
    df$year <- as.character(df$year)
    
    # Create study id column
    df$studyid <- paste(df$author, df$year, sep = "")
    
    #remove ineligible studies
    df <- df %>% filter(is.na(dont_include))
    
    #remove arms with less than 10 participants
    df <- df %>% filter(n>=10)
    
    #Drop specific studies
    df <- df %>% filter(studyid!='Buchanan1993') # Repackaging of Brier et al. 1992
    df <- df %>% filter(studyid!='Cannon2002') # Same cohort as Cannon 1998
    df <- df %>% filter(studyid!='Chakos1994')#Results reported in Lieberman 2001 and Bogerts 1990
    df <- df %>% filter(studyid!='Collin2012')#From two previously published cohorts (van Haren et al., 2008, Cahn et al., 2009) - originals included
    df <- df %>% filter(studyid!='Poletti2016')# effect sizes are erroneously large >3
    
    # Combine left and right
    regions <- str_sub(unique(colnames(select(df,contains("mn_l"))[0,])), end=-6)
    r = 0.7

    for (region in regions){
      for (study in unique(df$studyid)){
        single_study_df <- df %>% filter(studyid==study)
        if ((is.na(single_study_df[sprintf('%s_mn', region)][[1]]) & !is.na(single_study_df[sprintf('%s_mn_l', region)]))[[1]]){
          df <- calc_bilat(study, region, r, df)
        }
      }
    }
  
    #Rename Diagnoses
    df <- df %>% mutate(diagnosis = recode(diagnosis, control = "Control",  autism='ASD', trauma_control='Trauma Control', schizophrenia='Scz', schizoaffective='SczAff', mixed_anx_dep='Mixed Anx. + Dep.', depression='MDD', borderlinePD='BPD', psychotic_depres='Psychotic Dep.'))
    df$diagnosis_grouped <- ifelse(df$diagnosis %in% c('bipolar1_psychotic'), 'BPAD.P', df$diagnosis)
    df$diagnosis_grouped <- ifelse(df$diagnosis %in% c( 'mania', 'bipolar1_psychosis_unspecified','bipolar_1_or_2','bipolar2','bipolar1_nonpsychotic'), 'BPAD.NP', df$diagnosis_grouped)
    df$diagnosis_grouped <- ifelse(df$diagnosis %in% c('psychosis_mixed', 'psychosis_unspec', 'psychosis', 'affective_psychosis'), 'Psychosis', df$diagnosis_grouped)
    df$diagnosis_grouped <- ifelse(df$diagnosis %in% c('GAD', 'social_anx', 'panic'), 'Anxiety Disorder', df$diagnosis_grouped)
    df$diagnosis_grouped <- ifelse(df$diagnosis %in% c('CHR', 'schizotypalPD', 'UHR', "UHR - converted to psychosis", "UHR - did not convert to psychosis"), 'CHR', df$diagnosis_grouped)
    df$diagnosis_grouped <- df$diagnosis_grouped %>% plyr::revalue(c("Control" = "CON", "Trauma Control"="TCON", 
                                                                     "Anxiety Disorder"="ANX", "Scz"="SCZ",
                                                                     "Psychotic Dep."="PDEP", "Mixed Anx. + Dep."="A&D", "Psychosis"="PSY", 
                                                                     "SczAff"="SZAF"))
    #unlist pmid columnd
    df$pmid[sapply(df$pmid, is.null)] <- NA
    df$pmid <- unlist(df$pmid)
    df <- df %>% filter(diagnosis!='22q') 
    write_csv(df, '../data/cleaned_data.csv')
    return(df)
}

fix_data <- function(df, grouping='diagnosis_grouped', enigma){
  #Fixes studies with any doubling of diagnoses 
  #Identify studies with more than a single row for each diagnosis, then visually check
  to_fix=c()

  for (study in unique(df$studyid)){
    single_study_df <- df %>% filter(studyid==study)
    num_diagnoses = length(unique(single_study_df[[grouping]]))
    num_rows = length(single_study_df[[grouping]])
    if (num_rows != num_diagnoses){
      to_fix <- c(study, to_fix)
    }
  }
  
  #need to run 5 times for cases where there is more than one diagnosis with multiple rows
  for (i in (c(1:5))){
    for (study in to_fix){
      single_study_df <- df %>% filter(studyid==study)
      
      # Get rows where diagnosis doubled
      diagnosis_dbl <- names(which.max(table(single_study_df[[grouping]])))
      old_rows <- single_study_df %>% filter(eval(str2expression(grouping))==diagnosis_dbl)
      if(nrow(old_rows)>1){
        new_row <- old_rows[0:1,]
        new_row$n <- sum(old_rows$n)
        
        # Combine means
        weighted_avgs <- c('age', 'male_percent', 'black_ethnicity', 'white_ethnicity',	'other_ethnicity',	'age_onset',	'duration_illness_years', colnames(select(df,contains("mn"))[0,]))
        for (var in weighted_avgs){new_row <- combine_mean(var, new_row, old_rows)}                   
        
        # Combine sds
        weighted_avgs <- colnames(select(df,contains("mn"))[0,])
        for (var in weighted_avgs){new_row <- combine_sd(var, new_row, old_rows)}
        
        # Delete both rows in main df
        df <- df[!(df$studyid==study & df[[grouping]]==diagnosis_dbl),]
        
        # Put in combined row
        df <- rbind(df, new_row)}
    }
  }
  return(df)
}

get_net_data <- function(region,df, grouping='diagnosis_grouped'){
  net_data <- (matrix(ncol = 7, nrow = 0))
  colnames(net_data) <- c("TE", "seTE", "diagnosis1", "diagnosis2", "studyid", 'sd1', 'sd2')
  net_data <- data_frame(TE=double(), seTE=double(), diagnosis1=character(), diagnosis2=character(), studyid=character(), sd1=double(), sd2=double())
  
  region_df <- df %>% filter(!!as.symbol(sprintf('%s_mn', region))>0)
  for(study in unique(region_df$studyid)){
    #select individual study
    filtered_data<- region_df %>% filter(studyid==study)
    filtered_data <- filtered_data %>% drop_na(selected_grouping)
    i=1
    for(diagnosis1 in filtered_data[[grouping]]){
      trimmed_list <-  filtered_data[[grouping]][i:length(filtered_data[[grouping]])]
      i=i+1
      for(diagnosis2 in trimmed_list){
        if(diagnosis1!=diagnosis2){
          diagnosis1_n <- filtered_data %>% filter(eval(str2expression(grouping))==diagnosis1) %>% .$n
          diagnosis1_mean <- filtered_data %>% filter(eval(str2expression(grouping))==diagnosis1) %>% .[sprintf('%s_mn', region)]
          diagnosis1_sd <- filtered_data %>% filter(eval(str2expression(grouping))==diagnosis1) %>% .[sprintf('%s_sd', region)]
          diagnosis2_n <- filtered_data %>% filter(eval(str2expression(grouping))==diagnosis2) %>% .$n
          diagnosis2_mean <- filtered_data %>% filter(eval(str2expression(grouping))==diagnosis2) %>% .[sprintf('%s_mn', region)]
          diagnosis2_sd <- filtered_data %>% filter(eval(str2expression(grouping))==diagnosis2) %>% .[sprintf('%s_sd', region)]
          meta_result <- metacont(diagnosis1_n, diagnosis1_mean[[1]], diagnosis1_sd[[1]], diagnosis2_n, diagnosis2_mean[[1]], diagnosis2_sd[[1]], sm="SMD")
          # just to check if SEs have been used by mistake:
          sd2 <- diagnosis2_mean/diagnosis2_sd
          sd1 <- diagnosis1_mean/diagnosis1_sd
          net_data <- add_row(net_data, diagnosis1=diagnosis1, diagnosis2=diagnosis2, studyid=study, TE=meta_result$TE, seTE=meta_result$seTE, sd1=sd1[[1]], sd2=sd2[[1]])
        }
      }
    }
  }
  net_data <- net_data %>% mutate_at(vars(diagnosis1, diagnosis2, studyid), ~as.character(.)) 
  return(net_data)
}

combine_sd <- function(var, new_row, old_rows){
  # var should be the string for the mean value
  mean_var = var
  sd_var = str_replace(var, 'mn', 'sd')

  # Combine SDs
  a <- sum((old_rows$n-1)*(old_rows[sd_var])^2)
  b <- old_rows$n[1]*old_rows$n[2]
  c <- old_rows[mean_var][1,]^2+old_rows[mean_var][2,]^2-2*old_rows[mean_var][1,]*old_rows[mean_var][2,]
  d <- sum(old_rows$n)
  e <- (sum(old_rows$n)-1)
  f <- ((a+b*c/d)/e)^0.5
  new_row[sd_var] <- f[1,1]
  return(new_row)
}



remove_ageunmatched <- function(df){
  unmatched_groups = c()
  df$groupid <- paste0(df$studyid, df$diagnosis_grouped)
  for (study in unique(df$studyid)){
    single_study_df <- df %>% filter(studyid==study)
    age_diff <- diff(range(single_study_df$age))
    if(is.na(age_diff)){next}
    if (age_diff > 10){
      unmatched_group <- single_study_df$groupid[which.max(abs((single_study_df$age - mean(single_study_df$age))))]
      unmatched_groups <- c(unmatched_groups, unmatched_group)
      #print(paste0(study, age_diff))
    }
    if(nrow(single_study_df)==1){unmatched_groups <- c(unmatched_groups, single_study_df$groupid)}
  }
  
  df <- df %>% filter(groupid %nin% unmatched_groups)
  return(df)
}


remove_sexunmatched <- function(df){
  unmatched_groups = c()
  df$groupid <- paste0(df$studyid, df$diagnosis_grouped)
  for (study in unique(df$studyid)){
    single_study_df <- df %>% filter(studyid==study)
    sex_diff <- diff(range(single_study_df$male_percent))
    if(is.na(sex_diff)){next}
    if (sex_diff > .1){
      unmatched_group <- single_study_df$groupid[which.max(abs((single_study_df$male_percent - mean(single_study_df$male_percent))))]
      unmatched_groups <- c(unmatched_groups, unmatched_group)
      #print(paste0(study, sex_diff))
    }
    if(nrow(single_study_df)==1){unmatched_groups <- c(unmatched_groups, single_study_df$groupid)}
  }
  
  df <- df %>% filter(groupid %nin% unmatched_groups)
  return(df)
}




combine_mean <- function(var, new_row, old_rows){
  # var sould be the string for the mean value
  mean_var = var
  sd_var = str_replace(var, 'mn', 'sd')
  
  # Combine means
  new_row[mean_var] <- sum(old_rows$n*old_rows[mean_var])/new_row$n
  
  return(new_row)
}

calc_bilat <- function(study, region, r, df){
  #combines left and right for a given roi to fill the combined spaces
  single_study_df <- df %>% filter(studyid==study)
  mn_l = as.numeric(single_study_df[sprintf('%s_mn_l', region)][[1]])
  mn_r = as.numeric(single_study_df[sprintf('%s_mn_r', region)][[1]])
  sd_l = as.numeric(single_study_df[sprintf('%s_sd_l', region)][[1]])
  sd_r = as.numeric(single_study_df[sprintf('%s_sd_r', region)][[1]])
  single_study_df[sprintf('%s_mn', region)] = mn_l + mn_r
  single_study_df[sprintf('%s_sd', region)] = (sd_l^2 + sd_r^2 + 2*sd_l*sd_r*r)^0.5
  df[df['studyid']==study,] <- single_study_df
  return(df)
}


effect_heatmap <- function(mn0, title){
  #get effect size and ci
  league_df <- netleague(mn0,seq = netrank(mn0))$random
  league_df_es <- data.frame(lapply(league_df, function(y) gsub(" .*$", "", y)))
  league_df_ci <- data.frame(lapply(league_df, function(y) gsub(".*\\[ (.+) \\].*", "\\1", y)))
  
  #signficant squares
  l_ci <- data.frame(lapply(league_df, function(y) sub(".*\\[([^][]+);.*", "\\1", y)))
  l_ci <- data.frame(lapply(l_ci, function(y) as.numeric(as.character(y))))
  u_ci <- data.frame(lapply(league_df, function(y) sub(".*\\;([^][]+)].*", "\\1", y)))
  u_ci <- data.frame(lapply(u_ci, function(y) as.numeric(as.character(y))))
  league_df_signif <- data.frame((l_ci<0 & u_ci<0)|(l_ci>0 & u_ci>0))
  league_df_signif <- data.frame(lapply(league_df_signif, function(y) (ifelse(y, '*', ''))))
  league_df_signif %>% mutate_if(is.factor, as.character) -> league_df_signif
  
  # sort out data
  diags <- diag(as.matrix(league_df_es))
  league_df_es <- data.frame(lapply(league_df_es, function(y) as.numeric(as.character(y))))
  league_df_es <- league_df_es %>% rename_with(~ diags, colnames(league_df_es))
  league_df_signif <- league_df_signif %>% rename_with(~ diags, colnames(league_df_signif))
  league_df_es$diag <- factor(diags, levels=diags)
  league_df_signif$diag <- factor(diags, levels=diags)
  
  #plot
  league_df_es_long <- league_df_es %>% pivot_longer(!diag) 
  league_df_signif_long <- league_df_signif %>% pivot_longer(!diag) 
  league_df_es_long$signif <- league_df_signif_long$value
  league_df_es_long$name <- factor(league_df_es_long$name, levels=diags)
  league_df_es_long$g <- league_df_es_long$value
  league_df_es_long %>% ggplot(aes(x=diag,y=name, fill=g, label=signif))+geom_tile()+geom_text()+
    scale_fill_gradient2(low = "tomato", mid = 'white',high = "seagreen", limits=c(-0.6, 0.6), oob=squish, breaks= c(-0.6,0,0.6)) +xlab("")+ylab("")+ggtitle(title)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+# guides(fill=guide_legend(title='g'))+
    theme(plot.margin = margin(0.01,0.01,0.01,0.01, "cm"))
}

forestplot <- function(graph_results, title){
  ggplot(data=graph_results) +
    geom_pointrange( aes(x=label, y=estimate,  ymin=lb, ymax=ub), position = position_dodge(width = 0.9)) + 
    geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    ylab("") +
    xlab("")+
    ggtitle(title)+
    theme_minimal() + # use a white background
    theme(legend.position="none")+
    theme(panel.background = element_rect(fill = NA))+
    theme(strip.text.x=element_text(angle =180, hjust = 0.5))+
    theme(axis.text.x=element_text(angle = 0, hjust = 0.5))+
    theme(axis.text.y = element_text(angle = 0, hjust=0.5))+
    theme(strip.placement = 'outside')+coord_flip()+
    theme(plot.title = element_text(size=10, hjust=0.5))
}

amend_title <- function(region){
  title <- region
  if (region=='acc'){title <- 'Anterior Cingulate'}
  if (region=='frontal'){title <- 'Frontal Lobe'}
  if (region=='parahippocampal'){title <- 'Parahippocampal Gyrus'}
  if (region=='temporal'){title <- 'Temporal Lobe'}
  title <- str_to_title(str_replace(title, "_", " "))
  return(title)
}


get_changes<- function(x){
  evaluate <- c()
  for (i in 1:(length(x)-1)) {
    if(x[i] != x[i+1]){
      evaluate <- c(evaluate,i)
    }
  }
  return(evaluate)
}


