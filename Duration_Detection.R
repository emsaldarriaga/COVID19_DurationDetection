#Impetus: To estimate the mean duration of virus detection since syptom onset from Upper Respiraty Tract (URT), Lower Respiratory Tract (LRT), and/or Stool Samples.

#Sources:
##Data
# Epi Data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7323671/
#Epi data others: https://www.nature.com/articles/s41591-020-0869-5#Sec2; repo: https://github.com/ehylau/COVID-19

#Methods:
#1. Meta analysis in R: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/
#2. Convert mean to from median: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-135; first saw: https://www.researchgate.net/post/which_is_best_way_to_convert_median_to_mean_in_meta-analysis
#3. Medians vs Means: Go for Means. And consider: https://pubmed.ncbi.nlm.nih.gov/31553488/ for discussion
#Meta anlysis for simgle means in R: https://rdrr.io/cran/meta/man/metamean.html#:~:text=Fixed%20effect%20and%20random%20effects,(%20sm%20%3D%20%22MLN%22%20)
#4. Heterogeinity: https://www.statsdirect.com/help/meta_analysis/heterogeneity.htm; https://www.meta-analysis-workshops.com/download/common-mistakes1.pdf (more complete)


#Requiring####	
lbs = c("readr", "rvest","stringr","ggplot2","meta","dplyr","grid","gridExtra","cowplot")
sapply(c(lbs),require, character.only=T)

#Functions####
#Based on Wan et al, BMC 2014. Scenario 3. See methods 2.
convert_median = function (iq1, iq3, median, min, max, SampleSize,Method){ 
  if (Method== 1){ #1: MinMAx; 0:IQR
    mean = (min + 2*median +max)/4
    probs = (SampleSize - 0.375)/(SampleSize+0.25)
    sd = (max - min) / (2*qnorm(probs,0,1))
    vals = list(Mean = mean, SD = sd)
    return(vals)
  }else{
    mean = (iq1 +median +iq3)/3
    probs = (0.75*SampleSize - 0.125)/(SampleSize+0.25)
    sd = (iq3 - iq1) / (2*qnorm(probs,0,1))
    vals = list(Mean = mean, SD = sd)
    return(vals)  
  }
}
convert_median(min=0.02,max=0.04,median=0.03,SampleSize=100,Method=1) #validation. See excel downloaded from Methods 2. Research


convert_median_mean = function (iq1, iq3, median, min, max, SampleSize,Method){ 
  if (Method== 1){ #1: MinMAx; 0:IQR
    mean = (min + 2*median +max)/4
    return(mean)
  }else{
    mean = (iq1 +median +iq3)/3
    return(mean)  
  }
}

convert_median_sd = function (iq1, iq3, median, min, max, SampleSize,Method){ 
  if (Method== 1){ #1: MinMAx; 0:IQR
    probs = (as.numeric(SampleSize) - 0.375)/(as.numeric(SampleSize)+0.25)
    sd = (max - min) / (2*qnorm(probs,0,1))
    return(sd)
  }else{
    probs = (0.75*as.numeric(SampleSize) - 0.125)/(as.numeric(SampleSize)+0.25)
    sd = (iq3 - iq1) / (2*qnorm(probs,0,1))
    return(sd)  
  }
}

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} #inferior to lines 99:101 because gets only first occurance

numextract_multi <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
} #multiple occurances

pattern  = "[0-9]+| [0-9]+ |[0-9]+[.][0-9]+| [0-9]{1,2}[.][0-9]{0,2} |\\s[0-9]{1,2}[.][0-9]{0,2}"

extract_multi = function(string, pattern) { #finds numbers and extract them into data frame by order of finding. Quality assessment after apply it strongly recommended
  defaultW <- getOption("warn")
  options(warn = -1)
  string = data.frame(string)
  first = NULL
  second = NULL
  third = NULL
  fourth = NULL
  fifth = NULL
  for (i in 1:nrow(string)){
    ext  = unlist(str_extract_all(string[i,], pattern))
    first = rbind(first,as.numeric(ext[1]))
    second = rbind(second,as.numeric(ext[2]))
    third = rbind(third,as.numeric(ext[3]))
    fourth = rbind(fourth,as.numeric(ext[4]))
    fifth = rbind(fifth,as.numeric(ext[5]))
  }
  findings = data.frame(first, second, third, fourth, fifth)
  return(findings)
  options(warn = defaultW)
}

#Load data via web scrapping####
url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7323671/table/tbl0002/?report=objectonly"
epidata = url %>%
  read_html() %>%
  html_node(xpath = '//*[@id="tbl0002"]/div[2]/table') %>%
  html_table(fill = TRUE,header=TRUE)

urtdata = epidata[epidata$`Aggregate study-level duration of virus detection since symptom onset from URT samples*`!=epidata[2,5],-c(8,7,6)] 
#Selecting upper respiratory tract samples only
colnames(urtdata)[5] = "Values"
rownames(urtdata) = 1:nrow(urtdata)

#Clean data#### 
# Unit: seems all are days. Confirm
# Point estimates can be: Fixed, Mean, Median
# Unceratinty type: Range, range, IQR, 95% CI, SD

#Keep in mind that this are estimates of DURATION, meaning, the total period that the virus remains detectable. Meaning, %detectable=0 | duration > estimate+1day
#Once cleaned, evaluate conduct a meta analysis to fin Median, IQR/SD, and Max, and fit a Weibull or log-normal to the data

urtdata$unitid =  ifelse(str_detect(urtdata$Values,"day|days"),1,0)
urtdata[urtdata$unitid==0,]
urtdata = urtdata[-56,]
urtdata$unit =  ifelse(urtdata$unitid==1,"days","week")

urtdata$pestimate =  ifelse(str_detect(urtdata$Values,"Median|median"),"Median",
                            ifelse(str_detect(urtdata$Values,"Mean|mean"),"Mean","Fixed"))
View(urtdata[urtdata$pestimate=="Fixed",5])

urtdata$uncert =  ifelse(str_detect(urtdata$Values,"SD|sd|SD|standard deviation"),"SD",
                ifelse(str_detect(urtdata$Values,"IQR|Iqr|iqr"),"IQR",
                       ifelse(str_detect(urtdata$Values,"95% CI|95%|CI"),"95% CI",
                              ifelse(str_detect(urtdata$Values,"Range|range"),"Range","NA"))))
#By putting Range at the end we prioritize the other estimates if there's more than one, e.g. urtdata[2,5]
View(urtdata[urtdata$uncert=="NA",5])

#Both point estimates and uncertainty estimate match regarding their units. Now need to extract data

#Extract numbers: Start on population
numb <- "[0-9]+" #define a pattern
urtdata$popval = as.numeric(str_extract(urtdata$Population, numb))
summary(as.numeric(urtdata$popval))
ggplot(urtdata[as.numeric(urtdata$popval)<500,], 
       aes(x=as.numeric(popval))) + geom_histogram(bins=100)
length(urtdata[as.numeric(urtdata$popval)==1,1]) #21 out of 84 studies have only 1 observation

#Split data into dataframes by pestimate
#urtbypoint = urtdata %>% group_by(as.factor(pestimate))
urtby = split(urtdata,as.factor(urtdata$pestimate))
fixed = urtby$Fixed; median = urtby$Median; mean = urtby$Mean

#get the Mean
mean$Values = gsub(' '," ",mean$Values) #clean values

mean$pestimate_val = ifelse(is.na(str_extract(mean$Values, "[0-9]+[.][0-9]+")),
                            str_extract(mean$Values, "[0-9]+"),
                            str_extract(mean$Values, "[0-9]+[.][0-9]+")) 
#Done only for first match. Remember [0-9]+ equivalent [0-9]{m,n} m:min, n:max number of characters

mean$pestimate_val=NULL; mean$pestimate_val2=NULL; mean$pestimate_val3=NULL

pattern  = "[0-9]+| [0-9]+ |[0-9]+[.][0-9]+| [0-9]{1,2}[.][0-9]{0,2} |\\s[0-9]{1,2}[.][0-9]{0,2}" #\s replace any space character. need double \ to scape the first \

finds  = extract_multi(mean$Values,pattern)

mean$Mean = finds$first; mean$SD = finds$second
#quality control
rownames(mean) = 1:nrow(mean)
mean[1,]$SD = 6.6
mean[4,]$SD = 4.6
mean[7,]$SD = 6.7
mean[8,]$Mean = 9.71
mean[8,]$SD = (mean[8,]$popval)^0.5*(mean[8,]$Mean-8.21)/1.95
mean[9,]$SD = 10.7

#get the Median
median$Values = gsub(' '," ",median$Values) #clean values
finds  = extract_multi(median$Values,pattern)
finds = cbind(median[,c(5,9:10)],finds)
rownames(finds)= 1:nrow(finds)

#control
finds[1,c(5,6)]=(finds[1,c(6,7)])
finds[5,c(4,5,6)]=c(11,10,12)
finds[7,c(5,6)]=c(4.5,9.5)
finds[8,c(6)]=c(26.25)
finds[9,c(5,6)]=c(13.25,22)
finds[11,c(5,6)]=c(47.75,60.5)
finds[12,c(6)]=11
finds[15,c(5,6)]=c(3.5,13)
finds[21,c(5,6)]=c(11.5,16)


finds$Mean = ifelse(finds$uncert=="IQR",convert_median_mean(iq1=finds$second,
                                                            iq3 = finds$third,
                                                            median = finds$first,
                                                            SampleSize = finds$popval,
                                                            Method = 0),
                    ifelse(finds$uncert=="Range",convert_median_mean(min=finds$second,
                                                                      max = finds$third,
                                                                      median = finds$first,
                                                                      SampleSize = finds$popval,
                                                                      Method = 1), NA))

finds$SD = ifelse(finds$uncert=="IQR",convert_median_sd(iq1=finds$second,
                                                            iq3 = finds$third,
                                                            median = finds$first,
                                                            SampleSize = finds$popval,
                                                            Method = 0),
                    ifelse(finds$uncert=="Range",convert_median_sd(min=finds$second,
                                                                     max = finds$third,
                                                                     median = finds$first,
                                                                     SampleSize = finds$popval,
                                                                     Method = 1), NA))

finds[5,]$Mean = convert_median_mean(iq1=finds[5,]$second,
                                     iq3 = finds[5,]$third,
                                     median = finds[5,]$first,
                                     SampleSize = finds[5,]$popval,
                                     Method = 0)

finds[5,]$SD = convert_median_sd(iq1=finds[5,]$second,
                                     iq3 = finds[5,]$third,
                                     median = finds[5,]$first,
                                     SampleSize = finds[5,]$popval,
                                     Method = 0)

median$Mean = finds$Mean; median$SD = finds$SD

#Merge all bases:
#Fixed wasn't worked out because it doens't include a measure of uncertainty and shows poor assessment
fixed$Mean = NA; fixed$SD = NA

epidata_wo = rbind(mean, median, fixed)
epidata_wo$Subpop = with(epidata_wo,
                   ifelse(str_detect(Population,"and")==TRUE,"Children and Adults",
                          ifelse(str_detect(Population,"adult|adult|Adults|Adults")==TRUE,"Adults",
                                 ifelse(str_detect(Population,"children|Children")==TRUE,
                                        "Children","Unknown"))))
write.csv(epidata_wo, "epidata_upt.csv")

#Plot main effects####
maineffects = metagen(Mean,
                      SD,
                      data=epidata_wo[!is.na(epidata_wo$Mean),],
                      studlab=paste(`First Author`),
                      comb.fixed = T,
                      comb.random = T,
                      method.tau = "SJ", #different methods to estimate can tau lead to important differences in the final value and hence the prediction interval. Compare the default "DL" to "ML" and "SJ"
                      hakn = TRUE,
                      prediction=TRUE,
                      sm="MD") #SMD: standirez mean difference
forest(maineffects) 

png("DurationDetectionForest.png", width = 1490, height = 1000, res = 108)
forest(maineffects) 
dev.off()

maineffects_singlemean = metamean(n=as.numeric(popval),
                                  Mean,
                                  SD,
                                  studlab=paste(`First Author`),
                                  comb.fixed = T,
                                  comb.random = T,
                                  method.tau = "SJ",
                                  hakn = TRUE,
                                  prediction=TRUE,
                                  method.tau.ci = "QP",
                                  data = epidata_wo[!is.na(epidata_wo$Mean),])

forest(maineffects_singlemean)
png("DurationDetectionForestSM.png", width = 1490, height = 1000, res = 108)
forest(maineffects_singlemean) 
dev.off()

##use metagen to get the same results as metamean
maineffects_c = metagen(Mean,
                        SD/(as.numeric(popval))^0.5,
                        data=epidata_wo[!is.na(epidata_wo$Mean),],
                        studlab=paste(`First Author`),
                        comb.fixed = T,
                        comb.random = T,
                        method.tau = "SJ",
                        hakn = TRUE,
                        prediction=TRUE,
                        sm="MD")
forest(maineffects_c) 

png("DurationDetectionForest_corrected.png", width = 1490, height = 1000, res = 108)
forest(maineffects_c) 
dev.off()

# Analysis by sub-groups: alternatives Population (children, adult), Study Design (case series, cohort); Type of measure (mean, median)####
epidata = epidata_wo[!is.na(epidata_wo$Mean),]

epidata$Subpop = with(epidata,
                        ifelse(str_detect(Population,"and")==TRUE,"Children and Adults",
                        ifelse(str_detect(Population,"adult|adult|Adults|Adults")==TRUE,"Adults",
                        ifelse(str_detect(Population,"children|Children")==TRUE,
                                       "Children","Unknown"))))
table(epidata[,c(13)]) #17 adults, 4 children, 11 both

metamean_r = function(df){
  r = metamean(n=as.numeric(popval),
         Mean,
         SD,
         studlab=paste(`First Author`),
         comb.fixed = T,
         comb.random = T,
         method.tau = "SJ",
         hakn = TRUE,
         prediction=TRUE,
         data = df)
  return(r)
}

maineffects_subpop = update(maineffects_singlemean,byvar = epidata$Subpop)
forest(maineffects_subpop)
forest(maineffects_subpop, layout = "subgroup", test.subgroup = T)

maineffects_design = update(maineffects_singlemean,byvar = epidata$`Study design`)
forest(maineffects_design, layout = "subgroup", test.subgroup = T)

maineffects_estimate = update(maineffects_singlemean,byvar = epidata$pestimate)
forest(maineffects_estimate, layout = "subgroup", test.subgroup = T)

maineffects_country = update(maineffects_singlemean,byvar = epidata$Country)
forest(maineffects_country, layout = "subgroup", test.subgroup = T)

png("../CORONAVIRUS/Projects/DurationDetection/SM_by_subpopulation_all.png", width = 1490*1.3, height = 1000*1.3, res = 108)
forest(maineffects_subpop) 
dev.off()

plots = list("subpop", "design", "estimate","country")

for (i in plots) {
  namei = paste0("maineffects_",i)
  path = paste0("../CORONAVIRUS/Projects/DurationDetection/SM_by_",i,".png")
  png(path, width = 1490, height = 1000, res = 108)
  forest(get(namei), layout = "subgroup", test.subgroup = T)
  dev.off()
}  
