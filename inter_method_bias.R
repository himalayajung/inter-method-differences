# calculate and plot method biases for male and ASD subjects

library(reshape2)
library(ggplot2)
library(dplyr)
library(effsize)
library(gridExtra)
# -----
# calculate biases (in %) from the beta coefficients of the regression models
get_percent_diff = function (tissue, method, term, estimate, wrto="abs_mean"){
  tissue = as.character(tissue[1])
  estimate[term == "age"] = estimate[term == "age"]*sd(m$age)
  estimate[term == "I(age^2)"] = estimate[term == "I(age^2)"]*(sd(m$age))^2

  if (wrto=="abs_mean"){
  avg = mean(m_method[m_method["tissue"]==tissue, method], na.rm=TRUE)} 
    else if(wrto=="intercept"){
    avg=abs(estimate[term == "(Intercept)"])
  }
  percent_diff = estimate/avg*100
  
  if (method == "FS" & tissue == "CSF"){
    percent_diff = 0
  }
  print(percent_diff)
  return(percent_diff)
}

# data
m = read.csv("../../data/brain_volumes_long1.csv")

# changing factor levels to control the order in regression models and also in plots
m$control = factor(m$control, levels = c("TDC", "ASD"))
m$method = factor(m$method, levels = c("SPM", "FSL", "FS"))
m$tissue = factor(m$tissue, levels = c("TIV", "GM", "WM", "CSF"))
m_tissue = dcast(m, ...~tissue)
m_method = dcast(m, ...~method)

# SPM vs. FSL
spm_fsl = m_method %>% group_by(tissue) %>% 
  filter(!is.na(FIQ)) %>% # remove FIQ NAs
  mutate(value = SPM-FSL) %>%
  do(tidy(my.lm(value~control+sex+site+age+I(age^2), data = .))) %>%
  mutate(methods = "SPM-FSL", percent_diff=get_percent_diff(tissue, "FSL", term, estimate),
         percent_diff1 = get_percent_diff(tissue, "FSL", term, estimate, "intercept") )

# SPM vs. FS
spm_fs = m_method %>% group_by(tissue) %>% 
  filter(!is.na(FIQ)) %>% # remove FIQ NAs
  mutate(value = SPM-FS) %>%
  do(tidy(my.lm(value~control+sex+site+age+I(age^2), data = .))) %>%
  mutate(methods = "SPM-FS", percent_diff=get_percent_diff(tissue, "FS", term, estimate),
         percent_diff1 = get_percent_diff(tissue, "FS", term, estimate, "intercept") )

# FSL vs. FS
fsl_fs = m_method %>% group_by(tissue) %>% 
  filter(!is.na(FIQ)) %>% # remove FIQ NAs
  mutate(value = FSL-FS) %>%
  do(tidy(my.lm(value~control+sex+site+age+I(age^2), data = .))) %>%
  mutate(methods = "FSL-FS", percent_diff=get_percent_diff(tissue, "FS", term, estimate),
percent_diff1 = get_percent_diff(tissue, "FS", term, estimate, "intercept") )

im_bias = rbind(spm_fsl, spm_fs, fsl_fs)
im_bias = subset(im_bias, term %in% c("controlASD", "sexM"))

# multiple comparisons correction using FDR
im_bias$p_corrected = p.adjust(im_bias$p.value, method = "BH")

# some editing for plotting
im_bias[im_bias$term=="controlASD", "term"] = "ASD"
im_bias[im_bias$term=="sexM", "term"] = "Male"

write.csv(im_bias, "inter_method_bias.csv")

p = ggplot(subset(im_bias, p_corrected<0.05), aes(x=term, y=percent_diff1))+
  geom_bar(aes(fill=term), stat="identity")+
  geom_abline(intercept=0, slope=0)+
  geom_text(aes(label=term), angle=90, hjust=1)+
  ylab("% Differential Bias")+
  ggtitle("% Differential Bias for group and sex \n with resepect to \n the mean inter-method difference")+
  facet_grid(tissue~methods, scales = "free")+
  theme(legend.position="top")

q = ggplot(subset(im_bias, p_corrected<0.05), aes(x=term, y=percent_diff))+
  geom_bar(aes(fill=term), stat="identity")+
  geom_abline(intercept=0, slope=0)+
  geom_text(aes(label=term), angle=90, hjust=1)+
  ylab("% Differential Bias")+
  ggtitle("% Differential Bias for group and sex \n with resepect to \n the mean of the baseline method")+
  facet_grid(tissue~methods, scales = "free")+
  theme(legend.position="top")

pdf("inter_method_bias.pdf")
grid.arrange(p,q,nrow=1)
dev.off()
  

