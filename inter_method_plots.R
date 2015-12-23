# calculate and plot % inter-method differences 
# --
library(reshape2)
library(ggplot2)
library(dplyr)
library(effsize)
library(gridExtra)
library(RColorBrewer)

# calculate different statistics (Cohen's d, p-val, % ) of inter-method diffrences
inter_method_diff = function(df_method){
  # SPM vs. FSL & FS
  diff = df_method %>% group_by(tissue) %>% 
    summarise_each(funs(p = t.test.p.value(SPM, ., original=TRUE, pair = TRUE), 
                        cd = cohen.d(SPM, ., na.rm=TRUE)$estimate), vars=c(FSL, FS) )
  diff = melt(diff, id="tissue")
  diff = diff[order(as.character(diff$variable)), ]
  diff$metric = gsub(".*_", "", diff$variable)
  diff$variable = rep(c("SPM_FSL", "SPM_FS"), each=2*length(unique(diff$tissue)))
  
  # FSL vs. FS
  diff1 = df_method %>% group_by(tissue) %>% 
    summarise_each(funs(p = t.test.p.value(FSL, ., original=TRUE, pair = TRUE), 
                        cd = cohen.d(FSL, ., na.rm=TRUE)$estimate), vars=c(FS) )
  diff1 = melt(diff1, id="tissue")
  diff1$metric = gsub(".*_", "", diff1$variable)
  diff1$variable = rep(c("FSL_FS"), each=2*length(unique(diff1$tissue)))
  
  diff = rbind(diff, diff1)
  diff = dcast(diff, ...~metric)
  diff = diff[order(as.character(diff$variable)), ]
  # multiple comparisons correction FDR
  diff$p_corrected = p.adjust(diff$p, method="BH")
  return(diff)
}

# data
m = read.csv("../../data/brain_volumes_long1.csv")
# m = subset(m, site == "NYU") # if only NYU

m$control = factor(m$control, levels = c("TDC", "ASD"))
m$method = factor(m$method, levels = c("SPM", "FSL", "FS"))
m$tissue = factor(m$tissue, levels = c("TIV", "GM", "WM", "CSF"))
m_tissue = dcast(m, ...~tissue)
m_method = dcast(m, ...~method)

# boxplot for the distribution of brain volumes
m1 = m
m1[m1$method=="FS" & m1$tissue == "CSF", "value"] = 0
box_plot = ggplot(data = m1, aes(x=method, y=value)) +
  geom_boxplot(aes(fill=factor(method))) + 
  xlab("")+
  ylab("Volume (Liters)")+
  ggtitle("A) Distribution of the brain volumes estimated by SPM, FSL and FS")+
  scale_fill_manual(values=c("#CC79A7","#56B4E9","#E69F00"), name="Methods")+
  facet_wrap(~tissue, scales="free", nrow=1)+
  theme_bw(base_size=20)+
  theme(legend.position="right")+
  theme(plot.margin = unit(c(0,0,0.5, 0),"cm")) #top, right, bottom, left


# scatter plot to visualize inter-method distribution of brain volumes
m_method_cor = m_method %>% 
  group_by(tissue) %>%
  summarise(spm_fsl=paste0("r = ", round(cor.test(SPM, FSL, na.rm=TRUE)$estimate,2)), 
            spm_fs=paste0("r = ", round(cor.test(SPM, FS, na.rm=TRUE)$estimate,2)))
# to avoid plotting CSF correlation because CSF_FS is not the total CSF
m_method_cor[m_method_cor$tissue=="CSF", "spm_fs"]="" 

m_method1 = m_method
m_method1$SPM1 = m_method1$SPM
m_method1[m_method1$tissue=="CSF", "FS"] = NA # avoid plotting CSF_FS
for_scatter = melt(m_method1, measure.vars = c("SPM", "FSL", "FS"))
# to control the order of plotting
for_scatter$variable = factor(for_scatter$variable, levels=c("FSL", "FS", "SPM"))
# introducing dummy points to control the left bottom corner points
dummy = m_method1[1:4, ]
dummy["tissue"] = unique(m_method$tissue)
dummy[dummy$tissue=="TIV", c("SPM","FSL","FS", "SPM1")]=0.75
dummy[dummy$tissue=="GM", c("SPM","FSL","FS", "SPM1")]=0.4
dummy[dummy$tissue=="WM", c("SPM","FSL","FS", "SPM1")]=0.2
dummy[dummy$tissue=="CSF", c("SPM","FSL","FS", "SPM1")]=0.1
mdummy = melt(dummy, measure.vars = c("SPM", "FSL", "FS"))
mdummy$variable = factor(mdummy$variable, levels=c("FSL", "FS", "SPM"))


scatter_plot = ggplot(for_scatter, aes(x=SPM1, y=value)) + 
  geom_point(aes(color=variable, order=variable), alpha=0.7)+ 
  geom_text(data=m_method_cor, x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, 
            aes(label = spm_fsl), color="#56B4E9", size=9, family="Garamond")+
  geom_text(data=m_method_cor, x=-Inf, y=Inf, hjust=-0.2, vjust=2.3,
            aes(label = spm_fs), color="#E69F00", size=9, family="Garamond")+
  geom_point(data = mdummy, aes(color=variable), alpha=0 )+ 
  ylab("Volume (Liters) by FSL & FS")+
  xlab("Volume (Liters) by SPM")+
  ggtitle("B) Inter-method correlation of the brain volumes with respect to SPM")+
  # geom_abline(intercept=0, slope=1, size=1.5, aes(color="#CC79A7"))+
  scale_color_manual(limits = c("SPM", "FSL", "FS"),  values=c( "#CC79A7","#56B4E9","#E69F00"), name="Methods") +
  # increase the size of the points inside the legend key
  guides(colour = guide_legend(override.aes = list(size=10)))+
  facet_wrap(~tissue, scales="free", nrow = 1)+
  theme_bw(base_size=20)+
  theme(legend.position="right")+
  theme(plot.margin = unit(c(-0.2,0,1, 0),"cm")) #top, right, bottom, left

cairo_pdf('Figure1_method_diff_all.pdf', height=8, width=12)
grid.arrange(box_plot, scatter_plot)
dev.off()

# -----
# Absolute inter-method difference p-val and Cohen's d
pcd = inter_method_diff(m_method)


# -------------- this might be helpful too, just in case

# vars_to_test <- c("disp","hp","drat","wt","qsec")
# iv <- "vs"
# 
# mtcars %>%
#   group_by(am) %>%
#   summarise_each_(
#     funs_( 
#       sprintf("stats::t.test(.[%s == 0], .[%s == 1])$p.value",iv,iv)
#     ), 
#     vars = vars_to_test)
# 

