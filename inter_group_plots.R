# plots to show how inter-group differences are influences by inter-method differences
library(reshape2)
library(ggplot2)
library(dplyr)
library(effsize)
library(gridExtra)
library(broom) # to use tidy()

# custom boxplots
min.1sd.mean.1sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

min.1sd.median.1sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), median(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
# helper function to add stars to p-values
get_stars = function(x){
  if (x < 0.005){
    z = '**'
  } else if (x < 0.05){
    z = '* '
  } else{
    z = "  "
  }
  return(z)
}

# bar plots for methods 
diff_bar_plot = function(d, title, keep.legend=FALSE){
  d[d$method=="FS" & d$tissue == "CSF", "estimate"] = 0
  d[d$method=="FS" & d$tissue == "CSF", "p.value"] = 1
  
  d_intercept = subset(d, term=="(Intercept)")
  d = subset(d, term=="controlASD") # assuming ASD is second level in factors i.e. 1
  d$stars = sapply(d$p.value, get_stars)
  d$avg = d_intercept$estimate + 0.5*d$estimate #(average across all subjects)
  d$greater = factor(sign(d$estimate),  levels=c(1,-1))
  d$per_diff = round(d$estimate/d$avg*100,2)
  print(d[c("tissue","method", "term", "p.value", "per_diff")])
  
  p = ggplot(d, aes(x=factor(method), y=per_diff))+
    geom_bar(aes(fill=method), stat="identity", alpha=1) +
    # geom_text(aes(label=paste0(round(per_diff,1), stars)), vjust=0, size=6)+
    geom_text(aes(label=stars), vjust=0.55, size=10)+
    geom_text(aes(label=round(per_diff,1)), vjust=1, size=5)+
    geom_abline(intercept=0, slope=0, size=2)+
    xlab("")+
    ylab(" Mean ASD-TDC  \n Difference (%)")+
    ggtitle(title)+
    scale_fill_manual(values=c("#CC79A7",  "#56B4E9", "#E69F00"), name="Methods")+
    facet_wrap(~tissue, scales="free", nrow = 1)+
    theme_bw(base_size=20)+
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,-0.5, 0),"cm")) #top, right, bottom, left
  # print(p)
  return(p)
}


m = read.csv("../../data/brain_volumes_long1.csv")
# m = subset(m, site == "NYU") # if only NYU site

# changing factor levels to control the order in regression models and also in plots
m$control = factor(m$control, levels = c("TDC", "ASD"))
m$method = factor(m$method, levels = c("SPM", "FSL", "FS"))
m$tissue = factor(m$tissue, levels = c("TIV", "GM", "WM", "CSF"))


# inter-group box plot
m1 = m
m1[m1$method=="FS" & m1$tissue == "CSF", "value"] = 0 # to avoid plotting FS CSF
inter_group_box_plot = ggplot(m1, aes(x=factor(method), y=value))+
  stat_summary(aes(fill=factor(control), width=0.3), fun.data = min.1sd.mean.1sd.max, 
               geom = "boxplot", position=position_dodge(width=0.4))+
  # scale_fill_discrete(name="")+
  scale_fill_manual(name="", values = c("#009E73", "#D55E00"))+
  xlab("")+
  ylab("Volume (Liters)")+
  ggtitle("A(i) Inter-group distribution of raw brain volumes")+
  facet_wrap(~tissue, scales="free_y", nrow=1)+
  theme_bw(base_size=20)+
  theme(legend.position="right")+
  theme(plot.margin = unit(c(0,0,0, 0),"cm")) #top, right, bottom, left

# percent difference plots
m_stat = m %>% 
  filter(!is.na(value)) %>%
  group_by(method, tissue) %>% 
  do(tidy(lm(value~control, data = .)))
# multiple comparsion correction across method using FDR
m_stat = m_stat %>%
    group_by(tissue) %>%
    mutate(p.value = p.adjust(p.value, method="BH"))
raw = diff_bar_plot(m_stat, title='A(ii) Mean Inter-group difference in raw brain volumes')

# adjusting covariates
m_stat = m %>% 
  filter(!is.na(value)) %>% # remove FS CSF
  group_by(method, tissue) %>% 
  do(tidy(my.lm(value~control+sex+age+I(age^2)+site, data = .)))
# multiple comparsion correction across method using FDR
m_stat = m_stat %>%
    group_by(tissue) %>%
    mutate(p.value = p.adjust(p.value, method="BH"))
adjusted = diff_bar_plot(m_stat, title='B) Mean Inter-group volume difference after adjusting for covariates')

# adjusting covariates including FIQ
m_stat = m %>% 
  filter(!is.na(value)) %>% # remove FS CSF
  filter(!is.na(FIQ)) %>% # remove FIQ NAs
  group_by(method, tissue) %>% 
  do(tidy(my.lm(value~control+sex+age+I(age^2)+FIQ+site, data = .)))
# multiple comparsion correction across method using FDR
m_stat = m_stat %>%
    group_by(tissue) %>%
    mutate(p.value = p.adjust(p.value, method="BH"))
adjusted_FIQ = diff_bar_plot(m_stat, title='B) Mean Inter-group volume difference after adjusting for covariates')

# black line to put as border 
horizontal_line = ggplot()+
  geom_abline(aes(x=0,y=0), slope=0, intercept=0, size=1)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())+
  theme(plot.margin = unit(c(-0.2,0,0, -0.5),"cm")) #top, right, bottom, left

pdf('Figure2_group_diff.pdf', height=12, width=12)
grid.arrange(inter_group_box_plot, raw, nrow=2, heights=c(0.7, 1))
dev.off()

# with black line border
# pdf('Figure2_group_diff.pdf', height=12, width=12)
# grid.arrange(inter_group_box_plot, raw, horizontal_line,adjusted_FIQ, nrow=4, heights=c(0.7, 1, 0.2,1))
# dev.off()


# # ---------------sitewise analysis ----------
# diff_bar_plot_sitewise = function(d){
#   m_stat = subset(m_stat, term=="control1")
#   m_stat$stars = sapply(m_stat$p, get_stars)
#   m_stat$estimate = -m_stat$estimate # to make ASD-TDC
#   m_stat$greater = factor(sign(m_stat$estimate),  levels=c(1,-1))
#   # m_stat = rbind(m_stat, c("FS", "CSF", rep(NA, ncol(m_stat)-2)))
#   # m_stat[is.na(m_stat$estimate)] <- 0
#   p = ggplot(m_stat, aes(x=factor(method), y=estimate))+
#     geom_bar(aes(fill=greater), stat="identity") +
#     geom_text(aes(label=stars), size=10)+
#     xlab("Methods")+
#     ylab("ASD-TDC Volume Difference (Liters)")+
#     scale_fill_manual(values=c( "#D55E00", "#009E73"), name="", labels=c( "ASD>TDC", "ASD<TDC")) +
#     facet_grid(site~tissue, scales="free")+
#     theme(legend.position="top")
#   print(p)
#   return(p)
# }
# m_stat = m %>% 
#   filter(!is.na(value)) %>%
#   group_by(site, method, tissue) %>% 
#   do(tidy(lm(value~control, data = .))) 
# raw = diff_bar_plot_sitewise(m_stat)
# 
# m_stat = m %>% 
#   filter(!is.na(value)) %>% # remove FS CSF
#   filter(!is.na(FIQ)) %>% # remove FIQ NAs
#   group_by(site, method, tissue) %>% 
#   do(tidy(my.lm(value~control+age+sex+I(age^2)+FIQ, data = .))) 
# adjusted = diff_bar_plot_sitewise(m_stat)
# 
# pdf('group_diff_sitewise.pdf', height=14, width=8)
# grid.arrange(raw, adjusted, nrow=1)
# dev.off()

# without using broom

# m_stat = m %>% 
#   filter(!is.na(value)) %>%
#   group_by(site, method, tissue) %>% 
#   do(te=t.test(value~control, data=.)) %>% # with formula , summarise cannot be used directly
#   mutate(t=te$statistic, p=te$p.value)

# m_stat = m %>% 
#   filter(!is.na(value)) %>%
#   group_by(site, method, tissue) %>% 
#   do({tt=lm(value~control, data=.)
#   data.frame(t=tt$statistic, p=tt$p.value, delta=diff(tt$estimate)) 
#   })

# big boxplot
# box_plot = ggplot(data = m, aes(x=factor(control), y=value))+
#   # geom_boxplot(aes(fill=factor(control)),  notch = TRUE) + 
#   stat_summary(aes(fill=factor(method)), fun.data = min.1sd.median.1sd.max, geom = "boxplot",position=position_dodge(width=0.8))+
#   stat_summary(aes(fill=factor(method)), fun.data = min.1sd.mean.1sd.max, geom = "boxplot", position=position_dodge(width=0.8))+
#   ylab("Volume (Liters)")+
#   # scale_fill_manual(values=c("#D55E00", "#009E73"), name="ASD/TDC") +
#   facet_grid(tissue~site, scales="free")+
#   theme(legend.position="top")+
#   theme(plot.margin = unit(c(0,0,0, 0),"cm")) #top, right, bottom, left

