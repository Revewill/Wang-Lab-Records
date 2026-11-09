# Setup
PIPELINE_VERSION <- '3.0'
args <- commandArgs(trailingOnly=TRUE)
self <- grep('^--file=', commandArgs(), value=TRUE)
root <- if(length(args)) normalizePath(args[1]) else if(length(self)) dirname(normalizePath(sub('^--file=', '', self[1]))) else getwd()
raw <- file.path(root,'raw'); out <- file.path(root,'processed')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
read <- function(name) read.csv(file.path(raw,name),check.names=FALSE,stringsAsFactors=FALSE)

# Summary
summarise <- function(d,keys,value) {
  id <- interaction(d[keys],drop=TRUE,lex.order=TRUE)
  do.call(rbind,lapply(split(d,id),function(z){
    v<-z[[value]];v<-v[is.finite(v)]
    cbind(z[1,keys,drop=FALSE],n=length(v),mean=if(length(v))mean(v) else NA_real_,sd=if(length(v)>1)sd(v) else NA_real_)
  }))
}
annotate <- function(p) ifelse(is.na(p)|p>=.05,'',ifelse(p<.001,'***',ifelse(p<.01,'**','*')))

# HC3 contrasts
contrast_hc3 <- function(d,formula,contrast=NULL) {
 fit <- lm(formula,data=d); X<-model.matrix(fit)
 stopifnot(qr(X)$rank==ncol(X),nrow(X)>ncol(X))
 if(is.null(contrast)) contrast<-colMeans(X)
 h<-hatvalues(fit); bread<-solve(crossprod(X))
 V<-bread %*% crossprod(X,X*as.numeric(residuals(fit)^2/(1-h)^2)) %*% bread
 est<-sum(contrast*coef(fit)); se<-sqrt(drop(t(contrast)%*%V%*%contrast)); df<-df.residual(fit)
 tval<-est/se; p<-2*pt(-abs(tval),df)
 data.frame(estimate_log2=est,se_HC3=se,df=df,ci_low=est-qt(.975,df)*se,ci_high=est+qt(.975,df)*se,p_raw=p,n_colonies=nrow(d),n_trials=length(unique(d$run)),method='Fixed-effects OLS, HC3, approximate t inference',scope='Exploratory; conditional on observed trials and selected colonies')
}

# Inputs
counts<-read('colony_counts.csv')

# Lane mapping
make_lanes <- function(run, groups, plates, strengths, sizes) {
  n<-sum(sizes)
  data.frame(run=run,lane=seq_len(n),group=rep(groups,sizes),
             plate=rep(plates,sizes),strength=rep(strengths,sizes),
             colony_id=sprintf('%s_L%02d',run,seq_len(n)),review_flag='',
             stringsAsFactors=FALSE)
}
lanes<-rbind(
 make_lanes('R1',c('NT','NT','NT','NT','T','T','T','T'),
  c('Amp','Kan','Amp_Kan','Amp_Kan','Amp','Kan','Amp_Kan','Amp_Kan'),
  c(1,1,.5,1,1,1,.5,1),c(2,2,2,2,2,2,7,3)),
 make_lanes('R3',c('NT','NT','NT','T','T','T'),
  c('Amp','Kan','Amp_Kan','Amp','Kan','Amp_Kan'),rep(1,6),c(2,2,2,2,2,3)))
lanes$review_flag[(lanes$run=='R1' & lanes$lane %in% c(7,20)) |
                   (lanes$run=='R3' & lanes$lane==8)]<-'weak_signal_review'
stopifnot(nrow(counts)==20,!anyDuplicated(counts[c('run','group','plate','strength')]),all(counts$count>=0),all(counts$count==floor(counts$count)),!anyDuplicated(lanes[c('run','lane')]))

# Count normalization
normC<-aggregate(count~run,subset(counts,plate=='Kan' & strength==1),mean)
names(normC)[2]<-'denominator_mean_Kan_NT_T'
stopifnot(all(normC[[2]]>0))
counts<-merge(counts,normC,by='run',sort=FALSE)
counts$relative_count<-counts$count/counts$denominator_mean_Kan_NT_T
counts$zero_count<-counts$count==0
counts$log2_relative_count<-NA_real_
counts$log2_relative_count[!counts$zero_count]<-log2(counts$relative_count[!counts$zero_count])
counts$main_panel<-counts$strength==1
write(normC,'normalization_colony');write(counts,'colony_counts_normalized')

# Panel b
b<-subset(counts,main_panel);write(b,'panel_b_points')
write(summarise(b,c('plate','group'),'relative_count'),'panel_b_summary')
btests<-do.call(rbind,lapply(c('Amp','Kan','Amp_Kan'),function(pl){
 z<-subset(b,plate==pl); nt<-subset(z,group=='NT'); tt<-subset(z,group=='T');tt<-tt[match(nt$run,tt$run),]
 test<-t.test(tt$relative_count,nt$relative_count,paired=TRUE)
 data.frame(plate=pl,contrast='T minus NT',estimate=mean(tt$relative_count-nt$relative_count),ci_low=test$conf.int[1],ci_high=test$conf.int[2],p_raw=test$p.value,n_trials=nrow(nt),method='Paired two-sided t-test on linear relative counts; exploratory')
}))
btests$p_holm<-p.adjust(btests$p_raw,'holm');btests$annotation<-annotate(btests$p_holm);write(btests,'panel_b_tests')

# Panel c
wide<-reshape(b[c('run','group','plate','relative_count')],idvar=c('run','group'),timevar='plate',direction='wide')
A<-wide$relative_count.Amp;K<-wide$relative_count.Kan;D<-wide$relative_count.Amp_Kan
wide$co_selection_index<-D/(A+K-D)
wide$Kan_bias_index<-(K-A)/(A+K-2*D)
wide$inferred_Amp_only<-A-D;wide$inferred_Kan_only<-K-D
wide$population_model_valid<-(A>=D & K>=D & A+K>D)
wide$model_note<-ifelse(wide$population_model_valid,'Necessary count inequalities satisfied; model assumptions still required','Dual exceeds a single-selection yield; not a valid population fraction')
write(wide,'panel_c_population_indices')
eff<-do.call(rbind,lapply(unique(wide$run),function(r){
 z<-subset(wide,run==r);n<-subset(z,group=='NT');t<-subset(z,group=='T')
 rn<-n$relative_count.Amp_Kan/n$relative_count.Kan;rt<-t$relative_count.Amp_Kan/t$relative_count.Kan
 data.frame(run=r,NT_dual_over_Kan=rn,T_dual_over_Kan=rt,relative_dual_survival=rt/rn,apparent_editing_percent=100*(1-rt/rn))
}))
write(eff,'panel_c_apparent_editing')
cmean<-aggregate(cbind(co_selection_index,Kan_bias_index)~run,data=wide,FUN=mean)
cmean<-merge(cmean,eff[c('run','apparent_editing_percent')],by='run')
cpoints<-do.call(rbind,lapply(c('co_selection_index','Kan_bias_index','apparent_editing_percent'),function(v) {
 data.frame(run=cmean$run,metric=v,value=cmean[[v]],
 pooling=if(v=='apparent_editing_percent') 'Kan-corrected T/NT contrast per trial' else 'Arithmetic mean of NT and T indices within trial')
}))
write(cpoints,'panel_c_points')
write(summarise(cpoints,'metric','value'),'panel_c_summary')
write(data.frame(panel='c',test_status='not_applicable',reason='One descriptive pooled estimate per metric; no NT/T comparison'),'panel_c_test_status')

# Gel inputs
gels<-do.call(rbind,lapply(c('R1','R3'),function(r){
 file<-if(r=='R1') '20261105-Normalized lane signals.csv' else '20261106-Normalized lane signals.csv'
 z<-read(file);stopifnot(ncol(z)==3)
 names(z)<-c('lane','Amp_marker_normalized','Kan_marker_normalized');z$run<-r;z$source<-file;z
}))
stopifnot(nrow(gels)==35,all(gels$Amp_marker_normalized>0),all(gels$Kan_marker_normalized>0))
g<-merge(lanes,gels,by=c('run','lane'),sort=FALSE)
stopifnot(nrow(g)==35)

# Signal normalization
normG<-do.call(rbind,lapply(unique(g$run),function(r){
 ref<-subset(g,run==r & group=='NT' & plate=='Amp' & strength==1)
 stopifnot(nrow(ref)==2)
 a<-mean(ref$Amp_marker_normalized);k<-mean(ref$Kan_marker_normalized)
 data.frame(run=r,reference='NT Amp, both assays, arithmetic means then geometric combination',Amp_reference_mean=a,Kan_reference_mean=k,common_denominator=sqrt(a*k))
}))
g<-merge(g,normG[c('run','common_denominator')],by='run',sort=FALSE)
g$Amp_signal<-g$Amp_marker_normalized/g$common_denominator
g$Kan_signal<-g$Kan_marker_normalized/g$common_denominator
g$log2_Amp_over_Kan<-log2(g$Amp_signal/g$Kan_signal)
g$MA_A<-(log2(g$Amp_signal)+log2(g$Kan_signal))/2
g$MA_M<-g$log2_Amp_over_Kan;g$main_panel<-g$strength==1
stopifnot(max(abs(g$log2_Amp_over_Kan-log2(g$Amp_marker_normalized/g$Kan_marker_normalized)))<1e-12)
g<-g[order(g$run,g$lane),]
write(normG,'normalization_gel');write(g,'gel_pairs_normalized');write(subset(g,review_flag!=''),'gel_review_flags')

# Panel e
elong<-rbind(transform(g,assay='AmpR',signal=Amp_signal),transform(g,assay='KanR',signal=Kan_signal))
write(elong,'panel_e_points');write(summarise(elong,'assay','signal'),'panel_e_summary')
g$stratum<-interaction(g$group,g$plate,g$strength,drop=TRUE)
etest<-contrast_hc3(g,log2_Amp_over_Kan~factor(run)+stratum)
etest$contrast<-'Sample-weighted mean paired log2(AmpR/KanR); all 35 selected colonies'
etest$p_adjusted<-etest$p_raw;etest$annotation<-annotate(etest$p_adjusted);write(etest,'panel_e_test')

# Panel f
g$plate<-factor(g$plate,levels=c('Amp','Kan','Amp_Kan'))
write(g,'panel_f_MA_all');write(subset(g,main_panel),'panel_f_MA_main')

# Panel g
gmain<-subset(g,main_panel);write(gmain,'panel_g_points')
gtrial<-summarise(gmain,c('run','plate','group'),'log2_Amp_over_Kan');write(gtrial,'panel_g_trial_summary')
write(summarise(gtrial,c('plate','group'),'mean'),'panel_g_summary')
gtests<-do.call(rbind,lapply(c('Amp','Kan','Amp_Kan'),function(pl){
 z<-subset(gmain,plate==pl);z$group<-factor(z$group,levels=c('NT','T'))
 fit<-lm(log2_Amp_over_Kan~group+factor(run),z);con<-rep(0,length(coef(fit)));con[match('groupT',names(coef(fit)))]<-1
 res<-contrast_hc3(z,log2_Amp_over_Kan~group+factor(run),con);res$plate<-pl;res$contrast<-'T minus NT, log2(AmpR/KanR)';res
}))
gtests$p_holm<-p.adjust(gtests$p_raw,'holm');gtests$annotation<-annotate(gtests$p_holm);write(gtests,'panel_g_tests')

# Panel h
hc<-subset(counts,run=='R1' & plate=='Amp_Kan');write(hc,'panel_h_colony_points')
hw<-reshape(hc[c('group','strength','count','relative_count')],idvar='strength',timevar='group',direction='wide')
hw$T_over_NT<-hw$count.T/hw$count.NT;write(hw,'panel_h_colony_strengths')
write(data.frame(ratio_of_T_over_NT_half_to_full=hw$T_over_NT[hw$strength==.5]/hw$T_over_NT[hw$strength==1],interpretation='Descriptive agreement only; no equivalence test; no x2 concentration correction'),'panel_h_colony_comparison')
hg<-subset(g,run=='R1' & plate=='Amp_Kan');write(hg,'panel_h_gel_points')
for(v in c('Amp_signal','Kan_signal','log2_Amp_over_Kan'))write(summarise(hg,c('group','strength'),v),paste0('panel_h_',v,'_summary'))

# Test registry
registry<-rbind(
 data.frame(panel='b',comparison=btests$plate,p_raw=btests$p_raw,p_adjusted=btests$p_holm,annotation=btests$annotation),
 data.frame(panel='e',comparison='AmpR_vs_KanR',p_raw=etest$p_raw,p_adjusted=etest$p_adjusted,annotation=etest$annotation),
 data.frame(panel='g',comparison=gtests$plate,p_raw=gtests$p_raw,p_adjusted=gtests$p_holm,annotation=gtests$annotation))
registry$test_status<-ifelse(!is.finite(registry$p_adjusted),'not_estimable',ifelse(registry$p_adjusted<.05,'significant','not_significant'))
write(registry,'statistical_tests')
write(data.frame(version=PIPELINE_VERSION),'pipeline_version')

# Provenance
inputs<-file.path(raw,c('colony_counts.csv','20261105-Normalized lane signals.csv','20261106-Normalized lane signals.csv'))
write(data.frame(file=basename(inputs),bytes=file.info(inputs)$size,md5=unname(tools::md5sum(inputs))),'input_manifest')
write(data.frame(check=c('colony_rows','gel_pairs','gel_main_pairs','panel_e_points','real_zero_plates','invalid_population_model_rows'),value=c(nrow(counts),nrow(g),sum(g$main_panel),nrow(elong),sum(counts$zero_count),sum(!wide$population_model_valid))),'validation_counts')
writeLines(c('All p values are two-sided and exploratory.','b: 3 trial-paired tests; Holm correction across 3 plates.','e: all 35 selected colony pairs; trial/source fixed effects; HC3 conditional inference. Not 35 trial replicates.','g: 3 plate-specific conditional tests; Holm correction. Two observed gel trials only.','c: average NT/T indices within each trial, then mean and sample SD across three trial estimates. Editing uses one Kan-corrected estimate per trial. No NT/T test in c.','h: descriptive; no equivalence claim from absence of significance.','Endpoint-PCR signal ratios do not identify plasmid copy-number ratios.','Negative/undefined model outputs retained, never clipped. Zero counts have no finite log.','Marker-normalized signals divided by one common reference factor within each trial.','No p-value or significance direction is predetermined.'),file.path(out,'analysis_notes.txt'))
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
message('Completed: ',nrow(counts),' plate counts; ',nrow(g),' paired colonies; outputs in ',out)
