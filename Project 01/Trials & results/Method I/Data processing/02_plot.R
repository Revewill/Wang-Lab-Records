# Setup
args <- commandArgs(trailingOnly=TRUE)
self <- grep('^--file=', commandArgs(), value=TRUE)
root <- if(length(args)) normalizePath(args[1]) else if(length(self)) dirname(normalizePath(sub('^--file=', '', self[1]))) else getwd()
raw <- file.path(root,'raw')
processed <- file.path(root,'processed')
out <- file.path(root,'figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
packages <- c('ggplot2','cowplot','magick','svglite','ggbeeswarm')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install packages first: install.packages(c(',paste(sprintf('"%s"',missing),collapse=', '),'))')
suppressPackageStartupMessages(library(ggplot2))
FONT <- 'sans'
DPI <- 600
BAR <- .12
INK <- '#30343B'
PAL <- c(burgundy='#713E50',rose='#A56778',pink='#C997A3',pale='#E9D5DB',slate='#526879',grey='#8FA4B3',light='#C4CFD6',sand='#D5B78C')
GROUP <- c(NT=PAL['slate'],T=PAL['burgundy']);names(GROUP)<-c('NT','T')
ASSAY <- c(AmpR=PAL['slate'],KanR=PAL['burgundy']);names(ASSAY)<-c('AmpR','KanR')
TRIAL <- c(R1=PAL['burgundy'],R2=PAL['slate'],R3=PAL['sand']);names(TRIAL)<-c('R1','R2','R3')
SHAPES <- c(R1=21,R2=24,R3=22)
STRENGTH <- c('0.5'=PAL['grey'],'1'=PAL['burgundy']);names(STRENGTH)<-c('0.5','1')
PLATES <- c('Amp','Kan','Amp_Kan')
PLATE_LABELS <- c(Amp='Amp',Kan='Kan',Amp_Kan='Amp/Kan')
read <- function(name) read.csv(file.path(processed,paste0(name,'.csv')),stringsAsFactors=FALSE,check.names=FALSE)
if(!file.exists(file.path(processed,'pipeline_version.csv')) || as.numeric(read('pipeline_version')$version[1])!=3) stop('Run the updated 01_process.R first.')
blank <- function() cowplot::ggdraw()
theme_set(theme_classic(base_size=10,base_family=FONT)+theme(
 axis.line=element_line(colour=INK,linewidth=.35),axis.ticks=element_line(colour=INK,linewidth=.3),
 axis.text=element_text(colour=INK),axis.title=element_text(colour=INK,size=10),
 legend.title=element_blank(),legend.position='top',legend.text=element_text(size=8),
 legend.key.size=grid::unit(3,'mm'),legend.margin=margin(0,0,0,0),
 plot.title=element_text(size=10,face='bold',hjust=.5),plot.margin=margin(7,10,7,7),strip.background=element_blank(),strip.text=element_text(face='bold',size=9)))

# Helpers
save_plot <- function(plot,name,width,height) {
 ggsave(file.path(out,paste0(name,'.svg')),plot,device=svglite::svglite,width=width,height=height,bg='white',limitsize=FALSE)
 ggsave(file.path(out,paste0(name,'.png')),plot,device='png',width=width,height=height,dpi=DPI,bg='white',limitsize=FALSE)
}
tag <- function(p,letter) cowplot::ggdraw(p)+cowplot::draw_label(letter,x=.005,y=.995,hjust=0,vjust=1,fontfamily=FONT,fontface='bold',size=13,colour=INK)
star_rows <- function(d) {
 d$annotation[is.na(d$annotation)]<-''
 d[nzchar(d$annotation),,drop=FALSE]
}
vertical_stars <- function(p,d) {
 if(!nrow(d)) return(p)
 p+geom_segment(data=d,aes(x=x1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=.4,colour=INK)+
   geom_text(data=d,aes(x=(x1+x2)/2,y=y,label=annotation),inherit.aes=FALSE,vjust=-.35,size=4,colour=INK)
}
summary_frame <- function(d,keys,value) {
 do.call(rbind,lapply(split(d,interaction(d[keys],drop=TRUE)),function(z) {
  v<-z[[value]];cbind(z[1,keys,drop=FALSE],mean=mean(v),sd=if(length(v)>1)sd(v) else NA_real_)
 }))
}
offset <- function(group) ifelse(group=='NT',-.14,.14)
trial_jitter <- function(run) unname(c(R1=-.045,R2=0,R3=.045)[run])
fixed_jitter <- function(n,width=.045) if(n<=1) 0 else seq(-width,width,length.out=n)

# Image padding
plate_files <- c('Fig.1a-1.png','Fig.1a-2, 1~1.png','Fig.1a-3, 3~1.png')
gel_files <- c('Fig.1d-1-1.Tif','Fig.1d-1-2.Tif','Fig.1d-2-1.tif','Fig.1d-2-2.tif')
load_images <- function(files,block) {
 images<-lapply(files,function(f) magick::image_background(magick::image_read(file.path(raw,f))[1],color='black',flatten=TRUE))
 info<-do.call(rbind,lapply(images,function(im) magick::image_info(im)[1,c('width','height')]))
 W<-max(info$width);H<-max(info$height)
 padded<-lapply(images,function(im) magick::image_extent(im,geometry=paste0(W,'x',H),gravity='center',color='black'))
 manifest<-data.frame(block=block,file=files,width=info$width,height=info$height,canvas_width=W,canvas_height=H,
                      pad_left=floor((W-info$width)/2),pad_top=floor((H-info$height)/2))
 list(images=padded,width=W,height=H,manifest=manifest)
}
plates<-load_images(plate_files,'a');gels<-load_images(gel_files,'d')
write.csv(rbind(plates$manifest,gels$manifest),file.path(out,'image_layout.csv'),row.names=FALSE)
image_grob <- function(im) grid::rasterGrob(as.raster(im),width=grid::unit(1,'npc'),height=grid::unit(1,'npc'),interpolate=FALSE)
trial_key <- function(runs) {
 labels<-c(R1='Mass ratio 1:1',R2='Molar ratio 1:1',R3='Molar ratio 3:1')
 p<-cowplot::ggdraw()+cowplot::draw_label('Editor : target',x=.05,y=.87,hjust=0,size=9,fontface='bold',fontfamily=FONT,colour=INK)
 yy<-seq(.65,.65-.23*(length(runs)-1),length.out=length(runs))
 for(i in seq_along(runs)) {
  r<-runs[i]
  p<-p+cowplot::draw_grob(grid::rectGrob(gp=grid::gpar(fill=NA,col=TRIAL[r],lwd=2)),x=.05,y=yy[i]-.035,width=.10,height=.07)+
    cowplot::draw_label(labels[r],x=.20,y=yy[i],hjust=0,size=8.5,fontfamily=FONT,colour=INK)
 }
 p
}

# Panel a
ratio_labels <- c('Mass ratio 1:1','Molar ratio 1:1','Molar ratio 3:1')
plate_originals <- lapply(plate_files,function(f) {
 im<-magick::image_background(magick::image_read(file.path(raw,f))[1],color='black',flatten=TRUE)
 im
})
plate_columns <- c(4,3,3)
plate_scaled <- lapply(seq_along(plate_originals),function(i) magick::image_resize(plate_originals[[i]],paste0(plate_columns[i]*600-80)))
plate_height <- max(vapply(plate_scaled,function(im) magick::image_info(im)$height,numeric(1)))+80
plate_tiles<-lapply(seq_along(plate_scaled),function(i) {
 w<-plate_columns[i]*600
 im<-magick::image_extent(plate_scaled[[i]],paste0(w,'x',plate_height),gravity='center',color='black')
 grob<-grid::grobTree(image_grob(im),grid::rectGrob(width=.996,height=.996,gp=grid::gpar(fill=NA,col=TRIAL[i],lwd=2.5)))
 tile<-ggplot()+annotation_custom(grob,xmin=0,xmax=w,ymin=0,ymax=plate_height)+
  coord_fixed(xlim=c(0,w),ylim=c(0,plate_height),expand=FALSE)+theme_void()+theme(plot.margin=margin(2,4,2,4))
 label<-cowplot::ggdraw()+cowplot::draw_label(ratio_labels[i],fontfamily=FONT,size=9,colour=TRIAL[i])
 cowplot::plot_grid(label,tile,ncol=1,rel_heights=c(.09,1))
})
a_images<-cowplot::plot_grid(plotlist=plate_tiles,nrow=1,rel_widths=plate_columns)
p_a<-tag(a_images,'a')

# Panel b
b<-read('panel_b_points');bs<-read('panel_b_summary');bt<-star_rows(read('panel_b_tests'))
b$x<-match(b$plate,PLATES)+offset(b$group)+trial_jitter(b$run)
bs$x<-match(bs$plate,PLATES)+offset(bs$group)
positive<-b$relative_count[b$relative_count>0]
log_floor<-floor(log10(min(positive)))-.55
bs$low<-pmax(bs$mean-bs$sd,10^log_floor);bs$high<-bs$mean+bs$sd
bs$truncated<-(bs$mean-bs$sd)<=0
b$plot_y<-ifelse(b$relative_count>0,log10(pmax(b$relative_count,.Machine$double.xmin)),log_floor-.4)
if(nrow(bt)) {
 bt$x1<-match(bt$plate,PLATES)-.14;bt$x2<-match(bt$plate,PLATES)+.14
 bt$y<-vapply(bt$plate,function(pl) max(c(log10(bs$high[bs$plate==pl]),b$plot_y[b$plate==pl]))+.25,numeric(1))
}
ymax<-max(c(log10(bs$high),if(nrow(bt))bt$y+.4 else numeric()))
p_b<-ggplot()+
 geom_rect(data=bs,aes(xmin=x-BAR/2,xmax=x+BAR/2,ymin=log_floor,ymax=log10(mean),fill=group),colour=NA,alpha=.78)+
 geom_errorbar(data=bs,aes(x=x,ymin=log10(low),ymax=log10(high)),width=.075,linewidth=.4,colour=INK)+
 geom_point(data=subset(b,relative_count>0),aes(x=x,y=plot_y,shape=run,fill=group),size=2.2,stroke=.5,colour=INK)+
 geom_point(data=subset(b,relative_count==0),aes(x=x,y=plot_y,shape=run,fill=group),size=2.1,stroke=.5,colour=INK)+
 geom_hline(yintercept=log_floor-.2,linetype='dotted',linewidth=.3,colour=PAL['grey'])+
 geom_point(data=subset(bs,truncated),aes(x=x,y=log_floor),shape=25,fill='white',size=1.7,stroke=.35)+
 scale_fill_manual(values=GROUP,breaks=c('NT','T'))+scale_shape_manual(values=SHAPES,guide='none')+
 scale_x_continuous(breaks=1:3,labels=PLATE_LABELS,expand=expansion(add=.4))+
 scale_y_continuous(breaks=c(log_floor-.4,seq(ceiling(log_floor),ceiling(ymax))),labels=c('0',as.character(seq(ceiling(log_floor),ceiling(ymax)))),expand=expansion(mult=c(.03,.12)))+
 labs(x=NULL,y=expression(log[10]*' relative colony count'))
p_b<-tag(vertical_stars(p_b,bt),'b')

# Panel c
cp<-read('panel_c_points');cs<-read('panel_c_summary')
metric_plot <- function(metric,label,reference=NULL) {
 d<-cp[cp$metric==metric,,drop=FALSE];su<-cs[cs$metric==metric,,drop=FALSE]
 d$x<-1+trial_jitter(d$run)
 p<-ggplot()+geom_col(data=su,aes(x=1,y=mean),width=.12,fill=PAL['burgundy'],alpha=.72)+
  geom_errorbar(data=su,aes(x=1,ymin=mean-sd,ymax=mean+sd),width=.055,linewidth=.4)+
  geom_point(data=d,aes(x=x,y=value,shape=run),fill=PAL['burgundy'],colour=INK,stroke=.4,size=1.9)+
  scale_shape_manual(values=SHAPES,guide='none')+
  scale_x_continuous(breaks=1,labels='Pooled',limits=c(.55,1.45))+
  scale_y_continuous(expand=expansion(mult=c(.02,.16)))+labs(x=NULL,y=label)
 if(!is.null(reference))p<-p+geom_hline(yintercept=reference,linetype='dashed',linewidth=.3,colour=PAL['grey'])
 p
}
c1<-metric_plot('co_selection_index','Co-selection index')
c2<-metric_plot('Kan_bias_index','Kan bias index',1)
c3<-metric_plot('apparent_editing_percent','Apparent editing (%)')
p_c<-tag(cowplot::plot_grid(c1,c2,c3,nrow=1,align='h',axis='tb'),'c')

# Panel d
gel_row <- function(first,run) {
 gutter<-14
 canvas<-magick::image_blank(gels$width*2+gutter,gels$height,color='black')
 canvas<-magick::image_composite(canvas,gels$images[[first]],offset='+0+0')
 canvas<-magick::image_composite(canvas,gels$images[[first+1]],offset=paste0('+',gels$width+gutter,'+0'))
 grob<-grid::grobTree(image_grob(canvas),grid::rectGrob(width=.997,height=.995,gp=grid::gpar(fill=NA,col=TRIAL[run],lwd=3)))
 ggplot()+annotation_custom(grob,xmin=0,xmax=2*gels$width+gutter,ymin=0,ymax=gels$height)+
 coord_fixed(xlim=c(0,2*gels$width+gutter),ylim=c(0,gels$height),expand=FALSE)+theme_void()+theme(plot.margin=margin(3,3,3,3))
}
d_header<-cowplot::ggdraw()+cowplot::draw_label('AmpR',x=.25,y=.5,fontface='bold',fontfamily=FONT,size=10)+
 cowplot::draw_label('KanR',x=.75,y=.5,fontface='bold',fontfamily=FONT,size=10)
p_d<-tag(cowplot::plot_grid(d_header,gel_row(1,'R1'),gel_row(3,'R3'),ncol=1,rel_heights=c(.08,1,1)),'d')

# Panel e
e<-read('panel_e_points');es<-read('panel_e_summary');et<-star_rows(read('panel_e_test'))
e$x<-match(e$assay,c('AmpR','KanR'))

es$x<-match(es$assay,c('AmpR','KanR'))
p_e<-ggplot()+geom_col(data=es,aes(x=x,y=mean,fill=assay),width=.13,alpha=.72)+
 geom_errorbar(data=es,aes(x=x,ymin=mean-sd,ymax=mean+sd),width=.09,linewidth=.4)+
 ggbeeswarm::geom_quasirandom(data=e,aes(x=x,y=signal,fill=assay,shape=run,group=assay),width=.16,colour=INK,size=1.6,alpha=.8,stroke=.3)+
 scale_fill_manual(values=ASSAY,guide='none')+scale_shape_manual(values=SHAPES,guide='none')+
 scale_x_continuous(breaks=1:2,labels=c('AmpR','KanR'),limits=c(.55,2.45))+
 scale_y_continuous(expand=expansion(mult=c(.02,.18)))+labs(x=NULL,y='Relative PCR signal')
if(nrow(et)) {et$x1<-1;et$x2<-2;et$y<-max(c(e$signal,es$mean+es$sd))*1.07}
p_e<-tag(vertical_stars(p_e,et),'e')

# Panel f
f<-read('panel_f_MA_main');f$plate<-factor(f$plate,levels=PLATES)
f$group<-factor(f$group,levels=c('NT','T'))
p_f<-ggplot(f,aes(x=MA_A,y=MA_M))+geom_hline(yintercept=0,linetype='dashed',colour=PAL['grey'],linewidth=.4)+
 geom_point(aes(fill=plate,shape=group),size=2.6,colour=INK,stroke=.45,alpha=.9)+
 scale_fill_manual(values=c(Amp=PAL[['burgundy']],Kan=PAL[['slate']],Amp_Kan=PAL[['sand']]),breaks=PLATES,labels=PLATE_LABELS)+scale_shape_manual(values=c(NT=21,T=24))+
 labs(x=expression('Mean log'[2]*' signal'),y=expression(log[2]*'(AmpR/KanR)'))+
 guides(fill=guide_legend(order=1,override.aes=list(shape=21)),shape=guide_legend(order=2,override.aes=list(fill='white')))+
 theme(legend.box='vertical',legend.spacing.y=grid::unit(0,'mm'),legend.margin=margin(0,0,0,0))
p_f<-tag(p_f,'f')

# Panel g
g<-read('panel_g_points');gs<-read('panel_g_summary');gt<-star_rows(read('panel_g_tests'))
g$y<-4-match(g$plate,PLATES)+ifelse(g$group=='NT',.14,-.14)
for(pl in PLATES)for(gr in c('NT','T')) {
 ix<-which(g$plate==pl & g$group==gr);g$y[ix]<-g$y[ix]+fixed_jitter(length(ix),.045)
}
gs$y<-4-match(gs$plate,PLATES)+ifelse(gs$group=='NT',.14,-.14)
gs$low<-gs$mean-gs$sd;gs$high<-gs$mean+gs$sd
xmin<-min(c(0,g$log2_Amp_over_Kan,gs$low));xmax<-max(c(0,g$log2_Amp_over_Kan,gs$high));span<-xmax-xmin
p_g<-ggplot()+geom_vline(xintercept=0,linewidth=.35,colour=PAL['grey'])+
 geom_rect(data=gs,aes(xmin=pmin(0,mean),xmax=pmax(0,mean),ymin=y-BAR/2,ymax=y+BAR/2,fill=group),alpha=.45)+
 geom_segment(data=gs,aes(x=low,xend=high,y=y,yend=y),linewidth=.4,colour=INK)+
 geom_segment(data=gs,aes(x=low,xend=low,y=y-.035,yend=y+.035),linewidth=.4,colour=INK)+
 geom_segment(data=gs,aes(x=high,xend=high,y=y-.035,yend=y+.035),linewidth=.4,colour=INK)+
 geom_point(data=g,aes(x=log2_Amp_over_Kan,y=y,shape=run,fill=group),colour=INK,size=1.5,stroke=.3,alpha=.75)+
 scale_fill_manual(values=GROUP)+scale_shape_manual(values=SHAPES,guide='none')+
 scale_y_continuous(breaks=1:3,labels=rev(PLATE_LABELS),limits=c(.5,3.5))+
 scale_x_continuous(expand=expansion(mult=c(.04,.15)))+labs(x=expression(log[2]*'(AmpR/KanR)'),y=NULL)
if(nrow(gt)) {
 gt$y1<-4-match(gt$plate,PLATES)-.14;gt$y2<-gt$y1+.28
 gt$x<-xmax+.06*span
 p_g<-p_g+geom_segment(data=gt,aes(x=x,xend=x,y=y1,yend=y2),inherit.aes=FALSE,linewidth=.4)+
  geom_segment(data=gt,aes(x=x-.018*span,xend=x,y=y1,yend=y1),inherit.aes=FALSE,linewidth=.4)+
  geom_segment(data=gt,aes(x=x-.018*span,xend=x,y=y2,yend=y2),inherit.aes=FALSE,linewidth=.4)+
  geom_text(data=gt,aes(x=x+.025*span,y=(y1+y2)/2,label=annotation),inherit.aes=FALSE,hjust=0,size=3.8)
}
p_g<-tag(p_g,'g')

# Panel h
hc<-read('panel_h_colony_strengths');hc$x<-match(hc$strength,c(.5,1));hc$strength_key<-as.character(hc$strength)
h1<-ggplot(hc,aes(x=x,y=100*T_over_NT,fill=strength_key))+geom_col(width=.12,alpha=.78)+
 geom_point(shape=21,size=2,colour=INK,stroke=.4)+scale_fill_manual(values=STRENGTH,guide='none')+
 scale_x_continuous(breaks=1:2,labels=c('0.5x','1x'),limits=c(.5,2.5))+
 scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(x='Antibiotic strength',y='Dual-selection T/NT (%)')
hg<-read('panel_h_gel_points');hs<-read('panel_h_log2_Amp_over_Kan_summary')
hg$strength_key<-as.character(hg$strength);hs$strength_key<-as.character(hs$strength)
hg$x<-ifelse(hg$group=='NT',1,2)+ifelse(hg$strength==.5,-.14,.14)
for(gr in c('NT','T'))for(st in c(.5,1)) {
 ix<-which(hg$group==gr & hg$strength==st);hg$x[ix]<-hg$x[ix]+fixed_jitter(length(ix),.04)
}
hs$x<-ifelse(hs$group=='NT',1,2)+ifelse(hs$strength==.5,-.14,.14)
h2<-ggplot()+geom_hline(yintercept=0,linewidth=.3,colour=PAL['grey'])+
 geom_col(data=hs,aes(x=x,y=mean,fill=strength_key),width=BAR,alpha=.78)+
 geom_errorbar(data=hs,aes(x=x,ymin=mean-sd,ymax=mean+sd),width=.045,linewidth=.4)+
 geom_point(data=hg,aes(x=x,y=log2_Amp_over_Kan,fill=strength_key),shape=21,size=1.9,stroke=.4,colour=INK)+
 scale_fill_manual(values=STRENGTH,breaks=c('0.5','1'),labels=c('0.5x','1x'))+
 scale_x_continuous(breaks=1:2,labels=c('NT','T'),limits=c(.5,2.5))+
 scale_y_continuous(expand=expansion(mult=c(.1,.12)))+labs(x=NULL,y=expression(log[2]*'(AmpR/KanR)'))
p_h<-tag(cowplot::plot_grid(h1,h2,nrow=1,align='h',axis='tb'),'h')

# Exports
save_plot(p_b,'Fig.1b',4.2,3.3)
save_plot(p_c,'Fig.1c',5.3,3.3)
save_plot(c1,'Fig.1c-1',1.8,3.3)
save_plot(c2,'Fig.1c-2',1.8,3.3)
save_plot(c3,'Fig.1c-3',1.8,3.3)
save_plot(p_e,'Fig.1e',3.5,3.3)
save_plot(p_f,'Fig.1f',4.6,4.8)
save_plot(p_g,'Fig.1g',5.2,3.2)
save_plot(p_h,'Fig.1h',7.8,3.2)
save_plot(h1,'Fig.1h-1',3.9,3.2)
save_plot(h2,'Fig.1h-2',3.9,3.2)

# Assembly
trial_legend<-cowplot::ggdraw()
for(i in 1:3) {
 xx<-(i-1)/3+.06
 trial_legend<-trial_legend+cowplot::draw_grob(grid::pointsGrob(x=.5,y=.5,pch=SHAPES[i],size=grid::unit(2.2,'mm'),gp=grid::gpar(col=INK,fill='white')),x=xx,y=.25,width=.035,height=.5)+
 cowplot::draw_label(ratio_labels[i],x=xx+.045,y=.5,hjust=0,size=8,fontfamily=FONT,colour=INK)
}
row_bce<-cowplot::plot_grid(p_b,p_c,p_e,nrow=1,rel_widths=c(4.2,5.3,3.5))
row_bce<-cowplot::plot_grid(row_bce,trial_legend,ncol=1,rel_heights=c(1,.07))
row_df<-cowplot::plot_grid(p_d,p_f,nrow=1,rel_widths=c(8.4,4.6))
row_gh<-cowplot::plot_grid(p_g,p_h,nrow=1,rel_widths=c(5.2,7.8))
figure<-cowplot::plot_grid(p_a,row_bce,row_df,row_gh,ncol=1,rel_heights=c(3.1,3.55,4.8,3.2))
save_plot(figure,'Fig.1',13,14.65)
writeLines(c(
 'a: full original images retained; proportional scaling with at least 40 pixels of black padding on each side; no trimming. d: images padded black without cropping or intensity adjustment.',
 'Trial frames: burgundy, mass 1:1; slate, molar 1:1; sand, molar 3:1. Ratios are editor:target.',
 'b: linear relative-count means and sample SD displayed in log10 coordinates. The lower display floor is not zero; measured zeros appear in a separate band labelled 0; this band is not part of the logarithmic scale. Hollow downward triangles mark SD intervals extending to zero or below.',
 'b/c: points denote trials. Point shapes: circle R1; triangle R2; square R3.',
 'c: one bar per metric, with three trial-level points and sample SD. NT/T indices are averaged within each trial before pooling trials; editing remains a Kan-corrected T/NT contrast. No significance test in c. Model outputs are not clipped.',
 'e: 35 selected colony pairs, including the half-strength samples; bars show pooled means and sample SD. Stars test paired log2 signal contrast conditional on trial/source.',
 'f: 26 main-strength colony pairs. Colour denotes selection plate; shape denotes NT/T. No inferred clusters are imposed.',
 'g: six bars show means and SD of trial means; small points show individual colonies. Stars are conditional colony-level T/NT contrasts from script 01.',
 'h: trial02 only. Descriptive strength comparison; no formal equivalence claim.',
 'All stars are read from processed annotation columns; no raw p values or n.s. are drawn.',
 'PCR signal differences are not direct plasmid copy-number measurements.',
 'SVG contains vector chart elements and embedded original raster images; PNG is exported at 600 dpi.'
),file.path(out,'figure_notes.txt'))
message('Completed: individual panels and Fig.1.svg / Fig.1.png in ',out)
