### This script was written by Friederike Kessel,
### In the Workgroup of Prof. Hugo at the University Hospital Carl Gustav Carus
### Of the medical faculty of the University Dresden

setwd(choose.dir(caption="Directory with .lif-files", default=getwd()))
script_location<-dirname(sys.frame(1)$ofile) # directory with RScript (and ImageJ Macro)

inst_pack<-data.frame(installed.packages())$Package

# linear regression equation ----------------------------------------------
lm_eqn <- function(df){
  m <- lm(volume ~ frame, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 4),
                        r2 = format(summary(m)$r.squared, digits = 5)))
  as.character(as.expression(eq));
}

# ggplot2 -----------------------------------------------------------------
if("ggplot2"%in%inst_pack){
  library(ggplot2)
}else{
  warning("Package \"ggplot2\" is not installed")
  install.packages("ggplot2")
  library(ggplot2)
}

mytheme3 <- theme(legend.text = element_text(size = 15), 
                  axis.title = element_text(size = 15),
                  plot.subtitle = element_text(size = 15,hjust=0.5),
                  axis.text = element_text(size = 13), 
                  axis.line = element_line(size = 1,colour = "black"), 
                  axis.ticks = element_line(colour="black",size = rel(3)),
                  panel.background = element_rect(fill = "white", colour="black"), 
                  legend.key = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour="grey"),
                  panel.grid.minor = element_blank(),
                  legend.background = element_rect(fill="white",colour="black"),
                  legend.title = element_text(face = "bold", size = 15), 
                  plot.title = element_text(face = "bold",
                                            size = 18,hjust=0.5),
                  strip.text.x = element_text(size = 15))


# ImageJ Macro (only execute if files have not all been processed yet)------------------------------------------------------------
dir_fiji<-choose.files(caption="Select executable FIJI file")
system(paste0(dir_fiji," -macro ",script_location, "/cont-3D-snGFR.ijm ", getwd()))

data_process<-function(filename){
  df1<-read.table(paste0("Results/", filename), header=T, sep="\t") # intensity measurement in time series
  nframes<-max(df1$frame) # number of frames
  df2<-matrix(NA, ncol=nframes, nrow=nrow(df1)/nframes)
  
  colnames(df2)<-1:nframes
  rownames(df2)<-df1$length[1:nrow(df2)]
  
  # volume over length of tubulus-----
  if(file.size(paste0("Results/Volume_",filename))>3){
    dfvolume<-read.table(paste0("Results/Volume_",filename), header=T, sep="\t") #volume measurement
    if(max(dfvolume$Volume)>0){
      PT_length<-df1$length[df1$frame==1]
      PT_volume<-approx(dfvolume$Volume, n=length(PT_length))$y
      PT_difvolume<-c(0,diff(PT_volume))
      
      df1$volume<-rep(PT_volume, times=nframes)
      df1$difvolume<-rep(PT_difvolume, times=nframes)
      
      # crop table to where volume rendering vs. position is increasing
      df2<-df1[df1$frame==1,]
      df1<-subset(df1, difvolume>0.5*mean(df2$difvolume))
      df2<-df1[df1$frame==1,]
      
      # plot volume against position
      p<-ggplot(df2, aes(x=length, y=volume))+
        geom_point()+
        mytheme3+
        labs(title="PT Volume against length", subtitle=filename, x="[µm]", y="[µm³]")
      
      print(p)
      ggsave(paste0("Graphs_", Sys.Date(),"/Pos+Vol_", filename, ".png"), width=15, height=10, unit="cm")
      
      # smooth intensity
      realpos<-df1$volume[df1$frame==1]
      maxdifintens<-NA
      maxdifslope<-NA
      x1<-0
      for(x in realpos){
        x1<-x1+1
        df1$intensity[df1$volume==x]<-smooth.spline(df1$intensity[df1$volume==x], spar=0.5)$y
        df1$difintens[df1$volume==x]<-c(NA, diff(df1$intensity[df1$volume==x]))
        maxdifintens[x1]<-df1$intensity[df1$volume==x][which(diff(df1$intensity[df1$volume==x])==max(diff(df1$intensity[df1$volume==x])))]
        maxdifslope[x1]<-df1$difintens[df1$volume==x][which(diff(df1$intensity[df1$volume==x])==max(diff(df1$intensity[df1$volume==x])))]
      }
      
      # plot max intensities per frame and plot max intensity curve ------------------------------------------------
      dfstat<-aggregate(df1$intensity, by=list(df1$frame), FUN="mean")
      dfstat$std<-aggregate(df1$intensity, by=list(df1$frame), FUN="sd")$x
      dfstat$max<-aggregate(df1$intensity, by=list(df1$frame), FUN="max")$x
      dfstat$min<-aggregate(df1$intensity, by=list(df1$frame), FUN="min")$x
      dfstat$max<-smooth(dfstat$max)
      dfstat$diffmax<-c(0, diff(dfstat$max))
      
      maxdiffmax<-max(dfstat$diffmax) # maximum slope of intensity changes
      
      min_frame<-min(dfstat$Group.1[dfstat$diffmax>=0.5*maxdiffmax])-5
      max_frame<-max(dfstat$Group.1[dfstat$max>=0.8*max(dfstat$max)])+5
      
      TS_intens<-median(maxdifintens) # ts intensity in turning point of curve (median of max slope across positions)
      
      p<-ggplot(dfstat, aes(x=Group.1, y=max))+
        geom_jitter()+
        expand_limits(y=0)+
        labs(subtitle=filename, title="Maximum intensity per frame", y="intensity", x="frame")+
        mytheme3+
        geom_hline(yintercept=TS_intens)+
        scale_y_continuous()
      
      print(p)
      
      ggsave(paste0("Graphs_", Sys.Date(),"/Maximum_", filename, ".png"), width=15, height=10, unit="cm")
      
      # reduce dataframe: max intensity and slope -----------
      df2<-df1
      df2<-na.omit(df2)
      midx<-max(df2$volume)/2
      realpos<-df2$volume[df2$frame==min_frame]
      df2$difintens<-NA
      
      for(x in 1:length(realpos)){
        df2$intensity[df2$volume==realpos[x]]<-smooth.spline(df2$intensity[df2$volume==realpos[x]], spar=0.5)$y
        df2$difintens[df2$volume==realpos[x]]<-c(0, diff(df2$intensity[df2$volume==realpos[x]]))
        
        dfy<-subset(df2, volume==realpos[x])
        
        if(!is.na(maxdifslope[x])){
          if(maxdifslope[x]<0.8*median(maxdifslope, na.rm=T)){
            df2<-df2[-which(df2$volume==realpos[x]),]
          }
        }else{
          df2<-df2[-which(df2$volume==realpos[x]),]
        }
      }
      
      # reduce dataframe for printing (print only 20 positions)---------------
      realpos<-df2$volume[df2$frame==min_frame]
      df3<-df2
      nopos<-length(df2$volume[df2$frame==min_frame]) # number of different positions
      facpos<-floor(nopos/20) # which positions to take
      if(nopos>20){
        realpos<-df2$volume[df2$frame==min_frame][(1:20)*floor(nopos/20)]
      }else{
        realpos<-realpos
      }
      
      df3<-subset(df3, df3$volume%in%realpos)
      
      p<-ggplot(df3, aes(x=frame, y=intensity, fill=volume, colour=volume, group=volume))+
        scale_colour_gradient2(midpoint=midx, low="blue", mid="gray",
                               high="red", space ="Lab" )+
        scale_fill_gradient2(midpoint=midx, low="blue", mid="gray",
                             high="red", space ="Lab" )+
        geom_line(size=1)+
        labs(title="Intensity changes over time for every position", subtitle=filename, y="Intensity", x="Frame", colour="Volume")+
        mytheme3+
        geom_hline(yintercept=as.numeric(TS_intens), size=2, colour="green")+
        expand_limits(y=0)+
        scale_x_continuous(limits=c(min_frame, max_frame))+
        geom_vline(xintercept = min_frame:max_frame, colour="black")
      
      print(p)
      ggsave(paste0("Graphs_", Sys.Date(),"/Cleaner_Frames_", filename, ".png"), width=20, height=15, unit="cm")
      
      # reduce to max intensity/position ("wave")-------------------
      realpos<-df2$volume[df2$frame==1]
      df3<-df2
      for(x in realpos){
        dfx<-subset(df3, volume==x)
        fmax<-dfx$frame[which(dfx$intensity==max(dfx$intensity))[1]]
        if(fmax<max(df2$frame)){
          df3<-df3[-which(df3$volume==x&df3$frame>fmax),]
        }
        
      }
      df2<-df3
      
      # intercept of position with frame at ts-intensity----------------
      dfx<-data.frame("frame"=min_frame:max_frame)
      dfx$intensity<-NA
      dfx$length<-NA
      dfx$volume<-NA
      dfx$difintens<-NA
      
      for(x in 1:nrow(dfx)){
        vx<-dfx$frame[x]
        dfy<-subset(df2,frame==vx)
        
        if(nrow(dfy)>1){
          fitlength<-loess(intensity~length, data=dfy)
          fitvolume<-loess(length~volume, data=dfy)
          fitdifintens<-loess(intensity~difintens, data=dfy)
          
          dfx$length[x]<-approx(x=dfy$intensity, y=dfy$length, xout=TS_intens)$y
          dfx$volume[x]<-approx(x=dfy$intensity, y=dfy$volume, xout=TS_intens)$y
          
          dfx$intensity[x]<-TS_intens
          dfx$difintens[x]<-approx(x=dfy$intensity, y=dfy$difintens, xout=TS_intens)$y
        }
      }
      
      dfx<-na.omit(dfx)
      
      
      # check for closest volume near ts-intensity above/below given volume--------
      
      # below
      vx=min(dfx$frame)-1
      dfy<-subset(df2,frame==vx)
      if(nrow(dfy)>1){
        if(max(dfy$intensity)>=0.95*TS_intens){
          
          dfx<-rbind(c(vx, max(dfy$intensity), dfy$length[dfy$intensity==max(dfy$intensity)],
                       dfy$volume[dfy$intensity==max(dfy$intensity)],dfy$difintens[dfy$intensity==max(dfy$intensity)]), dfx)
        }
      }
      # above
      vx=max(dfx$frame)+1
      dfy<-subset(df2,frame==vx)
      if(nrow(dfy)>1){
        if(max(dfy$intensity)<=1.05*TS_intens){
          
          dfx<-rbind(c(vx, max(dfy$intensity), dfy$length[dfy$intensity==max(dfy$intensity)],
                       dfy$volume[dfy$intensity==max(dfy$intensity)],dfy$difintens[dfy$intensity==max(dfy$intensity)]), dfx)
        }
      }
      
      # adjust units (frames to seconds, µm³ to nl) -----------------------------
      dfx$frame<-dfx$frame/60/framerate
      dfx$volume<-dfx$volume/1000000
      dfx<-subset(dfx, difintens>0)
      
      # linear regression -------------------------------------------------------
      if(nrow(dfx)>0){
        midx<-max(dfx$volume)/2
        
        m <- lm(volume ~ frame, dfx)
        b = unname(coef(m)[2])
        r2=summary(m)$r.squared
        
        p<-ggplot(dfx, aes(x=frame, y=volume, colour=volume))+
          geom_point()+
          geom_smooth(method="lm", formula=y~x)+
          geom_line()+
          scale_colour_gradient2(midpoint=midx, low="blue", mid="gray",
                                 high="red", space ="Lab" )+
          geom_text(x = max(dfx$frame, na.rm=T), y = 1.01*min(dfx$volume, na.rm=T), label = lm_eqn(dfx), parse = TRUE, size=5, hjust=1, colour="black")+
          labs(title="snGFR", subtitle=filename, y="[nl]", x="[min]", colour="Volume")+
          mytheme3
        print(p)
        
        ggsave(paste0("Graphs_", Sys.Date(),"/Regression_", filename, ".png"), width=20, height=15, unit="cm")
        dfsum<-rbind(dfsum, c(filename, median(diff(dfx$volume)/diff(dfx$frame)), b, r2, max(dfvolume$Volume), max(df1$length), nrow(dfx)))
        assign("dfsum", dfsum, envir = .GlobalEnv)
        
      }
    }
  }
}

# R Script ----------------------------------------------------------------

if(!dir.exists(paste0("Graphs_", Sys.Date()))){
  dir.create(paste0("Graphs_", Sys.Date()))
}

framerate<-as.numeric(readline(prompt="Framerate in time series [fps]: "))

filelist<-list.files("Results", pattern=".txt")
dfsum<-c("Name", "Median slope", "snGFR", "R-squared", "Volume", "Length","Datapoints")

for(a in 1:length(filelist)){
  if(length(grep("Volume", filelist[a]))!=1){
    if(length(grep("Settings", filelist[a]))!=1){
      if(file.exists(paste0("Results/Volume_",filelist[a]))){
        if(file.exists(paste0("Results/",filelist[a]))){
          data_process(filelist[a])
        }
      }
    }
  }
}

write.table(dfsum, paste0(Sys.Date(),"-Result_summary.txt"), sep="\t", row.names=F, col.names=F)


