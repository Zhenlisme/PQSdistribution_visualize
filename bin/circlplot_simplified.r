library(circlize)

Drowing_circos=function(PQScoordf,PQStype,center_title,opfilename){
  
  Bed = read.csv2(PQScoordf,sep="\t",header = F,stringsAsFactors = F)
  colnames(Bed)=c("ncnumber","start","end","G4SCI","pqscount","genename","straincount")
  strain_count=Bed$straincount[1]
  align_anotation=paste("There are ",as.character(strain_count)," sequences"," aviable for this alignment.",sep = "")
  
  sondf3=Bed[,c(6,2,3,4,5)]
  sondf3$genename=factor(sondf3$genename,levels = unique(sondf3$genename),ordered = T)
  sondf3$G4SCI=as.numeric(sondf3$G4SCI)
  
  genecounts=length(unique(sondf3$genename))
  genomesize=max(sondf3$end)
  if(genecounts<=50){
      jpeg(file=opfilename,width = 2500,height = 2500)
      pointsize=4
      cexsize=2.6
      gapsize=2
      trackheight=0.12}
  else if(genecounts>50 & genecounts<=120){
      jpeg(file=opfilename,width = 3000,height = 3000)
      pointsize=3
      cexsize=2.6
      gapsize=2        
      trackheight=0.12}
  else if(genecounts>120 & genecounts<=150){
      jpeg(file=opfilename,width = 3200,height = 3200)
      pointsize=3
      cexsize=2.6
      gapsize=2        
      trackheight=0.12}
  else if(genecounts>150 & genecounts<=250){
      jpeg(file=opfilename,width = 3500,height = 3500)
      pointsize=2
      cexsize=2.5
      gapsize=1       
      trackheight=0.12}
  else if(genecounts>250 & genecounts<=300){
      jpeg(file=opfilename,width = 4000,height = 4000)
      pointsize=2
      cexsize=2.5
      gapsize=1       
      trackheight=0.12}
  else if(genecounts>300 & genecounts<=400){
      jpeg(file=opfilename,width = 4000,height = 4000)
      pointsize=2
      cexsize=2.5
      gapsize=0.8      
      trackheight=0.05}
  else if(genecounts>400 & genecounts<=500){
      jpeg(file=opfilename,width = 4000,height = 4000)
      pointsize=2
      cexsize=2.5
      gapsize=0.5     
      trackheight=0.05}
  else if(genecounts>500 & genecounts<=600){
      jpeg(file=opfilename,width = 4000,height = 4000)
      pointsize=2
      cexsize=1
      gapsize=0.5
      trackheight=0.05}
  else if(genecounts>600 & genecounts<=1000){
      jpeg(file=opfilename,width = 4000,height = 4000)
      pointsize=2
      cexsize=0.8
      gapsize=0.3
      trackheight=0.05}
  else{jpeg(file=opfilename,width = 4000,height = 4000)
      pointsize=2
      cexsize=0.8
      gapsize=0.1
      trackheight=0.05}
  circos.clear()
  circos.par(gap.after = gapsize, start.degree = 90, ADD = list(bg = NA, fg = NA),
               cell.padding = c(0.02, 0, 0.02, 0))

  
               
  circos.genomicInitialize(sondf3,axis.labels.cex = 3,tickLabelsStartFromZero = F,plotType = "axis",major.by=signif(signif(genomesize,1)%/%20,1))
  
  circo_colors=rep("grey",length(levels(sondf3$genename)))
  circo_colors[which(grepl('\\(\\+\\)',levels(sondf3$genename)))]="#F6F7D4"
  circo_colors[which(grepl('\\(-\\)',levels(sondf3$genename)))]="#D2F6C5"
  circos.track(ylim = c(0, 1),bg.col = circo_colors,bg.border = NA, track.height = 0.05)

    
  circos.track(ylim = c(0, 1),panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2],
                  CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE,
                  adj = c(1, 0.5), cex = cexsize,font = 2,col = 'blue')},
                  bg.border = NA,track.height = trackheight)
                                   
  
  
    
  circos.genomicTrack(sondf3,ylim=c(0,1),bg.col = NA,bg.border = NA,
                      panel.fun = function(region, value, ...) {                         
                          circos.genomicLines(region, value[[1]], type = "h",lwd=ifelse(value[[1]] != 0,pointsize+2,0),
                                              col = ifelse(value[[2]]> 0,"red",'black'))                                                                                                                      
                          circos.genomicPoints(region, abs(value[[1]]), 
                                             pch = 20,transparency = 0.8,cex=ifelse(value[[2]] != 0,pointsize,0),transparency = 0.5,
                                             col = ifelse(value[[2]]> 0,"red",'black'))                                            
                          circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 1, col = "black")
                          circos.yaxis(side = "left", at = seq(0,1,by = .2),
                                       sector.index = get.all.sector.index()[1], labels.cex = 2.5)})                             
  if(F){                                    
  circos.genomicTrack(sondf1, stack = TRUE,track.height = 0.05,bg.col=NA,bg.border = NA,
                      panel.fun = function(region, value, ...) {
                          circos.genomicPoints(region, value, pch = 16,transparency = 0.8,
                                               cex = ifelse(value[[1]] != 0,pointsize,0),
                                               col = ifelse(value[[1]] > 0,"red",ifelse(value[[1]]<0, "black","#C7E9C0")), ...)
                        })                                     
        }                              
    text(0, 0, center_title, cex = 4,font = 4,col = "blue")
    legend("bottom", legend = align_anotation,bty = 'n',cex=4,text.col = "blue",text.font = 4)
    mytitle=paste("The distribution of",PQStype,sep=" ")
    legend("top", legend = mytitle,bty = 'n',cex=4,text.col = "blue",text.font = 4)
    
    circos.clear()
    dev.off()
}

args=commandArgs(T)
Drowing_circos(PQScoordf=args[1],PQStype=args[2],center_title=args[3],opfilename=args[4])

