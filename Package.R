##### R package for SPORT analysis

message("Load Package SPORT....")


#Extern Functions:

##### Extern Functions#########
Barplot_SC=function(input,min_x=0,Text_cex=1, points=T, 
                    main="main", 
                    ylim=c(0,c(max(input,na.rm = TRUE)+3)),
                    mar=c(6,4,4,4), 
                    col_Schema=F,
                    col_Self=c("red", "grey"),
                    p_values=T,
                    AbstandxAxe=0.2,
                    AbstandBalken=0.1,
                    AbstandPWert=0.1){
  
  
  
  
  if(class(input)=="matrix"){len=ncol(input)}else{print("Input needs to be matrix")}
  par(mar=mar,xaxs = "i",yaxs = "i")
  maxx=0
  for(i in 1:ncol(input)){maxx=c(maxx,max(na.omit(input[,i])))}
  
  plot(x=c(1,len+1), y=c(c(min_x-c(max(maxx)/29)),max(maxx)), bty="n", type="n",ylab="", xlab="",xaxt="n", main=main,ylim=ylim)
  
  if (col_Schema==T){
    require(RColorBrewer); rf <- colorRampPalette(rev(brewer.pal(9,'Set1')))
    r <- rf(len)
    col=sample(r) 
  } else {col=col_Self}
  
  #plot
  for(i in 1: len){
    val=na.omit(input[,i])
    xx=c(i,i,i+0.8,i+0.8 )
    yy=c(min_x,mean(val), mean(val), min_x)
    polygon(x=xx,y=yy, col=col[i], border="black", cex=1)
  }
  for(i in 1: len){
    val=na.omit(input[,i])
    xx=c(i+0.4,i+0.4 )
    yy=c(mean(val), mean(val)+sd(val))
    polygon(x=xx,y=yy, col="black",cex=1)
    
    xx=c(i+0.3,i+0.5 )
    yy=c(mean(val)+sd(val), mean(val)+sd(val))
    polygon(x=xx,y=yy, col="black")
    
    
  }
  if(points==T){
    for(i in 1: len){
      val=na.omit(input[,i])
      points(x=rnorm(mean=i+0.4, sd=0.05, n=length(val)), y=as.numeric(val))
    }
  }
  
  if (p_values==T){
    if(ncol(input)>2){
      for(i in 1:c(ncol(input)-1)){
        p=t.test(as.numeric(na.omit(input[,i])),as.numeric(na.omit(input[,i+1])))$p.value
        text(x=c(i+0.8) ,y=c(max(input,na.rm = TRUE)+i/AbstandPWert), labels = paste("p=",round(p, digits = 3),cex=Text_cex, sep=""))
        polygon(x=c(i+0.4, i+1.4),y=c(max(input,na.rm = TRUE)+(i/AbstandBalken),max(input,na.rm = TRUE)+(i/AbstandBalken) ), col="black")
      }
      
    }
    else{
      p=t.test(as.numeric(na.omit(input[,1])),as.numeric(na.omit(input[,2])))$p.value
      text(x=1.8 ,y=c(max(input,na.rm = TRUE)+AbstandPWert), labels = paste("p=",round(p, digits = 3),cex=Text_cex, sep=""))
      polygon(x=c(1.4, 2.4),y=c(max(input,na.rm = TRUE)+AbstandBalken,max(input,na.rm = TRUE)+AbstandBalken ), col="black")
    }
  }
  
  
  #put Axisi
  polygon(x=c(0,len),y=c(0,0), col="black",cex=1)
  for(i in 1:len){
    polygon(y=c(c(-max(maxx)/35),max(maxx)/35),x=c(i+0.4,i+0.4), col="black",cex=1)
  }
  
  #input names
  
  text(seq(1,len,by=1)+0.4, par("usr")[3]-AbstandxAxe, srt = 60, adj= 1, xpd = TRUE,labels = colnames(input), cex=Text_cex)
}


vioplotDHH=function (x, range = 1.5,mar=c(6,4,4,4), h = NULL, ylim = NULL, names = colnames(x), 
                     horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                     lwd = 1, rectCol = "white", colMed = "white", pchMed = 19, 
                     at, add = FALSE, wex = 1, drawRect = TRUE, main="main") 
{
  par(mar=mar,las=1)
  datas = lapply(1:ncol(x), function(i){as.numeric(na.omit(x[,i]))})
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    
    if (!horizontal) {
      if (!add) {
        plot(NA,xlim = xlim, ylim = ylim, main=main, xaxt="n", bty="n", xlab="", ylab="")
        axis(2)
        
        text(1:length(names),par("usr")[3]-0.3, srt = 60, adj= 1, xpd = TRUE,labels = names, cex=0.8)
      }
      
      for (i in 1:n) {
        polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
                lty = lty, lwd = lwd)
        if (drawRect) {
          lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                lty = lty)
          rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
               q3[i], col = rectCol)
          points(at[i], med[i], pch = pchMed, col = colMed)
        }
      }
    }
  else {
    if (!add) {
      plot(NA,xlim = ylim, ylim = xlim,main=main, bty="n",xlab="", ylab="",xaxt="n")
      axis(1)
      
      text(1:length(names),par("usr")[3]-0.3, srt = 60, adj= 1, xpd = TRUE,labels = names, cex=0.8)
    }
    
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}
LineplotDHH= function(input,min_x=0,Text_cex=1, points=T,main=colnames(input)[i],
                      ylim=c(0.5,c(max(input,na.rm = TRUE)+0.3)),
                      lwd=1,
                      lty=1,
                      pAbstand=0.2,
                      pBalken=0.1,
                      mar=c(6,4,4,4),
                      col_Schema=F,
                      ParameterRücken=0.1,
                      col_Self=c("red", "grey"),
                      addtoPlot=F,
                      p_values=T){
  
  par(mar=mar,las=2)
  miny=min(input)-(c(min(input)/10))
  maxy=max(input)+(c(max(input)/10))
  if(addtoPlot==F){plot(x=c(1,ncol(input)), y=c(miny,maxy),type="n", xlab="", ylab="", xaxt="n", bty="n" )}
  
  for(i in 1:ncol(input)){points(x=i, y=mean(input[,i]), pch=19)
    sd=sd(input[,i])/sqrt(length(input[,i]))
    xx=c(i,i)
    yy=c(mean(input[,i])-sd, mean(input[,i])+sd)
    polygon(x=xx, y=yy)
    polygon(x=c(i-0.15, i+0.15), y=c( mean(input[,i])+sd, mean(input[,i])+sd))
    polygon(x=c(i-0.15, i+0.15), y=c( mean(input[,i])-sd, mean(input[,i])-sd))
  }
  points(colMeans(input), pch=19, type="l", lty=lty, lwd=lwd)
  if(addtoPlot==F){text(seq(1,ncol(input),by=1)+0.4, par("usr")[3]-ParameterRücken, srt = 60, adj= 1, xpd = TRUE,labels = colnames(input), cex=Text_cex)}
  
}
Heatmap_signifcance=function(input){
  sig_m=matrix(NA,ncol(input),ncol(input))
  for(i in 1:ncol(input)){
    data_sig1=as.numeric(na.omit(input[,i]))
    for(j in 1:ncol(input)){
      data_sig2=as.numeric(na.omit(input[,j])) 
      require(assertthat)
      if(are_equal(data_sig1,data_sig2)==F){sig=t.test(data_sig1,data_sig2)$p.value}else{sig=1}
      sig_m[i,j]=sig
    }
  }
  #multiple Testing
  for(i in 1: ncol(sig_m)){
    sig_m[,i]=p.adjust(sig_m[,i], n=length(sig_m[,i]))
  }
  #draw heatmap
  require(RColorBrewer)
  rf <- colorRampPalette((brewer.pal(9,'RdYlBu')))
  par(mar=c(6,6,6,6),xaxs = "i",yaxs = "i")
  image(sig_m,col=rf(11), axes=FALSE)
  axis(2, at=seq(0,1, length.out =ncol(input)), labels=colnames(input), lwd=0, pos=-0.25)
  text(seq(0,1, length.out =ncol(input)), par("usr")[1]-0.1, srt = 60, adj= 1, xpd = TRUE,labels = colnames(input), cex=1)
  x_y=seq(from=0, to=1, length.out = ncol(sig_m))
  
  for(i in 1: ncol(sig_m)){
    xxx=x_y[i]
    for(j in 1: ncol(sig_m)){
      yyy=x_y[j]
      text(x=xxx, y=yyy, labels = round(sig_m[i,j], digits = 3),cex=1)
      
    }
  }
}
map2color<-function(x,pal,limits=NULL){
  if(class(x)=="numeric"){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  }else{
    print(x[!duplicated(x)])
    da=data.frame(Terms=x[!duplicated(x)], Nr=seq(1,length.out = length(x[!duplicated(x)])))
    da$col=colorRampPalette(brewer.pal(6, "Set1"))(nrow(da))[da[,2]]
    daf=data.frame(x=x, col=1)
    for(i in 1:length(x)){
      daf[i, ]$col=da[da$Terms==daf[i, "x"], ]$col
      
    }
    
    return(list(daf$col, da))
  }
  
}






###Class and Getetter functions


SPORT=setClass("SPORT",slots=c(
  list_files="character", 
  coordinates="data.frame",
  clinical_data="data.frame",
  folder="character",
  Allign_list="data.frame",
  list_names_images="list", 
  list_Metabolite="list",
  list_id="list",
  TSNE_input="matrix",
  GroupMen="data.frame",
  phenodata="data.frame",
  TSNE_ID="data.frame",
  TSNE="data.frame",
  datout="matrix",
  ID_out="data.frame",
  PCA="data.frame",
  DimRed="data.frame",
  DE="data.frame"))


setMethod("initialize",signature = "SPORT", definition = function(.Object, list_files, clinical_data ){
  .Object@list_files <- list_files
  .Object@clinical_data <- clinical_data
  return(.Object)
})

setGeneric("Set_up_data", function(object,Image_Export=F, include_all=F) standardGeneric("Set_up_data"))
setMethod("Set_up_data",signature = "SPORT", definition = function(object, Image_Export=F, include_all=F){
  
  
  list_files=object@list_files
  print(paste("We found a total of ", length(list_files), " Files, all files will be processed"))
  
  print(paste("We start merging clinical data with uploaded files"))
  clindat=object@clinical_data
  
  #define a new var to define file names
  var_files=paste("A_",clindat$PIZ,".RDS", sep="")
  clindat=clindat[var_files %in%list_files, ]
  list_files=list_files[list_files %in% var_files]
  
  print(paste("We merged clinical data, process failed in ", 100-c(nrow(clindat)/length(list_files)*100), " %", sep=""))
  
  
  
  
  for(xx in 1:length(list_files)){
    
    print(paste("We start with ", list_files[xx], " , this is Nr", xx, "of ", length(list_files), " files in total"))
    Spx=readRDS(list_files[xx])
    
    # GET COORDINATES
    if(xx==1){
      
      Spx=readRDS(list_files[xx])
      object@coordinates <-
        Spx[[1]][1:2, ] %>% 
        t() %>% 
        as.data.frame() %>%  
        dplyr::mutate(RowCounts=as.numeric(1:nrow(.)) )
      
    }
 
    
    
    
    #head(Spx[[3]])
    
    Spx_files_raw=as.data.frame(Spx[[1]])
    #Format raw data into new format
    dim(Spx_files_raw)
    Metabolite=t(Spx_files_raw[6:39,])
    dim(Metabolite)
    
    
    #Seperate Different Regions and remove unknow regions
    
    if(include_all){
      
      message( "all files will be integrated, no filter for alligned voxel")
      
      Spx_files=as.data.frame(Spx[[3]]) %>% dplyr::mutate(ALL=1)
      id=Spx_files[, 1:2]
      names(id)[2]="ID"
      dim(id)
      rownames(Metabolite)=id$voxel_name
      id[,2]=c(NA)
      for(i in 2:(ncol(Spx_files))){
        index=Spx_files[Spx_files[, i]==1, ]$voxel_name
        id[id$voxel_name %in% index, 2]=names(Spx_files)[i]
      }
      id=na.omit(id)
      dim(id)
      object@list_id=c(object@list_id,list(id))
      print( paste( "We idendified ", nrow(id), " Identifier", sep=""))
    }
    else{
      
      Spx_files=as.data.frame(Spx[[3]]) #%>% dplyr::mutate(ALL=1)
      id=Spx_files[, 1:2]
      names(id)[2]="ID"
      dim(id)
      rownames(Metabolite)=id$voxel_name
      id[,2]=c(NA)
      for(i in 2:(ncol(Spx_files))){
        index=Spx_files[Spx_files[, i]==1, ]$voxel_name
        id[id$voxel_name %in% index, 2]=names(Spx_files)[i]
      }
      id=na.omit(id)
      dim(id)
      object@list_id=c(object@list_id,list(id))
      print( paste( "We idendified ", nrow(id), " Identifier", sep=""))
      
      
      
    }
    
    
    
    
    
    ###Output image
    #Input RDS file
    if(Image_Export==T){
      rds=Spx
      print( paste( "Images were Exported", sep=""))
      ### Change Folder for Image Output
      Images_dat=unlist(lapply(1:nrow(id), function(i){
        
        n=id[i,1]
        image_name=paste(list_files[2],"_",rownames(id)[i], "_Image_DeepL.png")
        #plot(rds[[2]][,n], type="l")
        
        #Example Spectra
        new_sp=data.frame(SP=rds[[2]][,n], row.names = rownames(rds[[2]]))
        
        #Convert 1D to 2D
        
        #plot(x=new_sp[,1],y=new_sp[,1])
        
        sqrt(nrow(new_sp))
        
        #Create a Matrix with 32x64 pixel
        
        image=matrix(1:nrow(new_sp),32,64)
        
        for(i in 1:nrow(new_sp)){
          image[image==i]=new_sp[i,1]
        }
        
        library(RColorBrewer)
        cols=(colorRampPalette(brewer.pal(9,"Greys"))(100))
        
        png(image_name)
        par(mai=c(0,0,0,0))
        image(1:nrow(image), 1:ncol(image), as.matrix(image), col=cols,axes=FALSE)
        dev.off()
        
        return(image_name)
      }))
      object@list_names_images=c(object@list_names_images,list(Images_dat))
    }
    
    
    
    
    
    #check metabolite 
    Metabolite=as.matrix((Metabolite[rownames(id), ]))
    dim(Metabolite)
    m <- mapply(Metabolite, FUN=as.numeric)
    dim(m)
    Metabolite1 <- matrix(data=m, ncol=ncol(Metabolite), nrow=nrow(Metabolite))
    
    colnames(Metabolite1)=colnames(Metabolite)
    rownames(Metabolite1)=rownames(Metabolite)
    
    dim(Metabolite1)
    
    object@list_Metabolite=c(object@list_Metabolite,list(Metabolite1))
    
    
    
  }
  
  #Allign list with spectra and origin
  infovox=object@list_id
  
  for(i in 1:length(infovox)){
    print(i)
    infovox[[i]]$file=list_files[i]
    infovox[[i]]$dataset=c(i)
    PIZ=gsub("A_","",list_files[i] )
    PIZ=gsub(".RDS","",PIZ )
    infovox[[i]]$PIZ=PIZ
    infovox[[i]]$Tumor=clindat[clindat$PIZ==PIZ, ]$Diag
    infovox[[i]]$WHO=clindat[clindat$PIZ==PIZ, ]$WHO.Grad
    infovox[[i]]$Geschlecht=clindat[clindat$PIZ==PIZ, ]$Geschlecht
  }
  infovox_new=as.data.frame(do.call(rbind, infovox))
  
  object@Allign_list=infovox_new
  
  return(object)
  
})

setGeneric("Check_Region", function(object,ylim=c(0,30),Index_Metabolite="Cr",Region=NA, Tumor=NA,listID="Region",plot="V", main="main") standardGeneric("Check_Region"))
setMethod("Check_Region",signature = "SPORT", definition = function(object,ylim=c(0,30),Index_Metabolite="Cr",Region=NA, Tumor=NA,listID="Region",plot="V", main="main"){
  
  Met=object@list_Metabolite
  ID_N=object@list_id
  
  
  if (is.null(object@Allign_list$ID_2)==T){object@Allign_list$ID_2=paste(object@Allign_list$PIZ,"_", object@Allign_list$voxel_name, sep="")}
  allign=object@Allign_list
  
  if(is.na(Region[1])==T){print("all regions");sepR=allign }else{
    sepR=allign[allign$ID %in% Region, ]
  }
  if(is.na(Tumor)[1]==T){print("all Tumor");sepT=sepR }else{
    sepT=sepR[sepR$Tumor %in% Tumor, ]
  }
  
  datasets=sepT[!duplicated(sepT$dataset), "dataset"]
  
  datout=as.data.frame(do.call(rbind, lapply(1:length(datasets), function(i){
    print(datasets[[i]])
    Met_1=as.data.frame(Met[[datasets[i]]][,Index_Metabolite])
    #get voxel nr
    nrow(Met_1)
    length(ID_N[[datasets[[i]]]]$voxel_name)
    
    
    set1=allign[allign$dataset==datasets[[i]], ]
    set1=set1[set1$voxel_name==rownames(Met_1), ]
    table(set1$voxel_name==rownames(Met_1))
    rownames(Met_1)=set1$ID_2
    
    names(Met_1)=Index_Metabolite
    Met_1$ID=ID_N[[datasets[i]]][,2]
    Met_1=Met_1[Met_1$ID %in% sepT[!duplicated(sepT$ID), "ID"], ]
    Met_1$Tumor=allign[allign$dataset==i, "Tumor" ][1]
    
    return(Met_1)
  })))
  
  
  if(listID=="Region"){
    #combine Matrix
    type=datout[!duplicated(datout$ID), "ID"]
    mat=matrix(NA, nrow(datout), length(type))
    colnames(mat)=type
    for(j in 1:length(type)){
      inp=datout[datout$ID==type[j], 1] 
      mat[, j]=c(inp,rep(NA, nrow(mat)-length(inp)))
    }
    
    #par(mfrow=c(1,2))
    if(plot=="Bar"){Barplot_SC(mat,ylim=ylim, main=main)}
    if(plot=="V"){vioplotDHH(mat,ylim=ylim, col=c("red", "grey", "white"),main=main)}
    Heatmap_signifcance(mat)
  }
  if(listID=="Index_Metabolite"){
    #combine Matrix
    type=Index_Metabolite
    mat=matrix(NA, nrow(datout), length(type))
    colnames(mat)=type
    for(j in 1:length(type)){
      inp=datout[, type[j]] 
      mat[, j]=c(inp,rep(NA, nrow(mat)-length(inp)))
    }
    
    #par(mfrow=c(1,2))
    if(plot=="Bar"){Barplot_SC(mat,ylim=ylim,main=main)}
    if(plot=="V"){vioplotDHH(mat,ylim=ylim, col=c("red", "grey", "white"),main=main)}
    Heatmap_signifcance(mat)
  }
  if(listID=="Tumor"){
    #combine Matrix
    type=Tumor
    mat=matrix(NA, nrow(datout), length(type))
    colnames(mat)=type
    for(j in 1:length(type)){
      inp=datout[datout$Tumor==type[j], 1] 
      mat[, j]=c(inp,rep(NA, nrow(mat)-length(inp)))
    }
    #par(mfrow=c(1,2))
    if(plot=="Bar"){Barplot_SC(mat,main=Index_Metabolite,ylim=ylim)}
    if(plot=="V"){vioplotDHH(mat,main=Index_Metabolite,ylim=ylim, col=c("red", "grey", "white"))}
    
    Heatmap_signifcance(mat)
  }
  if(listID=="Geschlecht"){
    #combine Matrix
    type=object@Allign_list[!duplicated(object@Allign_list$Geschlecht), "Geschlecht"]
    mat=matrix(NA, nrow(datout), length(type))
    colnames(mat)=type
    for(j in 1:length(type)){
      index=as.character(object@Allign_list[object@Allign_list$Geschlecht==type[j], "ID_2"])
      inp=as.numeric(na.omit(datout[index ,1]))
      print(paste("We detected", length(inp), "voxel"))
      #remove outlier
      remove_outliers <- function(x, na.rm = TRUE, ...) {
        qnt <- quantile(x, probs=c(.05, .95), na.rm = na.rm, ...)
        H <- 1.5 * IQR(x, na.rm = na.rm)
        y <- x
        y[x < (qnt[1] - H)] <- NA
        y[x > (qnt[2] + H)] <- NA
        y
      }
      inp1=remove_outliers(inp)
      print(paste("We removed", length(inp1[is.na(inp1)]), "voxel"))
      inp=inp1
      mat[, j]=c(inp,rep(NA, nrow(mat)-length(inp)))
      
      
    }
    
    
    
    
    print(main)
    #par(mfrow=c(1,2))
    if(plot=="Bar"){Barplot_SC(mat,ylim=ylim,main=main)}
    if(plot=="V"){vioplotDHH(mat,ylim=c(0, max(mat, na.rm=T)+0.5), col=c("red", "grey", "white"),main=main)}
    
    Heatmap_signifcance(mat)
  }
  
  return(object)
  
})

setGeneric("TSNE", function(object,basedat="M", number_of_k=4,RemoveOutlies=T, Region=NA,perplexity=30, colorSchema=c("red", "green", "blue", "purple", "orange", "brown", "black")) standardGeneric("TSNE"))
setMethod("TSNE",signature = "SPORT", definition = function(object,basedat="M", number_of_k=4,RemoveOutlies=T,Region=NA,perplexity=30,colorSchema=c("red", "green", "blue", "purple", "orange", "brown", "black")){
  
  Met=object@list_Metabolite
  ID_N=object@list_id
  alligned=object@Allign_list
  datout=as.matrix(do.call(rbind, lapply(1:length(Met), function(i){
    Met_1=as.data.frame(Met[[i]])
    return(Met_1)
  })))
  
  dim(datout)
  
  ID_out=as.data.frame(do.call(rbind, lapply(1:length(ID_N), function(i){
    ID_1=as.data.frame(ID_N[[i]])
    ID_1$col="white"
    if(nrow(ID_1[ID_1$ID=="tumor", ])>0 ){ID_1[ID_1$ID=="tumor", ]$col=colorSchema[1]}
    if(nrow(ID_1[ID_1$ID=="nam", ])>0 ){ID_1[ID_1$ID=="nam", ]$col=colorSchema[2]}
    if(nrow(ID_1[ID_1$ID=="flair_hyperintense", ])>0 ){ID_1[ID_1$ID=="flair_hyperintense", ]$col=colorSchema[3]}
    if(nrow(ID_1[ID_1$ID=="vent", ])>0 ){ID_1[ID_1$ID=="vent", ]$col=colorSchema[4]}
    if(nrow(ID_1[ID_1$ID=="contrast_agent_absorbing_tumor", ])>0 ){ID_1[ID_1$ID=="contrast_agent_absorbing_tumor", ]$col=colorSchema[5]}
    if(nrow(ID_1[ID_1$ID=="mixed", ])>0 ){ID_1[ID_1$ID=="mixed", ]$col=colorSchema[6]}
    if(nrow(ID_1[ID_1$ID=="necrosis", ])>0 ){ID_1[ID_1$ID=="necrosis", ]$col=colorSchema[7]}
    ID_1$Tumor=alligned[alligned$dataset==i, "Tumor"]
    ID_1$Sex=alligned[alligned$dataset==i, "Geschlecht"]
    ID_1$WHO=alligned[alligned$dataset==i, "WHO"]
    
    
    return(ID_1)
  })))
  
  
  
  ID_out$col_T="green"
  ID_out[ID_out$Tumor=="IDH_WT", ]$col_T="red"
  ID_out[ID_out$Tumor=="OGD", ]$col_T="blue"
  ID_out[ID_out$Tumor=="IDH_MUT", ]$col_T="purple"
  ID_out[ID_out$Tumor=="Lesion", ]$col_T="orange"
  ID_out[ID_out$ID=="nam", ]$col_T="green"
  
  ID_out$col_WHO="grey"
  ID_out[ID_out$WHO=="", ]$col_WHO="grey"
  ID_out[ID_out$WHO=="I", ]$col_WHO="black"
  ID_out[ID_out$WHO=="II", ]$col_WHO="yellow"
  ID_out[ID_out$WHO=="III", ]$col_WHO="orange"
  ID_out[ID_out$WHO=="IV", ]$col_WHO="red"
  ID_out[ID_out$ID=="nam", ]$col_WHO="green"
  
  ID_out$col_Sex="red"
  ID_out[ID_out$Sex=="m", ]$col_Sex="blue"
  
  print(ID_out[!duplicated(ID_out$ID), "ID"])
  
  pheno_dat=object@Allign_list
  pheno_dat$rownames=paste(pheno_dat$PIZ,"_",pheno_dat$voxel_name, sep="")
  nrow(pheno_dat)
  rownames(pheno_dat)=pheno_dat$rownames
  
  if(is.na(Region)[1]==T){print("No Region Specified")}else{
    index_R=ID_out$ID %in% Region
    datout=datout[index_R, ]
    ID_out=ID_out[index_R, ]
    pheno_dat=pheno_dat[index_R, ]
  }
  nrow(pheno_dat)
  print(table(ID_out[,2:3]))
  
  
  object@TSNE_ID=ID_out
  
  
  
  
  
  
  me=as.matrix(t(datout))
  dim(me)
  colnames(me)=pheno_dat$rownames
  
  dim(me)
  dim(pheno_dat)
  
  object@TSNE_input=t(me)
  
  
  library(AutoPipe)
  library(org.Hs.eg.db)
  
  
  if(RemoveOutlies==T){
    #tSNE MAP
    
    File_genes=Groups_Sup(me, me=me, number_of_k=number_of_k,TRw=c(0))
    
    groups_men=File_genes[[2]]
    me_x=File_genes[[1]]
    
    
    set.seed(99999)
    library(Rtsne)
    ana=Rtsne(t(me_x), check_duplicates = F,perplexity=perplexity, max_iter=500)
    
    data_plot=data.frame(ana$Y)
    rownames(data_plot)=colnames(me_x)
    object@phenodata=pheno_dat[rownames(data_plot), ]
    
    #plot(data_plot, bty="n", pch=".", col=ID_out$col)
    object@GroupMen=groups_men
    
  }else{
    #tSNE MAP
    #res<-AutoPipe::TopPAM(me,max_clusters = 12, TOP=34)
    #plot(res[[2]])
    number_of_k=number_of_k
    
    File_genes=Groups_Sup(me, me=me, number_of_k=number_of_k,TRw=c(-10))
    
    me=File_genes[[1]]
    groups_men=File_genes[[2]]
    
    
    #Sort for subcluster
    #C1
    subC=as.data.frame(do.call(rbind,lapply(1:5, function(i){
      
      
      SubC=me[, rownames(groups_men[groups_men$cluster==i, ])]
      #res<-AutoPipe::TopPAM(SubC,max_clusters = 4, TOP=34)
      #plot(res[[2]])
      number_of_k=3
      File_genes_SG=Groups_Sup(SubC, me=SubC, number_of_k=number_of_k,TRw=c(0))
      groups_men_SG=File_genes_SG[[2]]
      groups_men_SG$cluster=paste(i,"SG_",groups_men_SG$cluster, sep="")
      return(groups_men_SG)
    })))
    
    dim(subC)
    table(subC$cluster)
    
    meX=me[, rownames(subC)]
    
    set.seed(99999)
    library(Rtsne)
    ana=Rtsne(t(meX), check_duplicates = F,perplexity=perplexity, max_iter=500)
    data_plot=data.frame(ana$Y)
    rownames(data_plot)=colnames(meX)
    
    object@phenodata=pheno_dat[rownames(data_plot), ]
    plot(data_plot, bty="n", pch=".", col=ID_out$col)
    object@GroupMen=subC
  }
  
  
  object@TSNE=data_plot
  
  return(object)
  
})


setGeneric("plotTSNE", function(object,size=1, Select=F,Summary=F, Cluster=F,Regions=F, log=F, col_Tumor=F,col_WHO=F, col_Sex=F, Index_Metabolite=NA,bias=3,cex=0.5,size_box=3,sizex = 1,sizey = 1) standardGeneric("plotTSNE"))
setMethod("plotTSNE",signature = "SPORT", definition = function(object,size=1, Select=F,Summary=F,Cluster=F,Regions=F,log=F,col_Tumor=F,col_WHO=F, col_Sex=F,Index_Metabolite=NA,bias=3,cex=0.5,size_box=3,sizex = 1,sizey = 1){
  
  datout=object@TSNE_input
  ID_out=object@TSNE_ID
  print(paste("Possible Metabolite: ", colnames(datout)))
  
  data_plot= object@TSNE
  pheno=object@phenodata
  pheno$Tumor=as.character(pheno$Tumor)
  pheno$WHO=as.character(pheno$WHO)
  if(Select==F){
    pheno[pheno$ID=="nam", ]$Tumor=c("NAM")
    pheno[pheno$ID=="vent", ]$Tumor=c("NAM")
    pheno[pheno$ID=="nam", ]$WHO=c("NAM")
    pheno[pheno$ID=="vent", ]$WHO=c("NAM")
    
  }
  pheno$Cluster=object@GroupMen[rownames(pheno), ]$cluster
  
  if(nrow(pheno[is.na(pheno$Cluster), ])>1){pheno[is.na(pheno$Cluster), ]$Cluster=c("Outlier")}
  #correct phenodat
  data_plot=data_plot[rownames(pheno), ]
  datout=datout[rownames(pheno),]
  
  if(Summary==T){
    #pdf(paste("Data_Cluster_",Sys.time(),".pdf"),useDingbats=F)
    pheno$Cluster=object@GroupMen[rownames(pheno), ]$cluster
    ausw=table(pheno$Geschlecht,pheno$Cluster)
    colm=map2color(rownames(ausw), brewer.pal(9,"Set1"))
    barplot(ausw, col=colm[[1]], main="Sex")
    legend("bottomright",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
    
    ausw=table(pheno$dataset,pheno$Cluster)
    colm=map2color(as.numeric(rownames(ausw)), colorRampPalette(brewer.pal(9,"Set1"))(nrow(ausw))  )
    barplot(ausw, col=colm,main="Patients")
    
    
    ausw=table(pheno$ID,pheno$Cluster)
    colm=map2color(rownames(ausw), brewer.pal(9,"Set1"))
    barplot(ausw, col=colm[[1]], main="ID")
    legend("bottomright",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
    
    ausw=table(pheno$Tumor,pheno$Cluster)
    colm=map2color(rownames(ausw), brewer.pal(9,"Set1"))
    barplot(ausw, col=colm[[1]],main="Histo")
    legend("bottomright",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
    
    
    ausw=table(pheno$WHO,pheno$Cluster)
    colm=map2color(rownames(ausw), brewer.pal(9,"Set1"))
    barplot(ausw, col=colm[[1]],main="WHO")
    legend("bottomright",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
    
    dev.off()
  }
  
  
  
  par(mar=c(6,6,6,6))
  
  if (Regions==T){
    colm=map2color(pheno$ID, brewer.pal(9,'Set1'))
    plot(data_plot, xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
         bty="n",pch=16,cex=cex, col=colm[[1]], main=Index_Metabolite,
         axes=F, ylab="", xlab="")
    legend("bottomleft",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
    
  }else{
    if (col_Tumor==T){
      colm=map2color(pheno$Tumor, brewer.pal(9,'Set1'))
      plot(data_plot, xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
           bty="n",pch=16,cex=cex, col=colm[[1]], main=Index_Metabolite,
           axes=F, ylab="", xlab="")
      legend("bottomleft",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
      
    }else{
      if (col_WHO==T){
        colm=map2color(pheno$WHO, brewer.pal(9,'Set1'))
        plot(data_plot, xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
             bty="n",pch=16,cex=cex, col=colm[[1]], main=Index_Metabolite,
             axes=F, ylab="", xlab="")
        legend("bottomleft",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
        
      }else{
        if (col_Sex==T){
          colm=map2color(pheno$Geschlecht, brewer.pal(9,'Set1'))
          plot(data_plot, xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
               bty="n",pch=16,cex=cex, col=colm[[1]], main=Index_Metabolite,
               axes=F, ylab="", xlab="")
          legend("bottomleft",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
          
        }else{
          if (is.na(Index_Metabolite)==T){plot(data_plot,
                                               xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),
                                               ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
                                               bty="n",pch=16,cex=cex, col=ID_out$col, xlab="TSNE1", 
                                               ylab="TSNE2")}    
          if(Cluster==T){
            
            
            
            colm=map2color(as.character(pheno$Cluster), brewer.pal(9,'Set1'))
            
            plot(data_plot, xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
                 bty="n",pch=16,cex=cex, col=colm[[1]], main=Index_Metabolite,
                 axes=F, ylab="", xlab="")
            
            legend("bottomleft",legend=colm[[2]]$Terms, col=colm[[2]]$col, pch=19 )
          }
          else{
            mypal <- colorRampPalette(rev(brewer.pal(9,'RdYlBu')),bias=bias)(100)
            
            if(log==T){
              cozz=datout[,Index_Metabolite]
              cozz[cozz==0]=0.00001
              col_M=map2color(log2(cozz),mypal )}else{ col_M=map2color(datout[,Index_Metabolite],mypal )}
            
            
            
            zz=plot(data_plot, xlim=c(min(data_plot$X1)*1.3,max(data_plot$X1)*1.3),ylim=c(min(data_plot$X2)*1.3,max(data_plot$X2)*1.3), 
                    bty="n",pch=16,cex=cex, col=col_M, main=Index_Metabolite,
                    axes=F, ylab="", xlab="")
            
            ##add color legende
            library(plotrix)
            ypoints=min(data_plot$X2)+sizex
            xpoints=min(data_plot$X1)*2+sizey
            color.legend(ypoints,xpoints,c(max(data_plot$X2)/size),xpoints+size_box,
                         legend=c(round(seq(min(datout[,Index_Metabolite]), max(datout[,Index_Metabolite]), length.out = 8), digits=2) ), 
                         rect.col=mypal) 
            
            
            
            
          }
          
        }
      }
    }
  }
  
  
  
  print(nrow(data_plot))
  
  return(object)
  
})

setGeneric("realignment", function(object,PIZ, Gatename="frontal") standardGeneric("realignment"))
setMethod("realignment",signature = "SPORT", definition = function(object, PIZ, Gatename="frontal"){

  ##go into shiny app
  library(shiny)
  library(shinythemes)
  library(ggplot2)
  library(Cairo)   # For nicer ggplot2 output when deployed on Linux
  
  warning("This function will be manipulate the data within your object. It is recommended to save the object before!!!!")
  
  message("Create Gate")
  
  From_PIZ_to_df=function(object, PIZ){
    #get Metetabolote
    
    x=unique(object@phenodata[which(object@phenodata$PIZ==PIZ), ]$dataset)
    
    allignment_data=object@phenodata[which(object@phenodata$PIZ==PIZ), ]
    
    align_ob=object@Allign_list[which(object@Allign_list$PIZ==PIZ), ]
    
    
    
    df_use=data.frame(Voxel=align_ob$voxel_name, Region_Seg=align_ob$ID)
    
    return(list(df_use))
    
  }
  From_PIZ_to_df=From_PIZ_to_df(object, PIZ=PIZ)
  
  #create grid
  x=1:32
  y=32:1
  plot_df=as.data.frame(expand.grid(x,y));names(plot_df)=c("x", "y")
  
  #add features
  plot_df$feat=NA
  plot_df[rownames(plot_df) %in% From_PIZ_to_df[[1]]$Voxel, ]$feat=From_PIZ_to_df[[1]]$Region_Seg
  realign_shiny=function(object, plot_df){
    
    ui <- shiny::fluidPage(
      shiny::titlePanel("Re-Alignment of MRS"),
      
      # Sidebar layout with input and output definitions ----
      shiny::sidebarLayout(
        
        # Sidebar panel for inputs ----
        shiny::sidebarPanel(
          # done button
          shiny::actionButton("choose_toggle", "Choose/unchoose"),
          # clear button
          shiny::actionButton("reset", "Clear"),
          # done button
          shiny::actionButton("done", "Done"),
          
          
          #Inntroduction
          
          shiny::h3("Instructions:"),
          shiny::tags$ol(
            shiny::tags$li("Highlight Spots to change alignment"),
            shiny::tags$li("Click the 'Choose/unchoose' button."),
            shiny::tags$li("Click 'Done'.")
          ),
          
          shiny::h3("Selection:"),
          selectInput(inputId="Selection_feature", 
                      label="Feature:",
                      choices = "Region_Seg",
                      selected= "Region_Seg"
          ),
          selectInput(inputId="Selection_PIZ", 
                      label="PIZ:",
                      choices = c(object@clinical_data$PIZ),
                      selected= c(object@clinical_data$PIZ)[1]
          ),
          
        ),
        
        # Main panel for displaying outputs ----
        shiny::mainPanel(
          shiny::plotOutput("plot1", height = 350,
                            click = "plot1_click",
                            brush = shiny::brushOpts(id = "plot1_brush"))
        )
      )
    )
    
    
    
    
    
    
    
    
    server <- function(input, output, session) {
      
     
      
      
      print(names(plot_df))
      
      
      vals <- shiny::reactiveValues(
        keeprows = rep(TRUE, nrow(plot_df)))
      
      output$plot1 <- shiny::renderPlot({
        keep    <- plot_df[ vals$keeprows, , drop = FALSE]
        exclude <- plot_df[!vals$keeprows, , drop = FALSE]
        
        ggplot(keep, aes(x, y)) +
          geom_point(data = plot_df, aes(x=x+0.5, y = y+0.5, color=feat),
                     size = 6, alpha = 1) +
          geom_point(alpha = .7) +
          geom_point(data = exclude, shape = 21, fill = "red", color = "red") +
          labs(x="X-Grid", y="Y-Grid")+
          theme_void()+
          scale_color_viridis_d(na.value="white")+
          geom_hline(yintercept=1:33, linetype="dashed", color = "black", size=0.5)+
          geom_vline(xintercept=1:33, linetype="dashed", color = "black", size=0.5)
          
        
      })
      
      # Toggle points that are clicked
      shiny::observeEvent(input$plot1_click, {
        res <- shiny::nearPoints(plot_df,
                                 xvar = "x",
                                 yvar = "y",
                                 input$plot1_click,
                                 allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })
      
      # Toggle points that are brushed, when button is clicked
      shiny::observeEvent(input$choose_toggle, {
        res <- shiny::brushedPoints(plot_df, input$plot1_brush,
                                    xvar = "x",
                                    yvar = "y",
                                    allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })
      
      # Reset all points
      shiny::observeEvent(input$reset, {
        vals$keeprows <- rep(TRUE, nrow(plot_df))
      })
      
      
      #Done
      shiny::observeEvent(input$done, {
        shiny::stopApp(vals$keeprows)
      })
      
      
      
    }
    
    
    sel <- shiny::runApp(shiny::shinyApp(ui, server))
    
    return(rownames(plot_df)[!sel])
  }
  
  
  
  new=realign_shiny(object, plot_df)
  plot_df$new_allignment=NA
  plot_df[new, ]$new_allignment=Gatename
  
  message("Output new Gates")
  
  return(list(input_Matrix=plot_df,Gated_cells=new))
  
  

  
})

setGeneric("Gate_Scatter", function(object) standardGeneric("Gate_Scatter"))
setMethod("Gate_Scatter",signature = "SPORT", definition = function(object){
  
  
  ##go into shiny app
  library(shiny)
  library(shinythemes)
  library(ggplot2)
  library(Cairo)   # For nicer ggplot2 output when deployed on Linux
  
  mat_out=data.frame(object@TSNE); names(mat_out)=c("x", "y")
  
  
  Gate_points=function(mat_out){
    ui <- shiny::fluidPage(
      shiny::titlePanel("Choose your Cells of Interest for further analysis"),
      
      # Sidebar layout with input and output definitions ----
      shiny::sidebarLayout(
        
        # Sidebar panel for inputs ----
        shiny::sidebarPanel(
          # done button
          shiny::actionButton("choose_toggle", "Choose/unchoose"),
          # clear button
          shiny::actionButton("reset", "Clear"),
          # done button
          shiny::actionButton("done", "Done"),
          shiny::h3("Instructions:"),
          shiny::tags$ol(
            shiny::tags$li("Highlight Cells from your Input DimRed Plot by ",
                           "clicking and dragging."),
            shiny::tags$li("Click the 'Choose/unchoose' button."),
            
            shiny::tags$li("Click 'Done'.")
          ),
          shiny::h4("Details:"),
          shiny::tags$ul(
            shiny::tags$li(paste("Output of the function is the selected",
                                 "Cell names")),
            shiny::tags$li("To start over, click 'Clear'"),
            shiny::tags$li(paste("You can also choose/unchoose specific Cells",
                                 "by clicking on them directly"))
          )
          
        ),
        
        # Main panel for displaying outputs ----
        shiny::mainPanel(
          shiny::plotOutput("plot1", height = 350,
                            click = "plot1_click",
                            brush = shiny::brushOpts(id = "plot1_brush"))
        )
      )
    )
    
    
    
    
    
    
    
    
    server <- function(input, output, session) {
      
      ica_space_df=mat_out
      vals <- shiny::reactiveValues(
        keeprows = rep(TRUE, nrow(ica_space_df)))
      
      output$plot1 <- shiny::renderPlot({
        keep    <- ica_space_df[ vals$keeprows, , drop = FALSE]
        exclude <- ica_space_df[!vals$keeprows, , drop = FALSE]
        
        ggplot(keep, aes(x, y)) +
          geom_point(data = ica_space_df, aes(x=x, y = y),
                     size = .5, color = "gray", alpha = .3) +
          geom_point(alpha = .7) +
          geom_point(data = exclude, shape = 21, fill = "red", color = "red") +
          labs(x="Component 1", y="Component 2")+
          theme_void()
      }, height = function() {
        session$clientData$output_plot1_width
      })
      
      # Toggle points that are clicked
      shiny::observeEvent(input$plot1_click, {
        res <- shiny::nearPoints(ica_space_df,
                                 xvar = "x",
                                 yvar = "y",
                                 input$plot1_click,
                                 allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })
      
      # Toggle points that are brushed, when button is clicked
      shiny::observeEvent(input$choose_toggle, {
        res <- shiny::brushedPoints(ica_space_df, input$plot1_brush,
                                    xvar = "x",
                                    yvar = "y",
                                    allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })
      
      # Reset all points
      shiny::observeEvent(input$reset, {
        vals$keeprows <- rep(TRUE, nrow(ica_space_df))
      })
      
      shiny::observeEvent(input$done, {
        shiny::stopApp(vals$keeprows)
      })
      
      
      
    }
    sel <- shiny::runApp(shiny::shinyApp(ui, server))
    
    return(rownames(mat_out)[!sel])
  }
  
  
  
  Gate=Gate_points(mat_out)
  
  return(list(input_Matrix=mat_out,Gated_cells=Gate))
  
  
  
})


setGeneric("Create_feat_in_Allignment", function(object, new_feat) standardGeneric("Create_feat_in_Allignment"))
setMethod("Create_feat_in_Allignment",signature = "SPORT", definition = function(object,new_feat){
  newnames=c(names(object@Allign_list), new_feat)
  object@Allign_list$new=NA
  names(object@Allign_list)=newnames
  return(object)
})


setGeneric("Integrate_Gates_to_Allignment", function(object, new_feat, Gate_df, PIZ) standardGeneric("Integrate_Gates_to_Allignment"))
setMethod("Integrate_Gates_to_Allignment",signature = "SPORT", definition = function(object,new_feat, Gate_df, PIZ){
  
  #get_ID2
  #where?
  Gate_df$index=paste0(PIZ,"_",rownames(Gate_df))
  index=intersect(object@Allign_list$ID_2, Gate_df$index)
  Gate_df=Gate_df[Gate_df$index %in% index, ]
  object@Allign_list[object@Allign_list$ID_2 %in% index,  new_feat]=Gate_df$new_allignment
  return(object)
})


setGeneric("get_feature_Bar", function(object, feature) standardGeneric("get_feature_Bar"))
setMethod("get_feature_Bar",signature = "SPORT", definition = function(object, feature){
  
  barplot(table(object@Allign_list[,feature]))
  
})

setGeneric("Prepare_data_to_Analysis", function(object) standardGeneric("Prepare_data_to_Analysis"))
setMethod("Prepare_data_to_Analysis",signature = "SPORT", definition = function(object){
  
  Met=object@list_Metabolite
  ID_N=object@list_id
  alligned=object@Allign_list
  PIZ=alligned %>% dplyr::distinct(PIZ, .keep_all = TRUE) %>% dplyr::select(c(dataset, PIZ))
  
  
  
  if(nrow(PIZ)!=length(Met)) stop("Unequal nr of PIZ and Datasets")
  
  datout=as.matrix(do.call(rbind, lapply(1:length(Met), function(i){
    Met_1=as.data.frame(Met[[i]])
    index=paste0(PIZ[i,2], "_", rownames(Met_1))
    rownames(Met_1)=index
    return(Met_1)
  })))
  datout[1:10, 1:10]
  dim(datout)
  
  object@datout=scale(datout)
  
  index=paste0(alligned$PIZ,"_", alligned$voxel_name)
  ID_out=data.frame(alligned[,!(names(alligned) %in% c("file","voxel_name", "dataset", "ID_2" ))])
  rownames(ID_out)=index
  
  ID_out[!(ID_out$ID %in% c("flair_hyperintense", "tumor", "necrosis", "contrast_agent_absorbing_tumor")),  "Tumor"]="NAM"
  ID_out[!(ID_out$ID %in% c("flair_hyperintense", "tumor", "necrosis", "contrast_agent_absorbing_tumor")),  "WHO"]="NAM"
  
  table(rownames(ID_out) %in% rownames(datout))
  
  object@ID_out=ID_out
  
  message("We correctly prepared a dataset with: ", nrow(ID_out) , " Voxel and ", ncol(ID_out), " Features")
  
  return(object)
})

setGeneric("Get_Nr_PCA", function(object, PCA_max=20) standardGeneric("Get_Nr_PCA"))
setMethod("Get_Nr_PCA",signature = "SPORT", definition = function(object,PCA_max=20){
  
  if(is.null(object@datout)) stop("No Prepared Data")
  
  #Running Ranandom PCA Meethod to choose nr of relevant PCAs
  
  
  #check optimal number
  Real_data=(object@datout)
  Real_data[1:10, 1:10]
  #Create random data
  library(picante)
  Fake_data=randomizeMatrix(Real_data,null.model = "richness",iterations = 1000)
  
  pca_real=pcaMethods::pca(Real_data,method="svd", nPcs=PCA_max)
  pca_fake=pcaMethods::pca(Fake_data,method="svd", nPcs=PCA_max)
  
  
  par(mar=c(5,5,5,5))
  layout(matrix(1:4, 2,2))
  plot(pca_real@R2, type="l", ylim = c(0,0.3), lwd=2, bty="n")
  points(pca_fake@R2, type="l", col="red", lty=2)
  plot(pca_real@scores[,1:2], pch=16, bty="n");plot(pca_real@scores[,2:3],pch=16, bty="n");plot(pca_real@scores[,3:4],pch=16, bty="n")
  
  
  return(object)
  
  })


setGeneric("PCA", function(object, PCA_nr=5) standardGeneric("PCA"))
setMethod("PCA",signature = "SPORT", definition = function(object,PCA_nr=5){
  
  if(is.null(object@datout)) stop("No Prepared Data")
  
  #Running Ranandom PCA Meethod to choose nr of relevant PCAs
  #check optimal number
  Real_data=(object@datout)

  pca_real=pcaMethods::pca(Real_data,method="svd", nPcs=PCA_max)
  

  
  object@PCA=as.data.frame(pca_real@scores)
  
  return(object)
  
})

setGeneric("Dim_Red", function(object, type="TSNE", perplexity=10) standardGeneric("Dim_Red"))
setMethod("Dim_Red",signature = "SPORT", definition = function(object,type="TSNE", perplexity=10){
  
  if(is.null(object@datout)) stop("No Prepared Data")
  if(is.null(object@PCA)) stop("No Prepared PCA")
  
  PCA=object@PCA
  
  if(type=="TSNE"){
    message("running TSNE....")
    
    ts=Rtsne::Rtsne(PCA, perplexity=7)
    #plot(ts$Y)
    out=ts$Y
    rownames(out)=rownames(PCA)
    
    object@DimRed=as.data.frame(out)
    
  }
  
  if(type=="UMAP"){
    
    message("running UMAP....")
    umap=umap::umap(PCA)
    #plot(umap$layout)
    out=umap$layout
    rownames(out)=rownames(PCA)
    object@DimRed=as.data.frame(out)
    
  }
  
  return(object)
  
})

setGeneric("Lov_Cluster", function(object, num.nn=100) standardGeneric("Lov_Cluster"))
setMethod("Lov_Cluster",signature = "SPORT", definition = function(object, num.nn=100){
  
  if(is.null(object@datout)) stop("No Prepared Data")
  if(is.null(object@PCA)) stop("No Prepared PCA")
  

  message("Cluster Data.....")
  
  
  data.use=object@PCA
  library(RANN)
  library(reshape)
  nearest <- nn2(data.use, k = num.nn, treetype = "bd", 
                 searchtype = "priority")
  nearest$nn.idx <- nearest$nn.idx
  nearest$nn.dists <- nearest$nn.dists
  library(igraph)
  edges = melt(t(nearest$nn.idx[, 1:num.nn]))
  colnames(edges) = c("B", "A", "C")
  edges = edges[, c("A", "B", "C")]
  edges$B = edges$C
  edges$C = 1
  edges = unique(transform(edges, A = pmin(A, B), B = pmax(A, B)))
  names(edges) <- c("V1", "V2", "weight")
  edges$V1 <- rownames(data.use)[edges$V1]
  edges$V2 <- rownames(data.use)[edges$V2]
  g <- graph.data.frame(edges, directed = F)
  graph.out = cluster_louvain(g)
  
  clust.assign = factor(graph.out$membership, levels = sort(unique(graph.out$membership)))
  
  Cluster_out=data.frame(ID=rownames(data.use), Cluster=clust.assign)
  rownames(Cluster_out)=Cluster_out$ID
  
  out_ID_out=object@ID_out
  out_ID_out=out_ID_out[Cluster_out$ID,]
  out_ID_out$Cluster=Cluster_out$Cluster 
  rownames(out_ID_out)=Cluster_out$ID
  
  object@ID_out=out_ID_out
  
  
  return(object)
  
})

setGeneric("plot_DimRed", function(object, Feat="Cluster", Metabolite=NA, limit=NA) standardGeneric("plot_DimRed"))
setMethod("plot_DimRed",signature = "SPORT", definition = function(object, Feat="Cluster", Metabolite=NA, limit=NA){
  
  if(is.null(object@datout)) stop("No Prepared Data")
  if(is.null(object@PCA)) stop("No Prepared PCA")
  library(ggplot2)
  
  message(paste("Possible Features are:", names(object@ID_out),"             ", sep=" "))
  #message(paste("Possible Metabolites are:", colnames(object@datout),"             ", sep=" "))
  
  if(is.na(Metabolite)){
    
    
    index=rownames(object@DimRed)
    df_dim=data.frame(x=object@DimRed[,1],y=object@DimRed[,2], Feat=object@ID_out[index, Feat])
    
    p=ggplot() + 
      geom_point(data=df_dim, aes(x=x, y=y, color=as.factor(Feat)))+ 
      theme_classic()+ 
      scale_color_brewer(Feat,type='qual')+
      ggtitle(paste0("Dimension reduction: ", Feat))+
      theme(
        plot.margin = margin(t = 100, r = 100, b = 100, l = 100, unit = "pt"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black")
      )
  }else{
    
    if(length(Metabolite)==1){
      index=rownames(object@DimRed)
      df_dim=data.frame(x=object@DimRed[,1],y=object@DimRed[,2], Feat=object@datout[index, Metabolite])
      
      p=ggplot() + 
        geom_point(data=df_dim, aes(x=x, y=y, color=(Feat)))+ 
        theme_classic()+ 
        xlab("Dimension 1")+
        ylab("Dimension 2")+
        ggtitle(paste0("Dimension reduction: ", Metabolite))+
        theme(
          plot.margin = margin(t = 100, r = 100, b = 100, l = 100, unit = "pt"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black")
        )
      
      if(is.na(limit)){p=p+scale_color_viridis_c(Metabolite)}else{p=p+scale_color_viridis_c(Metabolite, limits=limit)}
      
      
    }else{
      index=rownames(object@DimRed)
      df_dim=data.frame(x=object@DimRed[,1],y=object@DimRed[,2], Feat=rowMeans(object@datout[index, Metabolite]) )
      
      p=ggplot() + 
        geom_point(data=df_dim, aes(x=x, y=y, color=(Feat)))+ 
        theme_classic()+ 
        xlab("Dimension 1")+
        ylab("Dimension 2")+
        ggtitle(paste0("Dimension reduction: Mean of Matabolite > 1"))+
        theme(
          plot.margin = margin(t = 100, r = 100, b = 100, l = 100, unit = "pt"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black")
        )
      
      if(is.na(limit)){p=p+scale_color_viridis_c("Mean")}else{p=p+scale_color_viridis_c("Mean", limits=limit)}
      
      
    }
    
    
    
    
    
  }
  
  
  #p
  
  return(p)
  
})

setGeneric("Get_Marker_Metabolites", function(object, Feat="Cluster") standardGeneric("Get_Marker_Metabolites"))
setMethod("Get_Marker_Metabolites",signature = "SPORT", definition = function(object, Feat="Cluster"){
  
  if(is.null(object@datout)) stop("No Prepared Data")
  if(is.null(object@PCA)) stop("No Prepared PCA")
  
  
  
  #Groups
  #Run loop for each metabolite 
  
  exp=t(object@datout)
  index=colnames(exp)
  group=object@ID_out[index, Feat]
  
  Met=rownames(exp)
  
  
  
  df_fit=data.frame(val=exp[3, ], group=group)
  cluster=unique(df_fit$group)
  #transform_Data
  for(i in 1:length(cluster)){
    df_fit$X=0
    names(df_fit)[2+i]=paste0("Cluster_",cluster[i])
    df_fit[which(df_fit$group==cluster[i]), paste0("Cluster_",cluster[i])]=1
    
  }
  
  library(pbapply)
  
  DE=as.data.frame(do.call(rbind,pblapply(1:nrow(exp), function(z){
  val=exp[z,]
  out_cluster=as.data.frame(do.call(rbind, lapply(3:ncol(df_fit), function(i){
    model=glm(val~ df_fit[,i])
    out=data.frame(coef(summary(model)))[2, ]
    rownames(out)=names(df_fit)[i]
    out$Metabolite=rownames(exp)[z]
    return(out)
  })))
  out_cluster$Cluster=rownames(out_cluster)
  
  
  return(out_cluster)
  })))
  rownames(DE)=NULL
  
  #adjustp
  DE$FDR=p.adjust(p=DE$Pr...t.., method="BH", nrow(DE))
  
  object@DE=DE
  
  return(object)
  
})




#Getter Functions:

setGeneric("get_meta", function(object) standardGeneric("get_meta"))
setMethod("get_meta",signature = "SPORT", definition = function(object){
  
  object@clinical_data
  
})

setGeneric("get_features", function(object) standardGeneric("get_features"))
setMethod("get_features",signature = "SPORT", definition = function(object){
  
  names(object@ID_out)
  
})



plot_DE_Feat=function(object, new, limit=c(-1,5), size_adp=1){
  #create grid 
  c=length(unique(new$Cluster))
  m=length(unique(new$Metabolite))
  
  met=unique(new$Metabolite)
  cl=unique(new$Cluster)
    
  cor.mat=matrix(0,c,m)
  rownames(cor.mat)=unique(new$Cluster)
  colnames(cor.mat)=unique(new$Metabolite)
  
  for(i in 1:ncol(cor.mat)){
    t=as.data.frame(new[which(new$Metabolite==colnames(cor.mat)[i]), ])
    
    for( z in 1:nrow(t)){
      cor.mat[t$Cluster[z], i]=t$Estimate[z]
    }
    
    }
  
  
  library(ggplot2)
  library(data.table)
  
  # compute corellation and p.values for all combinations of columns
  dt <- CJ(i=seq_len(c), j=seq_len(m))
  dt[, c("corr"):=({
    
    if(nrow(new[which(new$Cluster==cl[i] & new$Metabolite==met[j]), "Estimate"])>0){new[which(new$Cluster==cl[i] & new$Metabolite==met[j]), "Estimate"]}else{0}
      
    
    
    
    }), by=.(i,j)]
  dt[, c("lfdr"):=({
    
    if(nrow(new[which(new$Cluster==cl[i] & new$Metabolite==met[j]), "FDR"])>0){new[which(new$Cluster==cl[i] & new$Metabolite==met[j]), "FDR"]}else{0}
    
    
    
    
  }), by=.(i,j)]
  
  
  
  p=ggplot(dt, aes(x=i,y=j, col=corr,  size=(sqrt(1-lfdr)/size_adp) )) + 
    geom_point()+
    scale_color_viridis_c(direction=1, limits=limit,name="Estimation", na.value="white") +
    scale_color_viridis_c(direction=1, limits=limit,name="Estimation", na.value="white") +
    scale_x_continuous("variable", breaks = seq_len(c), labels = cl, na.value="white") +
    scale_y_continuous("variable", breaks = seq_len(m), labels = met, trans="reverse", na.value="white") +
    coord_fixed() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5),
          panel.background=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
    )

  return(p)

  
  }

Violin_DE=function(object, Metabolite=colnames(object@datout), Feat="Cluster", ylim=10){
  
  metabolite_show=Metabolite
  cluster=unique(object@ID_out[,Feat])
  
  
  index_df=data.frame(BC=rownames(object@ID_out), Cluster=object@ID_out[,Feat])
  
  exp=as.data.frame(object@datout[index_df$BC, ])
  exp$Cluster=index_df$Cluster
  
  plot_df=exp %>% dplyr::select(c(metabolite_show, "Cluster")) 
  plot_df_n=as.data.frame(do.call(rbind, lapply(1:length(metabolite_show), function(i){
    data.frame(Val=plot_df[,metabolite_show[i]],Metabolite_DF=metabolite_show[i], Cluster=plot_df$Cluster)
  })))

  p=ggplot(data=plot_df_n, aes(x=Cluster, y=Val))+
    geom_violin(aes(fill=Cluster))+
    facet_grid(rows=vars(Metabolite_DF), scales = "free_y")+
    theme_classic()+
    ylab("Metabolic Intensity")+
    scale_color_brewer(Feat,type='qual')+
    theme(
      plot.margin = margin(t = 50, r = 100, b = 50, l = 100, unit = "pt"),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black"),
      strip.text.y = element_text(angle = 0, face="italic", size=14),
      strip.placement = "outside",
      strip.background = element_rect(colour="white", fill="white"),
      panel.spacing.y = unit(10, "pt")
      
    )+
    ylim(0,ylim)
    
  
  return(p)
  
  
  }

Test_sig_plot=function(object, Metabolite=colnames(object@datout), Feat="Cluster", size=20){
  
  metabolite_show=Metabolite
  cluster=unique(object@ID_out[,Feat])
  
  
  index_df=data.frame(BC=rownames(object@ID_out), Cluster=object@ID_out[,Feat])
  
  exp=as.data.frame(object@datout[index_df$BC, metabolite_show])
  exp$Cluster=index_df$Cluster
  
  
  sig=exp %>%  group_split(Cluster) 
  pw_comp=matrix(NA,length(cluster), length(cluster));rownames(pw_comp)=colnames(pw_comp)=cluster
  for(i1 in 1:length(sig)){
    for(i2 in 1:length(sig)){
      pw_comp[i1, i2]=wilcox.test(as.numeric(data.frame(sig[[i1]])[,1]), as.numeric(data.frame(sig[[i2]])[,1]))$p.value
    }
  }
  library(reshape2)
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  plot_df=
    get_upper_tri(pw_comp) %>% melt(na.rm = F)
  plot_df$FDR=plot_df$test=p.adjust(plot_df$value, n=nrow(plot_df), method = "BH") 
  
  #plot_df$FDR=p.adjust(plot_df$value, n=nrow(na.omit(plot_df)), method = "BH")
  

  
  
  ggplot(data=plot_df, aes(x=Var1, y=Var2, color=FDR))+
    geom_point(size=size)+
    geom_text(aes(x=Var1, y=Var2, label = round(FDR, digits = 4)), color = "black", size = 4) +
    theme_classic()+ 
    scale_color_viridis_c(direction = -1,na.value="white")+
    theme(axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  
}


Compare_Barplot=function(object,Feat="Cluster", compare_to="Tumor"){
  plot_df=object@ID_out
  plot_df=plot_df[,c(Feat, compare_to)]
  names(plot_df)=c("Cx", "Tx")
  #plot_df$Cx=paste0("Cluster_",plot_df$Cx)
  
  p=ggplot(data=plot_df, aes(x=Tx, y=1, fill=Cx))+
    geom_bar(position="fill", stat="identity")+
    theme_classic()+
    xlab(compare_to)+
    scale_fill_brewer(Feat,type='qual')+
    ylab("Percentage")+
    theme(
        plot.margin = margin(t = 50, r = 100, b = 50, l = 100, unit = "pt"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black", angle = 75, vjust = .5)
        
      )
  
  
  
  
  
  
  
  return(p)
  
  
}

plotViolin_interactive=function(object){
  
  library(ggpubr)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(magrittr)
  library(shiny)
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("Choose Feature and Metabolites for plots"),
    
    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(
      
      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        # Select Feature
        shiny::selectInput(inputId="Feat", 
                           label="Select Sample that was used:",
                           choices = get_features(object),
                           selected= get_features(object)[6],
                           multiple = F
        ),
        
        # select Metabolite
        shiny::selectInput(inputId="Met", 
                           label="Select Sample that was used:",
                           choices = c(colnames(object@datout)),
                           selected= c(colnames(object@datout))[2],
                           multiple = T
        ),
        sliderInput("LIM", "LIM:", 0, 50, 10, step=0.51),
        
        # select Metabolite
        shiny::selectInput(inputId="Met_sig", 
                           label="Select Metabolite for pairwised comparison:",
                           choices = c(colnames(object@datout)),
                           selected= c(colnames(object@datout))[2],
                           multiple = F
        ),
        sliderInput("Size", "Size:", 0, 50, 10, step=0.51),
      ),
      
      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::plotOutput("plot1", height = 500),
        shiny::plotOutput("plot2", height = 300)
      )
    )
  )
  
  
  
  
  
  
  
  
  server <- function(input, output, session) {
    
    output$plot1<-renderPlot({
      
      if(input$Met=="ALL"){Met_s=colnames(object@datout)}else{Met_s=input$Met}
      
      Violin_DE(object, Metabolite=Met_s, Feat=input$Feat, ylim = input$LIM)
      
    })
    output$plot2<-renderPlot({
      
      Test_sig_plot(object, Metabolite=input$Met_sig, Feat=input$Feat, size=input$Size)
      
    })
    
    
  }
  shiny::runApp(shiny::shinyApp(ui, server))
  
  
}

plot_DE_Feat_interactive=function(object){
  
  library(ggpubr)
  library(ggplot2)
  library(dplyr)
  library(tidyverse)
  library(magrittr)
  library(shiny)
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("Choose Feature and Metabolites for plots"),
    
    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(
      
      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        # Select Feature
        shiny::selectInput(inputId="Feat", 
                           label="Select Sample that was used:",
                           choices = get_features(object),
                           selected= get_features(object)[5],
                           multiple = F
        ),
        
        
        sliderInput("Sig", "Significance Filter:", 0.001, 1, 0.05, step=0.01),
        sliderInput("Est", "Estimate Filter:", 0.1, 1, 0.5, step=0.01),
        actionButton("Plot", "Plot"),
        
      ),
      
      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::plotOutput("plot1", height = 500)
      )
    )
  )
  
  
  
  
  
  
  
  
  server <- function(input, output, session) {
    
    observeEvent(input$Plot , output$plot1<-renderPlot({
      
      object_n=Get_Marker_Metabolites(object, Feat=input$Feat)
      new= object_n@DE %>% filter(FDR<input$Sig & Estimate > input$Est) %>% group_by(Cluster)
      plot_DE_Feat(object_n, new, limit=c(1,2), size_adp=1)
      
    }))
    
    
    
  }
  shiny::runApp(shiny::shinyApp(ui, server))
  
  
}



MRS_to_SPATA=function(object){
  
  #Create SPATA object
  
  library(SPATA)
  library(magrittr)
  
  obj=new(Class="spata")
  

# Create Innputs ----------------------------------------------------------
  
  obj@data@norm_exp = t(object@datout)
  
  coordinates <-
    data.frame(barcodes=rownames(object@ID_out)) %>% 
    tidyr::separate(barcodes, sep="_", into=c("sample","RowCounts"), remove=F ) %>% 
    dplyr::left_join(. , object@coordinates, by="RowCounts") %>% 
    dplyr::select(barcodes, sample, X, Y)
  names(coordinates)[3:4]=c("x", "y")
  
  obj@coordinates=coordinates
  obj@fdata= object@ID_out %>% tibble::rownames_to_column("barcodes") %>% dplyr::mutate(sample=coordinates$sample)
  obj@fdata = obj@fdata %>% dplyr::select(barcodes, sample, ID , Geschlecht, Cluster)
  obj@fdata$sample=as.character(obj@fdata$sample)
  obj@coordinates$sample=as.character(obj@fdata$sample)
  
  obj@samples=unique(obj@coordinates$sample)
  
  
 return(obj)
  
  
}




# Split object ------------------------------------------------------------

setGeneric("Split_object", function(object, Regions="tumor") standardGeneric("Split_object"))
setMethod("Split_object",signature = "SPORT", definition = function(object, Regions="Tumor"){
  
  message("The output will be a object containing only a subset of voxel")
  seperate_index=object@ID_out %>% dplyr::filter(ID %in% Regions) %>% rownames()
  object@ID_out=object@ID_out[seperate_index, ]
  object@datout=object@datout[seperate_index, ]
  PCA_max=20
  object=PCA(object, PCA_nr=10)
  object=Lov_Cluster(object, num.nn=100)
  object=Dim_Red(object, type="UMAP")
  
  return(object)
  
})



setGeneric("Get_featureMRSprofile", function(object, Feat="Cluster", path_to_RDS) standardGeneric("Get_featureMRSprofile"))
setMethod("Get_featureMRSprofile",signature = "SPORT", definition = function(object, Feat="Cluster", path_to_RDS){
  
  message("Get the raw data ... take a coffee ")
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(ggplot2)
  
  #Collect data per cluster
  setwd(path_to_RDS)
  cluster_vox=object@ID_out %>% tibble::rownames_to_column("Voxel") 
  for(i in 1:ncol(cluster_vox)){cluster_vox[,i]=as.character(cluster_vox[,i])}
  
  
  get_data=lapply(1:n_distinct(cluster_vox %>%pull(!!sym(Feat))), function(i){
    
    Cluster_i=as.character(unique(cluster_vox %>%pull(!!sym(Feat)))[i])
    message(paste0("Start collect data for Cluster ", Cluster_i , " ... "))
   
    #Create list of PIZ and 
    files <- 
      cluster_vox %>% 
      filter(!!sym(Feat) == Cluster_i) %>% 
      tidyr::separate(Voxel, sep="_", c("PIZ", "Voxel")) %>% 
      dplyr::mutate(File_name=paste0("A_", PIZ, ".RDS")) 
    
    load=files %>% distinct(File_name) %>% pull(File_name)
    
    get_val=as.data.frame(do.call(cbind,lapply(1:length(load), function(ii){
      Spx=readRDS(load[ii])
      cols=as.numeric(files %>% filter(File_name==load[i]) %>% pull(Voxel))
      PIZ=unique(files %>% filter(File_name==load[i]) %>% pull(PIZ))
      rows=rownames(Spx[[2]])
      data=as.data.frame(Spx[[2]][,cols], row.names = rows)
      names(data)=paste0(PIZ, "_",cols)
      
      return(data)
      
    })))
    
    return(get_val)
    
    
  })
  
  message("Finish data ... summarize and create plots")
  plot_list=lapply(1:n_distinct(cluster_vox %>%pull(!!sym(Feat))), function(ix){
    
    plot_df =
      data.frame(ppm=as.numeric(rownames(get_data[[ix]])), 
                 val= rowMeans(get_data[[ix]]), 
                 val_median= apply(get_data[[ix]], 1, function(x) {median(x)} ),
                 sm= apply(get_data[[ix]], 1, function(x) {our=sd(x)/sqrt(length(x))} ),
                 sd= apply(get_data[[ix]], 1, function(x) {sd(x)} ),
                 IQR= apply(get_data[[ix]], 1, function(x) {IQR(x)} )
                 
      )
    
    p=ggplot(data=plot_df, aes(x=(ppm), y=(val_median)))+
      geom_line() +
      geom_ribbon(data=plot_df, aes(ymin = val_median- IQR, ymax = val_median + IQR), fill = "red", alpha=0.2)+
      scale_x_reverse()+
      theme_classic() + 
      xlim(5, 0) +
      ylab("Intensity") + xlab("PPM")
    
    
  })
  names(plot_list)=unique(cluster_vox %>%pull(!!sym(Feat)))
  
  
  
  
  return(plot_list)
  
})


#plots=Get_featureMRSprofile(object, Feat="Cluster", path_to_RDS)

####Done

message("Finish....")













