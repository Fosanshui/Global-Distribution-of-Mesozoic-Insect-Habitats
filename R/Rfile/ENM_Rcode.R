# All function  ##
remove_variable <- function(R = R,Coord = Coord, Name = Name, thresholds=0.8, Select = T){
  
  if(Select){
    #get coord climate "df"
    C1 <- NULL
    for(i in 1:nlayers(R)){
      C <- extract(R[[i]],Coord)
      C1 <- cbind(C1,C)
    }
    colnames(C1) <- Name
    df <- na.omit(C1)
  }else{
    df_c <- as.data.frame(R)
    df <- na.omit(df_c)
  }
  
  #cor process
  df1 <- cor(df)
  
  #PCA and get scores
  d_p <- princomp(df1,cor = T)
  M <- data.frame(value = apply(d_p$scores[,1:3], 1, crossprod))
  
  #remove variable
  data1 <- as.data.frame(which(abs(df1) > thresholds, arr.ind=TRUE))
  data2 <- data1[which(data1[,1]!=data1[,2]),]
  N <- data.frame(row = rownames(df1)[data2$row],col = rownames(df1)[data2$col])
  
  Char <- NULL
  for(i in 1:nrow(N)){
    char <- as.character(N[i,][which.min(M[as.character(N[i,]),])])
    Char <- c(Char,char)
  }
  return(unique(Char))
}

get_climate <- function(nc = nc, var= var, var2 = var2, 
                        p_mp = p_mp, p_np = p_np, path_ins = path_ins,
                        npc = 3, Select_point = T, thresholds = 0.8){
  timestamp()
  #model nc
  library(raster)
  modr<-raster(ncols=1800,nrows=900)
  
  clim_S <- raster(modr)
  for(i in 1:length(var)){
    setwd("E://CL/paper/climate")
    R <- rotate(stack(nc, varname = var[i]))
    Re <- overlay(R, fun = mean)
    clim_v <- resample(Re, modr, method = "ngb")
    clim_S <- stack(clim_S, clim_v)
  }
  print(paste0("/!\"process--------1/5--------@(","'_'",")@"))
  
  ma <- as.numeric(stringr::str_remove(nc, pattern = "ensmean_snapshots.nc"))
  
  path <- paste0(p_np,"climate_",ma,".nc")
  nc2 <- stack(path)
  
  for(i in 1:length(var2)){
    r1 <- rotate(stack(nc2[[i]]))
    clim_V <- resample(r1, modr, method = "ngb")
    clim_S <- stack(clim_S, clim_V)
  }
  
  names(clim_S) <- c(var,var2)
  print(paste0("/!\"process--------2/5--------@(","'_'",")@"))
  
  library(maptools)
  path <- paste0(p_mp,ma,paste0("/map",ma,".shp"))
  mp_f <- readShapePoly(path)
  library(dplyr)
  library(terra)
  C_M <- clim_S %>% mask(mp_f,na.rm=TRUE)
  print(paste0("/!\"process--------3/5--------@(","'_'",")@"))
  
  #path_ins <- paste0(p_inp,"/data_t",ma,".csv")
  data <- read.csv(path_ins)
  coord <- data.frame(x=data$Mlng,y=data$Mlat)
  Coord<-unique(coord[which(coord$x<300),])
  
  Remove <- remove_variable(R = C_M,Coord = Coord,Name = names(clim_S), thresholds=thresholds, Select = Select_point)
  
  C_M1 <- C_M[[which(!names(C_M)%in%Remove)]]
  
  library(RStoolbox)
  CM_PCA = rasterPCA(C_M1)
  print(paste0("/!\"process--------4/5--------@(","'_'",")@"))
  
  pc = npc
  modr<-raster(ncols=1800,nrows=900)
  clim_S <- raster(modr)
  for(i in 1:pc){
    R <- rotate(stack(CM_PCA$map[[i]]))
    clim_v <- resample(R, modr, method = "ngb")
    clim_S <- stack(clim_S, clim_v)
  }
  
  print(paste0("/!\"process--------5/5--------@(","'_'",")@"))
  timestamp()
  return(list(clim_S=clim_S,C_M = C_M, C_M1=C_M1, CM_PCA=CM_PCA))
}

distri <- function(filename = name,coord_1 = coord_1, clim2 = clim, models = models, eval = eval, clim1 = clim1, nb.rep = 1,
                   PA.nb.absences = 100, PA.nb.rep = 1, 
                   strategy = 'random', path = "E://CL/paper/model", maxent.path = 'E:\\CL\\paper\\model\\maxent\\maxent.jar'){
  
  
  A<- tryCatch({
    library(dismo)
    bio.ruber<-bioclim(clim2,coord_1)
  }, error = function(e){
    NULL
  })
  
  j=50
  b=1
  while(is.null(A)){
    b=b+1
    coord_1<-rbind(coord_1,data.frame(x=jitter(coord_1$x, j),y=jitter(coord_1$y,j)))
    A<- tryCatch({
      library(dismo)
      bio.ruber<-bioclim(clim2, coord_1)
    }, error = function(e){
      NULL
    })
    if(b>20){
      break
    }
  }
  
  set.seed(2)
  setwd(path)
  library(biomod2)
  myRespName <- filename
  myResp<-rep(1,nrow(coord_1))
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = clim2,
                                       resp.xy = coord_1, resp.name = myRespName,
                                       PA.nb.absences = PA.nb.absences,
                                       PA.nb.rep = PA.nb.rep, PA.strategy = strategy)
  myBiomodOption <- BIOMOD_ModelingOptions( 
    
    RF = list( do.classif = TRUE,
               ntree = 500,
               mtry = 'default',
               sampsize = NULL,
               nodesize = 5,
               maxnodes = NULL),
    
    MAXENT = list( path_to_maxent.jar = maxent.path, 
                   memory_allocated = 500,
                   initial_heap_size = NULL,
                   maximum_heap_size = NULL,
                   background_data_dir = 'default',
                   maximumbackground = 'default',
                   maximumiterations = 5000,
                   visible = FALSE,
                   linear = TRUE,
                   quadratic = TRUE,
                   product = TRUE,
                   threshold = TRUE,
                   hinge = TRUE,
                   lq2lqptthreshold = 80,
                   l2lqthreshold = 10,
                   hingethreshold = 15,
                   beta_threshold = -1,
                   beta_categorical = -1,
                   beta_lqp = -1,
                   beta_hinge = -1,
                   betamultiplier = 1,
                   defaultprevalence = 0.5)
  )
  myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                       models = models,
                                       modeling.id = 'AllModels',
                                       bm.options = myBiomodOption,
                                       metric.eval = eval,
                                       CV.strategy = strategy,
                                       CV.nb.rep = nb.rep,
                                       var.import=3,
                                       CV.perc = 0.6,
                                       seed.val = 12)
  
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    proj.name = 'Current',
                                    new.env = clim1,
                                    models.chosen = 'all')
  a<-raster(myBiomodProj@proj.out@link)
  return(list(out = a,out_models = myBiomodProj,myBiomodModelOut = myBiomodModelOut))
}

Distribution <- function(nc_list = nc_list, p_ins = p_ins, 
                         c.var = c.var, c.var2 = c.var2, c.p_mp = c.p_mp, c.p_np = c.p_np, Select_point = Select_point, thresholds = thresholds, c.npc = c.npc, 
                         rc = rc, nresamp = nresamp,
                         d.model = d.model, d.eval = d.eval, d.PA.nb.absences = d.PA.nb.absences, d.PA.nb.rep = d.PA.nb.rep, d.strategy = d.strategy, d.nb.rep = d.nb.rep, d.path = d.path, d.maxent.path = d.maxent.path,
                         outncpath = outncpath){
  library(raster)
  modr<-raster(ncols=1800, nrows=900)
  mean_dlist <- raster(modr)
  sd_dlist <- raster(modr)
  Cnameslist <- NULL
  dnameslist <- NULL
  Distribution <- NULL
  climatelist <- NULL
  for(k in 1:length(nc_list)){
    ma <- as.numeric(stringr::str_remove(nc_list[k], pattern = "ensmean_snapshots.nc"))
    climname <- paste0("clim_", ma)
    n <- substr(p_ins, 28,28)
    
    coord_th <- paste0(p_ins,"data_",n, ma, ".csv")
    coord <- read.csv(coord_th)
    coord1<-data.frame(x=coord$Mlng, y=coord$Mlat)
    Coord<-unique(coord1[which(coord1$x<300),])
    
    clim <- get_climate(nc = nc_list[k], var = c.var, var2 = c.var2, p_mp = c.p_mp, p_np = c.p_np, path_ins = coord_th, 
                        Select_point = Select_point, thresholds = thresholds, npc = c.npc)
    Cnames <- paste0("Clim_",ma)
    #assign(Cnames, clim)
    Cnameslist[[k]] <- Cnames
    climatelist[[k]] <- clim
    
    Distri <- NULL
    total_D <- raster(modr)
    for(l in 1:length(rc)){
      d_out <- raster(modr)
      d_model <- NULL
      ModelOut <- NULL
      for(m in 1:nresamp){
        set.seed(002)
        this.nrow <- sample(nrow(Coord),max(c(10,round(nrow(Coord)*rc[l]))),replace = T)
        d.filename <- paste0("insect_M_distribution",ma,n,rc[l]*10,"_",m)
        d1 <-tryCatch({distri(filename = d.filename, coord_1 = Coord[this.nrow,], 
                              clim2 = clim[[1]], models = d.model, eval = d.eval, 
                              clim1 = clim[[1]], PA.nb.absences = d.PA.nb.absences, 
                              PA.nb.rep = d.PA.nb.rep, strategy = d.strategy, nb.rep = d.nb.rep,
                              path = d.path, maxent.path = d.maxent.path)}, error = function(e){
                                NULL
                              })
        if(!is.null(d1)){
          d_out <- stack(d_out, d1$out)
          d_model[[m]] <- d1$out_models
          ModelOut[[m]] <- d1$myBiomodModelOut
        }
      }
      if(length(d_out[!is.na(d_out[[1]])])!=0){
        mean_d_out <- overlay(d_out, fun = mean)
        thislist <- list(mean_d_out = mean_d_out, d_out = d_out, d_model = d_model, ModelOut = ModelOut)
      }else{
        thislist <- NA
        mean_d_out <-NA
      }
      if(length(mean_d_out[!is.na(mean_d_out)])!=0){
        total_D <- stack(total_D,mean_d_out)
      }
      Distri[[l]] <- thislist
    }
    names(Distri) <- paste0(rc*100,"%","of","sample")
    
    mean_D <- overlay(total_D, fun = mean)
    sd_D <- overlay(total_D, fun = sd)
    
    dnames <- paste0("Distribution_",ma)
    #assign(dnames, Distri)
    dnameslist[[k]] <- dnames
    Distribution[[k]] <- Distri
    print(k) 
    
    setwd(outncpath)
    r_mean <- mean_D
    r_sd <- sd_D
    NCm <- paste0("mean",dnames,".nc")
    NCs <- paste0('sd',dnames,'.nc')
    writeRaster(r_mean,filename = NCm,overwrite=TRUE)
    writeRaster(r_sd,filename = NCs,overwrite=TRUE)
    
    mean_dlist <- stack(mean_dlist, mean_D)
    sd_dlist <- stack(sd_dlist, sd_D)
  }
  names(climatelist)<-as.character(Cnameslist)
  names(Distribution) <- as.character(dnameslist)
  
  timestamp()
  return(list(mean_dlist = mean_dlist, sd_dlist = sd_dlist,
              dnameslist = dnameslist, Cnameslist = Cnameslist,
              Distribution = Distribution, climatelist = climatelist))
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
setwd("E://CL/paper/climate")
nc_list <- list.files(pattern = ".nc")
#p_ins
rc<-c(0.1,0.3,0.5,0.7,0.9)
nresamp = 3

c.var = c("temp","salt","u","v","w","ustar","vstar","wstar","diff_cbt","hflx","sflx","taux",
          "tauy","convU","psiu","psiv","sff","sff_eta_smooth","eta","hmxl","hblt")
c.var2 = c("bio1","bio12","bio13","bio14","bio15","bio20","bio27","bio28","bio29","bio4","ws")
c.p_mp = "E:/CL/paper/paleomap/"
c.p_np = "E://CL/paper/paleoclimate/Combind/time1/"
Select_point = T
thresholds = 0.8
c.npc = 3

d.model = c("GAM", "CTA",  "ANN", "MAXENT","RF")
#d.model = c("GLM", "GBM", "GAM", "CTA", "ANN", "SRE", "FDA", "MARS", "RF", "MAXENT","MAXNET", "XGBOOST")
d.eval = c("TSS", "ROC","KAPPA","ETS")
d.PA.nb.absences = 100
d.PA.nb.rep = 2
d.strategy = 'random'
d.nb.rep = 2
d.path = "E://CL/paper/model"
d.maxent.path = 'E:\\CL\\paper\\model\\maxent\\maxent.jar'
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

p_ins = "E://CL/paper/insect data_M/total/"
outncpath = "E://CL/paper/out_nc/total"
path = "E://CL/paper/ModelT"
setwd("E://CL/paper/outRdata/total/")

#D <- Distribution(nc_list = nc_list, p_ins = p_ins, 
#                    c.var = c.var, c.var2 = c.var2, c.p_mp = c.p_mp, c.p_np = c.p_np, Select_point = Select_point, thresholds = thresholds, c.npc = c.npc, 
#                    rc = rc, nresamp = nresamp,
#                    d.model = d.model, d.eval = d.eval, d.PA.nb.absences = d.PA.nb.absences, d.PA.nb.rep = d.PA.nb.rep, d.strategy = d.strategy, d.nb.rep = d.nb.rep, d.path = d.path, d.maxent.path = d.maxent.path,
#                    outncpath = outncpath)
#setwd("E://CL/paper/outRdata/total/")
#save(D, file = paste0("E://CL/paper/outRdata/total/",names(D[[5]]),".RData"))
#

for(i in 1:19){
  D <- Distribution(nc_list = nc_list[i], p_ins = p_ins, 
                    c.var = c.var, c.var2 = c.var2, c.p_mp = c.p_mp, c.p_np = c.p_np, Select_point = Select_point, thresholds = thresholds, c.npc = c.npc, 
                    rc = rc, nresamp = nresamp,
                    d.model = d.model, d.eval = d.eval, d.PA.nb.absences = d.PA.nb.absences, d.PA.nb.rep = d.PA.nb.rep, d.strategy = d.strategy, d.nb.rep = d.nb.rep, d.path = d.path, d.maxent.path = d.maxent.path,
                    outncpath = outncpath)
  setwd("E://CL/paper/outRdata/total/")
  save(D, file = paste0("E://CL/paper/outRdata/total/",names(D[[5]]),".RData"))
}

