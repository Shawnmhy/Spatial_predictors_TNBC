require(ggplot2)
require(sp)
require(FSA)
## Loading required package: ggplot2
## Warning: package 'ggplot2' was built under R version 3.1.2
#set up the data frame.  The outer ring is defined in an anti-clockwise direction and the
#inner holes clockwise.  I have also labelled the holes.



holePolygon <- function(pos){
  
  
  pos[pos$hole == '0', 5] <- 'FALSE'
  pos[pos$hole == '1', 5] <- 'TRUE'
  
  
  
  closest_to_first_row <- function (xy){
    if (nrow(xy)==2){
      return (2)
    } else {
      return(1 + which.min(rowSums((xy[-1,] - xy[1,])^2)))
    }
  }
  
  bridge <- function(xy1, xy2){
    #which row in xy2 is closeset to mid point of xy1
    i.2 <- which.min(apply((colMeans(apply(xy1, 2,range))-xy2)^2,1,sum))
    i.1 <- which.min(apply((xy2[i.2,]-xy1)^2,1,sum)) #which row in xy1 is closest to this point in xy2
    i.2 <- which.min(apply((xy1[i.1,]-xy2)^2,1,sum)) #and again
    return(as.integer(c(i.1,i.2)))
  }
  
  require(data.table)
  ## Loading required package: data.table
  pos<-as.data.table(pos)
  
  #Add a new field to contain the original position as I need this for
  #drawing the orginal borders (which are correct)
  pos$draw=pos$ring
  
  while(TRUE){
    #summarise the information on each ring
    holes <- pos[hole==TRUE,.(mid.x=mean(range(x)), mid.y=mean(range(y))),by=.(id,ring)]
    holes <- holes[,num_holes:=.N ,by=id]
    if(max(holes[,num_holes])<=1) break  #Exit if only one hole per id remaining
    
    #get the first id that has more than one hole in it
    ii <- which(holes$num_holes>1)[1]
    h.id <- holes$id[ii]                  #which id are we dealing with           
    h1.ring <- holes$ring[ii]             #which hole are we dealing with first
    
    #get the hole ring which is closest to h1.ring  This ensures that the shortest path 
    #between h1.ring and h2.ring does not cross another hole ring.
    h2.ring <- holes[id==h.id][closest_to_first_row(as.matrix(holes[id==h.id, .(mid.x,mid.y)])),ring]
    cat('joining ring', h1.ring, 'to ring', h2.ring, '\n')
    
    #find the best bridging point
    h1.xy <- as.matrix(pos[id==h.id & ring==h1.ring, .(x, y)])             #xy matrix for ring1
    h2.xy <- as.matrix(pos[id==h.id & ring==h2.ring, .(x, y)])             #xy matrix for ring2
    h1.l  <- nrow(h1.xy)                                                   #number of points in ring1
    h2.l  <- nrow(h2.xy)                                                   #number of points in ring2
    h1.draw <- pos[id==h.id & ring==h1.ring, draw]                         #existing values for drawing border
    h2.draw <- pos[id==h.id & ring==h2.ring, draw]
    
    b <- bridge(h1.xy, h2.xy)   #b[1] is the row in h1 and b[2] is the row in h2 to bridge
    
    #reorder h2 values about the bridging point and insert into the bridge point in h1
    new.xy <-   rbind(
      h1.xy[seq(b[1]),]            #h1 points up to the bridge
      ,h2.xy[seq(b[2], h2.l-1),]    #h2 from over the bridge to one before the tail=head
      ,h2.xy[seq(1,b[2]),]          #h2 from the head to the bridge again
      ,h1.xy[seq(b[1], h1.l),]      #h1 from the bridge to the tail
    )
    new.draw <- c( h1.draw[seq(b[1])]            #arrange the 'draw' to line up with the orginal rings
                   ,h2.draw[seq(b[2], h2.l-1)]   #so can jump from one ring to another without drawing
                   ,h2.draw[seq(1,b[2])]         #a border over the jump
                   ,h1.draw[seq(b[1], h1.l)]  
    )
    
    #delete the old values and replace with the new values
    drop.rows <- which(pos$id==h.id & (pos$ring==h1.ring|pos$ring==h2.ring))
    
    #update the pos data frame by dropping the original and adding the new
    pos <- rbind(pos[-drop.rows,]
                 ,data.frame(id=h.id
                             ,ring=h1.ring
                             ,x=new.xy[,1]
                             ,y=new.xy[,2]
                             ,hole=TRUE
                             ,draw=new.draw)
    )
  }
  
  #reorder the pos data.frame according to the new rings (with the holes merged)
  pos<-pos[order(id,ring),]
  
  
  # return
  bdry_tr <- pos
}

Network <- function(cellPos, prox_thresh){
  
  #cellPos <- CellPos[, c('CellXPos', 'CellYPos', 'Phenotype', 'ExprPhenotype','CellID')]
  #prox_thresh <- 60
  r <- tri.mesh(cellPos$CellXPos, cellPos$CellYPos)
  delaunayTri <- tripack::triangles(r)
  
  Num_tri <- nrow(delaunayTri)
  # coord length of triangle
  a <- delaunayTri[,1] # node 1 index
  b <- delaunayTri[,2] # node 2 index
  c <- delaunayTri[,3] # node 3 index
  tri_sidelength <- cbind(sqrt( ( r$x[a] - r$x[b] ) ^ 2 + ( r$y[a] - r$y[b] ) ^ 2 ), sqrt( ( r$x[a] - r$x[c] ) ^ 2 + ( r$y[a] - r$y[c] ) ^ 2 ), sqrt( ( r$x[b] - r$x[c] ) ^ 2 + ( r$y[b] - r$y[c] ) ^ 2 ))
  
  
  # perimeter of triangles
  tri_perimeter <- rowSums(tri_sidelength)
  
  # area of triangles
  s <- 0.5 * tri_perimeter
  tri_area <- sqrt(s * (s - tri_sidelength[, 1]) * (s - tri_sidelength[, 2]) *
                     (s - tri_sidelength[, 3]))
  
  # clustering
  k <- seq_len( r$tlnew - 1 )
  i <- r$tlist[k]
  j <- r$tlist[r$tlptr[k]]
  keep <- i > 0
  
  
  i <- abs(i[keep])
  j <- abs(j[keep])
  
  
  distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
  
  i <- i[distances <= prox_thresh]
  j <- j[distances <= prox_thresh]
  #--------------------------------------------#
  
  Dual_NodeList <- cbind(seq(1, nrow(cellPos)),  cellPos) %>%
    `colnames<-` (c('nodes', 'x', 'y', 'Phenotype', 'ExprPhenotype', 'CellID'))
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2),
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  colnames(Dual_EdgeList) <- c('from', 'to')
  
  
  # return lists
  return(list(Dual_EdgeList, Dual_NodeList))
  
}


Network_val <- function(cellPos, prox_thresh){
  
  #cellPos <- coreData[, c('X', 'Y', 'L2ct.T','CID')]
  #prox_thresh <- 60
  r <- tri.mesh(cellPos$X, cellPos$Y)
  delaunayTri <- tripack::triangles(r)
  
  Num_tri <- nrow(delaunayTri)
  # coord length of triangle
  a <- delaunayTri[,1] # node 1 index
  b <- delaunayTri[,2] # node 2 index
  c <- delaunayTri[,3] # node 3 index
  tri_sidelength <- cbind(sqrt( ( r$x[a] - r$x[b] ) ^ 2 + ( r$y[a] - r$y[b] ) ^ 2 ), sqrt( ( r$x[a] - r$x[c] ) ^ 2 + ( r$y[a] - r$y[c] ) ^ 2 ), sqrt( ( r$x[b] - r$x[c] ) ^ 2 + ( r$y[b] - r$y[c] ) ^ 2 ))
  
  
  # perimeter of triangles
  tri_perimeter <- rowSums(tri_sidelength)
  
  # area of triangles
  s <- 0.5 * tri_perimeter
  tri_area <- sqrt(s * (s - tri_sidelength[, 1]) * (s - tri_sidelength[, 2]) *
                     (s - tri_sidelength[, 3]))
  
  # clustering
  k <- seq_len( r$tlnew - 1 )
  i <- r$tlist[k]
  j <- r$tlist[r$tlptr[k]]
  keep <- i > 0
  
  
  i <- abs(i[keep])
  j <- abs(j[keep])
  
  
  distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
  
  i <- i[distances <= prox_thresh]
  j <- j[distances <= prox_thresh]
  #--------------------------------------------#
  
  Dual_NodeList <- cbind(seq(1, nrow(cellPos)),  cellPos) %>%
    `colnames<-` (c('nodes', 'x', 'y', 'L2ct.T', 'CellID'))
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2),
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  colnames(Dual_EdgeList) <- c('from', 'to')
  
  # remove cells that do not connect any other cells
  
  #Dual_NodeList <- Dual_NodeList[Dual_NodeList$nodes %in% unique(c(
  #  as.vector(Dual_EdgeList$from), 
  #  as.vector(Dual_EdgeList$to))),]
  

  # return lists
  return(list(Dual_EdgeList, Dual_NodeList))
  
}





poisp <- function(pos_Target, nsim, spatial_polys){
  
  #pos_Target <- posTumor
  #nsim <- 500
  n <- nrow(pos_Target)
  #plot(NeighborDat)
  #simpp <- runifpoint(n = n, win = spatial_polys, nsim = nsim)
  
  simpp <- list()
  for(rounds in 1:nsim){
    #rounds <- 1
    # ensure each simulation is different
    set.seed(rounds)
    simpp_element <- spsample(spatial_polys, n = n, type = 'random', iter = 100)@coords %>%
      data.frame() 
    
    simpp <- c(simpp, list(simpp_element))
    names(simpp)[rounds] <- paste('Simulation', rounds)
    #ggplot() +
    #  geom_polygon(aes(Tissue_all[,1], Tissue_all[,2], group = Tissue_all[,3]), fill = NA, color = 'black') +
    #  geom_point(aes(simpp@coords[,1], simpp@coords[,2]))
    #geom_point(aes(simpp$`Simulation 32`$x, simpp$`Simulation 32`$y))
    
  }
  
  # return the result 
  return(simpp)
}


#-----------------------------#
# Find number of interactions #
#-----------------------------#

nnIntrxn <- function(sDF, simDF, radius){ # sDF: the coordiantes for the source point pattern
  
  #sDF <- sPos
  if(length(sDF) == 0 | length(simDF) == 0){
    simInt <- 0
  } else{
    dist <- nn2(simDF, query = sDF, k = nrow(simDF), treetype = 'kd', searchtype = 'radius', radius = radius)
    
    nn.idx <- data.frame(dist$nn.idx)
    
    nn.idx$source <- seq_len(nrow(sDF))
    
    nn.merge <- reshape2::melt(nn.idx, id.vars = 'source')
    
    # clear non-valid rows
    nn.merge <- nn.merge[nn.merge$value != 0, -2]
    
    
    # remove itself
    
    
    nn.merge$diff <- nn.merge$source - nn.merge$value
    nn.merge <- nn.merge[nn.merge$diff != 0, -3]
    
    simInt <- nrow(nn.merge)
  }
  
  
  # replace by real cell type
  return(simInt)
}


#----------- Check if the polygon is CLOCKWISE -------------#
clockwise <- function(x) {
  
  x.coords <- c(x[[1]], x[[1]][1])
  y.coords <- c(x[[2]], x[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 




bivarAnalysis.Kcross <- function(pts.type1, pts.type2, Region, neigbor_thresh, step_thresh){
  
  #Region <- Tissue_all[]
  
  #pts.type1 <- 'posCD163'
  #pts.type2 <- 'posCD8'
  # add id variable
  
  #buildings_list <- split(Region, Region$group)
  # only want lon-lats in the list, not the names
  #buildings_list <- lapply(buildings_list, function(x) rev(x[,1:2, drop = FALSE]))
  
  #ps <- lapply(buildings_list, sp::Polygon)
  
  #p1 <- lapply(seq_along(ps), function(i) Polygons(list(ps[[i]]), 
  #                                                 ID = names(buildings_list)[i]  ))
  
  # create SpatialPolygons object
  #my_spatial_polys <- SpatialPolygons(p1, proj4string = CRS("+proj=longlat +datum=WGS84") )  
  
  
  
  type1 <- get(eval(pts.type1)) %>%
    select('CellXPos', 'CellYPos')
  
  type2 <- get(eval(pts.type2)) %>%
    select('CellXPos', 'CellYPos')
  
  if(nrow(type1)*nrow(type2) != 0){
    
    # read pts dat
    #Region <- Region_HE
    
    type1$attr <- pts.type1
    
    type2$attr <- pts.type2
    
    # create multitype df
    pts_OI <- rbind(type1, type2)
    
    # define the type
    species <- factor(pts_OI$attr)
    
    # create multitype ppp
    #Region <- Region_CK56
    
    
    # check if empty  
    ppp1 <- ppp(type1$CellXPos, type1$CellYPos, owin(poly = Region))
    ppp2 <- ppp(type2$CellXPos, type2$CellYPos, owin(poly = Region))
    
    # prevent NA 
    
    if(is.empty(ppp1) == 'FALSE' & is.empty(ppp2) == 'FALSE'){
      
      
      
      
      multitype_ppp <- ppp(pts_OI$CellXPos, pts_OI$CellYPos, marks = species, owin(poly = Region))
      K.cross <- data.frame(Gcross(multitype_ppp, i = pts.type1, j = pts.type2, r = seq(0,neigbor_thresh,step_thresh), correction = 'border'))
      K.cross_ij <- K.cross[complete.cases(K.cross),]
      K.cross_ij$interval <- ifelse(K.cross_ij$r <= 20, 1, 2)
      K.cross_ij[K.cross_ij$r >= 40, 'interval'] <- 3 
      #plot(Gihc)
      # relocat DF
      
      
      
      
      # calculate the area (positive - negative )  
      #K.cross$km <- K.cross$km - K.cross$theo
      
      #i.to.j.diff.area <- trapz(K.cross$r, K.cross$km) 
      
      
      # j to i
      
      K.cross <- data.frame(Gcross(multitype_ppp, i = pts.type2, j = pts.type1, r = seq(0,neigbor_thresh,step_thresh), correction = 'border'))
      
      K.cross_ji <- K.cross[complete.cases(K.cross),]
      
      K.cross_ji$interval <- ifelse(K.cross_ji$r <= 20, 1, 2)
      K.cross_ji[K.cross_ji$r >= 40, 'interval'] <- 3 
      
      #K.cross$km <- K.cross$km - K.cross$theo
      #K.cross <- K.cross[complete.cases(K.cross),]
      #j.to.i.diff.area <- trapz(K.cross$r, K.cross$km) 
      
      
    }
  }
  return(list(K.cross_ij, K.cross_ji))
}






# area under the curve
bivarAnalysis.Kcross_AUC <- function(pts.type1, pts.type2, Region, neigbor_thresh){
  
  #Region <- Tissue_all[,1:2]
  
  #pts.type1 <- 'posTumor'
  #pts.type2 <- 'posCD163'
  # add id variable
  
  #buildings_list <- split(Region, Region$group)
  # only want lon-lats in the list, not the names
  #buildings_list <- lapply(buildings_list, function(x) rev(x[,1:2, drop = FALSE]))
  
  #ps <- lapply(buildings_list, sp::Polygon)
  
  #p1 <- lapply(seq_along(ps), function(i) Polygons(list(ps[[i]]), 
  #                                                 ID = names(buildings_list)[i]  ))
  
  # create SpatialPolygons object
  #my_spatial_polys <- SpatialPolygons(p1, proj4string = CRS("+proj=longlat +datum=WGS84") )  
  
  
  
  type1 <- get(eval(pts.type1)) %>%
    select('CellXPos', 'CellYPos')
  
  type2 <- get(eval(pts.type2)) %>%
    select('CellXPos', 'CellYPos')
  
  if(nrow(type1) > 10 & nrow(type2) > 10){
    
    # read pts dat
    #Region <- Region_HE
    
    type1$attr <- pts.type1
    
    type2$attr <- pts.type2
    
    # create multitype df
    pts_OI <- rbind(type1, type2)
    
    # define the type
    species <- factor(pts_OI$attr)
    
    # create multitype ppp
    #Region <- Region_CK56
    
    
    # check if empty  
    ppp1 <- ppp(type1$CellXPos, type1$CellYPos, owin(poly = Region))
    ppp2 <- ppp(type2$CellXPos, type2$CellYPos, owin(poly = Region))
    
    # prevent NA 
    
    if(is.empty(ppp1) == 'FALSE' & is.empty(ppp2) == 'FALSE'){
      
      
      
      #neigbor_thresh <- 60
      multitype_ppp <- ppp(pts_OI$CellXPos, pts_OI$CellYPos, marks = species, owin(poly = Region))
      K.cross <- data.frame(Gcross(multitype_ppp, i = pts.type1, j = pts.type2, r = seq(0,neigbor_thresh, 1), correction = 'border'))
      K.cross_ij <- K.cross[complete.cases(K.cross),]
      #plot(Gihc)
      # relocat DF
      
      # ------ to plot the exemplar diagram to show the AUC computation -------#
      #pt1 <- K.cross_ij[, c('r', 'km')] %>%
      #  `colnames<-` (c('r', 'y'))
      
      #pt2 <- K.cross_ij[, c('r', 'theo')] %>%
      #  `colnames<-` (c('r', 'y')) %>%
      #  mutate(r = rev(r), y = rev(y)) 
      
      #shade <- rbind.data.frame(pt1, pt2)
      
      #p <- ggplot() +
      #  theme_bw() +
      #  geom_polygon(data = shade, aes(r, y), fill ='grey50', alpha = 0.5) +
      #  geom_line(data = K.cross_ij, aes(r, theo)) +
      #  geom_line(data = K.cross_ij, aes(r, km), color = 'red') +
      #  theme(axis.title = element_text(size = 22),
      #        axis.text = element_text(size = 20),
      #        panel.grid.major = element_blank(),
      #        panel.grid.minor = element_blank()) +
      #  xlab(expression('r, ' ~ mu~'m')) +
      #  ylab('G(i, j)')
      #ggsave(p, file=paste0("./Figures/M9_Tumor_CD163.jpeg"), width = 8, height = 4, units = "in", dpi = 300)
      
      
      # calculate the area (positive - negative )  
      K.cross_ij$km <- K.cross_ij$km - K.cross_ij$theo
      
      #-------------- Plot the Kcross function -----------#
      
      i.to.j.diff.area <- trapz(K.cross_ij$r, K.cross_ij$km) 
      
      
      # j to i
      
      K.cross <- data.frame(Gcross(multitype_ppp, i = pts.type2, j = pts.type1, r = seq(0,neigbor_thresh,1), correction = 'border'))
      
      K.cross_ji <- K.cross[complete.cases(K.cross),]
      
      K.cross_ji$km <- K.cross_ji$km - K.cross_ji$theo
      j.to.i.diff.area <- trapz(K.cross_ji$r, K.cross_ji$km) 
      
      
    }
  }
  return(list(i.to.j.diff.area, j.to.i.diff.area))
}



#---------- Create a function to compute binned area -------------------#

binArea <- function(binStep, bdry_tr_){
  
  #
  #binStep <- 20
  # create the pixel matrix
  binMat <- expand.grid(seq(1, 2008), seq(1, 2008)) %>%
    `colnames<-` (c('CellXPos', 'CellYPos'))
  
  
  # point in polygon test
  pipTest <- bdry_tr %>%
    group_by(id) %>%
    group_map(~point.in.polygon(binMat$CellXPos, binMat$CellYPos,
                                .x$x, .x$y)) %>%
    data.frame() %>%
    as.matrix() %>%
    matrixStats::rowMaxs() %>%
    replace(. == 0, -1)
  
  
  
  dist2bdry <- nn2(data = bdry_tr_[, c('x', 'y')], query = binMat,
                   k = 1, treetype = 'kd')[['nn.dists']] %>%
    `colnames<-` ('distance') %>%
    data.frame() %>%
    mutate(pipTest = pipTest)
  
  
  binArea <- dist2bdry %>%
    create_bins(cutpoints = seq(from = 0, to = max(dist2bdry$distance), by = binStep)) %>%
    mutate(bin = distend / binStep,
           cutpoint = 0.5 * pipTest * (distbegin + distend) / 2) %>%
    group_by(cutpoint) %>%
    tally() %>%
    `colnames<-` (c('cutpoint', 'freq')) %>%
    mutate(Area = freq / 1000000)
  
    
  return(binArea)
}



#--------------------------------#
# Detect cell type milieu -------#
#--------------------------------#

milieu_detection <- function(cellData, label){
  
  require(ggvoronoi)
  require(ggforce)
  require(concaveman)
  require(sf)
  require(sp)
  #require(gissr)
  require(rgeos)
  #---------- filter the single-cell data file to get the cell data of interest -------#
  #cell_type <- 'Macrophage'

  #cellData <- nucleus_pos_core
  
  #label <- '64'
  milieu_cell_all <-data.frame(matrix(nrow = 0, ncol = 0))
  buffer_area_all <- data.frame(matrix(nrow = 0, ncol = 0))
  patch_area_all <- data.frame(matrix(nrow = 0, ncol = 0))
  
  tryCatch(
    expr = {
  
      cell_oi <- cellData %>%
        dplyr::filter(ExprPhenotype == label | ExprPhenotype == 68) %>%
        select(CellXPos, CellYPos, Phenotype, ExprPhenotype)
      
      # Assign new CellID
      cell_oi$CellID <- seq(1, nrow(cell_oi))
      
      
      # define ROI
      
      
      ggplot(cell_oi,aes(CellXPos, CellYPos, color = as.factor(ExprPhenotype))) +
        theme_void() +
        geom_point() +
      #  geom_voronoi(outline = rect, color = 'black', size = 0.3) +
      #  scale_fill_manual(values = c('0' = '#f8f8f8', '4' = '#82b5e5', '68' = '#82b5e5',
      #                               '64' = '#f8f8f8')) +
        theme(legend.position = 'NA') 
      
      
      # connect all cells of interest
      lists <- Network(cell_oi, 50)
      
      Dual_EdgeList <- lists[[1]]
      Dual_NodeList <- lists[[2]]
      
      
      Dual_EdgeList$from <- Dual_NodeList[Dual_EdgeList$from, 'nodes']
      Dual_EdgeList$to <- Dual_NodeList[Dual_EdgeList$to, 'nodes']
      
      Dual_EdgeList_types <- data.frame(from = cell_oi[match(Dual_EdgeList$from, cell_oi$CellID), 'ExprPhenotype'],
                                        to = cell_oi[match(Dual_EdgeList$to, cell_oi$CellID), 'ExprPhenotype'])
      
      #------------ Construct the igraph object to get connected components -----------#
      
      ig <- graph_from_data_frame(vertices = Dual_NodeList, d = Dual_EdgeList, directed = FALSE)
      
      
      
      milieu_id <- which(components(ig)$csize >= 10 & components(ig)$csize <= 100) # get the milieu id with at least 10 cells
      
      # get the node id associated with each selected milieu
      milieu_nodes <- components(ig)$membership %>%
        data.frame() %>%
        tibble::rownames_to_column("row_names") %>%
        'colnames<-' (c('nodeID', 'membership')) %>%
        dplyr::filter(membership %in% milieu_id)
      
      
    
      
      for(milieu in milieu_id){
        
        #milieu <- 45
        # get the nodeID from the current milieu
        nodeID <- milieu_nodes %>%
          dplyr::filter(membership == milieu) %>%
          select(nodeID) %>%
          as.matrix()
        
        # get the single cell data
        milieu_cell <- cell_oi %>%
          dplyr::filter(CellID %in% nodeID)
        
        
        # get the concave hull for the milieu 
        chull_id <- chull(as.matrix(milieu_cell[, c('CellXPos', 'CellYPos')]))
        
        milieu_concave <- concaveman(as.matrix(milieu_cell[, c('CellXPos', 'CellYPos')])) %>%
          as.data.frame() %>%
          'colnames<-' (c('x', 'y'))
        
        #milieu_convex <- milieu_cell[c(chull_id, chull_id[1]), c('CellXPos', 'CellYPos')] %>%
        #  'colnames<-' (c('x', 'y'))
        
        
        p = st_polygon(list(as.matrix(milieu_concave)))
        pbuf = st_buffer(p, 60)
        op1 <- pbuf[[1]] %>%
          as.data.frame() %>%
          'colnames<-' (c('x', 'y')) %>%
          rev()
        
        op1$id <- milieu
        
        op2 <- milieu_concave[rev(rownames(milieu_concave)),]
        
        op2$id <- milieu
        op2 <- rbind.data.frame(op2, op2[1,])
        
        
        
        # data frame for buffer area only
        
        buffer_area_all <- rbind.data.frame(buffer_area_all, rbind.data.frame(op1, op2))
        
        # data frame for patch area
        patch_area_all <- rbind.data.frame(patch_area_all,cbind.data.frame(milieu_concave, milieu))
        
        # data frame for milieu points
        milieu_cell_all <- rbind.data.frame(milieu_cell_all, cbind.data.frame(milieu_cell, milieu))
      }
      
      #p <- ggplot() +
      #  theme_void() +
      #  ylim(1200, 1600) +
      #  xlim(850, 1200) +
      #  geom_point(data = cellData, aes(CellXPos, CellYPos), color = 'grey', size = 6) +
      #  geom_path(data = patch_area_all[patch_area_all$milieu == 45,], aes(x, y), color = '#28adcf', size = 2) +
      #  geom_point(data = patch_area_all[patch_area_all$milieu == 45,], aes(x, y), color = 'red', size = 6) +
      #  geom_polygon(data = op1, aes(x, y), color = '#1815d5', fill = NA, size = 2)
      #p
      #ggsave(p, file=paste0("./Figures/Milieu_detection/Milieu_J8_45_exemp.jpeg"), width = 5, height = 6, units = "in", dpi = 300)  

      

    },
    error = function(e){ 
      # (Optional)
      # Do this if an error is caught...
      milieu_cell_all <- NA
    }
  )
  
  return(list(buffer_area_all, patch_area_all, milieu_cell_all))
  
}



ratio.sim <- function(ctype1, ctype2, spatial_polys, nsim){
  
  #nsim <- 500
  #IRPvars <- c("PD1.", "PDL1.", 'KI67.', 'EOMES.',  "IL10.", 'ICOS.')
  #ctype1 <-  posTumor
  #ctype2 <-  posFoxP3
  #nsim <- 500
  simpp_type1 <- poisp(ctype1, spatial_polys, nsim = nsim) # randomize cell type 1's location
  
  simpp_type2 <- poisp(ctype2, spatial_polys, nsim = nsim) # randomize cell type 2's location
  
  ratios_type1 <- data.frame(matrix(nrow = 0, ncol = 0 ))
  ratios_type2 <- data.frame(matrix(nrow = 0, ncol = 0 ))
  
  for(sim in seq_len(500)){
    
    #sim <- 4
    
    # current simulated CD4
    sim.type1 <- cbind.data.frame(simpp_type1[[sim]]$x, simpp_type1[[sim]]$y) %>%
      cbind.data.frame(ctype1$ExprPhenotype) %>%
      `colnames<-` (c('X', 'Y', 'ExprPhenotype'))
    

    sim.type2 <- cbind.data.frame(simpp_type2[[sim]]$x, simpp_type2[[sim]]$y) %>%
      cbind.data.frame(ctype2$ExprPhenotype) %>%
      `colnames<-` (c('X', 'Y', 'ExprPhenotype'))
    
    
    ggplot() +
      geom_point(data = sim.type1, aes(X, Y), color = 'red') +
      geom_point(data = sim.type2, aes(X, Y), color = 'green') 
    
    #set.seed(999)
    
    ### Case 1: CD163 to FoxP3
    
    distMatrix <- flexclust::dist2(sim.type1[, c('X', 'Y')], sim.type2[, c('X', 'Y')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.sim.type1.id <- which(apply(distMatrix, 1, FUN = min) < 60)
    
    # Get the single-cell data for these CD4 T cells
    adj.ctype1 <- sim.type1[adj.sim.type1.id, ]
    
    
    if(!(is.empty(adj.sim.type1.id))){
    
      
      ratios <- table(adj.ctype1$ExprPhenotype) %>%
        data.frame() %>%
        `colnames<-` (c('ExprPhenotype', 'count'))
      
      
      
      ratios[ratios$ExprPhenotype == '4', 'count'] <- ratios[ratios$ExprPhenotype == '4', 'count'] #+
        #ratios[ratios$ExprPhenotype == '68', 'count']
      
      ratios[ratios$ExprPhenotype == '64', 'count'] <- ratios[ratios$ExprPhenotype == '64', 'count'] #+
        #ratios[ratios$ExprPhenotype == '68', 'count']
      
      ratios <- ratios %>%
        mutate(ratio = count / nrow(adj.ctype1))
      
      ratios$nsim <- sim
      
      ratios_type1 <- rbind.data.frame(ratios_type1, ratios)
    }
    
    
    ### Case 2: FoxP3 to CD163
    
    distMatrix <- flexclust::dist2(sim.type2[, c('X', 'Y')], sim.type1[, c('X', 'Y')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.ctype2.id <- which(apply(distMatrix, 1, FUN = min) < 60)
    
    # Get the single-cell data for these CD4 T cells
    adj.ctype2 <- sim.type2[adj.ctype2.id, ]
    
    
    if(!(is.empty(adj.ctype2.id))){
      
      ratios <- table(adj.ctype1$ExprPhenotype) %>%
        data.frame() %>%
        `colnames<-` (c('ExprPhenotype', 'count'))
      
      
      
      ratios[ratios$ExprPhenotype == '4', 'count'] <- ratios[ratios$ExprPhenotype == '4', 'count'] #+
        #ratios[ratios$ExprPhenotype == '68', 'count']
      
      ratios[ratios$ExprPhenotype == '64', 'count'] <- ratios[ratios$ExprPhenotype == '64', 'count'] #+
        #ratios[ratios$ExprPhenotype == '68', 'count']
      
      ratios <- ratios %>%
        mutate(ratio = count / nrow(adj.ctype2))
      
      ratios$nsim <- sim
      
      ratios_type2 <- rbind.data.frame(ratios_type2, ratios)
    }
    
  }
  return(list(ratios_type1, ratios_type2))
  
}


#------- Shannon Entropy --------------#

ShannonE <- function(types, coreDat){
  
  #types <- ctype_tensor
  #coreDat <- metabricobj_imObj
  type.count <- length(types)
  
  Total <- nrow(coreDat) # total number of cells
  
  # all int/ext data
  ctype_stat_all <- data.frame(matrix(nrow =  0, ncol = 0))
  
  for (cseq in seq_len(type.count)) {   # outer loop, calcualte interior stats
    #cseq <- 1
    ctype_int <- 0
    ctype_ext <- 0
    p <- 0 # ratio
    
    # current cell type
    ctype <- types[cseq]
    
    # coordinates data for the current core, current cell type
    
    ctype.Dat <- coreDat[coreDat$Phenotype == ctype, c('Location_Center_X', 'Location_Center_Y')]
    
    p <- nrow(ctype.Dat)/Total
    
    # interior score
    ctype_int <- mean(as.matrix(dist(ctype.Dat[,c('Location_Center_X', 'Location_Center_Y')])))
    
    
    # other cell types
    ctype_other <- types[-cseq]
    
    ctype_stat <- cbind(p, ctype_int)
    
    for (cseq_other in ctype_other) {
      
      ctype_other.Dat <- coreDat[coreDat$Phenotype == cseq_other, c('Location_Center_X', 'Location_Center_Y')]
      
      # exterior score
      #test <- flexclust::dist2(ctype.Dat, ctype_other.Dat)
      ctype_ext <- ctype_ext +  mean(flexclust::dist2(ctype.Dat, ctype_other.Dat))
      
    }
    # number of computations = No. cell types - 1
    ctype_stat <- cbind(ctype_stat, ctype_ext / (type.count - 1))
    
    colnames(ctype_stat) <- c('p', 'int', 'ext')
    
    ctype_stat_all <- rbind(ctype_stat_all ,cbind(ctype_stat, ctype))
    
  }
  
  ShannonH <- 0
  # combine row data
  
  
  for (dat in seq_len(type.count)) {
    p <- as.numeric(as.character(ctype_stat_all[dat, 1]))
    d_int <- as.numeric(as.character(ctype_stat_all[dat, 2]))
    
    d_ext <- as.numeric(as.character(ctype_stat_all[dat, 3]))
    d_final <- d_int / d_ext
    
    if(isTRUE(d_int*d_ext == 0)){
      d_final <- 0
    }
    if(isTRUE(p != 0)){
      ShannonH <- -d_final*p*log2(p) + ShannonH
    }
  }
  
  return(ShannonH)
}

#---- modified function for Validation cohort -------#

ShannonE_val <- function(types, coreData){
  
  #types <- cell_types
  #coreDat <- coredata
  type.count <- length(types)
  
  Total <- nrow(coreData) # total number of cells
  
  # all int/ext data
  ctype_stat_all <- data.frame(matrix(nrow =  0, ncol = 0))
  
  for (cseq in seq_len(type.count)) {   # outer loop, calcualte interior stats
    #cseq <- 2
    ctype_int <- 0
    ctype_ext <- 0
    p <- 0 # ratio
    
    # current cell type
    ctype <- types[cseq]
    
    # coordinates data for the current core, current cell type
    #ctype <- 'CD163'
    # depending on the cell type:
    
    
    ctype.Dat <- coreData %>%
      select(X, Y, L2ct.T, CID) %>%
      dplyr::filter(L2ct.T == 'Tumor') %>%
      select(X, Y, CID)
    
    
    p <- nrow(ctype.Dat)/Total
    
    # interior score
    ctype_int <- mean(as.matrix(dist(ctype.Dat[,c('X', 'Y')])))
    
    
    # other cell types
    ctype_other <- types[-cseq]
    
    ctype_stat <- cbind(p, ctype_int)
    
    for (cseq_other in ctype_other) {
      
      ctype_other.Dat <- coreData[!(coreData$CID %in% ctype.Dat$CID), c('X', 'Y')]
      
      # exterior score
      #test <- flexclust::dist2(ctype.Dat, ctype_other.Dat)
      ctype_ext <- ctype_ext +  mean(flexclust::dist2(ctype.Dat[, c('X', 'Y')], ctype_other.Dat))
      
    }
    # number of computations = No. cell types - 1
    ctype_stat <- cbind(ctype_stat, ctype_ext / (type.count - 1))
    
    colnames(ctype_stat) <- c('p', 'int', 'ext')
    
    ctype_stat_all <- rbind(ctype_stat_all ,cbind(ctype_stat, ctype))
    
  }
  
  ShannonH <- 0
  # combine row data
  
  
  for (dat in seq_len(type.count)) {
    p <- as.numeric(as.character(ctype_stat_all[dat, 1]))
    d_int <- as.numeric(as.character(ctype_stat_all[dat, 2]))
    
    d_ext <- as.numeric(as.character(ctype_stat_all[dat, 3]))
    d_final <- d_int / d_ext
    
    if(isTRUE(d_int*d_ext == 0)){
      d_final <- 0
    }
    if(isTRUE(p != 0)){
      ShannonH <- -d_final*p*log2(p) + ShannonH
    }
  }
  
  return(ShannonH)
}


neighboRhood <- function(dat_cells, dat_relation, segment, n_perm, p_tresh){
  
  d = prepare_tables(dat_cells, dat_relation)
  
  
  dat_baseline = apply_labels(d[[1]], d[[2]]) %>%
    aggregate_histo()
  
  
  set.seed(12312)
  dat_perm = rbindlist(mclapply(1:n_perm, function(x){
    dat_labels = shuffle_labels(d[[1]])
    apply_labels(dat_labels, d[[2]]) %>%
      aggregate_histo()
  },mc.cores = 1
  ), idcol = 'run') 
  
  
  
  dat_p <- calc_p_vals(dat_baseline, dat_perm, n_perm = n_perm, p_tresh = p_tresh) 
  pmat = dcast(dat_p, 'FirstLabel ~ SecondLabel', value.var = 'sigval', fun.aggregate = sum,
               fill=0, drop=F) %>%
    data.frame()
  
  rname = pmat$FirstLabel
  
  
  pmat = pmat %>%
    select(-c('FirstLabel')) %>%
    data.frame()
  
  row.names(pmat) <- rname
  colnames(pmat) <- rname
  
  cols = rev(brewer.pal(11,'Spectral'))
  cmap = colorRampPalette(cols)
  
  #segment = 1
  pmat <- melt(setDT(pmat, keep.rownames = TRUE), "rn") %>%
    mutate(segment = segment) %>%
    data.frame() %>%
    `colnames<-` (c('col', 'row', 'n', 'segment'))
  
  return(pmat)
}




# 安装并加载sp包


# define a functio to determine if a point wihtin polygon
is_point_in_polygon <- function(polygon_points, point) {
  point.in.polygon(point[1], point[2], polygon_points$x, polygon_points$y) > 0
}

# Fucntion
# Parameters：
#   point: Points to be determined
#   polygons_df: a data frame contains multiple polygons
# Return：
#   TRUE: point locates in at least one polygon
#   FALSE: other cases

is_point_in_polygons <- function(point, polygons_df) {
  # dissect data with ID
  polygons_list <- polygons_df %>% group_by(id) %>% group_split()
  
  # determine if points in multiple polygonss
  in_polygons <- lapply(polygons_list, is_point_in_polygon, point = point)
  in_polygons <- unlist(in_polygons)
  
  # return
  any(in_polygons)
}



filter_vessel <- function(df, area_thresh) {
  # transfrom DataFrame to SpatialPolygonsDataFrame object
  #df <- pos_vessel
  df_split <- split(df, df$id)
  polygons <- lapply(df_split, function(df_group) {
    coords <- df_group[, c("x", "y")]
    p <- Polygon(coords)
    p <- Polygons(list(p), ID = df_group$id[1])
    return(p)
  })
  
  sp_polygons <- SpatialPolygons(polygons)
  areas <- gArea(sp_polygons, byid = TRUE)
  
  # filter potential artifacts
  area_filtered <- which(areas > area_thresh)
  
  #plot(df_test$x, df_test$y)
  
  return(area_filtered)
}


compute_polygons_area <- function(df, area_thresh) {
  # transfrom DataFrame to SpatialPolygonsDataFrame object
  #df <- pos_vessel
  df_split <- split(df, df$id)
  
  polygons <- lapply(df_split, function(df_group) {
    coords <- df_group[, c("x", "y")]
    p <- Polygon(coords)
    p <- Polygons(list(p), ID = df_group$id[1])
    return(p)
  })
  
  sp_polygons <- SpatialPolygons(polygons)
  areas <- gArea(sp_polygons, byid = TRUE)
  
  # filter potential artifacts
  area_filtered <- which(areas > area_thresh)
  
  return(sum(areas[area_filtered]/1000000))
}
