require('INLA')
require('raster')
require('plotly')
require('lattice')
require('grid')
require('gridExtra')
require('raster')
require('plotly')
require('dplyr')
require('tidyr')
require('tidyverse')
require('sp')
require('rgdal')

#setwd("/home/ljochems/PrelimHurdleModel/")

efb_data <- read.csv("FullVeg2019.csv")

plot(efb_data$hyd_bin~efb_data$wtr__)
plot(efb_data$EFB_cover~efb_data$wtr__)
plot(efb_data$hyd_bin~efb_data$NEAR_DIST)
plot(efb_data$EFB_cover~efb_data$NEAR_DIST)
plot(efb_data$hyd_bin~efb_data$typ_cover)
plot(efb_data$EFB_cover~efb_data$typ_cover)
plot(efb_data$hyd_bin~efb_data$MeanFetch)
plot(efb_data$EFB_cover~efb_data$MeanFetch)
#some of these relationships are not linear, nor necessarily quadratic...
#exponential? breakpoint?


### Inputs for the hurdle models
###################################################
######-----mesh building----#####
GLbuffer <- shapefile("Great_Lakes_Buffer.shp")
proj4string(GLbuffer) <- "+proj=longlat +datum=WGS84 +no_defs"
GL_utm <- spTransform(GLbuffer,
                      CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(GL_utm)
points(efb_data[,c(44,45)], col = "red", pch = 16)
#lake polys overlap all sites

#check for replicates
nrow(unique(efb_data[,c(44,45)])) #none

# detect coordinates duplicated
dup <- efb_data[duplicated(efb_data[,c(44,45)]),c(44,45)]
dup

# distance between points
dis <- as.matrix(dist(coordinates(efb_data[, c(44,45)])))
diag(dis) <- NA
dis2 <- dis[lower.tri(dis, diag = FALSE)] #returns lower triangle of matrix distance values 
range(dis2) #very small min distance (transects ~0.1-0.2m apart)
hist(dis2) #in degrees, convert to km but from conversion factor at equator, need to adjust for great lakes region 
cutoff <- 9.558920e-01 

# build the mesh
mesh1 <- inla.mesh.2d(boundary=GL_utm,loc=efb_data[c(44,45)],cutoff=cutoff, max.edge=c(222000,888000),
                      crs = crs(GL_utm)) # may need to convert to m, max.edge = c(2,8)
# plot the mesh and sampled locations
plot(mesh1);points(efb_data[,c(44,45)], col = "red", pch = 16)
#looks decent with more or less all points on vertices 
#do I need a outer boundary? 
summary(mesh1)
print("Number of nodes");mesh1$n

# SPDE model, spatial random effect considered as gaussian field  and GMRF to get cov matrices 
#with many zeros based via Stochastic partial differential equation 
spde <- inla.spde2.matern(mesh1)

# to assist with keeping track of which elements relate to what effects, create an index, represents mesh
s.index <- inla.spde.make.index("s.index_mY", n.spde = spde$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh1, loc = as.matrix(efb_data[,c(44,45)]))
# number of observations x number of nodes
A_dim <- dim(A_matrix)
print("dim of A matrix"); A_dim

# how many points are on nodes: 1s because there is only 1 value diferent from zero
node_pts <- table(apply(A_matrix,1,nnzero))
print("Observations on nodes - within triangles"); node_pts
#all obs on nodes! (all but 1 obs for 2011-21)

# how many columns (nodes) have no points (either on or within the triangles)
nodes_nopts <- table(apply(A_matrix, 2, sum) > 0)
print("Nodes not associated with any observation - associated"); nodes_nopts
# FALSE because 0 > 0 is FALSE

pts_assoc <- table(apply(A_matrix, 2, nnzero))
print("How nodes are associated with observations"); pts_assoc

nzero_wt <- which(apply(A_matrix, 2, sum)>0)[1]
print("The first node with a nonzero weight"); nzero_wt


#####-----create inla stacks-----##### 
stack.PA <- inla.stack(data = list(z = efb_data$hyd_bin),
                       A = list(A_matrix,1),
                       effects = list(s.index_mY = spde$n.spde,
                                      list(b0 = rep(1, nrow(efb_data)),
                                           data.frame(depth=scale(efb_data$wtr__)),data.frame(typha=scale(efb_data$typ_cover)),
                                           data.frame(boats=scale(efb_data$NEAR_DIST), data.frame(fetch=scale(efb_data$MeanFetch))))),
                       tag = "pred")


formula_all <- z ~ 0 + b0 +
  depth + typha + boats + fetch +
  f(s.index_mY, model=spde) 
#consider adding human mod binomial predictor (model='iid')
#as juanmi did in his paper 

# #model to explore spatial error deviations from intercept of response variable
EFB.hurdlemodel.inla <- inla(formula_all,
                             data = inla.stack.data(stack.PA),
                             control.predictor = list(A = inla.stack.A(stack.PA), compute = TRUE),
                             control.compute = list(dic = T, waic = T, config = T,
                                                    hyperpar = T, return.marginals=T), 
                             family = "binomial",
                             control.family = list(list(link = 'logit')), verbose = T)
# control.family = list(beta.censor.value = cens) seems confusing for joint hurdlemodel distribution
#instead just took off a little bit from 1's
#control.family(list(control.link('logit') or loglog cloglog
#cloglog gives same results... trying other links
#control.compute=list(return.marginals.predictor=TRUE)
#control.family(list(control.link('logit')))
#control.inla = list(strategy="gaussian")?
#EFB.hurdlemodel.inla <- summary(EFB.hurdlemodel.inla)
#saveRDS(EFB.hurdlemodel.inla, "PrelimHurdleModel_scale1.rds")


#trying paula morega's ch 9 tutorial to get mappped predictions 
head(efb_data)

library(dplyr)
index <- inla.stack.index(stack = stack.PA , tag = "pred")$data

mod_mean <- EFB.hurdlemodel.inla$summary.fitted.values[index,"0.975quant"]

library(viridis)
library(leaflet)

sps <- SpatialPoints(efb_data[,c("UTM_X","UTM_Y")],
                     proj4string = CRS("+proj=utm +zone=16,17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84"))
efb_data[,c("long","lat")] <- coordinates(spst)


tag.map.title <- htmltools::tags$style(htmltools::HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 28px;
  }
"))

title <- htmltools::tags$div(
  tag.map.title, htmltools::HTML("Statewide EFB Risk Map (2019)")
)  


pal <- colorNumeric("viridis", c(0.1, 0.8), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(lng = efb_data[,49], lat = efb_data[,50],
                   radius = 1,
                   color= "black") %>%
  addCircles(
    lng = efb_data[,49], lat = efb_data[,50],
    color = pal(mod_mean), 
    weight = 5
  ) %>%
  addLegend("bottomright",
            pal = pal, values = mod_mean,
            title = "Probability (%) of EFB Habitat Suitability"
  ) %>%
  addScaleBar(position = c("bottomleft")) %>%
  addControl(title,position = "topright", className="map-title")

title <- htmltools::tags$div(
  tag.map.title, htmltools::HTML("EFB Risk at Burnt Cabin Point Wetland (2019)")
)   

#for inset 
leaflet() %>%
  addProviderTiles('Esri.WorldImagery') %>%
  # addCircleMarkers(lng = efb_data[,49], lat = efb_data[,50],
  #                  radius = 1,
  #                  color= "black") %>%
  addCircles(
    lng = efb_data[,49], lat = efb_data[,50],
    color = pal(mod_mean), 
    weight = 5,
    opacity = 1
  ) %>%
  # addLegend("bottomright",
  #           pal = pal, values = mod_mean,
  #           title = "Probability (%) of EFB Habitat Suitability"
  # ) %>%
  addScaleBar(position = c("bottomleft"),
              options = scaleBarOptions(maxWidth=500)) %>%
  addControl(title,position = "topright", className="map-title")



# #one way to save model object, but need to specify certain findings to troubleshoot
EFB.hurdlemodel.inla.complete <- list(summary.fixed = EFB.hurdlemodel.inla$summary.fixed,
                                      summary.hyperpar = EFB.hurdlemodel.inla$summary.hyperpar,
                                      summary.fitted.values = EFB.hurdlemodel.inla$summary.fitted.values,
                                      summary.random = EFB.hurdlemodel.inla$summary.random,
                                      marginals.fixed = EFB.hurdlemodel.inla$marginals.fixed,
                                      marginals.random = EFB.hurdlemodel.inla$marginals.random,
                                      marginals.hyperpar = EFB.hurdlemodel.inla$marginals.hyperpar,
                                      internal.marginals.hyperpar = EFB.hurdlemodel.inla$internal.marginals.hyperpar,
                                      marginals.spde2.blc = EFB.hurdlemodel.inla$marginals.spde2.blc,
                                      marginals.spde3.blc = EFB.hurdlemodel.inla$marginals.spde3.blc,
                                      internal.marginals.hyperpar = EFB.hurdlemodel.inla$internal.marginals.hyperpar)
save(EFB.hurdlemodel.inla.complete, file="HurdleModel2019.RData")
# #not sure why this isn't saving...
#
# #to obtain range of spatial error terms across the nodes
length(EFB.hurdlemodel.inla$summary.random$s.index._mY$mean)
spatial_error <- range(EFB.hurdlemodel.inla$summary.random$s.index._mY$mean)
print("Range of Spatial Errors"); spatial_error
# #
# # # project the spatial random effect
gproj <- inla.mesh.projector(mesh1)
g.mean <- inla.mesh.project(gproj, EFB.hurdlemodel.inla$summary.random$s.index_mY$mean)
g.sd <- inla.mesh.project(gproj, EFB.hurdlemodel.inla$summary.random$s.index_mY$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))
#
# #not sure what's happening?
#
# # get the spatial parametres of the spatial random effect
spatial.parameters <- inla.spde2.result(inla = EFB.hurdlemodel.inla, name = "s.index_mY",
                                        spde = spde,
                                        do.transform = T)
# # nominal variance (the sill)
sill <- inla.emarginal(function(x) x, spatial.parameters$marginals.variance.nominal[[1]])
sill
# # plot posterior distribution
plot(spatial.parameters$marginals.variance.nominal[[1]],type='l',main='Sill')
#
# # range
range <- inla.emarginal(function(x) x, spatial.parameters$marginals.range.nominal[[1]])
range
# # plot posterior distribution
plot(spatial.parameters$marginals.range.nominal[[1]],type='l',main='Range')
#
# # # nugget, doesn't seem to work
nugget <- inla.emarginal(function(x) 1/x, EFB.hurdlemodel.inla$marginals.hyperpar$`Precision for the Gaussian observations`)
nugget
# # # plot posterior distribution
plot(inla.tmarginal(function(x) 1/x, EFB.hurdlemodel.inla$marginals.hyperpar$`Precision for the Gaussian observations`),
     type='l', main='Nugget')
#

# # plot model residuals
fitted <- EFB.hurdlemodel.inla$"summary.fitted.values"[,1][1:length(efb_data$hyd_bin)]
residuals <- (fitted - efb_data$hyd_bin)
plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')

#need to figure out how to calculate r2 for hurdle model
observed <- efb_data$hyd_bin
meanOb <- mean(efb_data$hyd_bin)
numerator<- sum((fitted - meanOb)^2)
denominator <- sum((observed - meanOb)^2)
r2<- (numerator) / (denominator)
print('R2');r2


