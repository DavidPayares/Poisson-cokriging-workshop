
# Setting Up

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
path.folder = "C:/Users/PayaresGarciaDE/Desktop/Workshop"
setwd(path.folder)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Libraries
packages <- c("tigris", "sf" , "gstat", "proxy", "pbapply", "utils", "leaflet", "leafsync", "htmltools")
# Install packages
#install.packages(packages)
# Load packages
lapply(packages, require, character.only = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Functions
source('./PCK-Functions.R')


# Data

## ----message=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
#Data
lung.cancer = read.csv('./Data/Lung_Mortality_Indiana.csv', sep = ";")
diabetes  = read.csv('./Data/Diabetes_Incidence_Indiana.csv', sep = ";")
comorbidity  <-read.csv('./Data/Lung_Diabetes_Comorbidity_Indiana.csv', sep = ";")

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
class(lung.cancer)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
names(lung.cancer)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(lung.cancer)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(diabetes)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(comorbidity)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load counties and reproject
counties = counties(state = 'Indiana') %>%  st_transform(26915)

# Remove unwanted attributes
counties = counties[,c('GEOID')]
counties$FIPS = as.numeric(counties$GEOID)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(counties$geometry)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
lung.diabetes = merge(lung.cancer, diabetes, by = "FIPS")
names(lung.diabetes) = c("FIPS", "County", "Yl", "N", "Yd")
head(lung.diabetes)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Observed risk lung cancer
lung.diabetes$Rl <- lung.diabetes$Yl / lung.diabetes$N * 100000
# Observed risk diabetes
lung.diabetes$Rd <- lung.diabetes$Yd / lung.diabetes$N * 100000
# Data
head(lung.diabetes)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Observed cases comordbidity
lung.diabetes$Yld <- comorbidity$Counts


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Merge with spatial dataframe
ld.sp = merge(counties, lung.diabetes, by = "FIPS")
class(ld.sp)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Color palette
pal.l <- colorNumeric(palette = "magma", domain = lung.diabetes$Rl)
# Labels with information per county
labels.l <- sprintf("<strong>%s</strong><br/>Cases: %s <br/>Risk: %s",
                  ld.sp$County, ld.sp$Yl,  round(ld.sp$Rl, 2)) %>% lapply(htmltools::HTML)

# Leaflet map
map.l <- leaflet(ld.sp %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.l(Rl), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.l, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.l, values = ~Rl, opacity = 0.7, title = "Risk", position = "bottomright")
# Map
map.l


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Color palette
pal.d <- colorNumeric(palette = "magma", domain = lung.diabetes$Rd)
# Labels with information per county
labels.d <- sprintf("<strong>%s</strong><br/>Cases: %s <br/>Risk: %s",
                  ld.sp$County, ld.sp$Yd,  round(ld.sp$Rd, 2)) %>% lapply(htmltools::HTML)

# Leaflet map
map.d <- leaflet(ld.sp %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.d(Rd), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.d, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.d, values = ~Rd, opacity = 0.7, title = "Risk", position = "bottomright")
# Map
map.d


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
leafsync::sync(map.l, map.d)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.l = pck.variogram(data = ld.sp, cases = "Yl", pop = "N", pop.rate = 100000, maxd = 150, nbins = 15)
plot(var.l)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.d = pck.variogram(data = ld.sp, cases = "Yd", pop = "N", pop.rate = 100000, maxd = 150, nbins = 15)
plot(var.d)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
var.ld = pck.crossvariogram(data = ld.sp, cases.a = "Yl", cases.b = "Yd", cases.ab = "Yld", pop = "N", pop.rate = 100000, maxd = 150, nbins = 12)
plot(var.ld)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Inital semivariogram values 
baseline.var =  vgm(5,"Exp", 100, add.to = vgm(20, "Sph", 50))
# LMC fitting
lmc.ld <- lmc.poisson.cokrige(var.l ,var.d , var.ld, ld.sp, cases.a = "Yl", cases.b = "Yd",pop = "N"  , baseline.var)


## ----------------, results='hide'----------------------------------------------------------------------------------------------------------------------------------
# Locations for predicition
prediction.loc <- st_coordinates(st_centroid(ld.sp)) / 1000


## ----------------, results='hide'----------------------------------------------------------------------------------------------------------------------------------
# Perform smoothing
smooth <- poisson.cokrige(data.a = ld.sp, data.b = ld.sp, coords.pred = prediction.loc, cases.a = 'Yl',cases.b = 'Yd', cases.ab = 'Yld', pop = 'N', lmc = lmc.ld, pop.rate = 100000)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(smooth)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Compare smoothed rates against observed rates
plot(smooth$rate.pred, smooth$observed, xlab = 'smoothed risk', ylab = 'observed risk', pch = 20)
abline(lm(smooth$rate.pred ~ smooth$observed), col = "black", lty = 2)
rho.s =round(cor(smooth$rate.pred, smooth$observed), 3)
text(40,60, bquote(rho==.(rho.s)))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Split data-frame
set.seed(257)
sample = sample.int(n = nrow(ld.sp), size = floor(.30 * nrow(ld.sp)), replace = F)
ld.train =  ld.sp[sample,]
ld.test =  ld.sp[-sample,]


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Locations for predicition
prediction.loc <- st_coordinates(st_centroid(ld.test)) / 1000


## ----------------, results='hide'----------------------------------------------------------------------------------------------------------------------------------
# Predict
predictions <- poisson.cokrige(data.a = ld.train, data.b = ld.sp, coords.pred = prediction.loc, cases.a = 'Yl',cases.b = 'Yd', cases.ab = 'Yld', pop = 'N', lmc = lmc.ld, pop.rate = 100000)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(predictions)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Compute observed risk 
observed.rl = ld.test$Yl/ld.test$N * 100000

# Compare predicted risk against observed risk
plot(predictions$rate.pred, observed.rl, xlab = 'predicted risk', ylab = 'observed risk', pch = 20)
abline(lm(predictions$rate.pred ~ observed.rl), col = "black", lty = 2)
rho.p =round(cor(predictions$rate.pred, observed.rl), 3)
text(40,60, bquote(rho==.(rho.p)))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
ld.sp$Smooth.Rl <- smooth$rate.pred
ld.sp$var <- smooth$rate.var


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Color palettes
pal.l <- colorNumeric(palette = "magma", domain = ld.sp$Smooth.Rl)
pal.v <- colorNumeric(palette = "viridis", domain = ld.sp$var)
# Labels with information per county
labels.l <- sprintf("<strong>%s</strong><br/>Cases: %s <br/>Smoothed Risk: %s", ld.sp$County, ld.sp$Yl,  round(ld.sp$Smooth.Rl, 2)) %>% lapply(htmltools::HTML)
labels.v <- sprintf("<strong>%s</strong> <br/>Variance: %s", ld.sp$County, round(ld.sp$var, 2)) %>% lapply(htmltools::HTML)

# Leaflet maps
map.l <- leaflet(ld.sp %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.l(Smooth.Rl), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.l, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.l, values = ~Smooth.Rl, opacity = 0.7, title = "Smooth Risk", position = "bottomright")

map.v <- leaflet(ld.sp %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.v(var), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.v, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.v, values = ~var, opacity = 0.7, title = "Variance", position = "bottomright")

leafsync::sync(map.l, map.v)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
ld.test$risk.est <- predictions$rate.pred
ld.test$var <- predictions$rate.var


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Color palettes
pal.l <- colorNumeric(palette = "magma", domain = ld.test$Rl)
pal.p <- colorNumeric(palette = "magma", domain = ld.test$risk.est)
pal.v <- colorNumeric(palette = "viridis", domain = ld.test$var)


# Labels with information per county
labels.l <- sprintf("<strong>%s</strong><br/>Cases: %s <br/>Smoothed Risk: %s", ld.test$County, ld.test$Yl,  round(ld.test$Rl, 2)) %>% lapply(htmltools::HTML)
labels.p <- sprintf("<strong>%s</strong><br/>Cases: %s <br/>Smoothed Risk: %s", ld.test$County, ld.test$Yl,  round(ld.test$risk.est, 2)) %>% lapply(htmltools::HTML)
labels.v <- sprintf("<strong>%s</strong> <br/>Variance: %s", ld.test$County, round(ld.test$var, 2)) %>% lapply(htmltools::HTML)

# Leaflet maps
map.l <- leaflet(ld.test %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.l(Rl), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.l, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.l, values = ~Rl, opacity = 0.7, title = "Observed Risk", position = "bottomright")

map.p <- leaflet(ld.test %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.p(risk.est), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.p, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.p, values = ~risk.est, opacity = 0.7, title = "Predicted Risk", position = "bottomright")

map.v <- leaflet(ld.test %>% st_transform('+proj=longlat +datum=WGS84')) %>% 
         addTiles() %>% addProviderTiles(providers$CartoDB.Positron) %>%
         addPolygons(color = "grey", weight = 1, fillColor = ~pal.v(var), fillOpacity = 0.7, 
                     highlightOptions = highlightOptions(weight = 4), 
                     label = labels.v, labelOptions = labelOptions(
                       style = list("font-weight" = "normal", padding = "3px 8px"),
                       textsize = "15px", direction = "auto")) %>% 
         addLegend(pal = pal.v, values = ~var, opacity = 0.7, title = "Variance", position = "bottomright")


leafsync::sync(map.l, map.p, map.v, ncol = 3)

