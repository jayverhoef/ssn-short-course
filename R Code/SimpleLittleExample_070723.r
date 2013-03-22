library(maptools)
library(lattice)
source("C:\\00NMML\\StreamNetworks\\StreamNetworksCourse\\R Code\\sglm.stream_functions070723.r")

setwd("C:\\00NMML\\StreamNetworks\\StreamNetworksCourse\\datasets\\statdata")

data1 <- read.table("sitedata.csv", header = TRUE, sep = ",")
data1 <- as.data.frame(data1)
N.all <- length(data1[,1])

dist.junc <- read.table("hydrodist.txt",
		header = TRUE, sep = ",")
dist.junc <- as.matrix(dist.junc[1:N.all,2:(N.all+1)])
dist.junc <- dist.junc/1000
dist.hydro <- dist.junc + t(dist.junc)
		
flow.matrix <- read.table("pimatrix.txt",
		header = TRUE, sep = ",")
flow.matrix <- as.matrix(flow.matrix[1:N.all,2:(N.all+1)])

obs.data <- data1[!is.na(data1[,"MTLSCULP"]),]
xy <- obs.data[,c("albersX",  "albersY")]
pred.data <- data1[is.na(data1[,"MTLSCULP"]),]
pxy <- pred.data[,c("albersX",  "albersY")]
plot(pxy, pch = 19, col = "blue", cex = 1.8)
points(xy, col = "green", pch = 19, cex = 1.8)

emp.variogram.flow.connect(data1, response.col = "MTLSCULP",
  flow.matrix, dist.junc, transformation = "sqrt",
  nlag = 4, maxlag = 15,
  nlag.cutoff = 10)

# -------------------------------------------------------------
#            MODEL COMPARISONS
# -------------------------------------------------------------

MTLSCULP.out1 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "LinearSill.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out2 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Spherical.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out3 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Mariah.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out4 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out5 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.tailup", CorModel2 = "LinearSill.taildown",
  EstMeth = "ML")
MTLSCULP.out6 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out7 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Mariah.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out8 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Spherical.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out9 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "LinearSill.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out10 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "LinearSill.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out11 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Spherical.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out12 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Mariah.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out13 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.tailup", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out14 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.tailup", CorModel2 = "LinearSill.taildown",
  EstMeth = "ML")
MTLSCULP.out15 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out16 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Mariah.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out17 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Spherical.taildown", CorModel2 = NULL, EstMeth = "ML")
MTLSCULP.out18 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "LinearSill.taildown", CorModel2 = NULL, EstMeth = "ML")

# compare model-fitting criteria
MTLSCULP.IC.out <- Info.crit.compare(
  list(
    MTLSCULP.out1, MTLSCULP.out2, MTLSCULP.out3, MTLSCULP.out4, MTLSCULP.out5,
    MTLSCULP.out6, MTLSCULP.out7, MTLSCULP.out8, MTLSCULP.out9, MTLSCULP.out10,
    MTLSCULP.out11, MTLSCULP.out12, MTLSCULP.out13, MTLSCULP.out14,
    MTLSCULP.out15, MTLSCULP.out16, MTLSCULP.out17, MTLSCULP.out18
  )
)
MTLSCULP.IC.out

# look at variance components
parsil.up <- MTLSCULP.out5$Covariance.Parameters1[["parsil"]]
parsil.dn <- MTLSCULP.out5$Covariance.Parameters2[["parsil"]]
nugget <- MTLSCULP.out5$nugget
v.comp <- c(parsil.up, parsil.dn, nugget)/(parsil.up + parsil.dn + nugget)
names(v.comp) <- c("tail up", "tail down", "nugget")
barplot(v.comp, main = "Variance Components", col = "blue", cex.axis = 2,
  cex = 2, cex.main = 2)
MTLSCULP.out1$Covariance.Parameters1[["parsil"]]
MTLSCULP.out1$nugget

# fixed effects table, use REML after selecting model
MTLSCULP.out1 <- slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix,  trans.power = .5,
  CorModel1 = "LinearSill.tailup", CorModel2 = NULL, EstMeth = "REML")
MTLSCULP.out1$fixed.effects.estimates
MTLSCULP.out14 <- slm.stream(MTLSCULP ~ perAgri, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix, trans.power = .5,
  CorModel1 = "Exponential.tailup", CorModel2 = "LinearSill.taildown",
  EstMeth = "REML")
MTLSCULP.out14$fixed.effects.estimates

# model diagnostics
resids <- MTLSCULP.out1$fit
boxplot(resids[,"resid.student"], col = "blue",
  main = "Studentized Residuals", cex.main = 2, cex.axis = 2)
qqnorm(resids[,"cv.stndr"], col = "blue",
  main = "Standardized Cross-validation Residuals", cex.main = 1.5,
  pch = 19, cex = 3,
  cex.axis = 1.2, cex.lab = 1.5)
ind <- !is.na(data1[, "MTLSCULP"])
emp.variogram.flow.connect(resids, response.col = "resid",
  flow.matrix[ind, ind], dist.junc[ind, ind],
  nlag = 4, maxlag = 15,
  nlag.cutoff = 10)

# predictions
preds <- MTLSCULP.out1$data
pred.graph.data <- preds[!is.na(preds[,"predictions"]),
  c("predictions", "pred.std.err")]
pred.graph.data <- cbind(pxy, pred.graph.data)
obs.graph.data <- preds[is.na(preds[,"predictions"]),
  "MTLSCULP"]
win.graph()
plot(pred.graph.data[,c("albersX","albersY")], pch = 1,
  cex = (pred.graph.data[,"predictions"]+ 200)/200)
points(xy, pch = 19, cex = (obs.graph.data + 200)/200)
win.graph()
plot(pred.graph.data[,c("albersX","albersY")], pch = 1,
  cex = (pred.graph.data[,"pred.std.err"] + 200)/100)
points(xy, pch = 19, cex = (obs.graph.data + 200)/300)

# block kriging predictions
win.graph()
plot(pred.graph.data[,c("albersX","albersY")], pch = 1,
  cex = (pred.graph.data[,"predictions"]+ 200)/200)
points(xy, pch = 19, cex = (obs.graph.data + 200)/200)
segmentA <- c(240:250)
points(pxy[segmentA,], pch = 19, col = "red")
segmentB <- c(300:310)
points(pxy[segmentB,], pch = 19, col = "green", cex = 2)
slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix,
  CorModel1 = "LinearSill.tailup", CorModel2 = NULL,
  EstMeth = "REML", prediction.indices = segmentA)$block.krige
slm.stream(MTLSCULP ~ 1, data = data1, dist.junc = dist.junc,
  flow.matrix = flow.matrix,
  CorModel1 = "LinearSill.tailup", CorModel2 = NULL,
  EstMeth = "REML", prediction.indices = segmentB)$block.krige

