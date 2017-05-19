################# TEST BIAS CORRECTION #################
source("code_biasCorrection.R")

### Test 0
res_1_iris <- test.real.data(dataset="Data/1_iris.csv", target_nx.=c(25,30,20), nTest=10000)
res_2_iono <- test.real.data(dataset="Data/2_ionosphere.csv", target_nx.=c(63,75), nTest=10000)
res_3_segment <- test.real.data(dataset="Data/3_segment.csv", target_nx.=c(120,220,120,220,120,220,120), nTest=10000)
res_4_ohscal <- test.real.data(dataset="Data/4_ohscal.csv", target_nx.=rep(400,10), nTest=10000)
res_5_wave <- test.real.data(dataset="Data/5_wave.csv", target_nx.=c(1092,753,455), nTest=10000)
res_6_chess <- test.real.data(dataset="Data/6_chess.csv", target_nx.=c(669,1027), nTest=10000)

res_boxplot(res_1_iris)
res_boxplot_error(res_1_iris)
res_boxplot(res_2_iono)
res_boxplot_error(res_2_iono)
res_boxplot(res_3_segment)
res_boxplot_error(res_3_segment)
res_boxplot(res_4_ohscal)
res_boxplot_error(res_4_ohscal)
res_boxplot(res_5_wave)
res_boxplot_error(res_5_wave)
res_boxplot(res_6_chess)
res_boxplot_error(res_6_chess)

saveRDS(res_1_iris, "Result_Real/res_1_iris.RDS")
saveRDS(res_2_iono, "Result_Real/res_2_iono.RDS")
saveRDS(res_3_segment, "Result_Real/res_3_segment.RDS")
saveRDS(res_4_ohscal, "Result_Real/res_4_ohscal.RDS")
saveRDS(res_5_wave, "Result_Real/res_5_wave.RDS")
saveRDS(res_6_chess, "Result_Real/res_6_chess.RDS")

res_1_iris <- readRDS("Result_Real/res_1_iris.RDS")
res_2_iono <- readRDS("Result_Real/res_2_iono.RDS")
res_3_segment <- readRDS("Result_Real/res_3_segment.RDS")
res_4_ohscal <- readRDS("Result_Real/res_4_ohscal.RDS")
res_5_wave <- readRDS("Result_Real/res_5_wave.RDS")
res_6_chess <- readRDS("Result_Real/res_6_chess.RDS")

################# TEST DETERMINANT ################# 
source("code_maxDeterminant.R")

### Test 1
det_1_iris <- test.det.var(dataset="Data/1_iris.csv", target_nx.=c(20,20,15), test_nx.=c(25,25,25))
det_2_iono <- test.det.var(dataset="Data/2_ionosphere.csv", target_nx.=c(50,50), test_nx.=c(50,100))
det_3_segment <- test.det.var(dataset="Data/3_segment.csv", target_nx.=c(rep(100,7)), test_nx.=c(100,200,100,200,100,200,100))
det_4_ohscal <- test.det.var(dataset="Data/4_ohscal.csv", target_nx.=rep(400,10), test_nx.=c(rep(100,5), rep(200,5)))
det_5_wave <- test.det.var(dataset="Data/5_wave.csv", target_nx.=rep(300,3), test_nx.=c(300,600,900))
det_6_chess <- test.det.var(dataset="Data/6_chess.csv", target_nx.=c(300,500), test_nx.=c(1000,500))

### Test 2 - Test > Target
det_1_iris <- test.det.var(dataset="Data/1_iris.csv", target_nx.=c(10,10,15), test_nx.=c(25,25,25))
det_2_iono <- test.det.var(dataset="Data/2_ionosphere.csv", target_nx.=c(30,30), test_nx.=c(50,100))
det_3_segment <- test.det.var(dataset="Data/3_segment.csv", target_nx.=c(rep(50,7)), test_nx.=c(100,200,100,200,100,200,100))
det_4_ohscal <- test.det.var(dataset="Data/4_ohscal.csv", target_nx.=rep(200,10), test_nx.=c(rep(100,5), rep(200,5)))
det_5_wave <- test.det.var(dataset="Data/5_wave.csv", target_nx.=rep(200,3), test_nx.=c(300,600,900))
det_6_chess <- test.det.var(dataset="Data/6_chess.csv", target_nx.=c(200,300), test_nx.=c(1000,500))

### Test 3 - Test < Target
det_1_iris <- test.det.var(dataset="Data/1_iris.csv", test_nx.=c(10,10,15), target_nx.=c(25,25,25))
det_2_iono <- test.det.var(dataset="Data/2_ionosphere.csv", test_nx.=c(30,30), target_nx.=c(50,100))
det_3_segment <- test.det.var(dataset="Data/3_segment.csv", test_nx.=c(rep(50,7)), target_nx.=c(100,200,100,200,100,200,100))
det_4_ohscal <- test.det.var(dataset="Data/4_ohscal.csv", test_nx.=rep(200,10), target_nx.=c(rep(100,5), rep(200,5)))
det_5_wave <- test.det.var(dataset="Data/5_wave.csv", test_nx.=rep(200,3), target_nx.=c(300,600,900))
det_6_chess <- test.det.var(dataset="Data/6_chess.csv", test_nx.=c(200,300), target_nx.=c(1000,500))

res_varplot(det_1_iris)
res_varplot(det_2_iono)
res_varplot(det_3_segment)
res_varplot(det_4_ohscal)
res_varplot(det_5_wave)
res_varplot(det_6_chess)

res_cordet(det_1_iris)
res_cordet(det_2_iono)
res_cordet(det_3_segment)
res_cordet(det_4_ohscal)
res_cordet(det_5_wave)
res_cordet(det_6_chess)

saveRDS(det_1_iris, "Result_Det/det_1_iris.RDS")
saveRDS(det_2_iono, "Result_Det/det_2_iono.RDS")
saveRDS(det_3_segment, "Result_Det/det_3_segment.RDS")
saveRDS(det_4_ohscal, "Result_Det/det_4_ohscal.RDS")
saveRDS(det_5_wave, "Result_Det/det_5_wave.RDS")
saveRDS(det_6_chess, "Result_Det/det_6_chess.RDS")

det_1_iris <- readRDS("Result_Det/det_1_iris.RDS")
det_2_iono <- readRDS("Result_Det/det_2_iono.RDS")
det_3_segment <- readRDS("Result_Det/det_3_segment.RDS")
det_4_ohscal <- readRDS("Result_Det/det_4_ohscal.RDS")
det_5_wave <- readRDS("Result_Det/det_5_wave.RDS")
det_6_chess <- readRDS("Result_Det/det_6_chess.RDS")

