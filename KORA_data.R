require(RODBC)
db = odbcConnect("KORA_remote")

sqlQuery(db, "show databases")
sqlQuery(db, "use KORA")
sqlQuery(db, "show tables")

tmp = sqlQuery(db, "select ltgluk2a from KORA.s4f4")
S4 = sqlQuery(db, "select * from s4f4, s4_bioc where s4_bioc.s4metabo_zz=s4f4.zz_nr_s4_bio")

rst = NULL
for(m in 197:362){
  S4$m = scale(log(S4[,m]))
  model = lm(ltgluk2a~m
             + ltalteru + as.factor(lcsex)
             + ltbmi + ltsysmm + ll_hdla
             , subset = which(S4$lp_diab_who06 == 2|S4$lp_diab_who06==0 &S4$ltnuecht==1)
             , S4)
  rst = rbind(rst,summary(model)$coef[2,])
}
rownames(rst) = colnames(S4)[197:362]
write.csv(rst, file = "2h glucose and metabolites_multivariate model_S4.csv")

