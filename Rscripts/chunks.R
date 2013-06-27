## @knitr Demographics
demographics.df <- read.csv(file="data/Demographics.csv", sep=";")

## @knitr DataSetup
rec.data <- read.csv(file="data/Rec.csv", sep=";")
dists.list <- emp.dists(rec.data)
fit.data <- do.fitting(dists.list, MPT1HTM(), MPT1HTM(TRUE), MPT1HTM(FALSE, TRUE), MPT1HTM(FALSE, TRUE, TRUE), EVSDT(), EVSDT(TRUE), EVSDT(FALSE, TRUE), EVSDT(FALSE, TRUE, TRUE))
dists.list <- fit.data$data
model.param <- fit.data$fit.par
model.res <- fit.data$model.res

## @knitr AllPlots
#plot.vps(dists.list, model.param)

## @knitr LoadExp1
res.ses1 <- model.res[[c(1,1)]]
res.ses2 <- model.res[[c(1,2)]]

## @knitr LoadExp2
res.ses1 <- model.res[[c(2,1)]]
res.ses2 <- model.res[[c(2,2)]]

## @knitr Res01
n.part.ses1 <- nrow(res.ses1[[1]]$data$observed$individual)
n.part.ses2 <- nrow(res.ses2[[1]]$data$observed$individual)
select.ses1 <- select.mpt(res.ses1)
select.ses2 <- select.mpt(res.ses2)
comp1mpt <- gsq.test("MPT1HTM", "MPT1HTM(res.d.scale)", res.ses1)
comp1sdt <- gsq.test("EVSDT", "EVSDT(res.mu.scale)", res.ses1)
comp2mpt <- gsq.test("MPT1HTM(res.d.scale)", "MPT1HTM(res.d.scale,res.map)", res.ses1)
comp2sdt <- gsq.test("EVSDT(res.mu.scale)", "EVSDT(res.mu.scale,res.crit)", res.ses1)
n.part.ses1
n.part.ses2
select.ses1
select.ses2
comp1mpt
comp1sdt
comp2mpt
comp2sdt

## @knitr Res02
select.ses1 <- select.mpt(res.ses1[c("MPT1HTM(res.d.scale)", "EVSDT(res.mu.scale)")], output="full")
select.ses2 <- select.mpt(res.ses2[c("MPT1HTM(res.d.scale)", "EVSDT(res.mu.scale)")], output="full")
diff.ses1 <- g.sq.diff("MPT1HTM(res.d.scale)","EVSDT(res.mu.scale)", res.ses1)
sum.diff.ses1 <- summary(diff.ses1[[1]][1:diff.ses1[[2]][1,3],4])
diff.ses2 <- g.sq.diff("MPT1HTM(res.d.scale)","EVSDT(res.mu.scale)", res.ses2)
sum.diff.ses2 <- summary(diff.ses2[[1]][1:diff.ses2[[2]][1,3],4])
select.ses1
select.ses2
diff.ses1
sum.diff.ses1
diff.ses2
sum.diff.ses2

## @knitr Res02p1
plot.g.sq.diff(diff.ses1)

## @knitr Res02p2
plot.g.sq.diff(diff.ses2)

## @knitr Res03
stab.diff <- g.sq.diff.stability(diff.ses1, diff.ses2)
stab.diff[[2]]["n.stable"]

## @knitr unused
print.p(1)
gsq <- exctractModelCriteria(model.param)
ggplot(data=gsq[order(gsq[, "group"],gsq[, "session"],gsq[, "EVSDT_G.Squared"]-gsq[, "MPT1HTM2g_G.Squared"]),], aes(x=seq_along(code), y=EVSDT_G.Squared-MPT1HTM2g_G.Squared, group=paste(group, session), colour=factor(paste("group:", ifelse(group==1, "'fully encoded',", "'forced guessing',"), "session:", session))))+geom_line()+geom_point()+scale_y_continuous(expression(G[EVSDT]^2 - G[MPT1HTM2g]^2))+scale_x_continuous(expression(participant %.% session))+geom_hline(aes(yintercept=0))+scale_color_discrete("", guide=guide_legend(nrow=2))+theme(legend.position="top")
