## @knitr dataSetup
rec.data <- read.csv(file="data/Rec.csv", sep=";") 
dists.list <- emp.dists(rec.data)
fit.data <- do.fitting(dists.list, EVSDT, MPT1HTM2g, UVSDT, MPT2HTM)
dists.list <- fit.data$data
model.param <- fit.data$fit.par
model.comp <- fit.data$model.comp
gsq <- exctractModelCriteria(model.param)
ggplot(data=gsq[order(gsq[, "group"],gsq[, "session"],gsq[, "EVSDT_G.Squared"]-gsq[, "MPT1HTM2g_G.Squared"]),], aes(x=seq_along(code), y=EVSDT_G.Squared-MPT1HTM2g_G.Squared, group=paste(group, session), colour=factor(paste("group:", ifelse(group==1, "'fully encoded',", "'forced guessing',"), "session:", session))))+geom_line()+geom_point()+scale_y_continuous(expression(G^2 (EVSDT) - G^2 (MPT1HTM2g)))+scale_x_continuous("participant * session")+geom_hline(aes(yintercept=0))+scale_color_discrete("", guide=guide_legend(nrow=2))+theme(legend.position="top")

## @knitr AllPlots
plot.vps(dists.list, model.param)