## input: nmds, factors to group by, ordisurf grid data

ggplot.NMDS.2f<-function(XX, df_2f, COLORS){
        library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2

NMDS<-data.frame(MDS1,MDS2,df_2f)
#print(names(NMDS))


veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(unique(NMDS[, 4]))){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[unique(NMDS[, 4])==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
}


X1<-ggplot() +
geom_point(data = NMDS, aes_string(x="MDS1", y="MDS2", fill = colnames(NMDS)[3], shape = colnames(NMDS)[4]) ,size=3,alpha=0.75) +geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=2, linetype=5)+ scale_shape_manual(values=c(21:25)) +
scale_fill_manual(values=c(COLORS), guide = guide_legend(override.aes = list(shape = 23)))+
theme_bw() +
theme(aspect.ratio=1) +
theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20)) +
theme(legend.title=element_text(size=15),legend.text=element_text(size=15))

X1    
}
