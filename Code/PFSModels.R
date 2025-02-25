rm(list=ls())
set.seed(3534)
## ichorCNA

RECIST <- read.table("TableS3.txt", header=T, sep="\t")
ichorCNA <- read.table("TableS4.txt", header=T, sep="\t")
NGTAS <- read.table("TableS5.txt", header=T, sep="\t")
CA15.3 <- read.table("TableS6.txt", header=T, sep="\t")
Other <- read.table("TableS7.txt", header=T, sep="\t")
colnames(Other)[4] <- "z.score"
Patients <- read.table("TableS2.txt", header=T, sep="\t")
Patients <- split(Patients, Patients$Public.ID)
Patients <- lapply(Patients, function(x)
    data.frame(Public.ID=x$Public.ID[1], Group=x$Group[1],
               Earliest.Treatment=min(x$Start.Day),
               Latest.Treatment=max(x$End.Day),
               Date.of.brain.mets=x$Day.Brain.Mets[1],
               Date.RIP=x$Day.Death[1]))
Patients <-do.call("rbind", Patients)
Patients$Scans <- NA
Patients$Earliest.Scan <- NA
Patients$Latest.Scan <- NA
Patients$ichorCNA <- NA
Patients$Earliest.ichorCNA <- NA
Patients$Latest.ichorCNA <- NA
Patients$CA15.3 <- NA
Patients$Earliest.CA15.3 <- NA
Patients$Latest.CA15.3 <- NA
for (i in 1:nrow(Patients)) {
    ids <- which(RECIST$Public.ID==Patients$Public.ID[i])
    Patients$Scans[i] <- length(ids)
    Patients$Earliest.Scan[i] <- min(RECIST$Day[ids])
    Patients$Latest.Scan[i] <- max(RECIST$Day[ids])
    ids <- which(ichorCNA$Public.ID==Patients$Public.ID[i])
    Patients$ichorCNA[i] <- length(ids)
    Patients$Earliest.ichorCNA[i] <- min(ichorCNA$Day[ids])
    Patients$Latest.ichorCNA[i] <- max(ichorCNA$Day[ids])
    ids <- which(CA15.3$Public.ID==Patients$Public.ID[i])
    Patients$CA15.3[i] <- length(ids)
    if (length(ids)>0) {
        Patients$Earliest.CA15.3[i] <- min(CA15.3$Day[ids])
        Patients$Latest.CA15.3[i] <- max(CA15.3$Day[ids])
    }

}
Patients$Entry.Date <- pmin(Patients$Earliest.Treatment,
                           Patients$Earliest.Scan,
                           Patients$Earliest.ichorCNA,
                           Patients$Earliest.CA15.3, na.rm=T)
Patients$Exit.Date <- pmax(Patients$Latest.Treatment,
                           Patients$Latest.Scan,
                           Patients$Latest.ichorCNA,
                           Patients$Latest.CA15.3, na.rm=T)

############################################################
############################################################
############################################################
## Comparison of other measures
############################################################
############################################################
############################################################

     panel.hist <- function(x, ...)
     {
         usr <- par("usr")
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col = "olivedrab", ...)
     }

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
     {
         par(usr = c(0, 1, 0, 1))
         r <- cor(x, y)
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         txt <- paste("r=", txt, sep="")
         if(missing(cex.cor)) cex.cor <- 0.3/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * r)
     }
pdf("FigureS2.pdf")
pairs(Other[,c(3:5)], pch=19, lower.panel=panel.smooth,
      upper.panel=panel.cor, diag.panel=panel.hist)
dev.off()


X <- NULL
for (i in 1:nrow(Patients)) {
    tmp1 <- subset(RECIST, Public.ID==Patients$Public.ID[i])
    if (nrow(tmp1)>0) {
        tmp2 <- subset(Other, Public.ID==Patients$Public.ID[i])
        possible.dates <- data.frame(Date=c(tmp1$Day, tmp2$Day,
                                        Patients$Exit.Date[i]),
                                 Event=c(tmp1$Progression,
                                         rep("ichor",
                                             length(tmp2$Day)),
                                         "exit"))
        possible.dates <- possible.dates[order(possible.dates$Date),]
        possible.dates <- possible.dates[which(possible.dates$Event!="NO"),]
        relapse <- which(possible.dates$Event=="YES")
        if (length(relapse)>0) possible.dates <- possible.dates[1:relapse[1],]
        for (j in 1:nrow(possible.dates)) {
        if (possible.dates$Event[j]=="ichor") {
            start <- as.numeric(possible.dates$Date[j] -
                                Patients$Entry.Date[i])
            end <- as.numeric(possible.dates$Date[j+1] - Patients$Entry.Date[i])
            status <- ifelse(possible.dates$Event[j+1]=="YES", 1, 0)
            tmp <- data.frame(Public.ID=Patients$Public.ID[i],
                              start=start, stop=end, status=status,
                              ichorCNA=tmp2$ichorCNA[which(tmp2$Day==possible.dates$Date[j])],
                              z.score=tmp2$z.score[which(tmp2$Day==possible.dates$Date[j])],
                              t.Mad=tmp2$t.Mad[which(tmp2$Day==possible.dates$Date[j])],
                              Date.Start=possible.dates$Date[j],
                              Date.Stop=possible.dates$Date[j+1])
            X <- rbind(X, tmp)
            }
        }
    }
}

library(survival)
X$StoppingRule1 <- 1 * (X$ichorCNA>=7)
m1 <- coxph(Surv(start, stop, status) ~ ichorCNA, data=X)
m2 <- coxph(Surv(start, stop, status) ~ z.score, data=X)
m3 <- coxph(Surv(start, stop, status) ~ t.Mad, data=X)
m4 <- coxph(Surv(start, stop, status) ~ StoppingRule1, data=X)

c.score <- rbind(summary(m1)$concordance,
summary(m2)$concordance,
summary(m3)$concordance,
summary(m4)$concordance)
rownames(c.score) <- c("ichorCNA", "z.score", "t.Mad","ichorCNA>=7")


## save samples for alluvial plot
alluvial.other <- X[,c('Public.ID', 'ichorCNA', 'Date.Start')]
save(alluvial.other, file="Alluvial.Other.RData")

############################################################
############################################################
############################################################
## Progression-free survival
############################################################
############################################################
############################################################


res <- NULL
for (i in 1:nrow(Patients)) {
    tmp1 <- subset(RECIST, Public.ID==Patients$Public.ID[i])
    if (nrow(tmp1)>0) {
        tmp2 <- subset(ichorCNA, Public.ID==Patients$Public.ID[i])
        possible.dates <- data.frame(Date=c(tmp1$Day, tmp2$Day,
                                        Patients$Exit.Date[i]),
                                 Event=c(tmp1$Progression,
                                         rep("ichor",
                                             length(tmp2$Day)),
                                         "exit"))
        possible.dates <- possible.dates[order(possible.dates$Date),]
        possible.dates <- possible.dates[which(possible.dates$Event!="NO"),]
        relapse <- which(possible.dates$Event=="YES")
        if (length(relapse)>0) possible.dates <- possible.dates[1:relapse[1],]
        for (j in 1:nrow(possible.dates)) {
            if (possible.dates$Event[j]=="ichor") {
                start <- as.numeric(possible.dates$Date[j] -
                                    Patients$Entry.Date[i])
                end <- as.numeric(possible.dates$Date[j+1] - Patients$Entry.Date[i])
                status <- ifelse(possible.dates$Event[j+1]=="YES", 1, 0)
                tmp <- data.frame(Public.ID=Patients$Public.ID[i],
                                  start=start, stop=end, status=status,
                                  ichorCNA=tmp2$ichorCNA[which(tmp2$Day==possible.dates$Date[j])],
                                  Date.Start=possible.dates$Date[j],
                                  Date.Stop=possible.dates$Date[j+1])
                res <- rbind(res, tmp)
            }
        }
    }
}


res <- merge(res, Patients[,1:2], all.x=T)
ichorCNA.res <- res
## Estimate the change rate

library(segmented)
library(rms)
dd <- datadist(res)
options(datadist='dd')
o <- NULL
mod <- list()
mod.cont <- list()
for (i in 1:50) {
    sel <- sample(unique(res$Public.ID), 70)
    m <- cph(Surv(start, stop, status) ~ rcs(ichorCNA), data=subset(res, Public.ID %in% sel),
             iter.max = 500000)
    preds <- Predict(m, ichorCNA=NA)
    lin.mod <- lm(yhat ~ ichorCNA, data=preds)
    o[i] <- segmented(lin.mod, seg.Z=~ichorCNA, npsi=2,
                      control=seg.control(seed=NULL))$psi[1,'Est.']
    res$StoppingRule1 <- 1 * (res$ichorCNA>=round(o[i]))
    res$StoppingRule1 <- factor(res$StoppingRule1,
                                levels=c(0,1),
                            labels=c("ichorCNA<7",  "ichorCNA>=7"))
    mod[[i]] <- coxph(Surv(I(start/30), I(stop/30), status) ~ StoppingRule1,  data=subset(res, !Public.ID %in% sel),
         iter.max = 500000)
    mod.cont[[i]] <- coxph(Surv(I(start/30), I(stop/30), status) ~ ichorCNA,  data=subset(res, !Public.ID %in% sel),
         iter.max = 500000)

}
runif(1)

pdf("FigureS3.pdf", width=10, height=10)
par(mfrow=c(3,2))
plot(density(o), xlab="ichorCNA", cex.axis=1.3, cex.lab=1.5,
     cex.main=1.5, main="a)")
plot(density(sapply(mod, function(x) summary(x)$concordance[1])), xlab="c-index", cex.axis=1.3, cex.lab=1.5,
     cex.main=1.5, main="b)")
plot(density(sapply(mod, function(x) log10(summary(x)$coefficients[1,5]))),
     xlab="Estimated log10 p-value", cex.axis=1.3, cex.lab=1.5,
     cex.main=1.5, main="c)")
m <- cph(Surv(start, stop, status) ~ rcs(ichorCNA), data=res,
             iter.max = 500000)
p <- Predict(m)
plot(yhat ~ ichorCNA, data=p, ylab="Effect of ichorCNA", xlab="ichorCNA",
     type="l", cex.axis=1.3, cex.lab=1.5,
     cex.main=1.5, main="d)")
abline(v=7, col=2, lty=2)
m <- cph(Surv(start, stop, status) ~ rcs(ichorCNA), data=subset(res, Group=="ER+/Her2-"),
             iter.max = 500000)
p <- Predict(m)
plot(yhat ~ ichorCNA, data=p, ylab="Effect of ichorCNA", xlab="ichorCNA",
     type="l", cex.axis=1.3, cex.lab=1.5,
     cex.main=1.5, main="d) ER+/Her2-")
m <- cph(Surv(start, stop, status) ~ rcs(ichorCNA), data=subset(res, Group=="Her2+"),
             iter.max = 500000)
p <- Predict(m)
plot(yhat ~ ichorCNA, data=p, ylab="Effect of ichorCNA", xlab="ichorCNA",
     type="l", cex.axis=1.3, cex.lab=1.5,
     cex.main=1.5, main="e) Her2-")
dev.off()

summary(o)

summary(sapply(mod, coef))
summary(sapply(mod, function(x) summary(x)$concordance[1]))
summary(sapply(mod.cont, function(x) summary(x)$concordance[1]))

res$StoppingRule1 <- 1 * (res$ichorCNA>=7)
res$StoppingRule1 <- factor(res$StoppingRule1,
                            levels=c(0,1),
                            labels=c("ichorCNA<7",  "ichorCNA>=7"))

## Model
dd <- datadist(res)
options(datadist='dd')
summary(cph(Surv(I(start/30), I(stop/30), status) ~ StoppingRule1,  data=res,
         iter.max = 500000, x=TRUE, y=TRUE))



pdf("Figure3.pdf", width=12, height=6)
par(mfrow=c(1,2))
dd <- datadist(res)
options(datadist='dd')
m <- cph(Surv(I(start/30.25), I(stop/30.25), status) ~ StoppingRule1,
         data=subset(res, Group=="ER+/Her2-"),
         iter.max = 500000, x=TRUE, y=TRUE)
m2 <- coxph(Surv(I(start/30.25), I(stop/30.25), status) ~ StoppingRule1,
         data=subset(res, Group=="ER+/Her2-"),
         iter.max = 500000)
summary(m)
survplot(m, StoppingRule1=c("ichorCNA<7",
                            "ichorCNA>=7"),
         col=c("blue", "red"), xlim=c(0,90),lwd=1.5,lty=c(1,1),
         ylab="Progression-free survival",
         xlab="Months since plasma sample")
pval <- round(coef(summary(m2))[1,5], 3)
if (pval < 0.001) pval <- "<0.001)"
tt <- paste("Hazard-ratio:", round(exp(coef(m)), 2), "(p-value", pval)
text(40, 0.9, tt)
## Median time to event

feo <- survest(m, newdata=data.frame(Group=c(rep("ER+/Her2-", 2), rep("Her2+", 2)), StoppingRule1=rep(c("ichorCNA<7", "ichorCNA>=7"), 2)), times=seq(from=0, to=150, by=0.1))

feo$time[apply(feo$surv, 1, function(x) which.min(abs(0.5-x)))]

data.frame(Group=c(rep("ER+/Her2-", 2), rep("Her2+", 2)), StoppingRule1=rep(c("ichorCNA<7", "ichorCNA>=7"), 2),
           Average.Time.Progr=feo$time[apply(feo$surv, 1, function(x) which.min(abs(0.5-x)))])


title("ER+/Her2-")
dd <- datadist(res)
options(datadist='dd')
m <- cph(Surv(I(start/30.25), I(stop/30.25), status) ~ StoppingRule1,
         data=subset(res, Group=="Her2+"),
         iter.max = 500000, x=TRUE, y=TRUE)
m2 <- coxph(Surv(I(start/30.25), I(stop/30.25), status) ~ StoppingRule1,
         data=subset(res, Group=="Her2+"),
         iter.max = 500000)
summary(m)
survplot(m, StoppingRule1=c("ichorCNA<7",
                            "ichorCNA>=7"),
         col=c("blue", "red"), xlim=c(0,90), lwd=1.5,lty=c(1,1),
         ylab="Progression-free survival",
         xlab="Months since plasma sample")
pval <- round(coef(summary(m2))[1,5], 3)
if (pval < 0.001) pval <- "<0.001)"
tt <- paste("Hazard-ratio:", round(exp(coef(m)), 2), "(p-value", pval)
text(40, 0.9, tt)
title("Her2+")
## Median time to event

feo <- survest(m, newdata=data.frame(Group=c(rep("ER+/Her2-", 2), rep("Her2+", 2)), StoppingRule1=rep(c("ichorCNA<7", "ichorCNA>=7"), 2)), times=seq(from=0, to=150, by=0.1))

feo$time[apply(feo$surv, 1, function(x) which.min(abs(0.5-x)))]

data.frame(Group=c(rep("ER+/Her2-", 2), rep("Her2+", 2)), StoppingRule1=rep(c("ichorCNA<7", "ichorCNA>=7"), 2),
           Average.Time.Progr=feo$time[apply(feo$surv, 1, function(x) which.min(abs(0.5-x)))])
dev.off()

alluvial.PFS <- res[,c('Public.ID', 'ichorCNA', 'Date.Start')]
save(alluvial.PFS, file="Alluvial.PFS.RData")


############################################################
############################################################
############################################################
## Survival
############################################################
############################################################
############################################################
OS <- Patients[,c(1, 6)]
OS$OS <- 0
OS$OS[which(!is.na(OS$Date.RIP))] <- 1
OS <- subset(OS, OS==1)
OS <- na.omit(OS)
ichorCNA <- merge(ichorCNA, OS, all.x=T)
ichorCNA$OS[which(is.na(ichorCNA$OS))] <- 0
ichorCNA <- ichorCNA[order(ichorCNA$Public.ID, ichorCNA$Date.RIP),]


ichor.OS <- NULL
for (i in unique(Patients$Public.ID)) {
    min.date <- Patients$Entry.Date[which(Patients$Public.ID==i)]
    max.date <- Patients$Exit.Date[which(Patients$Public.ID==i)]
    sub.ichor <- subset(ichorCNA, Public.ID==i)
    sub.ichor <- sub.ichor[order(sub.ichor$Day),]
    if (nrow(sub.ichor)>1) {
        for (j in 1:(nrow(sub.ichor)-1)) {
            tmp1 <- data.frame(Public.ID=i, start=sub.ichor$Day[j]-min.date,
                               stop=sub.ichor$Day[j+1]-min.date, OS=0,
                               ichorCNA=sub.ichor$ichorCNA[j],
                               date.start=sub.ichor$Day[j])
            ichor.OS <- rbind(ichor.OS, tmp1)
        }
        j <- nrow(sub.ichor)
        if(i %in% OS$Public.ID)  {
            max.date <- OS$Date.RIP[which(OS$Public.ID==i)]
            tmp1 <- data.frame(Public.ID=i, start=sub.ichor$Day[j]-min.date,
                               stop=max.date-min.date, OS=1,
                               ichorCNA=sub.ichor$ichorCNA[j],
                               date.start=sub.ichor$Day[j])
        } else {
            tmp1 <- data.frame(Public.ID=i, start=sub.ichor$Day[j]-min.date,
                               stop=max.date-min.date, OS=sub.ichor$OS[j],
                               ichorCNA=sub.ichor$ichorCNA[j],
                               date.start=sub.ichor$Day[j])
            }
        ichor.OS <- rbind(ichor.OS, tmp1)
    }
}

ichor.OS <- merge(ichor.OS, Patients[,c('Public.ID', 'Group')], all.x=T)
ichor.OS$stop[which(ichor.OS$start==ichor.OS$stop)] <- ichor.OS$stop[which(ichor.OS$start==ichor.OS$stop)] + 0.1
ichor.OS$StoppingRule1 <- 1 * (ichor.OS$ichorCNA>=7)
ichor.OS$StoppingRule1 <- factor(ichor.OS$StoppingRule1,
                            levels=c(0,1),
                            labels=c("ichorCNA<7",  "ichorCNA>=7"))


m <- coxph(Surv(I(start/30), I(stop/30), OS) ~ Group + ichorCNA,  data=ichor.OS)
summary(m)
m <- coxph(Surv(I(start/30), I(stop/30), OS) ~ Group + StoppingRule1,  data=ichor.OS)
summary(m)
alluvial.OS <- ichor.OS[,c('Public.ID', 'ichorCNA', 'date.start')]

save(alluvial.OS, file="Alluvial.OS.RData")

##
pdf("~/Desktop/FigureS4.pdf", height=3, width=8)
dd <- datadist(ichor.OS)
options(datadist='dd')
par(mfrow=c(1,3))
for (i in c("ER+/Her2-", "Her2+")) {
m <- cph(Surv(I(start/30), I(stop/30), OS) ~ Group + StoppingRule1,
           data=ichor.OS, surv=T, iter.max=5000)
survplot(m, StoppingRule1, Group=i,
         col=c("blue", "red"), xlim=c(0,90),
         ylab="Overall survival",
         xlab="Months since plasma sample", cex.ylab=1.2, cex.xlab=1.2, lwd=1.5, lty=c(1,1), label.curves=list(cex=1.3),
         mark.time=T)
title(paste(i, "\n(n=", length(unique(subset(ichor.OS, Group==i)[,1])), ")", sep=""))
}
dev.off()



##################################################
##################################################
## Now NGTAS
##################################################
##################################################
ALLNGTAS <- merge(NGTAS, ichorCNA)
ALLNGTAS <- ALLNGTAS[order(ALLNGTAS$Public.ID, ALLNGTAS$Day),]
ALLNGTAS$Visit <- NA
ALLNGTAS$VisitProp <- NA
for (i in unique(ALLNGTAS$Public.ID)) {
    ids <- which(ALLNGTAS$Public.ID==i)
    ALLNGTAS$Visit[ids] <- as.numeric(as.factor(ALLNGTAS$Day[ids]))
    ALLNGTAS$VisitProp[ids] <- (ALLNGTAS$Visit[ids] -
                                max(ALLNGTAS$Visit[ids]))/ max(ALLNGTAS$Visit[ids])
    }
id <- unique(ALLNGTAS$Public.ID)
id <- sub("P", "", id)
id <- sort(as.numeric(id))
ALLNGTAS$Public.ID <- factor(ALLNGTAS$Public.ID,
                             levels=paste0("P",id))
library(latticeExtra)
library(RColorBrewer)
cols <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
cols <- cols[c(1:9, 11,15)]
p <- xyplot(VAF ~ Visit|Public.ID,
            groups=Gene, type="b", pch=19, col=cols,
            ichorCNA=ichorCNA, data=ALLNGTAS,
            par.settings = list(superpose.symbol =
                                    list(pch=19, lwd=2,
                                         col=cols),
                                superpose.line = list(
                                    col = cols,
                                    lwd = 2,
                                    lty = 1)),
            auto.key=T, cex=0.6)
p <- p + as.layer(xyplot(I(ichorCNA/100) ~ Visit|Public.ID, type="b", pch="+", ichorCNA=ichorCNA, cex=2, col="black",
            data=ALLNGTAS))
pdf("NGTAS-ichorCNA-Legend.pdf")
print(p)
dev.off()

p <- xyplot(VAF ~ Visit|Public.ID,
            groups=Gene, type="b", pch=19,
            ichorCNA=ichorCNA, data=ALLNGTAS,
            par.settings = list(superpose.symbol =
                                    list(pch=19, lwd=2,
                                         col=cols),
                                superpose.line = list(
                                    col = cols,
                                    lwd = 2,
                                    lty = 1)),
            auto.key=F, cex=0.6,
            ylab="Tumour fraction estimation",
            scales = list(x=list(at=seq(from=1, to=13, by=2),
                                 relation="free")))
p <- p + as.layer(xyplot(I(ichorCNA/100) ~ Visit|Public.ID, type="b", pch="+", ichorCNA=ichorCNA, cex=1.7, col="grey",
            data=ALLNGTAS))
pdf("NGTAS-ichorCNA-V1.pdf", width=12/1.4, height=8/1.4)
print(p)
dev.off()
p <- xyplot(VAF ~ Visit|Public.ID,
            groups=Gene, type="b", pch=19,
            ichorCNA=ichorCNA, data=ALLNGTAS,
            par.settings = list(superpose.symbol =
                                    list(pch=19, lwd=2,
                                         col=cols),
                                superpose.line = list(
                                    col = cols,
                                    lwd = 2,
                                    lty = 1)),
            auto.key=F, cex=0.6,
            ylab="Tumour fraction estimation")
p <- p + as.layer(xyplot(I(ichorCNA/100) ~ Visit|Public.ID, type="b", pch="+", ichorCNA=ichorCNA, cex=1.7, col="grey",
            data=ALLNGTAS))
pdf("NGTAS-ichorCNA-V2.pdf", width=12/1.4, height=8/1.4)
print(p)
dev.off()

p <- xyplot(VAF ~ VisitProp|Public.ID,
            groups=Gene, type="b", pch=19,
            ichorCNA=ichorCNA, data=ALLNGTAS,
            par.settings = list(superpose.symbol =
                                    list(pch=19, lwd=2,
                                         col=cols),
                                superpose.line = list(
                                    col = cols,
                                    lwd = 2,
                                    lty = 1)),
            auto.key=F, cex=0.6,
            ylab="Tumour fraction estimation", xlab="Plasma Sample",
            scales=list(x=list(at=seq(from=0, to=1, by=0.25), labels=rep("", 5))))
p <- p + as.layer(xyplot(I(ichorCNA/100) ~ VisitProp|Public.ID, type="b", pch="+", ichorCNA=ichorCNA, cex=1.7, col="grey",
            data=ALLNGTAS))
pdf("NGTAS-ichorCNA-V3.pdf", width=12/1.4, height=8/1.4)
print(p)
dev.off()

NGTASMax <- aggregate(ALLNGTAS[,c('ichorCNA', 'VAF')], by=list(ALLNGTAS$Public.ID, ALLNGTAS$Day), FUN=max)
NGTASMean <- aggregate(ALLNGTAS[,c('ichorCNA', 'VAF')], by=list(ALLNGTAS$Public.ID, ALLNGTAS$Day), FUN=mean)
NGTASSD <- aggregate(ALLNGTAS[,c('ichorCNA', 'VAF')], by=list(ALLNGTAS$Public.ID, ALLNGTAS$Day), FUN=sd)

par(mfrow=c(1,3))
plot(VAF ~ ichorCNA, data=NGTASMax)
plot(VAF ~ ichorCNA, data=NGTASMean)
plot(NGTASMax$ichorCNA, NGTASSD$VAF)

colnames(NGTASMax) <- c("Public.ID", "Days", 'ichorCNA', "VAF")
NGTASMax <- NGTASMax[order(NGTASMax$Public.ID, NGTASMax$Day),]

res <- NULL
for (i in 1:nrow(Patients)) {
    tmp1 <- subset(RECIST, Public.ID==Patients$Public.ID[i])
    if (nrow(tmp1)>0) {
        tmp2 <- subset(NGTASMax, Public.ID==Patients$Public.ID[i])
        possible.dates <- data.frame(Date=c(tmp1$Day, tmp2$Day,
                                        Patients$Exit.Date[i]),
                                 Event=c(tmp1$Progression,
                                         rep("ichor",
                                             length(tmp2$Day)),
                                         "exit"))
        possible.dates <- possible.dates[order(possible.dates$Date),]
        possible.dates <- possible.dates[which(possible.dates$Event!="NO"),]
        relapse <- which(possible.dates$Event=="YES")
        if (length(relapse)>0) possible.dates <- possible.dates[1:relapse[1],]
        for (j in 1:nrow(possible.dates)) {
            if (possible.dates$Event[j]=="ichor") {
                start <- as.numeric(possible.dates$Date[j] -
                                    Patients$Entry.Date[i])
                end <- as.numeric(possible.dates$Date[j+1] - Patients$Entry.Date[i])
                status <- ifelse(possible.dates$Event[j+1]=="YES", 1, 0)
                tmp <- data.frame(Public.ID=Patients$Public.ID[i],
                                  start=start, stop=end, status=status,
                                  ichorCNA=tmp2$ichorCNA[which(tmp2$Day==possible.dates$Date[j])],
                                  NGTAS=tmp2$VAF[which(tmp2$Day==possible.dates$Date[j])],
                                  Date.Start=possible.dates$Date[j],
                                  Date.Stop=possible.dates$Date[j+1])
                res <- rbind(res, tmp)
            }
        }
    }
}



res <- merge(res, Patients[,c(1,2)], all.x=T)

sub.ichorPFS <- res

m.ichor <- coxph(Surv(start, stop, status) ~ ichorCNA, data=sub.ichorPFS)
m.NGTAS <- coxph(Surv(start, stop, status) ~ sqrt(NGTAS), data=sub.ichorPFS)
m.both <- coxph(Surv(start, stop, status) ~ ichorCNA + sqrt(NGTAS), data=sub.ichorPFS)

anova(m.ichor, m.both)
alluvial.NGTAS <- unique(ALLNGTAS[,c(1, 2)])
alluvial.NGTAS.PFS <- sub.ichorPFS[,c('Public.ID', 'NGTAS', 'Date.Start')]
save(alluvial.NGTAS.PFS, alluvial.NGTAS, file="Alluvial.NGTAS.RData")


##################################################
##################################################
## Now CA15-3
##################################################
##################################################
RECIST <- RECIST[order(RECIST$Public.ID, RECIST$Day),]
CA15.3 <- CA15.3[order(CA15.3$Public.ID, CA15.3$Day),]
res <- NULL
for (i in 1:nrow(Patients)) {
    tmp1 <- subset(RECIST, Public.ID==Patients$Public.ID[i])
    if (nrow(tmp1)>0) {
        tmp2 <- subset(CA15.3, Public.ID==Patients$Public.ID[i])
        if (nrow(tmp2)>0) {
        possible.dates <- data.frame(Day=c(tmp1$Day, tmp2$Day,
                                        Patients$Exit.Date[i]),
                                 Event=c(tmp1$Progression,
                                         rep("ichor",
                                             length(tmp2$Day)),
                                         "exit"))
        possible.dates <- possible.dates[order(possible.dates$Day),]
        possible.dates <- possible.dates[which(possible.dates$Event!="NO"),]
        relapse <- which(possible.dates$Event=="YES")
        if (length(relapse)>0) possible.dates <- possible.dates[1:relapse[1],]
        for (j in 1:nrow(possible.dates)) {
        if (possible.dates$Event[j]=="ichor") {
            start <- as.numeric(possible.dates$Day[j] -
                                Patients$Entry.Date[i])
            end <- as.numeric(possible.dates$Day[j+1] - Patients$Entry.Date[i])
            status <- ifelse(possible.dates$Event[j+1]=="YES", 1, 0)
            tmp <- data.frame(Public.ID=Patients$Public.ID[i],
                              start=start, stop=end, status=status,
                              CA15.3=tmp2$Units.CA[which(tmp2$Day==possible.dates$Day[j])],
                              Day.Start=possible.dates$Day[j],
                              Day.Stop=possible.dates$Day[j+1])
            res <- rbind(res, tmp)
            }
        }
        }
    }
}



res <- merge(res, Patients[,c(1,2)], all.x=T)

CA153.PFS <- res

ss <- unique(intersect(CA153.PFS$Public.ID, ichorCNA.res$Public.ID))

sub.CA153.PFS <- subset(CA153.PFS, Public.ID %in% ss)
summary(coxph(Surv(start, stop, status) ~ log(CA15.3), data=sub.CA153.PFS))

## Compare with a model with ichorCNA

sub.ichor <- ichorCNA.res[which(ichorCNA.res$Public.ID %in% ss),]
summary(coxph(Surv(start, stop, status) ~ ichorCNA, data=sub.ichor))


## Plot CA15.3

ALLCA15.3 <- CA15.3
ALLCA15.3$Type <- "sqrt(CA15-3)"
colnames(ALLCA15.3)[3] <- "Score"
ALLCA15.3$Score <- sqrt(ALLCA15.3$Score)
ALLCA15.3 <- rbind(ALLCA15.3, data.frame(Public.ID=ichorCNA$Public.ID, Day=ichorCNA$Day, Score=ichorCNA$ichorCNA, Type="ichorCNA"))
ss <- unique(intersect(CA15.3$Public.ID, ichorCNA$Public.ID))
ALLCA15.3 <- ALLCA15.3[order(ALLCA15.3$Public.ID, ALLCA15.3$Day),]
ALLCA15.3 <- subset(ALLCA15.3, Public.ID %in% ss)

ALLCA15.3$Visit <- NA
ALLCA15.3$VisitProp <- NA
for (i in unique(ALLCA15.3$Public.ID)) {
    ids <- which(ALLCA15.3$Public.ID==i)
    ALLCA15.3$Visit[ids] <- as.numeric(as.factor(ALLCA15.3$Day[ids]))
    ALLCA15.3$VisitProp[ids] <- (ALLCA15.3$Visit[ids] -
                                max(ALLCA15.3$Visit[ids]))/ max(ALLCA15.3$Visit[ids])
    }
id <- unique(ALLCA15.3$Public.ID)
id <- sub("P", "", id)
id <- sort(as.numeric(id))
ALLCA15.3$Public.ID <- factor(ALLCA15.3$Public.ID,
                             levels=paste0("P",id))

pdf("CA15.3-ichorCNA.pdf", width=12/1.4, height=8/1.4)
xyplot(Score ~ VisitProp|Public.ID, groups=Type, data=ALLCA15.3, type="b", pch=19, auto.key=F, cex=0.4,
                   par.settings = list(superpose.symbol =
                                    list(pch=19, lwd=2,
                                         col=cols),
                                superpose.line = list(
                                    col = cols,
                                    lwd = 2,
                                    lty = 1)),

                   ylab="Tumour fraction estimation", xlab="Plasma Sample",
            scales=list(x=list(at=seq(from=0, to=1, by=0.25), labels=rep("", 5))))
dev.off()
pdf("CA15.3-ichorCNALegend.pdf", width=12/1.4, height=8/1.4)
xyplot(Score ~ VisitProp|Public.ID, groups=Type, data=ALLCA15.3, type="b", pch=19, auto.key=T, cex=0.4,
                   par.settings = list(superpose.symbol =
                                    list(pch=19, lwd=2,
                                         col=cols),
                                superpose.line = list(
                                    col = cols,
                                    lwd = 2,
                                    lty = 1)),

                   ylab="Tumour fraction estimation", xlab="Plasma Sample",
            scales=list(x=list(at=seq(from=0, to=1, by=0.25), labels=rep("", 5))))
dev.off()

## Correlation between both using interpolation
cors.spearman <- NULL
cors.linear <- NULL
for (i in unique(ALLCA15.3$Public.ID)) {
    tmp <- subset(ALLCA15.3, Public.ID==i)
    tmp1 <- subset(tmp, Type=="ichorCNA")
    tmp2 <- subset(tmp, Type=="sqrt(CA15-3)")
    dat <- c(max(min(tmp1$Day), min(tmp2$Day)), min(max(tmp1$Day), max(tmp2$Day)))
    xout <- sort(unique(c(tmp1$Day, tmp2$Day)))
    xout <- xout[which(xout >= dat[1] & xout<=dat[2])]
    if (length(xout) > 1) {
        inter.ichorCNA <- approx(tmp1$Day, tmp1$Score, xout=xout, method="linear")$y
        inter.CA15.3 <- approx(tmp2$Day, tmp2$Score, xout=xout, method="linear")$y
        cors.linear <- c(cors.linear, cor(inter.ichorCNA, inter.CA15.3))
        cors.spearman <- c(cors.spearman, cor(inter.ichorCNA, inter.CA15.3, method="spearman"))
    } else {
        cors.spearman <- c(cors.spearman, NA)
        cors.linear <- c(cors.linear, NA)
        }
}
cors.obt <- data.frame(Public.ID=unique(ALLCA15.3$Public.ID), cor=cors.linear, cor.spearman=cors.spearman)



## Agreement based on threshold
colnames(NGTASMax)[2] <- "Day"
ichorCNA <- merge(ichorCNA, NGTASMax[,c(1, 2, 4)], all.x=T)

ichorCNA$CA15.3 <- NA
ichorCNA$CA15.3Day <- NA
ichorCNA$Day.CA <- NA
for (i in 1:nrow(ichorCNA)) {
    tmp <- subset(CA15.3, Public.ID==ichorCNA$Public.ID[i])
    tmp <- tmp[which.min(abs(as.numeric(tmp$Day-ichorCNA$Day[i]))),]
    if (nrow(tmp)>0) {
        ichorCNA$CA15.3[i] <- tmp$Units.CA
        ichorCNA$CA15.3Day[i] <- ichorCNA$Day[i] - tmp$Day
        ichorCNA$Day.CA[i] <- tmp$Day
        }
}

ichorCNA$CA15.3[which(abs(ichorCNA$CA15.3Day)>15)] <- NA

## we should make sure that we only use CA15.3 once
for (i in unique(ichorCNA$Public.ID)) {
    for (j in subset(ichorCNA, Public.ID==i)$Day.CA) {
        ids <- which(ichorCNA$Public.ID==i & ichorCNA$Day.CA==j)
        if(sum(!is.na(ichorCNA$CA15.3[ids]))>1) {
            ichorCNA$CA15.3[ids[-which.min(abs(ichorCNA$CA15.3Day[ids]))]] <- NA
            }

        }
    }
ichorCNA$logCA15 <- log(ichorCNA$CA15.3)
pairs(ichorCNA[,c('ichorCNA', 'VAF', 'logCA15')])

ichorCNA <- merge(ichorCNA, NGTAS[,c(1, 2, 4)], all.x=T)

## Compare pred.ability
library(caTools)
ichorCNA$ichorRule <- 1 * (ichorCNA$ichorCNA>6.5)
ichorCNA$ichorRule <- factor(ichorCNA$ichorRule, levels=c(0,1), labels=c("Continue", "Stop"))
ichorCNA$VAFRule <- 1 * (ichorCNA$VAF>0.025)
ichorCNA$VAFRule <- factor(ichorCNA$VAFRule,  levels=c(0,1), labels=c("Continue", "Stop"))
ichorCNA$CA15.3Rule <- 1 * (ichorCNA$CA15.3>=31)
ichorCNA$CA15.3Rule <- factor(ichorCNA$CA15.3Rule,  levels=c(0,1), labels=c("Continue", "Stop"))

round(sum(diag(table(ichorCNA$ichorRule, ichorCNA$VAFRule))) / sum((table(ichorCNA$ichorRule, ichorCNA$VAFRule))), 2)
round(sum(diag(table(ichorCNA$ichorRule, ichorCNA$CA15.3Rule))) / sum((table(ichorCNA$ichorRule, ichorCNA$CA15.3Rule))), 2)


##################################################
##################################################
## Survival with CA15.3
##################################################
##################################################
CA15.3 <- CA15.3[order(CA15.3$Public.ID, CA15.3$Day),]
CAOS <- merge(CA15.3, OS, all.x=T)
CAOS$OS[which(is.na(CAOS$OS))] <- 0
CAOS <- CAOS[order(CAOS$Public.ID, CAOS$Date.RIP),]


CAOS.OS <- NULL
for (i in unique(Patients$Public.ID)) {
    min.date <- Patients$Entry.Date[which(Patients$Public.ID==i)]
    max.date <- Patients$Exit.Date[which(Patients$Public.ID==i)]
    sub.CAOS <- subset(CAOS, Public.ID==i)
    sub.CAOS <- sub.CAOS[order(sub.CAOS$Day),]
    if (nrow(sub.CAOS)>1) {
        for (j in 1:(nrow(sub.CAOS)-1)) {
            tmp1 <- data.frame(Public.ID=i, start=sub.CAOS$Day[j]-min.date,
                               stop=sub.CAOS$Day[j+1]-min.date, OS=0,
                               Units.CA=sub.CAOS$Units.CA[j],
                               date.start=sub.CAOS$Day[j])
            CAOS.OS <- rbind(CAOS.OS, tmp1)
        }
        j <- nrow(sub.CAOS)
        if(i %in% OS$Public.ID)  {
            max.date <- OS$Date.RIP[which(OS$Public.ID==i)]
            tmp1 <- data.frame(Public.ID=i, start=sub.CAOS$Day[j]-min.date,
                               stop=max.date-min.date, OS=1,
                               Units.CA=sub.CAOS$Units.CA[j],
                               date.start=sub.CAOS$Day[j])
        } else {
            tmp1 <- data.frame(Public.ID=i, start=sub.CAOS$Day[j]-min.date,
                               stop=max.date-min.date, OS=sub.CAOS$OS[j],
                               Units.CA=sub.CAOS$Units.CA[j],
                               date.start=sub.CAOS$Day[j])
            }
        CAOS.OS <- rbind(CAOS.OS, tmp1)
    }
}

CAOS.OS <- merge(CAOS.OS, Patients[,c('Public.ID', 'Group')], all.x=T)
CAOS.OS$stop[which(CAOS.OS$start==CAOS.OS$stop)] <- CAOS.OS$stop[which(CAOS.OS$start==CAOS.OS$stop)] + 0.1
m <- coxph(Surv(I(start/30), I(stop/30), OS) ~ Group + sqrt(Units.CA),  data=CAOS.OS)
summary(m)

summary(coxph(Surv(I(start/30), I(stop/30), OS) ~ Group + ichorCNA,  data=ichor.OS[which(ichor.OS$Public.ID %in% CAOS.OS$Public.ID),]))
##############################################################
##############################################################
##############################################################

discr <- ichorCNA[which((ichorCNA$ichorCNA < 6.5 & ichorCNA$VAF>0.05 ) |
                        (ichorCNA$ichorCNA < 6.5 & ichorCNA$CA15.3>=31) |
                        (ichorCNA$ichorCNA > 6.5 & ichorCNA$VAF<0.05 ) |
                        (ichorCNA$ichorCNA > 6.5 & ichorCNA$CA15.3<31)
                        ),]



discr <- merge(discr, Patients[,c(1,2)])

## Closest scan
discr$Progression <- NA
discr$Day.SCAN <- NA
discr$Treatment <- NA
discr$Start.Treatment <- NA
discr$End.Treatment <- NA

Treatments <- read.table("TableS1.txt", header=T, sep="\t")
for (i in 1:nrow(discr)) {
    tmp1 <- subset(Treatments, Public.ID==discr$Public.ID[i])
    tmp1 <- tmp1[which(tmp1$Start.Day < discr$Day[i] &
                       tmp1$End.Day >= discr$Day[i]),]
    if (nrow(tmp1)>1) stop()
    if (nrow(tmp1)>0) {
        discr$Treatment[i] <- tmp1$Treatment
        discr$Start.Treatment[i] <- tmp1$Start.Day
        discr$End.Treatment[i] <- tmp1$End.Day

        }
    tmp2 <- subset(RECIST, Public.ID==discr$Public.ID[i])
    tmp2 <- tmp2[which.min(abs(tmp2$Day-discr$Day[i])),]
    if (nrow(tmp2)>0) {
        discr$Progression[i] <- tmp2$Progression
        discr$Day.SCAN[i] <- tmp2$Day
    }


    }



discr$T.Scan <- as.numeric(discr$Day.SCAN - discr$Day)
discr$T.StartTreat <- as.numeric(discr$Start.Treatment - discr$Day)
discr$T.EndTreat <- as.numeric(discr$End.Treatment - discr$Day)

discr <- discr[which(abs(discr$T.Scan) <=90),]
pdf("discrepancies.pdf", width=12, height=10)
par(mfrow=c(1,3), oma= rep(3,4), mar=rep(2,4))
for (g in unique(discr$Group)) {
    tmp <- subset(discr, Group==g)
    tmp <- tmp[order(tmp$Progression, tmp$ichorRule),]
    tmp$startmark <- 0
    tmp$startmark[which(tmp$T.StartTreat<min(tmp$T.Scan))] <- 1
    tmp$endmark <- 0
    tmp$endmark[which(tmp$T.EndTreat>max(tmp$T.Scan))] <- 1
    tmp$T.StartTreat[which(tmp$T.StartTreat<min(tmp$T.Scan))] <- min(tmp$T.Scan)
    tmp$T.EndTreat[which(tmp$T.EndTreat>max(tmp$T.Scan))] <- max(tmp$T.Scan)
    plot(0,0, type="l",ylim=c(-0.1, 0.1 + max(table(discr$Group))),
         xlim=c(min(-5, min(tmp$T.Scan)), max(5, max(tmp$T.Scan))), axes=F,
     xlab="Days", main=g, cex.main=2)
    for (i in 1:nrow(tmp)) {
        if (tmp$startmark[i]==1 & tmp$endmark[i]==0) {
            polygon(x=c(tmp$T.StartTreat[i]+3, tmp$T.EndTreat[i],
                    tmp$T.EndTreat[i],
                    tmp$T.StartTreat[i]+3, tmp$T.StartTreat[i],tmp$T.StartTreat[i]+3),
                y=c(max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i +0.4,
                    max(table(discr$Group)) - nrow(tmp) + i+0.4,
                    max(table(discr$Group)) - nrow(tmp) + i,
                    max(table(discr$Group)) - nrow(tmp) + i-0.4),
                lwd=0.5)

        } else if (tmp$startmark[i]==1 & tmp$endmark[i]==1) {
            polygon(x=c(tmp$T.StartTreat[i]+3, tmp$T.EndTreat[i]-3,
                    tmp$T.EndTreat[i], tmp$T.EndTreat[i]-3,
                    tmp$T.StartTreat[i]+3, tmp$T.StartTreat[i],tmp$T.StartTreat[i]+3),
                y=c(max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i,
                    max(table(discr$Group)) - nrow(tmp) + i +0.4,
                    max(table(discr$Group)) - nrow(tmp) + i+0.4,
                    max(table(discr$Group)) - nrow(tmp) + i,
                    max(table(discr$Group)) - nrow(tmp) + i-0.4),
                lwd=0.5)
    } else if (tmp$startmark[i]==0 & tmp$endmark[i]==1) {
            polygon(x=c(tmp$T.StartTreat[i], tmp$T.EndTreat[i]-3,
                    tmp$T.EndTreat[i], tmp$T.EndTreat[i]-3,
                    tmp$T.StartTreat[i], tmp$T.StartTreat[i]),
                    y=c(max(table(discr$Group)) - nrow(tmp) + i-0.4,
                        max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i,
                    max(table(discr$Group)) - nrow(tmp) + i +0.4,
                    max(table(discr$Group)) - nrow(tmp) + i+0.4,
                    max(table(discr$Group)) - nrow(tmp) + i),

                lwd=0.5)

        } else {
            polygon(x=c(tmp$T.StartTreat[i], tmp$T.EndTreat[i],
                    tmp$T.EndTreat[i],
                    tmp$T.StartTreat[i], tmp$T.StartTreat[i]),
                y=c(max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i-0.4,
                    max(table(discr$Group)) - nrow(tmp) + i +0.4,
                    max(table(discr$Group)) - nrow(tmp) + i+0.4,
                    max(table(discr$Group)) - nrow(tmp) + i-0.4),
                lwd=0.5)
            }
}
for (i in 1:nrow(tmp)) {
    points(tmp$T.Scan[i], max(table(discr$Group)) - nrow(tmp) + i, pch=16, cex=1.5, col=ifelse(tmp$Progression[i]=="YES", "red", "blue"))
    col <- "grey"
    if(!is.na(tmp$VAF[i])) {
        col <- ifelse(tmp$VAF[i]>=0.01, "red", "blue")
        }
    points(0, max(table(discr$Group)) - nrow(tmp) + i, pch=ifelse(tmp$ichorCNA[i]>6.5, 24, 25),
           col=col, bg=col, cex=1.5)
    if(!is.na(tmp$CA15.3[i])) {
        points(tmp$CA15.3Day[i], max(table(discr$Group)) - nrow(tmp) + i+0.35, pch=ifelse(tmp$CA15.3[i]>=30, "+", "-"), cex=2)
        }
}
    pos.axis <- c(-3, -2, -5)
    names(pos.axis) <- unique(discr$Group)
    axis(1, line=nrow(tmp) - max(table(discr$Group))+pos.axis[g], at=round(seq(from=min(tmp$T.Scan)+5, to=max(tmp$T.Scan)-5, length=5)))
axis(2, at=(max(table(discr$Group)) - nrow(tmp)+1):max(table(discr$Group)),
     labels=tmp$Public.ID, las=2, cex.axis=1.2)
mtext("Days", side=1, line=nrow(tmp) - max(table(discr$Group))+pos.axis[g] +3, cex=1.1, at=median(round(seq(from=min(tmp$T.Scan)+5, to=max(tmp$T.Scan)-5, length=5))))
}

legend(-8, 40, pch=c(16, 16), col=c("blue", "red"), legend=c("Stable dis./Partial resp.", "Progressive disease"), cex=2,
       bty="n")

legend(-9, 35, pch=24, legend="", bty="n", col="grey", cex=2, pt.bg="grey")
legend(-6, 35, pch=24, legend="", bty="n", cex=2, col="red", pt.bg="red")
legend(-3, 35, pch=24, legend="", bty="n", cex=2, col="blue", pt.bg="blue")
legend(0, 35,
       legend="ichorCNA>=7",
       bty="n", cex=2)
legend(-9, 33, pch=25, legend="", bty="n", cex=2, col="grey", pt.bg="grey")
legend(-6, 33, pch=25, legend="", bty="n", cex=2, col="red", pt.bg="red")
legend(-3, 33, pch=25, legend="", bty="n", cex=2, col="blue", pt.bg="blue")
legend(0, 33,
       legend="ichorCNA<=7",
       bty="n", cex=2)


legend(-9, 30, pch=24, legend="", bty="n", col="red", cex=2, pt.bg="red")
legend(-6, 30, pch=25, legend="", bty="n", cex=2, col="red", pt.bg="red")
legend(-7, 30,
       legend="NGTAS positive",
       bty="n", cex=2)
legend(-9, 27, pch=24, legend="", bty="n", cex=2, pt.bg="blue")
legend(-6, 27, pch=25, legend="", bty="n", cex=2, pt.bg="blue")
legend(-7, 27,
       legend="NGTAS negative",
       bty="n", cex=2)
legend(-15, 25,
       legend="+/- CA15-3 positive/negative",
       bty="n", cex=2)
dev.off()

