## Threshold-based model
rm(list=ls())
RECIST <- read.table("TableS3.txt", header=T, sep="\t")
ichorCNA <- read.table("TableS4.txt", header=T, sep="\t")
Treatments <- read.table("TableS2.txt", header=T, sep="\t")

res <- NULL
for (i in unique(Treatments$Public.ID)) {
    tmp1 <- subset(Treatments, Public.ID==i & !is.na(Treatments$End.Day))
    for (j in 1:nrow(tmp1)) {
        baseline.ichor <- subset(ichorCNA, Public.ID==i & Day<tmp1$Start.Day[j])
        if (nrow(baseline.ichor)>0) {
        baseline.ichor.day <- max(baseline.ichor$Day)
        baseline.ichor <- baseline.ichor$ichorCNA[which.max(baseline.ichor$Day)]
        } else {
            baseline.ichor.day <- NA
            baseline.ichor <- NA
            }

        tmp2 <- subset(RECIST, Public.ID==i & Day > tmp1$Start.Day[j] &
                               Day <= tmp1$End.Day[j] & !is.na(Progression))
        tmp3 <- subset(ichorCNA, Day > tmp1$Start.Day[j] &
                               Day <= tmp1$End.Day[j] & Public.ID==i)
        ## Now get all scans and closest ichor/CA15.3
        if (nrow(tmp2)>0 | nrow(tmp3)>0) {
            if (nrow(tmp3)==0) {
                res <- rbind(res,
                             data.frame(tmp1[j,],
                                        baseline.ichor=baseline.ichor,
                                        baseline.ichor.day=baseline.ichor.day,
                                        Day.ichor=NA,
                                             ichorCNA=NA,
                           Day.SCAN=tmp2$Day,
                           Progression=tmp2$Progression))
            }
            if (nrow(tmp2)==0) {
                res <- rbind(res,
                             data.frame(tmp1[j,],
                                        baseline.ichor=baseline.ichor,
                                        baseline.ichor.day=baseline.ichor.day,
                                        Day.ichor=NA,
                                        ichorCNA=tmp3$ichorCNA,
                                        Day.SCAN=NA,
                                        Progression=NA))
                baseline.ichor=tmp3$ichorCNA
                baseline.ichor.day=NA
            }
            if (nrow(tmp2)>0 & nrow(tmp3)>0) {
                for (r in 1:nrow(tmp2)) {
                    ids <- which(tmp3$Day <= tmp2$Day[r])
                    if (length(ids)>0) {
                        for (k in 1:length(ids)) {
                            bas <- subset(tmp3, Day < tmp3$Day[k])
                            if (nrow(bas)==0) {
                                baseline.ichor <- NA
                                baseline.ichor.day <- NA
                            } else {
                                baseline.ichor <- bas$ichorCNA[which.max(bas$Day)]
                                baseline.ichor.day <- max(bas$Day)
                                }
                            res <- rbind(res,
                                         data.frame(tmp1[j,],
                                                    baseline.ichor=baseline.ichor,
                                                    baseline.ichor.day=baseline.ichor.day,
                                                    Day.ichor=tmp3$Day[k],
                                                    ichorCNA=tmp3$ichor[k],
                                                    Day.SCAN=tmp2$Day[r],
                                                    Progression=tmp2$Progression[r]))
                            baseline.ichor=tmp3$ichorCNA[k]
                            baseline.ichor.day=tmp3$Day[k]
                        }
                    } else {
                        res <- rbind(res,
                                     data.frame(tmp1[j,],
                                                baseline.ichor=baseline.ichor,
                                                baseline.ichor.day=baseline.ichor.day,
                                                Day.ichor=NA,
                                                ichorCNA=NA,
                                                Day.SCAN=tmp2$Day[r],
                                                Progression=tmp2$Progression[r]))
                    }
                }
            }
        }
    }
}

res$t.ichor <- as.numeric(res$Day.ichor - res$Start.Day)
res$t.scan <- as.numeric(res$Day.SCAN - res$Start.Day)
res$t.ichor.baseline <- as.numeric(res$Start.Day-res$baseline.ichor.day)


res$Treat.ID <- paste(res$Treatment, res$Start.Day, sep="/")

res$Public.ID <- as.character(res$Public.ID)
res$Treat.ID <- as.character(res$Treat.ID)

res <- res[order(res$Public.ID, res$Day.ichor),]

ALLRES <- res
res <- res[which(res$t.scan <= 30*6 & res$t.scan > 21 & res$t.ichor <= 30*6
                 & res$t.ichor > 21),]

res <- res[which(!is.na(res$Group)),]

res$StoppingRule <- NA
res$StoppingRule[which(res$ichorCNA < 7)] <- 0
res$StoppingRule[which((res$ichorCNA >= 7) & (res$ichorCNA >= res$baseline.ichor))] <- 1
res$StoppingRule[which((res$ichorCNA >= 7) & (res$ichorCNA < res$baseline.ichor))] <- 0

res <- res[which(!is.na(res$StoppingRule)),]

res$Diff.scan <- res$t.scan - res$t.ichor
res$far <- 1 * (res$Diff.scan > 90)

res$Action <- factor(res$StoppingRule, levels=c(0,1),
                     labels=c("Continue", "Stop"))

res1 <- subset(res, far==0)

alluvial.ichor.Respos <- unique(res1[,c(1,12)])

table(res1$Action, res1$Progression)

library(caret)
obs <- factor(res1$Action, levels=c("Continue", "Stop"), labels=c("NO", "YES"))
reference <- factor(res1$Progression, levels=c("NO", "YES"), labels=c("NO", "YES"))
round(sensitivity(obs, reference, positive = "YES"), 2)
round(specificity(obs, reference, negative = "NO"), 2)
round(posPredValue(obs, reference, positive = "YES"), 2)
round(negPredValue(obs, reference, negative = "NO"), 2)

aggregate(I(res1$t.scan - res1$t.ichor), by=list(Action=res1$Action, Progression=res1$Progression), mean)

res1$time.before <- res1$t.scan - res1$t.ichor
t.test(time.before ~ Progression,
       data=subset(res1, Action=="Stop"))
t.test(time.before ~ Progression,
       data=subset(res1, Action=="Continue"))

write.table(res, file="ThresholdPredictions.txt", row.names=F, sep="\t", quote=F)
write.table(res1, file="ThresholdPredictions_90days.txt", row.names=F, sep="\t", quote=F)



BAYML <- read.table("BAYML_Preds.txt", header=T, sep="\t")
th <- res1[,c('Public.ID', 'Day.ichor', 'Day.SCAN', 'Action', 'Progression')]
th$T <- th$Day.SCAN - th$Day.ichor
th$ID <- 1:nrow(th)
BAYML <- unique(BAYML)
lastCT <- merge(BAYML, th)

Good <- subset(lastCT, mean>=0.23 & Progression=="YES")
Good <- Good[order(Good$T),]


pdf("SwimmerPlot_1.pdf")
par(mfrow=c(1,1), oma= rep(2,4), mar=rep(2,4))
plot(0,0, type="l",ylim=c(0.8, nrow(Good)), xlim=c(-4, 100), axes=F,
     xlab="Days", main="Concordant progressions dynamic model")
Good$Treatment_new_final <- factor(Good$Treatment_new_final)
Cols <- c(brewer.pal(6, "Set1"), brewer.pal(length(levels(Good$Treatment_new_final))-6, "Set2"))
Cols.adj <- adjustcolor(Cols, alpha=0.3)
for (i in 1:nrow(Good)) {
    polygon(x=c(-3, Good$T[i], Good$T[i], -3, -3),
            y=c(i-0.4, i-0.4, i+0.4, i+0.4, i-0.4),
            col="lightgrey",
            lwd=0.5)
}
for (i in 1:nrow(Good)) {
    tmp <- Good[i,]
    col <- "red"
    if (tmp$Action=="Continue") col="blue"
    points(0, i, pch=17, col=col, cex=1.5)
    points(tmp$T, i, pch=19, cex=1.5, col="red")
    segments(0, i, tmp$T, i)
        }
axis(1)
axis(2, at=1:nrow(Good),
     labels=Good$Public.ID, las=2, cex.lab=0.7, cex.axis=0.7)
mtext("Days", side=1, line=2)
box()
legend(50, 23.5, pch=c(17, 19), legend=c("ichorCNA", "CT Scan"), bty="n", col=c("grey", "red"), ncol=2)
legend(50, 22.0, pch=c(17, 17), col=c("red", "blue"),
       legend=c("Threshold model concordant", "Threshold model discordant"),
       bty="n")
dev.off()


## Ca15.5 + NGTAS

NGTAS <- read.table("TableS4.txt", header=T, sep="\t")
colnames(NGTAS)[2] <- "Day.ichor"
lastCT <- merge(lastCT, NGTAS, all.x=T)

CA15.3 <- read.table("TableS5.txt", header=T, sep="\t")
colnames(CA15.3)[2] <- "Date.CA"

lastCT$CA15.3 <- NA
lastCT$CA15.3Days <- NA

for (i in 1:nrow(lastCT)) {
    tmp <- subset(CA15.3, Public.ID==lastCT$Public.ID[i])
    tmp <- tmp[which(tmp$Date.CA<=lastCT$Day.SCAN[i]),]
    tmp <- tmp[which.min(abs(as.numeric(tmp$Date.CA-lastCT$Day.ichor[i]))),]
    if (nrow(tmp)>0) {
        lastCT$CA15.3[i] <- tmp$Units.CA
        lastCT$CA15.3Days[i] <- lastCT$Day.ichor[i] - tmp$Date.CA
        }
}

Wrong <- subset(lastCT, (mean<0.23 & Progression=="YES"))


Wrong$CA15.3[which(abs(Wrong$CA15.3Days)>15)] <- NA

Wrong <- Wrong[order(Wrong$T),]

pdf("SwimmerPlot_2.pdf")
par(mfrow=c(1,1), oma= rep(2,4), mar=rep(2,4))

plot(0,0, type="l",ylim=c(0.7, 0.3 + nrow(Wrong)), xlim=c(-4,110), axes=F,
     xlab="Days", main="Discordant progressions dynamic model")
Wrong$Treatment_new_final <- factor(Wrong$Treatment_new_final)
Cols <- c(brewer.pal(3, "Set1"), brewer.pal(length(levels(Wrong$Treatment_new_final))-3, "Set2"))
Cols <- adjustcolor(Cols, alpha=0.3)
for (i in 1:nrow(Wrong)) {
    polygon(x=c(-3, Wrong$T[i], Wrong$T[i], -3, -3),
            y=c(i-0.4, i-0.4, i+0.4, i+0.4, i-0.4), col="lightgrey",
            lwd=0.5)
}
for (i in 1:nrow(Wrong)) {
    tmp <- Wrong[i,]
    col <- "red"
    if (tmp$Action=="Continue") col="blue"
    points(tmp$T, i, pch=19, cex=1.5, col="red")
    points(0, i, pch=17, col=col, cex=1.5)
    segments(0, i, tmp$T, i)
    if(!is.na(tmp$CA15.3)) {
        points(tmp$CA15.3Days, i+0.35, pch=ifelse(tmp$CA15.3>=30, "+", "-"), cex=2)
        }
    if(!is.na(tmp$VAF)) {
        points(0, i-0.35, pch=ifelse(tmp$VAF>=0.01, "+", "-"), cex=2)
        }

        }
axis(1)
axis(2, at=1:nrow(Wrong),
     labels=Wrong$Public.ID, las=2, cex.axis=1.2)
mtext("Days", side=1, line=2, cex=1.1)
box()
legend(46, 6.8, pch=c(17, 19), legend=c("ichorCNA", "CT Scan"), bty="n", ncol=2, col=c("grey", "red"))
legend(46, 6.1, pch=c(17, 17), col=c("red", "blue"),
       legend=c("Threshold model concordant", "Threshold model discordant"),
       bty="n")
legend(41, 4.55,
       legend="+/- CA15-3 positive/negative",
       bty="n")
legend(46, 4.25, pch=17, legend="", bty="n")
legend(46, 3.25, pch=17, legend="", bty="n")
legend(41, 3.0,
       legend="+/- NGTAS positive/negative",
       bty="n")
dev.off()

## Decisions according to treatment

types <- read.table("tmp_treats.txt", header=T, sep="\t")

lastCT$Treat <- sapply(strsplit(lastCT$Treatment, "/"), function(x) x[1])
lastCT <- merge(lastCT, types, by.y="Treatment", by.x="Treat", all.x=T)

lastCT$Group <- interaction(lastCT$ER.status, lastCT$Her2.status)
lastCT$Group <- factor(lastCT$Group, levels=levels(lastCT$Group),
                       labels=c("ER-/Her2-",
                                "ER+/Her2-", "ER-/Her2+",
                                "ER+/Her2+"))

pdf("Figure6.pdf", width=10, height=8)
par(mfrow=c(2,4), mar=c(2, 3, 4, 2), oma=c(4, 4, 2, 3))
K <- subset(lastCT, Type=="chemotherapy ")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Type=="chemotherapy "))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Type=="chemotherapy "), col=c("olivedrab", "grey"),
        main="Chemotherapy",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)
K <- subset(lastCT, Type=="endocrine based ")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Type=="endocrine based "))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Type=="endocrine based "), col=c("olivedrab", "grey"), main="Endocrine-based",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)
K <- subset(lastCT, Type=="targeted therapy ")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Type=="targeted therapy "))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Type=="targeted therapy "), col=c("olivedrab", "grey"), main="Targeted therapy",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)
K <- subset(lastCT, Type=="chemotherapy + targeted")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Type=="chemotherapy + targeted"))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Type=="chemotherapy + targeted"), col=c("olivedrab", "grey"),
        main="Targeted+Chemo",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)

K <- subset(lastCT, Group=="ER-/Her2-")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Group=="ER-/Her2-"))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Group=="ER-/Her2-"), col=c("lightblue", "indianred"), main="ER-/Her2-",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)
K <- subset(lastCT, Group=="ER-/Her2+")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Group=="ER-/Her2+"))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Group=="ER-/Her2+"), col=c("lightblue", "indianred"), main="ER-/Her2+",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
K <- subset(lastCT, Group=="ER+/Her2-")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Group=="ER+/Her2-"))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Group=="ER+/Her2-"), col=c("lightblue", "indianred"), main="ER+/Her2-",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)
K <- subset(lastCT, Group=="ER+/Her2+")
n1 <- length(unique(K$Public.ID))
n2 <- nrow(K)
pval <- round(t.test(mean ~ Progression, data=subset(lastCT, Group=="ER+/Her2+"))$p.val, 3)
pval <- ifelse (pval < 0.001, "p-value<0.001", paste0("p-value=", pval))
boxplot(mean ~ Progression, data=subset(lastCT, Group=="ER+/Her2+"), col=c("lightblue", "indianred"), main="ER+/Her2+",
        xlab="CT Scan Progression", ylab="Probability of progression", cex.main=1.8, cex.lab=1.2, cex.axis=1.2)
mtext(paste0("(n=",n1, "/", n2, ")"), side=1, line=2)
mtext(pval, side=3, line=0)
mtext("Probability of progression", side=2, outer=T, line=1)
mtext("CT Scan progression", side=1, outer=T, line=2)
dev.off()
