library(GGally)
library(bnlearn)
library(Rgraphviz)
library(penalized)
library(visNetwork)
library(dplyr)

"plot_network" = function(dag,strength_df=NULL,undirected=FALSE,
                          group=NA,title=NULL,height=NULL,width=NULL)
    {
    edge_size = ifelse(is.null(strength_df),NA,
                   round(right_join(strength_df, data.frame(dag$arcs[,c(1,2)]))$strength))
    
    nodes = names(dag$nodes)
    nodes = data.frame(id   = nodes,
                       label= nodes,
                       size = 16,
                       font.size= 18,
                       shadow   = TRUE,
                       group    = group)
    
    edges = data.frame(from   = dag$arcs[,1],
                       to     = dag$arcs[,2],
                       value  = edge_size,
                       arrows = list(to=list(enabled=TRUE,scaleFactor=.5)),
                       shadow = TRUE)
    
    if(is.na(group[1]))     nodes = nodes[,-6] # without group
    if(is.na(edge_size)) edges = edges[,-3] # without edge_size
    if(undirected)       edges$arrows.to.enabled=FALSE
    
    network=visNetwork(nodes,edges,main=title,height=height, width=width)%>% 
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
    return(network)
}

#----------------------------------------------------------
# VISUALIZATION 

diff=read.csv("dental.csv") %>% 
    mutate_if(is.integer, as.numeric)
# CoGo, dGoPg and dCoA & Treatment, dANB and dCoA
diff_delta = cbind(Treatment=diff$Treatment,
                   Growth=diff$Growth,
                   sapply(diff[, 1:6], function(x) x / diff$dT))

ggcorr(diff_delta,palette="RdBu",label=TRUE,label_round=3)
# ggsave("images/dental-cor.png",width=7,height=7)

# pairs
ggpairs(diff[,c("dANB", "dPPPM", "dIMPA", "dCoA", "dGoPg", "dCoGo")])
# ggsave("images/dental-pairs.png",width=7,height=7)

rho = cor(data.frame(diff_delta))
undirected_graph = empty.graph(colnames(rho))
amat(undirected_graph) = (rho > 0.4) + 0L - diag(1L, nrow(rho))


par(mfrow=c(1,1))
graphviz.plot(undirected_graph, layout = "fdp", shape = "ellipse")
plot_network(undirected_graph,title="상관관계",undirected=TRUE)

#----------------------------------------------------------
# static Bayesian network as a difference model
#----------------------------------------------------------
# step1: black list와 white list 정하기
# DAG에 포함되지 말아야할 엣지들
black_list = tiers2blacklist(list("dT", "Treatment", "Growth",
                                  c("dANB", "dPPPM", "dIMPA", "dCoA", "dGoPg", "dCoGo")))
black_list = rbind(black_list, c("dT", "Treatment"))
black_list

# DAG에 포함되어야할 엣지들
white_list=matrix(c("dANB", "dIMPA","dPPPM", "dIMPA","dT", "Growth"),
                  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))
white_list

#----------------------------------------------------------
# step2: DAG 만들기
# 첫번째 DAG
#  hill-climbing (HC) search for finding the structure

DAG = hc(diff, whitelist = white_list, blacklist = black_list)
DAG2 = gs(diff, whitelist = white_list, blacklist = black_list)

group = ifelse(names(DAG$nodes)%in%c("Treatment","Growth"),2,1)
plot_network(DAG,group=group,title="Hill-Climbing")
plot_network(DAG2,group=group,title="Markov-Blanket")

#----------------------------------------------------------
# step3: DAG 튜닝하기
# Bootstrap을 통한 DAG 조정: 하는 이유? 효과?
strength_diff = boot.strength(diff, R = 200, algorithm = "hc",
                              algorithm.args = list(whitelist = white_list, blacklist = black_list))

strength_diff2 = boot.strength(diff, R = 200, algorithm = "gs",
                              algorithm.args = list(whitelist = white_list, blacklist = black_list))
head(strength_diff)
attr(strength_diff, "threshold")
average_diff = averaged.network(strength_diff)
average_diff2 = averaged.network(strength_diff2)

strength.plot(average_diff, strength_diff, shape = "ellipse", highlight = list(arcs = white_list))

plot_network(average_diff,strength_diff, group=group) 
plot_network(average_diff2,strength_diff2, group=group) 

# 하나 했을 때랑 비교
par(mfrow = c(1, 2))
graphviz.compare(average_diff, DAG, shape = "ellipse", main = c("averaged DAG", "single DAG"))
graphviz.compare(average_diff2, DAG2, shape = "ellipse", main = c("averaged DAG", "single DAG"))
compare(average_diff, DAG,arcs=TRUE)

# cpdag: find the equivalent class
undirected.arcs(cpdag(average_diff, wlbl = TRUE))
undirected.arcs(cpdag(DAG, wlbl = TRUE))
compare(cpdag(average_diff, wlbl = TRUE), cpdag(DAG, wlbl = TRUE))

#----------------------------------------------------------
# Step4. Pruning
plot(strength_diff)
plot(strength_diff2)
abline(v = 0.75, col = "tomato", lty = 2, lwd = 2)
abline(v = 0.85, col = "steelblue", lty = 2, lwd = 2)


nrow(strength_diff[strength_diff$strength > attr(strength_diff, "threshold") &strength_diff$direction > 0.5, ])
nrow(strength_diff[strength_diff$strength > 0.75 & strength_diff$direction > 0.5, ])
nrow(strength_diff[strength_diff$strength > 0.85 & strength_diff$direction > 0.5, ])

average_diff_simple = averaged.network(strength_diff, threshold = 0.85)
plot_network(average_diff_simple,strength_diff, group=group)
strength.plot(average_diff_simple, strength_diff, shape = "ellipse", highlight = list(arcs = white_list))


#----------------------------------------------------------
# Step 5: learning parameters
fit = bn.fit(average_diff_simple, diff)

# same result
fit$dANB
lm(dANB~Growth+Treatment, diff)

# 만약 coefficient를 조정하고 싶다면 penalize를 할 수 있음
fitted_new$dANB = penalized(diff$dANB, penalized = diff[, parents(average_diff_simple, "dANB")],
                            lambda2 = 20, model = "linear", trace = FALSE)
fitted_new$dANB

#----------------------------------------------------------
# Step 6: Model validation
# loss: Posterior Predictive Correlation 

# 10-fold cross validation, # of runs = 10 총 100개
xval = bn.cv(diff, bn = "hc", 
             algorithm.args = list(blacklist = black_list, whitelist = white_list),
             loss = "cor-lw", 
             loss.args = list(target = "Growth", n = 500), 
             runs = 10)
err = numeric(10)
for (i in 1:10){
    tt = table(unlist(sapply(xval[[i]], '[[', "observed")),
               unlist(sapply(xval[[i]], '[[', "predicted")) > 0.50)
    err[i] = (sum(tt) - sum(diag(tt))) / sum(tt)
}

err

xval[[1]] # first run
xval[[1]][[1]] # first run, first fold
xval[[1]][[1]]$observed
xval[[1]][[1]]$predicted # probability


round(err,3)

# for numeric variable
predcor = structure(numeric(6),names = c("dCoGo", "dGoPg", "dIMPA", "dCoA", "dPPPM", "dANB"))

for (var in names(predcor)) {
    xval = bn.cv(diff, bn = "hc", algorithm.args = list(blacklist = black_list, whitelist = white_list),
                 loss = "cor-lw", loss.args = list(target = var, n = 200), runs = 10)
    predcor[var] = mean(sapply(xval, function(x) attr(x, "mean")))
}

round(predcor, digits = 3)
mean(predcor)

#----------------------------------------------------------
# Knowledge
# cpdist generates random samples conditional on 
# the evidence using the method specified in the method argument.

#1. "An excessive growth of CoGo should induce a reduction in PPPM."
sim = cpdist(fit, nodes = c("dCoGo", "dPPPM"), n = 10^4,evidence = (Treatment <0.5))
nrow(sim)
plot(sim, col = "grey")
abline(v = 0, col = 2, lty = 2, lwd = 2)
abline(h = 0, col = 2, lty = 2, lwd = 2)
abline(coef(lm(dPPPM ~ dCoGo, data = sim)), lwd = 2)

#2. "A small growth of CoGo should induce an increase in PPPM."
nrow(sim[(sim$dCoGo <= 0) & (sim$dPPPM > 0), ]) / nrow(sim[(sim$dCoGo <= 0), ])

# a small growth of CoGo (dCoGo ∈ [0, 2]) unfortunately dPPPM ⋜ 0 with probability ≈ 0.50
nrow(sim[(sim$dCoGo > 0) & (sim$dCoGo < 2) & (sim$dPPPM > 0), ]) /
    nrow(sim[(sim$dCoGo) > 0 & (sim$dCoGo < 2),  ])

# 3. "If ANB decreases, IMPA decreases to compensate."
sim = cpdist(fit, nodes = c("dIMPA", "dANB"), n = 10^4,
             evidence = (Treatment < 0.5))
plot(sim, col = "grey")
abline(v = 0, col = 2, lty = 2, lwd = 2)
abline(h = 0, col = 2, lty = 2, lwd = 2)
abline(coef(lm(dIMPA ~ dANB, data = sim)), lwd = 2)

#4. "If `GoPg` increases strongly, then both ANB and IMPA decrease." 
# If we simulate dGoPg, dANB and dIMPA from the BN 
# assuming dGoPg > 5 (i.e. GoPg is increasing) 
# we estimate the probability that dANB > 0 (i.e. ANB is increasing) at ≈ 0.70 and 
# that dIMPA < 0 at only ≈ 0.58.
sim = cpdist(fit, nodes = c("dGoPg", "dANB", "dIMPA"), n = 10^4,
             evidence = (dGoPg > 5) & (Treatment < 0.5))
nrow(sim[(sim$dGoPg > 5) & (sim$dANB < 0), ]) / nrow(sim[(sim$dGoPg > 5), ])
nrow(sim[(sim$dGoPg > 5) & (sim$dIMPA < 0), ]) / nrow(sim[(sim$dGoPg > 5), ])

# 5. "Therapy attempts to stop the decrease of ANB. 
# If we fix ANB is there any difference treated and untreated patients?"
sim = cpdist(fit, nodes = c("Treatment", "Growth"), n = 5 * 10^4,
         evidence = abs(dANB) < 0.01)
tab = table(TREATMENT = sim$Treatment < 0.5, GOOD.GROWTH = sim$Growth > 0.5)
round(prop.table(tab, margin = 1), 2)

#6. If we fix GoPg, is there any difference in dANB of treated and untreated patients?
sim6 = cpdist(fit, nodes = c("Treatment", "dANB", "dGoPg"),
              n=10^5,
              evidence = abs(dGoPg) < 0.1)

sim6$Treatment = ifelse(sim6$Treatment>0.5, 1, 0)
boxplot(dANB ~ Treatment, data = sim6)
