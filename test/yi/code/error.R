library(ggplot2)
for(j in 1:length(allstepsizes)){
  L <- gradient_descent(Y,r,alpha_bnd,allstepsizes[j],1,100,matrix(rep(1,nrow(Y)*ncol(Y)),nrow=nrow(Y)))[1]
  alle1[j] <- gradient_descent(Y,r,alpha_bnd,allstepsizes[j],1,100,matrix(rep(1,nrow(Y)*ncol(Y)),nrow=nrow(Y)))[2]
  Lp <- gradient_descent( Y1,r,alpha_bnd,allstepsizes[j]/p,1,100,sparsity)[1]
  alle2[j] <- gradient_descent( Y1,r,alpha_bnd,allstepsizes[j]/p,1,100,sparsity)[2]
}



ee1=alle1

xaixs <- c(1:length(ee1)-1)
i <- 1
η_0.05 <- ee1[i][2:length(ee1[i])]
i=i+1
η_0.1 <- ee1[i][2:length(ee1[i])]
i <- i+1
η_0.2 <- ee1[i][2:length(ee1[i])]
i <- i+1
η_0.4 <- ee1[i][2:length(ee1[i])]
i <- i+1
η_0.7 <- ee1[i][2:length(ee1[i])]
i <- i+1
η_1 <- ee1[i][2:length(ee1[i])]
i <- i+1
η_1.5 <- ee1[i][2:length(ee1[i])]
i <- i+1
η_2.5 <- ee1[i][2:length(ee1[i])]
df <- data.frame(xaixs,η0.05,η0.1,η0.2,η0.4,η0.7,η1,η1.5,η2.5)
df2 <- melt(data=df,id.var="xaixs")

ggplot(data = df2, aes(x = xaxis, y = value, colour = variable)) + geom_line()


h_legend <- legend('\eta=0.05','\eta=0.1','\eta=0.2','\eta=0.4','\eta=0.7','\eta=1','\eta=1.5','\eta=2.5','Location','NorthEast')
set(h_legend,'FontSize',12)
title(xlab="Iterations", ylab="Error")


