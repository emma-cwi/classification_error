library(ggplot2)
library(stringr)

############### MAKE DATA #####################

threshold_base <- c((1:10)/10) %>% as.character
class_base <- c('C0','C1')
output_base <- c('selected','discarded')
tp=c((10:1)*10)
fn=c(100-(10:1)*10)*-1
tn=c((1:10)*20) *-1
fp=c(200-(1:10)*20)

nbin <- length(threshold_base)

threshold <- class <- output <- nitem <- c()
for(i in 1:nbin){
  for(c in class_base){
    for(o in output_base){
      threshold <- c(threshold,threshold_base[i])
      class <- c(class, c)
      if(c=='C0'){
        if(o=='selected') {
          output <- c(output, 'FP')
          nitem <- c(nitem,fp[i])
        }
        else {
          output <- c(output, 'TN')
          nitem <- c(nitem,tn[i])
        }
      } else {
        if(o=='selected') {
          output <- c(output, 'TP')
          nitem <- c(nitem,tp[i])
        }
        else {
          output <- c(output, 'FN')
          nitem <- c(nitem,fn[i])
        }
      }  
    }    
  }
}

data <- data.frame(cbind(threshold=as.character(threshold), class, output, nitem), stringsAsFactors=FALSE)
data$nitem <- as.numeric(data$nitem)

ggplot(data, aes(x=threshold, y=nitem, fill=output, group=class)) + 
  geom_bar(stat='identity', position='dodge', color='#ffffff') +
  scale_fill_manual(values=c("#999999", "#078ac6","#c46362", "#333333")) +
  theme_minimal() 


############### OLD #####################
threshold_base <- c((1:10)/10) %>% as.character
class_base <- c('C0','C1')
tp=c((10:1)*10)
fn=c(100-(10:1)*10)*-1
tn=c((1:10)*20) *-1
fp=c(200-(1:10)*20)

nbin <- length(threshold_base)

threshold <- class <- selected <- discarded <- c()
for(i in 1:nbin){
  for(c in class_base){
    threshold <- c(threshold,threshold_base[i])
    class <- c(class, c)
    if(c=='C0'){
      selected <- c(selected,fp[i])
      discarded <- c(discarded,tn[i])
    } else {
      selected <- c(selected,tp[i])
      discarded <- c(discarded,fn[i])
    }      
  }
}

data <- data.frame(cbind(threshold=as.character(threshold), class, selected, discarded), stringsAsFactors=FALSE)
data$selected <- as.numeric(data$selected)
data$discarded <- as.numeric(data$discarded)

up <- ggplot(data, aes(x=threshold, y=selected, fill=class)) + 
  geom_bar(stat='identity',position='dodge') +
  scale_fill_manual(values=c("#333333", "#078ac6")) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0))

down <- ggplot(data, aes(x=threshold, y=discarded, fill=class)) + 
  geom_bar(stat='identity',position='dodge') +
  scale_fill_manual(values=c("#cccccc", "#c46362")) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0))



############### OLD #####################

cm <- matrix(c(90,10,10,
               5,85,10,
               5,5,80),3,3)

data_tbl <- data.frame(cbind( threshold=c((1:10)/10),
                              tp=c((10:1)*10),
                              fn=c(100-(10:1)*10),
                              tn=c((1:10)*20),
                              fp=c(200-(1:10)*20) ))


data$fn <- -1 * data$fn
data$tn <- -1 * data$tn                   

barplot(t(cbind(data$fn,data$tp)))


ggplot(data, aes(x=threshold, y=fn)) + geom_bar(stat="identity")
