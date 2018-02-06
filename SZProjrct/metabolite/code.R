### input data
library(randomForest)
setwd("C:/Users/zoumanli/Desktop")
data_pos = read.csv("m_pos_pick.csv")
data_pos = data_pos[,-1]
data_pos$state2 <- as.factor(datamrf$state2)

### Random Forest
set.seed(1234)
metabolism_rf <- randomForest(state2~., data=data_pos, ntree= 2000, proximity=TRUE)
table(predict(metabolism_rf),data_pos$state2)
print(metabolism_rf)



### input data
require("vegan")
data_neg = read.csv("m_pos_pick.csv")
data_posv = data_pos[ ,3:ncol(data_pos)]
data_pos$state2 = as.numeric(data_pos$state2)
ncol(data_posv)

### permanova
set.seed(1234)
adonis(data_posv~state2, data=data_pos, permutations=1000,methods="bray")
