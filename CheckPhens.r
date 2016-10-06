nrs = c(0, 1, 2, 3)

par(mfrow = c(1, length(nrs)))
for(i in nrs){
    file = paste("pop", i, ".phens.txt", sep = "")
    data = read.delim(file, header = F, sep = " ")
    data = data[, !is.na(colSums(data))]
    avg = colMeans(data)
    avg = matrix(avg, 7, 7)
    image(avg[, 7:1], col = rev(grey.colors(100)))
    }

    
    file = paste("pop", 3, ".phens.txt", sep = "")
    data = read.delim(file, header = F, sep = " ")
    data = data[, !is.na(colSums(data))]
    par(mfrow = c(3, 3))
    for(i in 1:25){
        avg = as.numeric(data[2,])
        avg = matrix(avg, 7, 7)
        image(avg[, 7:1], col = rev(grey.colors(100)))
        }
        
F1s = read.delim("pop3_F1genotypes.txt", header = F, sep = "")
BCs = read.delim("pop3_BCgenotypes.txt", header = F, sep = "")     
    

dim(BCs)    


## list all existing alleles in parental genomes (*looking at F1s)
allgenes1 = NULL
allgenes2 = NULL
for(i in 1:100){
    allgenes1 = rbind(allgenes1, matrix(unlist(F1s[i, ]), 7, 7))
    allgenes2 = rbind(allgenes2, matrix(unlist(F1s[i + 100, ]), 7, 7))
    }
allgenes1 = unique(allgenes1)
allgenes2 = unique(allgenes2)

## assign parental origin in BCs genomes
ddd = function(x, y) abs(sum(x-y))

getdelta(x, mat){
    apply(apply(df1, 2 ,`==`, df2), 1, any)
    }
    
    
orig1 = NULL
orig2 = NULL
for(i in 1:100){
    w1_BC = matrix(unlist(BCs[1, ]), 7, 7)
    w2_BC = matrix(unlist(BCs[101, ]), 7, 7)

    apply(w1_BC, 1, )
    }

test = NULL
for(loc in 1:7){
    a = max(rowSums(unlist(w1_BC[loc,]) == allgenes1))
    b = max(rowSums(unlist(w1_BC[loc,]) == allgenes2))
    test = rbind(test, c(a, b))
    }
    
    
    