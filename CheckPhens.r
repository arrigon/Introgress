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
