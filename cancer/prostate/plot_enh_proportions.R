prostate.lncapRT <- data.frame(quintile=1:5

                               , gains=c(11568, 24645, 33907, 42118, 47822)

                               , losses=c(2926, 6784, 11805, 17631, 29666)

                               , unchanged=c(675, 2723, 6428, 11457, 19572))


prostate.prec <- data.frame(quintile=1:5

                               , gains=c(17472, 27396, 33938, 40482, 40911)

                               , losses=c(586, 4002, 10503, 18680, 35043)

                               , unchanged=c(215, 1726, 5482, 11606, 21796))
                               
                               Fisher.rt.quintiles <- function(x, test_enh){

  # li <- list()

  for(i in 1:4){

    m <- as.matrix(x[i:(i+1), c(test_enh, "unchanged")])

    # print(m)

    out <- fisher.test(m, alternative = "greater")

    print(paste("test:", paste0("Q", i), "vs", paste0("Q", i+1)))

    print(paste("odds.ratio", out$estimate))

    print(paste("p.value", out$p.value))

  }

}


#PROSTATE CANCER - LNCAP REPLICATION TIME

#gains

Fisher.rt.quintiles(prostate.lncapRT, "gains")


# [1] "test: Q1 vs Q2"

# [1] "odds.ratio 1.89350586187734"

# [1] "p.value 7.6884899766226e-52"

# [1] "test: Q2 vs Q3"

# [1] "odds.ratio 1.71578467766024"

# [1] "p.value 1.786992221817e-114"

# [1] "test: Q3 vs Q4"

# [1] "odds.ratio 1.4348711141922"

# [1] "p.value 7.20080428488543e-100"

# [1] "test: Q4 vs Q5"

# [1] "odds.ratio 1.50452217699273"

# [1] "p.value 7.74954217799775e-204"


#losses

Fisher.rt.quintiles(prostate.lncapRT, "losses")


# [1] "test: Q1 vs Q2"

# [1] "odds.ratio 1.73985929105627"

# [1] "p.value 2.4579099953357e-32"

# [1] "test: Q2 vs Q3"

# [1] "odds.ratio 1.35656823283948"

# [1] "p.value 2.78348269983041e-29"

# [1] "test: Q3 vs Q4"

# [1] "odds.ratio 1.19338758043738"

# [1] "p.value 8.59180006147813e-20"

# [1] "test: Q4 vs Q5"

# [1] "odds.ratio 1.01527349977971"

# [1] "p.value 0.160008929854532"


#PROSTATE CANCER - PREC REPLICATION TIME


#gains

Fisher.rt.quintiles(prostate.prec, "gains")

# [1] "test: Q1 vs Q2"

# [1] "odds.ratio 5.11970777303141"

# [1] "p.value 7.62916359150769e-161"

# [1] "test: Q2 vs Q3"

# [1] "odds.ratio 2.56368304622898"

# [1] "p.value 1.05716124545566e-263"

# [1] "test: Q3 vs Q4"

# [1] "odds.ratio 1.77485040167123"

# [1] "p.value 6.27556539921969e-233"

# [1] "test: Q4 vs Q5"

# [1] "odds.ratio 1.85831445658451"

# [1] "p.value 0"

#losses

Fisher.rt.quintiles(prostate.prec, "losses")

# [1] "test: Q1 vs Q2"

# [1] "odds.ratio 1.17550004001553"

# [1] "p.value 0.0301374452635817"

# [1] "test: Q2 vs Q3"

# [1] "odds.ratio 1.21020200957253"

# [1] "p.value 4.39397035402631e-09"

# [1] "test: Q3 vs Q4"

# [1] "odds.ratio 1.19035834003938"

# [1] "p.value 6.38310325494893e-18"

# [1] "test: Q4 vs Q5"

# [1] "odds.ratio 1.00108239719068"

# [1] "p.value 0.4735124353159"


