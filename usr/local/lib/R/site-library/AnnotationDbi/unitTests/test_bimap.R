require(org.Hs.eg.db)
require(RUnit)
## map is just to represent a classic Bimap
map  <- org.Hs.egSYMBOL
## map2 represents an AnnotationDbMap mapping made by some other process for BC
## map2 <- new("AnnotationDbMap", AnnotDb=org.Hs.eg.db, columns="SYMBOL")
map3 <- AnnotationDbi:::makeFlatBimapUsingSelect(org.Hs.eg.db, col="SYMBOL")
## And another map object expected by most other tests.
map2 <- AnnotationDbi:::makeFlatBimapUsingSelect(org.Hs.eg.db, col="ONTOLOGY")

##map3 <- AnnotationDbi:::flatten(map)
# molecular function, cellular component, biological process
bioconcepts <- c("MF", "CC", "BP")


## test ls 
test_ls <- function(){
  res <- ls(map)
  checkTrue(is.character(res))
  checkTrue(length(res) > 0)
  checkEquals(c("1","10","100"), head(res,n=3))
  res2 <- ls(map3)
  checkTrue(is.character(res2))
  checkTrue(length(res2) > 0)
  checkEquals(c("1","10","100"), head(res2,n=3))
}

## test Lkeys and Rkeys
test_Lkeys <- function(){
lres <- Lkeys(map)
checkTrue(length(lres) >0)
checkTrue(!any(is.na(lres)))
checkTrue(length(lres) == length(unique(lres)))

lres2 <- Lkeys(map2)  
checkTrue(length(lres2) >0)
checkTrue(!any(is.na(lres2)))
checkTrue(length(lres2) == length(unique(lres2)))
}

test_Rkeys <- function(){
rres <- Rkeys(map)
checkTrue(length(rres) >0)
checkTrue(!any(is.na(rres)))
checkTrue(length(rres) == length(unique(rres)))

rres2 <- Rkeys(map2)
checkTrue(length(rres2) >0)
checkTrue(!any(is.na(rres2)))
checkTrue(length(rres2) == length(unique(rres2)))
}

## test revmap (add a test now that it seems to work...
test_revmap <- function(){
  rmap <- revmap(map)
  checkTrue(rmap@direction == -1)
  rmap2 <- revmap(map2)
  checkTrue(rmap2@direction == -1)
}

## test mget
test_mget <- function(){
  k <- c("1","2")
  res <- mget(k, map)
  checkEquals(names(res), k)
  checkEquals(res[[1]], "A1BG")
  checkTrue(length(res)==length(k))
  
  res2 <- mget(k, map2) 
  checkEquals(names(res2), k)
  checkTrue(all(res2[[1]] %in% bioconcepts))
  checkTrue(length(res2)==length(k))

  ## reverse test 
  kr <- c("CC","MF")
  res3 <- mget(kr, revmap(map2))
  checkEquals(names(res3), kr)
  checkEquals(res3[[1]][1], "1")
  checkTrue(length(res3)==length(kr))
}

## test as.list
test_as.list <- function(){
  res <- as.list(map)
  checkEquals(names(res)[1], "1")
  checkEquals(res[[1]][1], "A1BG")
  checkTrue(length(res)>1000)
  
  res2 <- as.list(map2)
  checkEquals(names(res2)[[1]], "1")
  checkTrue(all(res2[[1]] %in% bioconcepts))
  checkTrue(length(res2)>1000)

  ## reverse test  
  res3 <- as.list(revmap(map2))
  checkTrue(all(names(res3) %in% bioconcepts))
  checkEquals(res3[[1]][1], "1")
  checkIdentical(length(res3), 3L)
}

## test as.character
test_as.character <- function(){
  res <- as.character(map)
  checkEquals(names(res)[1], "1")
  checkEquals(res[[1]][1], "A1BG")

  res2 <- as.character(map2)       
  checkEquals(names(res2[1:3]), rep("1", 3))
  checkTrue(all(res2[names(res2) == "1"] %in% bioconcepts))
  
  ## reverse test
  res3 <- as.character(revmap(map2)) 
  checkTrue(all(names(res3[res3 == "1"]) %in% bioconcepts))
  checkIdentical(unname(res3[1:3]), rep("1", 3))
}

## test eapply
test_eapply <- function(){
  res <- eapply(map, length)
  checkEquals(names(res)[1], "1")
  checkTrue(res[[1]][1] == 1)

  res2 <- eapply(map2, length)
  checkEquals(names(res2)[1], "1")
  checkTrue(res2[[1]][1] == 3)
}

## test get
test_get <- function(){
  k <- "1"
  res <- get(k, map)
  checkTrue(res == "A1BG")
  
  res2 <- get(k, map2)
  checkTrue(all(res2 %in% bioconcepts))

  ## reverse test 
  kr <- "CC"
  res3 <- get(kr, revmap(map2))
  checkTrue(res3[[1]][1] == "1")
}

## test exists
test_exists <- function(){
  checkTrue(exists("2", map) == TRUE)     
  checkTrue(exists("titi", map) == FALSE)
  
  checkTrue(exists("9", map2) == TRUE)  
  checkTrue(exists("titi", map2) == FALSE)  
}

## test "[["
test_dblBrackets <- function(){
  res <- map[["1"]]
  checkTrue(res == "A1BG")
  res2 <- map2[["1"]]
  checkTrue(all(res2 %in% bioconcepts))
}

test_head <- function(){
  res <- head(map, n=3)
  checkTrue( class(res) == "AnnDbBimap" )

  res2 <- head(map2, n=3)  ## implement Lkeys and Rkeys
  checkTrue( class(res2) == "data.frame" )
  checkTrue( dim(res2)[1] == 3 )
  checkTrue( dim(res2)[2] == 2 )
}

test_tail <- function(){
  res <- tail(map, n=3)
  checkTrue( class(res) == "AnnDbBimap" )

  res2 <- tail(map2, n=3)
  checkTrue( class(res2) == "data.frame" )
  checkTrue( dim(res2)[1] == 3 )
  checkTrue( dim(res2)[2] == 2 )
}

## test "$"
test_Dollar <- function(){
  res <- map$"1"
  checkTrue(res == "A1BG")
  res2 <- map2$"1"
  checkTrue(all(res2 %in% bioconcepts))
}

## test toTable as.data.frame
test_toTable <- function(){
  res <- toTable(map)
  resdf <- as.data.frame(map)
  checkEquals(res, resdf)
  checkEquals(colnames(res), c("gene_id","symbol"))
  checkTrue(res[1,1]==1)
  checkTrue(res[1,2]=="A1BG")
  
  ## So one potential issue I have is that I get the "wrong" sort of headings?
  ## this is largely a cosmetic issue though...
  res2 <- toTable(map2)
  resdf2 <- as.data.frame(map2)
  checkEquals(res2, resdf2)
  checkEquals(colnames(res2), c("ENTREZID","ONTOLOGY"))
  checkTrue(res2[1,1]==1)
  checkTrue(
      all(res2[res2[["ENTREZID"]]==1, "ONTOLOGY"] %in% bioconcepts)
  )
}

test_sample <- function(){
  res <- sample(map,size=2)
  checkTrue(length(res)==2)
  checkTrue(class(res)=="list")
  
  res2 <- sample(map2,size=2)
  checkTrue(length(res2)==2)
  checkTrue(class(res2)=="list")
}

