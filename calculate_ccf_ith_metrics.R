library(stringr)
library(GenomicRanges)
library(data.table)
library(caTools)
library(KernSmooth)
library(bezier)
library(dplyr)

path2files <- "example_files/" #directory with files for analysis
records = read.delim(paste0(path2files,"region_annotation.txt")) #dataframe with annotation of each region
full_data <- read.delim(paste0(path2files,"SRC331_SNV_info.txt")) #dataframe with sequencing information all SNVs
clust = read.delim(paste0(path2files,"TitanCNA_purities.txt")) #dataframe with purity for each region
first = T

weightAFs <- function(af.data, mafs, depths, minAF=0) { 
  keep = vector()
  for (i in 1:length(mafs)) {
    cmaf = mafs[i]
    cmafi = match( cmaf, colnames(af.data) )
    if (length(keep) == 0) {
      keep = af.data[,cmafi] > minAF
    } else {
      keep = keep | (af.data[,cmafi] > minAF)
    }
  }
  af.data = af.data[keep,]
  tdepth = rep(0,dim(af.data)[1])
  talt = rep(0,dim(af.data)[1])
  for (i in 1:length(mafs)) {
    cmaf = mafs[i]
    cmafi = match(cmaf,colnames(af.data))
    cdep = depths[i]
    cdepi = match(cdep,colnames(af.data))
    tdepth = tdepth + af.data[,cdepi]
    talt = talt + af.data[,cdepi]*af.data[,cmafi]
  }
  tmaf = talt/tdepth
  return(tmaf)
}

rAUC <- function(data, lower, mafs, depths) {
  if (length(mafs) <= 8) {
    lower = round((0.10/length(mafs)),2)
  } else {
    lower = max(0.005, round((0.10/length(mafs)),3))
  }
  weightAF = weightAFs(data, mafs, depths)
  weightAF = weightAF[which(weightAF >= lower & weightAF <= 0.39)]
  vafs = seq(lower,0.39,by=0.01)
  nstep = length(vafs)
  counts = vector()
  ncounts = ((1/vafs)-(1/0.39))/((1/lower)-(1/0.39))
  
  for ( i in 1:length(vafs) ) {
    counts = append(counts, length(which(weightAF > vafs[i]))-length(which(weightAF > vafs[nstep])))
  }
  counts = counts/counts[1]
  to_bei = cbind(vafs, counts)
  t = seq(0,1,length=length(counts))
  bei = bezier(t,to_bei)
  AUC = trapz(bei[,1],bei[,2])
  nAUC = trapz(vafs,ncounts)
  rAUC = AUC/nAUC
  rAUC = round(rAUC,8)
  return(list(rAUC=rAUC,weightAF=weightAF))
}

fst.hudson <- function(af, minAF=0.08) {  #combine for multi SSNV sites using ratio of averages
  mafis = which(grepl("maf", colnames(af)))
  keep = as.vector(apply(af, 1, function(x, mafis) {
    maxmaf = max(as.numeric(x[mafis]))
    if (maxmaf > minAF) {
      TRUE
    } else {
      FALSE
    }
  }, mafis=mafis))     #filter data
  af = af[keep,]
  
  Ns = c()
  Ds = c()
  for(k in 1:nrow(af)) {
    n1 = af$depth1[k]
    n2 = af$depth2[k]
    p1 = af$maf1[k]
    p2 = af$maf2[k]
    N = (p1-p2)^2-(p1*(1-p1))/(n1-1)-(p2*(1-p2))/(n2-1)     # covariance
    D = p1*(1-p2)+p2*(1-p1)                                 # standard deviations
    Ns = c(Ns, N)
    Ds = c(Ds, D)
  }
  Fst.h = mean(Ns)/mean(Ds)
  return(Fst.h)
}

subclonalMut <- function(sampAB, regions, minAF=0.08, statsAF=0.08, highAF=0.2, ratio=1)  {      #determinine subclonal mutations
  fhSubv = c()
  fHrsv = c()
  KSDv = c()
  FSTv = c()
  m = 0
  fHrs_noregion = 0
  rAUCv = c()
  k = length(regions)
  n = k*(k-1)/2
  pobi = match("pob", colnames(sampAB))
  for (i in 1:k){
    #calculate fhSub based upon thresholds for fhSub 
    sn = regions[i]
    mafaAi = match(paste(sn, "mafa", sep=""), colnames(sampAB))
    subAi = which(grepl("private", sampAB[,pobi]) & sampAB[,mafaAi] > minAF)
    highAi = which( grepl("private", sampAB[,pobi]) & sampAB[,mafaAi] > highAF)
    fhSubv = append(fhSubv, length(highAi)/(length(subAi)))
  }
  
  for (i in 1:length(regions)){
    snA = regions[i]
    mafaAi = match(paste(snA, "mafa", sep=""), colnames(sampAB))
    nbAi = match(paste(snA, "nb", sep=""), colnames(sampAB))
    depthAi = match(paste(snA, "depth", sep=""), colnames(sampAB))
    if (i<k){
      for (j in (i+1):k){
        #pairwise metrics
        m = m+1
        snB = regions[j]
        mafaBi = match(paste(snB, "mafa", sep=""), colnames(sampAB))
        depthBi = match(paste(snB, "depth", sep=""), colnames(sampAB))
        nbBi = match(paste(snB, "nb", sep=""), colnames(sampAB))
        sampAB2 = sampAB[which((sampAB[,nbAi] > 0 & sampAB[,nbBi] > 0) | (sampAB[,nbAi] == 0 & sampAB[,nbBi] == 0)),]
        mafa1Index = match(paste(snA,"mafa",sep=""), colnames(sampAB))    #for mafa
        mafa2Index = match(paste(snB,"mafa",sep=""), colnames(sampAB))    #for mafa
        subAi = which( grepl("private", sampAB2[,pobi]) & sampAB2[,mafaAi] > minAF & ((sampAB2[,mafaBi] == 0 & (sampAB2[,nbBi] != 0 | sampAB2[,nbAi] == 0)) | sampAB2[,mafaBi] != 0) )
        subBi = which( grepl("private", sampAB2[,pobi]) & sampAB2[,mafaBi] > minAF & ((sampAB2[,mafaAi] == 0 & (sampAB2[,nbAi] != 0 | sampAB2[,nbBi] == 0)) | sampAB2[,mafaAi] != 0) )
        allSubRows = union(subAi,subBi)
        mafs = paste(c(snA,snB), "mafa", sep="")
        depths = paste(c(snA,snB), "depth", sep="")
        rAUCout = rAUC(sampAB2[allSubRows,], 0.04, mafs, depths)     #need to add for lower bound
        rAUCv = append(rAUCv,round(rAUCout$rAUC,8))
        #auc calculations
        depthPowerKeep <- as.vector(apply(sampAB2, 1, function(x,mafa1i,mafa2i,dp1i,dp2i) {
          if(as.numeric(x[mafa1i]) == 0){vaf = as.numeric(x[mafa2i])
          if (vaf > 1 | vaf < 0){FALSE}  else if (pbinom(0,as.numeric(x[dp1i]),vaf) < 0.05){TRUE} else {FALSE}}
          else if (as.numeric(x[mafa2i]) == 0){vaf = as.numeric(x[mafa1i])
          if (vaf > 1 | vaf < 0){FALSE}  else if (pbinom(0,as.numeric(x[dp2i]),vaf) < 0.05){TRUE} else {FALSE}}
          else {TRUE}
        }, mafa1i=mafa1Index,mafa2i=mafa2Index,dp1i=depthAi,dp2i=depthBi))
        sampAB3 = sampAB2[depthPowerKeep,]
        subAi = which( grepl("private", sampAB3[,pobi]) & sampAB3[,mafaAi] > minAF & ((sampAB3[,mafaBi] == 0 & (sampAB3[,nbBi] != 0 | sampAB3[,nbAi] == 0)) | sampAB3[,mafaBi] != 0) )
        mutsA = sampAB3[subAi,mafaAi]/ratio
        subBi = which( grepl("private", sampAB3[,pobi]) & sampAB3[,mafaBi] > minAF & ((sampAB3[,mafaAi] == 0 & (sampAB3[,nbAi] != 0 | sampAB3[,nbBi] == 0)) | sampAB3[,mafaAi] != 0) )
        mutsB = sampAB3[subBi,mafaBi]/ratio
        allSubRows = union(subAi,subBi)
        mutsASp2 = sampAB3[intersect(subAi, which( sampAB3[,mafaAi] > statsAF & sampAB3[,mafaBi] == 0)), mafaAi]
        mutsASph2 = sampAB3[intersect(subAi, which( sampAB3[,mafaAi] > highAF & sampAB3[,mafaBi] == 0)), mafaAi]
        mutsBSp2 = sampAB3[intersect(subBi, which( sampAB3[,mafaBi] > statsAF & sampAB3[,mafaAi] == 0)), mafaBi]
        mutsBSph2 = sampAB3[intersect(subBi, which( sampAB3[,mafaBi] > highAF & sampAB3[,mafaAi] == 0)), mafaBi]
        ratioHighA = length(mutsASph2)/length(mutsASp2)
        ratioHighB = length(mutsBSph2)/length(mutsBSp2)
        mutsA2 = sampAB3[intersect(subAi, which( sampAB3[,mafaAi] > statsAF )), mafaAi]
        mutsAh2 = sampAB3[intersect(subAi, which( sampAB3[,mafaAi] > highAF )), mafaAi]
        mutsB2 = sampAB3[intersect(subBi, which( sampAB3[,mafaBi] > statsAF )), mafaBi]
        mutsBh2 = sampAB3[intersect(subBi, which( sampAB3[,mafaBi] > highAF )), mafaBi]
        ratioHighA2 = length(mutsAh2)/length(mutsA2)
        ratioHighB2 = length(mutsBh2)/length(mutsB2)
        mutsSub = sampAB3[allSubRows,]
        mutsSub = data.frame( maf1 = mutsSub[,mafaAi], depth1=mutsSub[,depthAi], maf2 = mutsSub[,mafaBi], depth2=mutsSub[,depthBi] )
        R1 = nrow(sampAB[which(sampAB[,mafaAi] > .39 & sampAB[,mafaBi] < .05),])
        R2 = nrow(sampAB[which(sampAB[,mafaBi] > .39 & sampAB[,mafaAi] < .05),])
        C = nrow(sampAB[which(sampAB[,mafaAi] > .39 & sampAB[,mafaBi] > .39),])
        if (!(is.na(ratioHighA))&&!(is.na(ratioHighB))){
          fHrsv = append(fHrsv,(ratioHighA+ratioHighB)/2)
        }
        else if (!(is.na(ratioHighA))){
          fHrsv = append(fHrsv,(ratioHighA))
        }
        else if (!(is.na(ratioHighB))){
          fHrsv = append(fHrsv,(ratioHighB))
        }
        #remove any NAs prior to calculating fHrs
        if (length(mutsA[which(mutsA > statsAF)])!=0&&length(mutsB[which(mutsB > statsAF)]!=0)){
          kst = as.numeric(ks.test( mutsA[which(mutsA > statsAF)], mutsB[which(mutsB > statsAF)] )$statistic)
          KSDv = append(KSDv,round(kst,8))
          #calculate KSD
        }
        fstVal = round(fst.hudson(mutsSub),8)
        FSTv = append(FSTv, fstVal)
        #calculate FST
      }
    }
  }
  fhSub = mean(fhSubv, na.rm=TRUE)
  fHrs = mean(fHrsv, na.rm=TRUE)
  KSD = mean(KSDv, na.rm=TRUE)
  FST = mean(FSTv, na.rm=TRUE)
  rAUC = mean(rAUCv, na.rm=TRUE)
  metrics_df = data.frame(fhSub = fhSub, fHrs = fHrs, KSD=KSD, FST=FST, rAUC=rAUC)
  return (metrics_df)
}

smallerdf <- function(sampAB){
  #compresses dataframe into format for easy calculation of metrics
  unique = data.frame(chrpos = unique(sampAB$chrpos))
  regions = unique(sampAB$sn)
  for (i in regions){
    subsAB = sampAB[which(sampAB$sn==i),]
    headings = data.frame(subsAB$chrpos, subsAB$mafa, subsAB$nb, subsAB$TOTAL_DEPTH, subsAB$cpsType)
    colnames(headings)=c('chrpos', paste0(i,'mafa'), paste0(i,'nb'), paste0(i,'depth'), paste0('pob_',i))
    unique <- merge(unique, headings, by='chrpos', all=T)
  }
  pubs <- colnames(unique[ grep('pob', colnames(unique))])
  uniquesubs <- unique[,pubs]
  pob <- coalesce(!!!uniquesubs)
  unique = unique[,!(names(unique) %in% pubs)]
  unique$pob = pob
  mafas <- colnames(unique[ grep('mafa', colnames(unique))])
  nbs <- colnames(unique[ grep('nb', colnames(unique))])
  setnafill(unique[,mafas],fill=0)
  setnafill(unique[,nbs],fill=1)
  depths <- colnames(unique[ grep('depth', colnames(unique))])
  if(length(depths)>1){
    setnafill(unique[,depths],fill=30)
  }
  else{
    unique <- nafill(unique[,depths],fill=30)
  }
  if (0 %in% rowMeans(unique[,grep("mafa",colnames(unique))])){
    unique =  unique[-which(rowMeans(unique[,grep("mafa",colnames(unique))])==0),]
  }
  return(unique)
}

pubOrSub.prob.binom <- function(A, S, pu, sAGP, nt, nb, probpub, pa=1, prior="unknown") {      #using binomial test
  N = A + S
  nc = nt * sAGP + 2 * (1 - sAGP)
  f.abs <- 0.02
  f.pub <- (pu * (nt - nb))/nc
  
  f.pub = ifelse(nb == 0 & nt == 0, 0.5*pu,
                 ifelse(nb == 0 & prior == "pub", max((pu - sAGP)/nc, 0),
                        ifelse(nt >= 2, pu/nc, (pu * (nt - nb))/nc)))
  
  f.pub = min(1,f.pub)
  p.pub = pbinom(S, N, f.pub)       #p for pub (under pub vaf, the prob reaching S)
  p.abs = pbeta(f.abs, S+1, A+1)
  return(c(probpub, p.abs))
  # use pbeta / pbinom distributions to calculate probabilities of being public or absent
}

pubOrSub.Calc <- function(mutVector, regions, samples, maxPa = vector(), maxsAGP = vector()) {
  aboveContri = length(which(as.numeric(mutVector$ccf) +2.58*as.numeric(mutVector$ccfsd) >= 1))
  if (nrow(mutVector)<length(samples)){
    CCFaboveOne = F
    VAFaboveQua = F
  }
  #see if VAF is above vaf cutoff or CCF + 2.58 sd is above cutoff
  else{
    CCFaboveOne = all(as.numeric(mutVector$ccf) +   2.58*as.numeric(mutVector$ccfsd) >= 1)
    VAFaboveQua = all(as.numeric(mutVector$mafa) >= 0.39)
  }
  cpop = sapply(regions, function(x, regions, samples, mutVector, aboveContri, maxPa, maxsAGP) {
    index = which(regions==x)
    ccf = mutVector$ccf[mutVector$sn==x]
    ccfsd = mutVector$ccfsd[mutVector$sn==x]
    mafai = mutVector$mafa[mutVector$sn==x]
    refci = mutVector$ref[mutVector$sn==x]
    altci = mutVector$alt[mutVector$sn==x]
    pui = mutVector$pu[mutVector$sn==x]
    pai = mutVector$pa[mutVector$sn==x]
    sAGPi = mutVector$sAGP[mutVector$sn==x]
    nti = mutVector$nt[mutVector$sn==x]
    nbi = mutVector$nb[mutVector$sn==x]
    depthi = mutVector$TumorDepth[mutVector$sn==x]
    cPa=as.numeric(pai)
    probpubi = mutVector$probpub[mutVector$sn==x]
    prior="unknown"
    if (cPa == maxPa[index] || cPa == 0){
      cPa == 1
    }
    if (aboveContri >= 2 | (aboveContri >= 1 & length(samples) == 2)) {                    #looks like public, change prior
      prior = "pub"
    }
    cpop.probs = pubOrSub.prob.binom(as.integer(refci), as.integer(altci), as.numeric(maxsAGP[index]), as.numeric(sAGPi), as.integer(nti), as.integer(nbi), probpubi, pa=cPa, prior=prior)
    cpop.probs
  }, regions=regions, samples=samples, mutVector=mutVector, aboveContri=aboveContri, maxPa=maxPa, maxsAGP=maxsAGP)
  cpp = prod(cpop[1,])
  cpa = prod(cpop[2,])
  #calculate cumulative probabilities of being public / absent
  return(c(cpp, cpa, as.numeric(CCFaboveOne | VAFaboveQua)))
}

pubOrSub <- function(sampAB, minAF=0.05) {
  #classify mutations as public (clonal) or private (subclonal)
  maxPa = vector()
  maxsAGP = vector()
  mutations = as.character(unique(sampAB$chrpos))
  samples = as.numeric(unique(sampAB$sn))
  minDepTotal=length(samples)*5
  output = data.frame()
  for (i in 1:length(samples)) {
    sn = samples[i]
    pai = sampAB$pa[which(sampAB$sn==sn)]
    sAGPi = sampAB$sAGP[which(sampAB$sn==sn)]
    maxPa = c(maxPa, max(pai))
    maxsAGP = c(maxsAGP, max(sAGPi))
  }
  for (mut in mutations){
    mutVector = sampAB[sampAB$chrpos==mut,]
    regions = as.numeric(unique(mutVector$sn))
    cppres = pubOrSub.Calc(mutVector, regions, samples, maxPa, maxsAGP)
    cpstype = "unknown"
    totalDepth = sum(as.numeric(mutVector$TOTAL_DEPTH))
    #
    if (totalDepth < minDepTotal) {
      cpstype = "unknown"
    } else if (cppres[1]> 0.05 | cppres[3] == 1) {
      cpstype = "public"
    } else if (cppres[2] >= 0.05) { 
      #accept absent if binomial probability of seeing number of reads in the case of an absent mut > 0
      cpstype = "absent"
    } else {                                           #subclonal
      foundSamples = ""
      for (i in 1:length(regions)) {                 #get subclonal info for each sample
        sn = regions[i]
        mafai = mutVector$mafa[mutVector$sn==sn]
        nbi = mutVector$nb[mutVector$sn==sn]
        depthi = mutVector$TumorDepth[mutVector$sn==sn]
        if (as.numeric(mafai) > minAF) {
          foundSamples = paste(foundSamples, sn, ",", sep="")
        }
      }
      if (foundSamples != ""){
        foundSamples = gsub(",$","", foundSamples)
        cpstype = paste("private", foundSamples, sep="=")
      }
    }
    output = rbind(output, data.frame(chrpos=mut, sn=regions, cpsType=cpstype))
  }
  sampAB=merge(sampAB,output,by=c('chrpos','sn'),all=T)
  return (sampAB)
}

computeSD <- function(N, S, f, cc=seq(0.02, 1, by = 0.01)) {
  M1list <- c()
  M2list <- c()
  MLElist <- c()
  if (length(f)!=length(cc)){
    print('bad sd call')
    print(length(f))
    print(length(cc))
  }
  for (ii in 1:length(N)) {
    PF <- sum(dbinom(S[ii], N[ii], f), na.rm = TRUE)
    M1 <- sum(dbinom(S[ii], N[ii], f) * cc, na.rm = TRUE)/PF
    M2 <- sum(dbinom(S[ii], N[ii], f) * cc^2, na.rm = TRUE)/PF
    M1list <- c(M1list, M1)
    M2list <- c(M2list, M2)
    MLElist <- c(MLElist, cc[which.max(dbinom(S[ii], N[ii], f))])
  }
  SD = sqrt(M2list - M1list^2)
  if (is.nan(SD)){
    SD = .0002
  }
  return(list(M1 = MLElist, SD = sqrt(M2list - M1list^2)))
  # find standard deviation for ccf
}

partialRound <- function(x) {
  r = x
  for (i in 1:length(x)) {
    r[i] = round(x[i])
    if (abs(r[i]-x[i]) > 0.3 & x[i] > 1) {
      r[i] = round(x[i],1)
    }
    if (r[i] < 0) {
      r[i] = 0
    }
  }
  return(r)
}

computeCCF <- function(f, A, S, pu, pa, sAGP, nt, nb, prior="unknown", overadj=1.6, sigTh=0.90, sn) {
  na = nt-nb
  ccf = 0
  ccf2 = 0
  sd = 0
  cc <- seq(0.02, 1, by = 0.01)
  evoType = "A1/A2/B/C"
  N = A + S
  nc = nt * pa + 2 * (1 - pa)
  nc2 = nt * sAGP + 2 * (1 - sAGP)
  toohigh = F
  adjust = F
  otherwisebad = F
  mixA = F
  mixB = F
  if (nb == 1 & nt == 2) {   #normal diploid
    ccf = 2*(f/pu)
    ff = pu*cc/2
    Ms = computeSD(N, S, ff)
    ccf2 <- Ms$M1
    sd <- Ms$SD
    probpub = pbinom(S, N, max(ff)) 
    # no CNA
  }
  else if (nt == 1) {        #het deletion
    ff.C <- pu*cc/nc
    Ms.C <- computeSD(N, S, ff.C)
    ccf2 <- Ms.C$M1
    sd <- Ms.C$SD
    probpub = pbinom(S, N, max(ff.C)) 
    fh.ea = na*sAGP/nc2+(1-pa)*pu/nc2 
    fl.ea <- (sAGP * (nt - nb))/nc2
    fh.t <- sAGP/nc2
    fh.e <- (1 - pa)*pu/nc2
    pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1) #beta binom of early
    pLate <- pbeta(fh.t, S+1, A+1) #beta binomial probability of late mutation
    pEuploid <- pbeta(fh.e, S+1, A+1) #beta binomial probability of euploid mutation
    Ptot <- pEarly.a + pLate + pEuploid
    if (Ptot == 0){
      Ptot = 10^-323.6
      pLate = 10^-323.6
    }    # prevent division by zero errors
    cp.A <- pEarly.a/Ptot
    cp.CD <- 1 - cp.A
    cp.C <- pLate/Ptot
    cp.D <- pEuploid/Ptot
    cp.AC <- 1 - cp.D
    cp.AD <- 1 - cp.C
    if (cp.A >= sigTh | cp.A > cp.CD) {
      evoType <- "A1"
      ccf = (f*nc2-na*sAGP)/pu+pa
      ff.A = na*sAGP/nc2+(cc[cc>=(pa)]-pa)*pu/nc2
      Ms.A <- computeSD(N, S, ff.A,cc=cc[cc>=pa])
      ccf2 <- Ms.A$M1
      sd <- Ms.A$SD
      probpub = pbinom(S, N, max(ff.A)) 
      
    } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh) {
      evoType <- "B/C"
      ccf = (f/pu)*nc2 
      ff.BC = pu*cc/nc2
      Ms.BC = computeSD(N, S, ff.BC)
      probpub = pbinom(S, N, max(ff.BC)) 
      ccf2 <- Ms.BC$M1
      sd <- Ms.BC$SD
    } else if (cp.C >= sigTh) {
      evoType <- "B"
      ccf = (f/pu)*nc2 
      ff.BC = pu*cc/nc2
      Ms.BC = computeSD(N, S, ff.BC)
      ccf2 <- Ms.BC$M1
      sd <- Ms.BC$SD
      probpub = pbinom(S, N, max(ff.BC)) 
    } else if (cp.D >= sigTh) {
      evoType <- "C"
      ccf = (f/pu)*nc2 
      ff.BC = pu*cc/nc2
      Ms.BC = computeSD(N, S, ff.BC)
      ccf2 <- Ms.BC$M1
      sd <- Ms.BC$SD
      probpub = pbinom(S, N, max(ff.BC)) 
    } 
    else{
      evoType <- "B/C"
      ccf = (f/pu)*nc2 
      ff.BC = pu*cc/nc2
      Ms.BC = computeSD(N, S, ff.BC)
      ccf2 <- Ms.BC$M1
      sd <- Ms.BC$SD
      probpub = pbinom(S, N, max(ff.BC)) 
      
    }
    
    # infer evolution type and fit to it
  } else if (nb == 0 & nt == 0) {        #homozygous deletiohead
    evoType <- "C"
    # if you see mutations it has to be euploid
    pu_corr = (1-pa*pu)/(1-sAGP)
    
    if (pa!=1){
      ccf = f/pu*nc2
      ff.C = cc[cc<=(1-pa)]*pu/nc2
      Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=(1-pa)])
      ccf2 <- (Ms.C$M1)                  #dbinom
      sd <- Ms.C$SD
      probpub = pbinom(S, N, max(ff.C)) 
    }
    else{
      ccf = 2*(f/pu)
      ff = pu*cc/2
      Ms = computeSD(N, S, ff)
      ccf2 <- Ms$M1
      sd <- Ms$SD
      evoType = "A1/A2/B/C"
      probpub = pbinom(S, N, max(ff)) 
      
    }
  } else if (nb == 0 | nt == 2 * nb) {   #NLOH or other balanced CNAs
    fh.ea = na*sAGP/nc2+(1-pa)*pu/nc2
    fl.ea <- (sAGP * (nt - nb))/nc2
    fh.t <- sAGP/nc2
    fh.e <- (1 - pa)*pu/nc2
    pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
    pLate <- pbeta(fh.t, S+1, A+1)
    pEuploid <- pbeta(fh.e, S+1, A+1)
    Ptot <- pEarly.a + pLate + pEuploid
    cpEarly.a <- pEarly.a/Ptot
    cpLate.eup <- 1 - cpEarly.a
    cpLate <- pLate/Ptot
    cpEup <- pEuploid/Ptot
    if (Ptot > 0) {
      if (cpEarly.a >= sigTh){
        evoType <- "A1"
      } else if (cpLate.eup >= sigTh){
        evoType <- "B/C"
      } else if (cpLate >= sigTh){
        evoType <- "B"
      } else if (cpEup >= sigTh){
        evoType <- "C"
      }
      else if (cpEarly.a>.5){
        evoType <- "A1"
      }
      else {
        evoType <- "B/C"
      }
    }
    allprobs = c(pEarly.a, pLate, pEuploid)
    names(allprobs) = c("pEarly.a", "pLate", "pEuploid")
    maxType = names(allprobs[match(max(allprobs),allprobs)])
    if (evoType == "A1" & prior != "late") {
      ccf = (f*nc2-na*sAGP)/pu+pa
      ff.A = na*sAGP/nc2+((cc[cc>=(pa)]-pa)*pu/nc2)
      probpub = pbinom(S, N, max(ff.A)) 
      
      Ms.A <- computeSD(N, S, ff.A,cc=cc[cc>=pa])               #dbinom
      ccf2 <- Ms.A$M1  
      sd <- Ms.A$SD
      if (pa > ccf){ #assignment operates on the assumption that
        ccf = nc2 * f / (pu * na) 
        ff = cc * na*pu/nc2
        Ms.A <- computeSD(N, S, ff,cc=cc)
        probpub = pbinom(S, N, min(max(ff),1)) #calculate probability of finding mutation
        print (c(max(ff), N, S, probpub))
        ccf2 <- Ms.A$M1                             #dbinom
        sd <- Ms.A$SD
        evoType = "A1.5"
      }
      if ((abs(ccf-ccf2)>.02 & ccf<1)){ #correct any other problematic overflows (should not occur)
        ccf = (f*pu*(nt*pa+(1-pa)*2+2*(1-pu)))/(pu*(na*pa+1-pa))
        ff.A2 = (cc*(pu*(na*pa+1-pa)))/(pu*(nt*pa+(1-pa)*2+2*(1-pu)))
        probpub = pbinom(S, N, max(ff.A2)) 
        evoType = "A1.5"
        Ms.A <- computeSD(N, S, ff.A2, cc)               #dbinom
        ccf2 <- Ms.A$M1 
      }
    }
    else {
      ccf = (f/pu)*nc2
      ff.C <- pu*cc/nc2                        #dbinom
      Ms.C <- computeSD(N, S, ff.C)           #dbinom
      ccf2 <- Ms.C$M1                         #dbinom
      sd <- Ms.C$SD
      probpub = pbinom(S, N, max(ff.C)) 
      
    }
  } else if (nb >= 1 & nt > 2) {
    fh.ea = na*sAGP/nc2+(1-pa)*pu/nc2
    fl.ea <- na*sAGP/nc2
    fh.eb = nb*sAGP/nc2+(1-pa)*pu/nc2
    fl.eb <- nb*sAGP/nc2
    fh.t <- sAGP/nc2
    fh.e <- pu*(1-pa)/nc2
    pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
    pEarly.b <- pbeta(fh.eb, S+1, A+1) - pbeta(fl.eb, S+1, A+1)
    pLate <- pbeta(fh.t,S+1, A+1)
    pEuploid <- pbeta(fh.e, S+1, A+1)
    Ptot <- pEarly.a + pEarly.b + pLate + pEuploid
    if (Ptot == 0){
      Ptot = 10^-323.6
      pLate = 10^-323.6
    }
    cp.A <- pEarly.a/Ptot
    cp.B <- pEarly.b/Ptot
    cp.C <- pLate/Ptot
    cp.D <- pEuploid/Ptot
    cp.AB <- 1 - cp.C - cp.D
    cp.AC <- 1 - cp.B - cp.D
    cp.AD <- 1 - cp.B - cp.D
    cp.BC <- 1 - cp.A - cp.D
    cp.BD <- 1 - cp.A - cp.C
    cp.CD <- 1 - cp.A - cp.B
    cp.ABC <- 1 - cp.D
    cp.ABD <- 1 - cp.C
    cp.ACD <- 1 - cp.B
    cp.BCD <- 1 - cp.A
    if (Ptot > 0) {
      if (cp.A >= sigTh) {                   # estimate type of tumor evolution (notation)
        evoType = "A1"
      } else if (cp.B >= sigTh){
        evoType <- "A2"
      } else if (cp.C >= sigTh){
        evoType <- "B"
      } else if (cp.D >= sigTh){
        evoType <- "C"
      } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh){
        evoType <- "B/C"
      } else if (cp.AB >= sigTh & cp.A < sigTh & cp.B < sigTh){
        evoType <- "A1/A2"
      } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh){
        evoType <- "A1/B"
      } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh){
        evoType <- "A1/C"
      } else if (cp.BC >= sigTh & cp.B < sigTh & cp.C < sigTh){
        evoType <- "A2/B"
      } else if (cp.BD >= sigTh & cp.B < sigTh & cp.D < sigTh){
        evoType <- "A2/C"
      } else if (cp.BCD >= sigTh & cp.BC < sigTh & cp.BD < sigTh & cp.CD < sigTh & cp.B < sigTh & cp.C < sigTh & cp.D < sigTh){
        evoType <- "A2/B/C"
      } else if (cp.ABC >= sigTh & cp.BC < sigTh & cp.AB < sigTh & cp.AC < sigTh & cp.B < sigTh & cp.C < sigTh & cp.A < sigTh){
        evoType <- "A1/A2/B"
      } else if (cp.ABD >= sigTh & cp.AB < sigTh & cp.AD < sigTh & cp.BD < sigTh & cp.B < sigTh & cp.D < sigTh & cp.A < sigTh){
        evoType <- "A1/A2/C"
      } else if (cp.ACD >= sigTh & cp.AC < sigTh & cp.AD < sigTh & cp.CD < sigTh & cp.A < sigTh & cp.D < sigTh & cp.C < sigTh){
        evoType <- "A1/B/C"
      }
    }
    allprobs = c(pEarly.a, pEarly.b, pLate, pEuploid)
    names(allprobs) = c("pEarly.a", "pEarly.b", "pLate", "pEuploid")
    maxType = names(allprobs[match(max(allprobs),allprobs)])
    if (max(cp.A,cp.B,cp.CD)==cp.A & prior != "late") {
      ccf = (f*nc2-na*sAGP)/pu+pa
      ff.A = na*sAGP/nc2+((cc[cc>=pa]-pa)*pu/nc2)

      Ms.A <- computeSD(N, S, ff.A, cc=cc[cc>=pa])               #dbinom
      ccf2 <- Ms.A$M1
      sd <- Ms.A$M1
      probpub = pbinom(S, N, max(ff.A)) 
      
      if (pa > ccf){
        # check the alternate condition of the maximum
        ccf = nc2 * f / (pu * na)
        ff = cc * na*pu/nc2
        Ms.A <- computeSD(N, S, ff,cc=cc)
        ccf2 <- Ms.A$M1                             #dbinom
        sd <- Ms.A$SD
        evoType = "A1.5"
        probpub = pbinom(S, N, max(ff)) 
        
      }
      
      if ((abs(ccf-ccf2)>.02 & ccf<1)){
        ccf = mean(c(ccf,ccf2))
        ccf2 = round(ccf,2)
        ccf = (f*pu*(nt*pa+(1-pa)*2+2*(1-pu)))/(pu*(na*pa+1-pa))
        ff.A2 = (cc*(pu*(na*pa+1-pa)))/(pu*(nt*pa+(1-pa)*2+2*(1-pu)))
        evoType = "A1.5"
        Ms.A2 <- computeSD(N, S, ff.A2, cc)
        probpub = pbinom(S, N, max(ff.A2)) 
        #dbinom
        ccf2 <- Ms.A2$M1
        sd <- Ms.A2$SD
        print (c(evoType, nt, nb))
      }
      #dbinom
      sd <- Ms.A$SD
    } else if (max(cp.A,cp.B,cp.CD)==cp.B & prior != "late") {
      
      ccf = (f*nc2-nb*sAGP)/pu+pa
      #early A2
      ff.B <- nb*sAGP/nc2+(cc[cc>=(pa)]-pa)*pu/nc2  #dbinom
      Ms.B <- computeSD(N, S, ff.B, cc=cc[cc>=pa]) 
      probpub = pbinom(S, N, max(ff.B)) 
      ccf2 <- Ms.B$M1 
      sd <- Ms.B$SD
      if (pa > ccf){
        ccf = nc2 * f / (pu * nb)
        ff = cc * nb*pu/nc2
        Ms.A <- computeSD(N, S, ff,cc=cc)
        ccf2 <- Ms.A$M1                             #dbinom
        sd <- Ms.A$SD
        probpub = pbinom(S, N, max(ff)) 
        
        evoType = "B1.5"
      }
      
      if (abs(ccf-ccf2)>.02 & ccf<1){

        ccf = (f*pu*(nt*pa+(1-pa)*2+2*(1-pu)))/(pu*(nb*pa+1-pa))
        ff.B2 = (cc*(pu*(nb*pa+1-pa)))/(pu*(nt*pa+(1-pa)*2+2*(1-pu)))
        evoType = "B1.5"
        probpub = pbinom(S, N, max(ff.B2)) 
        
        mixB = T
        
        Ms.A <- computeSD(N, S, ff.B2, cc) 
        #dbinom
        ccf2 = Ms.A$M1
        sd <- Ms.A$SD

      }
    } 
    else {
      ccf = (f/pu)*nc2
      cc <- seq(0.02, 1, by = 0.01)#other
      ff.C <- pu*cc/nc2                       #dbinom
      Ms.C <- computeSD(N, S, ff.C)   
      probpub = pbinom(S, N, max(ff.C)) 
      #dbinom
      ccf2 <- Ms.C$M1                        #dbinom
      sd <- Ms.C$SD
    }
  }
  if (f > .1 & ccf >= 1.5 & pa==1 & nt>1){
    total_probs = c()
    pLate = pbeta(fh.t, S+1, A+1)
    for (val in (2:na)){
      prob = val*sAGP/nc2+(1-pa)*pu/nc2
      peclonal = pbeta(prob, S+1, A+1) - pbeta(prob*.95, S+1, A+1)
      total_probs = append(total_probs, peclonal)
    }
    if (max(total_probs)>pLate){
      ncopies = which(total_probs==max(total_probs))+1
      ccf = (f*nc2-ncopies*sAGP)/pu+pa
      evoType = "Aclonal"
      ff = ncopies*sAGP/nc2+((cc-pa)*pu/nc2)
      probpub = pbinom(S, N, max(ff)) 
      
      Ms.A <- computeSD(N, S, ff,cc=cc)               #dbinom
      ccf2 <- Ms.A$M1                             #dbinom
      sd <- Ms.A$SD
      if (ccf > 4){
        print ('pu too low')
        print (pu)
      }
      else if (1 > ccf){
        ccf = nc2 * f / (sAGP * ncopies)
        ff = cc * ncopies*sAGP/nc2
        probpub = pbinom(S, N, max(ff)) 
        
        Ms.A <- computeSD(N, S, ff,cc=cc)
        ccf2 <- Ms.A$M1                             #dbinom
        sd <- Ms.A$SD
      }
    }
  }  else if ( f > 0.1 & ccf >= overadj & pu >.05) {
    otherwisebad = T
  } 
  if (f > .1 & ccf > 1.5 & nt>1 & pa!=0 & pa!=1){
    total_probs = c()
    fh.t <- sAGP/nc2
    pLate = pbeta(fh.t, S+1, A+1)
    for (val in (2:na)){
      prob = val*sAGP/nc2+(1-pa)*pu/nc2 
      prob2 = val*sAGP/nc2
      peclonal = pbeta(prob, S+1, A+1) - pbeta(prob2, S+1, A+1)
      total_probs = append(total_probs, peclonal)
    }
    if (max(total_probs)>pLate){
      ncopies = which(total_probs==max(total_probs))+1
      ccf = (f*nc2-ncopies*sAGP)/pu+pa
      evoType = "Aclonal"
      ff = ncopies*sAGP/nc2+((cc[cc>=pa]-pa)*pu/nc2)
      probpub = pbinom(S, N, max(ff)) 
      
      Ms.A <- computeSD(N, S, ff,cc=cc[cc>=pa])               #dbinom
      ccf2 <- Ms.A$M1                             #dbinom
      sd <- Ms.A$SD
      if (ccf > 4){
        print ('pu too low')
        print (pu)
      }
      else if (pa > ccf){
        ccf = nc2 * f  / (pu * ncopies)
        ff = cc * ncopies*pu/nc2
        probpub = pbinom(S, N, max(ff)) 
        
        Ms.A <- computeSD(N, S, ff,cc=cc)
        ccf2 <- Ms.A$M1                             #dbinom
        sd <- Ms.A$SD
      }
    }
  }
  if ( f > 0.7 & ccf2 < 0.05 & is.nan(sd)) {  #un-identified LOH
    ccf = f*2
    ccf2 = 1
    sd = 0.01
  }
  if (length(ccf2)==2){
    print("duplicates")
    print (ccf2)
  }
  if (!(is.na(ccf))){
    if (ccf >= 1.5){
      ccf = 1.5
    }
  }
  if (ccf > 1.8 & is.nan(sd)){
    ccf = 1.5
    ccf2= 1.5
    sd = .002
  }
  return(c(ccf, evoType, ccf2, sd, toohigh, adjust, otherwisebad, mixA, mixB, probpub))
}

adjustccfmulti <- function(sampAB, tumor, sn1, clustall, ploidall, t =.02, overadj = 1.6, correctColname=FALSE, sigTh=.9){
  #function which combines SNV data with CNA data, estimates evolutionary type and calculates ccf (cancer cell fraction)
  files = list.files(path2files)
  samples2 = files[grepl('.ichor.segfull.txt',files)]
  samples = samples2[grepl(tumor, samples2)]
  clustnum <- clust$region
  num <- as.numeric(str_split_fixed(samples, "[_]", 3)[,2])
  # finding the value for purity
  num = append(num, sn1)
  sampAB$pu = NA
  sampAB$pa = NA
  sampAB$nt = NA
  sampAB$nb = NA
  sampAB$sAGP = NA
  sampAB$seg = NA
  # begin extracting values from file
  sampAB$ref = sampAB$TOTAL_DEPTH-sampAB$TUMOR_DEPTH
  sampAB$alt = sampAB$TUMOR_DEPTH
  sampAB$maf = sampAB$TUMOR_DEPTH/sampAB$TOTAL_DEPTH
  # done extracting data from file
  to_remove = which(!(sampAB$sn %in% num))
  if (length(to_remove)>0){
    sampAB = sampAB[-c(to_remove),]
  }
  unique = unique(sampAB$sn)
  for (i in unique){
    print(paste0(tumor,"_",i,"_"))
    tag = grep(paste0(tumor,"_",i,"_"),samples)
    if (!(i %in% clustnum)){
      # catching error with no purity data properly provided
      sampAB = sampAB[(sampAB$sn!=i),]
      print(paste("no purity data for",i))
    }
    else if (length(tag)>0 | i %in% sn1){
      n = which(clustnum==i)
      subsAB = sampAB[which(sampAB$sn==i),]
      n = which(clustnum==i)
      subsAB = sampAB[which(sampAB$sn==i),]
        # picking the purity from file
        cnv.inputA = paste0(path2files,samples[tag])
        cnvA = read.delim(cnv.inputA)
        missing = cnvA[which(cnvA$Start.snp=="Fill_In"),]
        con1 = as.numeric(clust$norm[n]) 
        if (length(missing)>0){
          cnvA = cnvA[-which(cnvA$Start.snp=="Fill_In"),]
        }
        cnvA = cnvA[which(!is.na(cnvA$Cellular_Prevalence)),]                #skip NA
      cnvA$nt = cnvA$Corrected_Copy_Number
      cnvA$nb = cnvA$Corrected_MinorCN
      cnvSeqNames = cnvA$Chromosome
      if (!grepl("chr",cnvA$Chromosome[1])){
        cnvSeqNames = paste("chr",cnvA$Chromosome,sep="")
      }
      snvSeqNames = subsAB$CHR
      if (!grepl("chr",subsAB$CHR[1])){
        snvSeqNames = paste("chr",subsAB$CHR,sep="")
      }
      # merging the cnv / snv data
      cnvRangeA = GRanges(seqnames = cnvSeqNames, ranges = IRanges(cnvA$Start, end=cnvA$End), strand=rep('+',dim(cnvA)[1]))
      snvRange = GRanges(seqnames = snvSeqNames, ranges = IRanges(subsAB$POSITION, end=subsAB$POSITION), strand=rep('+',dim(subsAB)[1]))
      foA = findOverlaps(snvRange, cnvRangeA)
      message("info table building")
      queHits = queryHits(foA)
      subHits = subjectHits(foA)
      infos = t(sapply(1:dim(subsAB)[1], function(x, queHits, subHits) {
        if (x %in% queHits){
          queMIndex = match(x, queHits)
          subMIndex = subHits[queMIndex]
          c(cnvA$Cellular_Prevalence[subMIndex],     #pa
            cnvA$nt[subMIndex],                     #nt
            cnvA$nb[subMIndex],                     #nb
            subMIndex)                              #seg
        } else {
          c(0,2,1,0)
        }}, queHits = queHits, subHits = subHits))
      infos = data.frame(infos)
      message("info table built")
      subsAB$pa = as.numeric(infos[,1])
      subsAB$nt = as.numeric(infos[,2])
      subsAB$nb = as.numeric(infos[,3])
      subsAB$seg = as.numeric(infos[,4])
      # finish merging the two data points
      pu_uncorrected = 1 - (con1 + (1-max(as.numeric(subsAB$pa)))*(1-con1))
      subsAB$pu = pu_uncorrected
      subsAB$sAGP = subsAB$pa*subsAB$pu
      region = i
      pu1 = subsAB$pu
      oldSamp = sampAB[which(sampAB$sn!=i),]
      sampAB = rbind(oldSamp,subsAB)
    }
  }
  sampAB$chrpos <- paste(sampAB$CHR, sampAB$POSITION, sep='_')
  mutations <- unique(sampAB$chrpos)
  print(paste("number of mutations is:", length(mutations)))
  sampAB$foundSites=0
  if (min(sampAB$pu)<.1){
    sampAB = sampAB[which(sampAB$pu>=.1),]
  }
  sampAB= na.omit(sampAB)
  unique = unique(sampAB$sn)
  print ('unique')
  mut_df <- data.frame()
  mutations <- unique(sampAB$chrpos)
  vector_out = data.frame(name = paste0("diagnostic_metrics/", tumor,"_region", min(sampAB$sn),"-",max(sampAB$sn)), lowpurity = 0, adjustment = 0, unclear_overall = 0, evoAunclear = 0 , evoBunclear = 0)
  name = paste0("diagnostic_metrics/", tumor,"_region", min(sampAB$sn))
  individual_outputs = data.frame()
  #using combination of CNV and SNV data to estimate tumor evolution type
  for (mut in mutations){
    n = which(mutations == mut)
    if (n%%1000==0){
      print(paste(n, "mutations analyzed out of", length(mutations)))
    }
    foundSites = 0             # count how many sites found
    for (j in unique){
      presentregions = unique(sampAB$sn[which(sampAB$chrpos==mut)])
      tum1 <- sampAB[which(sampAB$chrpos==mut & sampAB$sn==j),]
      if (nrow(tum1>0)){
        maf1 <- as.numeric(tum1$maf)
        if (maf1 >t){
          foundSites=foundSites+1
        }
      }
      if (j == presentregions[length(presentregions)]){
        sampAB$foundSites[which(sampAB$chrpos==mut)] <- foundSites
      }
    }
    for (j in unique){
      tum1 <- sampAB[which(sampAB$chrpos==mut & sampAB$sn==j),]
      if (nrow(tum1)>0 & max(sampAB$maf[which(sampAB$chrpos==mut)]>=.025)){
        #only analyze regions which have a >.05 vaf in ANY of the 
        maf1 <- as.numeric(tum1$maf)
        pa1 <- as.numeric(tum1$pa)
        nt1 <- as.numeric(tum1$nt)
        nb1 <- as.numeric(tum1$nb)
        refc1 <- as.numeric(tum1$ref)
        altc1 <- as.numeric(tum1$alt)
        pu1 <- as.numeric(tum1$pu)
        sAGP <- as.numeric(tum1$sAGP)
        foundSites <- as.numeric(tum1$foundSites)
        #input in all the known values
        if (maf1 > t && foundSites/length(unique(sampAB$sn))>=0){
          #no threshold for a presumed "late" mutation based upon # of reads in region due to heterogeneity
          CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"unknown",overadj=1.6,sigTh=.9, j)
          #calculate ccf values
          ccf1 <- as.numeric(CCF1[3])
          ccfsd1 <- as.numeric(CCF1[4])
          mafa1 <- as.numeric(CCF1[1])/2
          evoType1 <- as.character(CCF1[2])
          prob_pub <- as.numeric(CCF1[10])
          #output diagnostic metrics
          if (CCF1[5]==T){
            vector_out$lowpurity = vector_out$lowpurity + 1
            ind_supplement = data.frame(name = name, chrpos = mut, region = j, error = "too low purity")
            individual_outputs = rbind(individual_outputs, ind_supplement)
          }
          if (CCF1[6]==T){
            vector_out$adjustment = vector_out$adjustment + 1
            ind_supplement = data.frame(name = name, chrpos = mut, region = j,error = "adjustment")
            individual_outputs = rbind(individual_outputs, ind_supplement)
          }
          if (CCF1[7]==T){
            vector_out$unclear_overall = vector_out$unclear_overall + 1
            ind_supplement = data.frame(name = name, chrpos = mut, region = j,error = "evaluation is uncertain")
            individual_outputs = rbind(individual_outputs, ind_supplement)
          }
          if (CCF1[8]==T){
            vector_out$evoAunclear = vector_out$evoAunclear + 1
            ind_supplement = data.frame(name = name, chrpos = mut, region = j,error = "unclear if it is truly early major")
            individual_outputs = rbind(individual_outputs, ind_supplement)
          }
          if (CCF1[9]==T){
            vector_out$evoBunclear = vector_out$evoBunclear + 1
            ind_supplement = data.frame(name = name, chrpos = mut, region = j, error = "unclear if it is truly early minor")
            individual_outputs = rbind(individual_outputs, ind_supplement)
            
          }
        }
        else{
          CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"late",overadj=1.6,sigTh=0.9, j)
          ccf1 = as.numeric(CCF1[3])
          ccfsd1 = as.numeric(CCF1[4])
          mafa1 = as.numeric(CCF1[1])/2
          evoType1 <- as.character(CCF1[2])
          prob_pub <- as.numeric(CCF1[10])
        }
        mut_df <- rbind(mut_df, data.frame(chrpos=mut, sn=j, mafa=mafa1, ccf=ccf1, ccfsd=ccfsd1, e=evoType1, probpub = prob_pub))
      }
    }
  }
  
  #print(head(mut_df))
  print('output metrics')
  overall = merge(sampAB, mut_df, by=c('chrpos','sn'))
  #create new file with ccfs
  return(overall)
}

combineCCF = function(path, region=c(), not_full=F,location=path2files){
  #function which combines data for snvs across all regions for a single patient, specifically designed to take in preRT / postRT input separately (two distinct passes)
  first = TRUE
  temp = unique(full_data$SAMPLE)
  num = str_split_fixed(temp, "[-]", 2)[,2]
  samples = unique(full_data$SAMPLE)
  if (not_full==T){
    print(num)
    num = num[which(num %in% paste0("T",region))]
    print(num)
  }
  for (i in num){
    if (first == T){
      right_sample = samples[match(paste0(path,"-",i),samples)]
      region1 = full_data[which(full_data$SAMPLE==right_sample),]
      region1 = region1[,c("TUMOR_DEPTH", "TOTAL_DEPTH", "POSITION", "CHR")]
      duplicates = duplicated(region1[,1:4])
      position = which(duplicates==T)
      if (length(position)>0){
        print(paste("duplicates found in region", i))
        region1 <- region1[-c(position),]
      }
      region1$sn = substr(i,2,nchar(i))
      first = F
    }
    else {
      right_sample = samples[grep(paste0(path,"-",i),samples)]
      region = full_data[which(full_data$SAMPLE==right_sample),] 
      region = region[,c("TUMOR_DEPTH", "TOTAL_DEPTH", "POSITION", "CHR")]
      duplicates = duplicated(region[,1:4])
      position = which(duplicates==T)
      if (length(position)>0){
        print(paste("duplicates found in region", i))
        region <- region[-c(position),]
      }
      region$sn = substr(i,2,nchar(i))
      region1 = rbind(region1, region)
    }
  }
  return (region1)
}

for (sample in unique(records$sample)){ #loops though all data and runs functions in tandem
timepoints = c(0,1)
region_tp0 = length(unique(records$region[which(records$sample==sample & records$timepoint==0)]))
region_tp1 = length(unique(records$region[which(records$sample==sample & records$timepoint==1)]))
min_regions = min(region_tp0, region_tp1)
region_tp10 = length(unique(records$region[which(records$sample==sample & records$timepoint==0)]))
region_tp11 = length(unique(records$region[which(records$sample==sample & records$timepoint==1)]))
min_total_regions = min(region_tp10, region_tp11)
for (tp in timepoints){
  print(tp)
  max_timepoint = max(records$timepoint[records$sample==sample])
  region_all=records$region[which(records$sample==sample & records$timepoint==tp)]
  region_all2=records$region[which(records$sample==sample & records$timepoint==tp)]
  category=records$category[which(records$sample==sample & records$timepoint==tp)]
    region=region_all
    sampAB = combineCCF(sample, region=region_all, not_full=T)
    sampAB2 = sampAB
    sampAB2$chrpos=paste0(sampAB$CHR,'_',sampAB$POSITION)
    for (mutation in unique(sampAB2$chrpos)){
      max1 = max(sampAB2$TUMOR_DEPTH[which(sampAB2$chrpos==mutation)]/sampAB2$TOTAL_DEPTH[which(sampAB2$chrpos==mutation)])
    }
    sampAB=sampAB2[which(sampAB2$sn %in% region),]
    output = adjustccfmulti(sampAB, sample, c(), "filler")
    if (min(output$pu)<.1){
      output = output[which(output$pu>=.1),]
    }
    if (max(output$pu)<.2){
      output = output[which(output$pu>=.2),]
    }
    if (length(unique(output$sn))>1){
      pubcalced = pubOrSub(output)
      compressed = smallerdf(pubcalced)
      regions = unique(pubcalced$sn)
      metrics = subclonalMut(compressed, regions)
      metrics$sample=sample
      metrics$timepoint=tp
      metrics$region=length(regions)
      metrics$category=unique(category)
      if (first == T){
        metrics_final=metrics
        first = F
      }  else {
        metrics_final = rbind(metrics_final, metrics)
        write.table(metrics_final, file=paste0(path2files,"metrics.tsv"), sep='\t', col.names=T, row.names=F)
      }
    }
}
}
