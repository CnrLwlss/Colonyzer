# Get output files from Colonyzer2
outfiles=list.files(pattern="*\\.dat")
outfiles=outfiles[outfiles!="Output.dat"]
col2=data.frame()
for(out in outfiles){
	print(out)
	data=read.delim(out,header=FALSE,sep="\t",stringsAsFactors=FALSE)
	col2=rbind(col2,data)
}

colonyzercolumns=c("FILENAME","ROW","COLUMN","TOPLEFTX","TOPLEFTY","WHITEAREA","TRIMMED","THRESHOLD","INTENSITY","EDGEPIXELS","COLR","COLG","COLB","BKR","BKG","BKB","EDGELEN","XDIM","YDIM")
ordercols=sort(colonyzercolumns)
colnames(col2)=ordercols
col2$BARCODE=substr(col2$FILENAME,1,11)
col2$DATETIME=substr(col2$FILENAME,13,nchar(col2$FILENAME[1]))

# Get original Colonyzer equivalents
outfiles=list.files(path="Output_Data/",pattern="*\\.dat")
col1=data.frame()
for(out in outfiles){
	print(out)
	data=read.delim(file.path("Output_Data",out),header=FALSE,sep="\t",stringsAsFactors=FALSE)
	col1=rbind(col1,data)
}

colonyzercolumns=c("FILENAME","ROW","COLUMN","TOPLEFTX","TOPLEFTY","WHITEAREA","TRIMMED","THRESHOLD","INTENSITY","EDGEPIXELS","COLR","COLG","COLB","BKR","BKG","BKB","EDGELEN","XDIM","YDIM")
colnames(col1)=colonyzercolumns
col1$BARCODE=substr(col1$FILENAME,1,11)
col1$DATETIME=substr(col1$FILENAME,13,nchar(col1$FILENAME[1]))

col1=col1[col1$BARCODE%in%unique(col2$BARCODE),]
col1=col1[order(col1$ROW,col1$COLUMN,col1$DATETIME),]
col2=col2[order(col2$ROW,col2$COLUMN,col2$DATETIME),]

# Add a time-since-inoculation column
fmt="%Y-%m-%d_%H-%M-%S"
t1<-as.POSIXlt(as.character(col1$DATETIME),format=fmt)
t0=min(t1)
col1$TIME=as.numeric(difftime(t1,t0,units="days"))

t1<-as.POSIXlt(as.character(col2$DATETIME),format=fmt)
t0=min(t1)
col2$TIME=as.numeric(difftime(t1,t0,units="days"))

plot(col1$INTENSITY,col2$INTENSITY,xlab="Colonyzer Intensity (AU)",ylab="Colonyzer2 Fractional Intensity",cex.lab=1.45,pch=16,cex=0.1)
plot(col1$TRIMMED,col2$TRIMMED,xlab="Colonyzer Trimmed Intensity (AU)",ylab="Colonyzer2 Fractional Trimmed Intensity (AU)",cex.lab=1.45,pch=16,cex=0.1)
plot(col1$WHITEAREA,col2$WHITEAREA,xlab="Colonyzer Area (px)",ylab="Colonyzer2 Fractional Area",cex.lab=1.45,pch=16,cex=0.1)

col=10
row=10

dat1=col1[(col1$ROW==row)&(col1$COLUMN==col),]
dat2=col2[(col2$ROW==row)&(col2$COLUMN==col),]

op<-par(mfrow=c(2,3),oma=c(0,2,0,2))
clab=1.5
plot(dat1$TIME,dat1$INTENSITY,type="b",xlab="Time (d)",ylab="Intensity (AU)",cex.lab=clab)
plot(dat1$TIME,dat1$TRIMMED,type="b",xlab="Time (d)",ylab="Trimmed Intensity (AU)",main="Colonyzer",cex.lab=clab)
plot(dat1$TIME,dat1$WHITEAREA,type="b",xlab="Time (d)",ylab="Area (px)",cex.lab=clab)

plot(dat2$TIME,dat2$INTENSITY,type="b",xlab="Time (d)",ylab="Fractional Intensity",cex.lab=clab)
plot(dat2$TIME,dat2$TRIMMED,type="b",xlab="Time (d)",ylab="Fractional Trimmed Intensity",main="Colonyzer2",cex.lab=clab)
plot(dat2$TIME,dat2$WHITEAREA,type="b",xlab="Time (d)",ylab="Fractional Area",cex.lab=clab)
par(op)
