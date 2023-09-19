import os,sys
##Usage: python vcf-variantQC.py in.vcf field cutoff > op.vcf
ivcf=sys.argv[1] ## input VCF
ifield=sys.argv[2] ## FORMAT Field or INFO Field
ival=sys.argv[3] ## cutoff or filter value

op=ivcf.replace(".vcf","")
#print(op)

def formatfield(vformat,field): ##Find the position of the format field in the format column
	tmp=-999
	vf=vformat.split(":")
	for n,v in enumerate(vf):
		if v==field:
			tmp=int(n+1)
			break
		else:
			continue
	return tmp

def singlefilter(vardict,pos,val,sams): ## format filters which has single value for Alt
	qcvar={}
	badsam=[]
	for nk,k in enumerate(sams):
		if nk<=8:
			qcvar[k]=vardict[k]
		elif nk>8:
			fvararr=str(vardict[k]).split(":")
			if fvararr[0]=="./.":
				qcvar[k]=str(":".join(fvararr))
			else:
				try:
					try:
						if int(fvararr[pos-1])<int(val):
							fvararr[0]="0/0"
							qcvar[k]=str(":".join(fvararr))
							badsam.append(k)
						else:
							qcvar[k]=str(":".join(fvararr))
					except ValueError:
						qcvar[k]=str(":".join(fvararr))
				except IndexError:
					qcvar[k]=str(":".join(fvararr))
	bsams=";".join(badsam)
	return (qcvar,bsams)

def multifilter(vardict,pos,val,sams,field):
	qcvar={}
	badsam=[]
	for nk,k in enumerate(sams):
		if nk<=8:
			qcvar[k]=vardict[k]
		elif nk>8:
			fvararr=str(vardict[k]).split(":")
			if fvararr[0]=="./.":
				qcvar[k]=str(":".join(fvararr))
				#print(fvararr)
			else:
				refgeno=str(fvararr[0].split("/")[0])
				altgeno=str(fvararr[0].split("/")[1])
				if str(fvararr[pos-1]).split(",")[0]=="." or str(fvararr[pos-1]).split(",")[1]==".":
					fvararr[0]="0/0"
					qcvar[k]=str(":".join(fvararr))
				else:
					refad=int(str(fvararr[pos-1]).split(",")[0])
					altad=int(str(fvararr[pos-1]).split(",")[1])
				#print(refad,altad)
					try:
						ab=float(altad)/float(refad+altad)
					except ZeroDivisionError:
						ab=0
					#print(refgeno,altgeno,refad,altad,ab,val,float(float(val)/100.0))
					if field=="AB":
						if refgeno==altgeno: ##Homozygous variantsi
							#print(refgeno,altgeno,refad,altad,ab,val,float(float(val)/100.0))
							if refgeno=="0" and altgeno=="0":
								qcvar[k]=str(":".join(fvararr))
							else:
								#print(refgeno,altgeno,refad,altad,ab,val,float(float(val)/100.0))
								if float(ab) < float(float(val)/100.0):
									fvararr[0]="0/0"
                        		                        	qcvar[k]=str(":".join(fvararr))
									badsam.append(k)
								else:
									qcvar[k]=str(":".join(fvararr))
						else:
							qcvar[k]=str(":".join(fvararr))
					elif field=="ABHet":
						if refgeno!=altgeno: #Heterozygous
							#print(refgeno,altgeno,refad,altad,ab,val,float(float(val)/100.0))
							if float(ab) > float(float(1)-float(float(val)/100.0)): ## changing to alt if AB < 100-20 the genotype is changed to the alternative variant
								fvararr[0]=str(altgeno)+"/"+str(altgeno)
								qcvar[k]=str(":".join(fvararr))
								badsam.append(k)
							elif float(ab) < float(float(val)/100.0): ## changing to alt if AB < 20 the genotype is changed to reference variant
								fvararr[0]=str(refgeno)+"/"+str(refgeno)
								qcvar[k]=str(":".join(fvararr))
								badsam.append(k)
							else:
								qcvar[k]=str(":".join(fvararr))
						else:
							qcvar[k]=str(":".join(fvararr))
	bsams=";".join(badsam)
	return (qcvar,bsams)
	

def acanaf_snpindel(vardict,sams): ## Recalculate AC, AN, AF for Snps and Indels
	ac=0
	an=0
	af=0
	#print(vardict['INFO'].split(";"))
	varinfodict={}
	for nk,k in enumerate(sams):
                if nk>8:
			#print(nk,k,vardict[k])
			if "0/0" in vardict[k].split(":")[0] or "0|0" in vardict[k].split(":")[0]:
				ac=ac+0
				an=an+2
			elif "0/1" in vardict[k].split(":")[0] or "0|1" in vardict[k].split(":")[0]:
				ac=ac+1
				an=an+2
			elif "1/1" in vardict[k].split(":")[0] or "1|1" in vardict[k].split(":")[0]:
				ac=ac+2
				an=an+2
			elif "./." in vardict[k].split(":")[0] or ".|." in vardict[k].split(":")[0]:
				ac=ac+0
				an=an+0
			elif "./1" in vardict[k].split(":")[0] or ".|1" in vardict[k].split(":")[0]:
				#print(vardict[k])
				ac=ac+1
				an=an+1
			elif "./0" in vardict[k].split(":")[0] or ".|0" in vardict[k].split(":")[0]:
                                #print(vardict[k])
				ac=ac+0
                                an=an+1
			elif "1/." in vardict[k].split(":")[0] or "1|." in vardict[k].split(":")[0]:
                                #print(vardict[k])
				ac=ac+1
                                an=an+1
			elif "0/." in vardict[k].split(":")[0] or "0|." in vardict[k].split(":")[0]:
                                #print(vardict[k])
				ac=ac+0
                                an=an+1
			else:
				#print(vardict[k].split(":")[0])
				ac=ac+0
				an=an+0
	if ac==an==0:
		af=0.0
	else:
		af=float(ac)/float(an)
	varinfo=vardict['INFO'].split(";")
	tt=[]
	for varinf in varinfo:
		if "=" in varinf:
			varinfodict[varinf.split("=")[0]]=varinf.split("=")[1]
		else:
			tt.append(varinf)
	#print(varinfodict['AC'],varinfodict['AN'],varinfodict['AF'])
	varinfodict['AC']=ac
	varinfodict['AN']=an
	varinfodict['AF']=af
	#print(varinfodict['AC'],varinfodict['AN'],varinfodict['AF'])
	newinfoarr=[]
	for kk in varinfodict.keys():
		newinfoarr.append(str(kk)+"="+str(varinfodict[kk]))
	if tt==[]:
		newvarinfo=";".join(newinfoarr)
	else:
		newvarinfo=";".join(newinfoarr)+";"+";".join(tt)
	vardict['INFO']=newvarinfo
	return (vardict)

def acanaf_mav(vardict,sams): ## Recalculate AC, AN, AF for MAVs
	#print(vardict)
	ac=[]
	an=0
	af=[]
	varinfodict={}
	altlen=len(vardict['ALT'].split(","))
	#print(vardict['ALT'])
	for altl in range(0,altlen+1,1):
		ac.append(int(0))
	for nk,k in enumerate(sams):
		if nk>8:
			for altl in range(0,altlen+1,1):
				if str(altl)==str(str(vardict[k]).split(":")[0].split("/")[0]):
					#print(str(altl),str(str(vardict[k]).split(":")[0].split("/")))
					ac[altl]=int(ac[altl])+1
					#an=an+1
				elif str(altl)==str(str(vardict[k]).split(":")[0].split("/")[1]):
					#print(str(altl),str(str(vardict[k]).split(":")[0].split("/")))
					ac[altl]=int(ac[altl])+1
					#an=an+1
			if "./." in str(vardict[k]).split(":")[0] or ".|." in  str(vardict[k]).split(":")[0]:
				#print(str(vardict[k]).split(":")[0])
				g=0
			else:
				an=an+2
	varinfo=vardict['INFO'].split(";")
        tt=[]
        for varinf in varinfo:
                if "=" in varinf:
                        varinfodict[varinf.split("=")[0]]=varinf.split("=")[1]
                else:
                        tt.append(varinf)
	for nx,x in enumerate(ac):
		af.append(str(float(x)/float(an)))
		ac[nx]=str(x)
        varinfodict['AC']=",".join(ac[1:])
        varinfodict['AN']=an
        varinfodict['AF']=",".join(af[1:])
        newinfoarr=[]
	for kk in varinfodict.keys():
                newinfoarr.append(str(kk)+"="+str(varinfodict[kk]))
	if tt==[]:
		newvarinfo=";".join(newinfoarr)
	else:
		newvarinfo=";".join(newinfoarr)+";"+";".join(tt)
        vardict['INFO']=newvarinfo
        return (vardict)


def dictvartoline(vardict,sams):
	varlinearr=[]
	for nk,k in enumerate(sams):
		varlinearr.append(vardict[k])
	varline="\t".join(varlinearr)
	return (varline)
		

samples=[]
ovcf=open(os.path.expanduser(str(op)+"_"+str(ifield)+str(ival)+".vcf"),'w')
olog=open(os.path.expanduser(str(op)+"_"+str(ifield)+str(ival)+".log"),'w')
for line in open(os.path.expanduser(str(ivcf)),'r'):
	#print(line)
	line=line.strip()
	if "##" in line:
		a=0
		ovcf.write(line+"\n")
		olog.write(line+"\n")
	elif "#CHROM" in line:
		samples=line.split("\t")
		ovcf.write(line+"\n")
		olog.write(line+"\n")
	else:
		arr=line.split("\t")
		#print(arr)
		'''variable declaration'''
		chrom=str(arr[0])
		pos=int(arr[1])
		varid=str(arr[2])
		ref=str(arr[3])
		altarr=str(arr[4]).split(",")
		varqual=float(arr[5])
		varfilter=str(arr[6])
		varinfo=str(arr[7])
		varformat=str(arr[8])
		'''convert variant line into dictionary '''
		vararr={}
		qcvararr={}
		badqcsams=""
		newqcvar={}
		tmp=''
		for ng,geno in enumerate(arr):
			vararr[samples[ng]]=geno
		'''single field variant filter '''
		if ifield=="DP" or ifield=="GQ":
			formatpos=formatfield(varformat,ifield) # '''get the pos of filter field'''
			if formatpos==-999:
				qcvararr=vararr
				badqcsams=""
			else:
				qcres=singlefilter(vararr,formatpos,ival,samples)
				qcvararr=qcres[0]
				badqcsams=qcres[1]
		elif ifield=="AB" or ifield=="ABHet":
			formatpos=formatfield(varformat,"AD") # '''get the pos of filter field'''
			if formatpos==-999:
                                qcvararr=vararr
                                badqcsams=""
                        else:
	                        qcres=multifilter(vararr,formatpos,ival,samples,ifield)
                        	qcvararr=qcres[0]
                        	badqcsams=qcres[1]
		''' recalculate AC, AN, AF'''
		if len(altarr)>1:
			g=0
			newqcvararr=acanaf_mav(qcvararr,samples)
			ovcf.write(str(dictvartoline(newqcvararr,samples))+"\n")
			olog.write(str(chrom)+"\t"+str(pos)+"\t"+str(varid)+"\t"+str(ref)+"\t"+str(",".join(altarr))+"\t"+str(varqual)+"\t"+str(varfilter)+"\t"+str(varinfo)+"\t"+str(varformat)+"\t"+str(badqcsams)+"\n")
		else:
			g=0
			newqcvararr=acanaf_snpindel(qcvararr,samples)
			ovcf.write(str(dictvartoline(newqcvararr,samples))+"\n")
			olog.write(str(chrom)+"\t"+str(pos)+"\t"+str(varid)+"\t"+str(ref)+"\t"+str(",".join(altarr))+"\t"+str(varqual)+"\t"+str(varfilter)+"\t"+str(varinfo)+"\t"+str(varformat)+"\t"+str(badqcsams)+"\n")
