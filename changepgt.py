import os,sys

for line in open(os.path.expanduser(str(sys.argv[1])),'r'):
	line=line.strip()
	if "##" in line:
		print(line)
	elif "#CHROM" in line:
		print(line)
	else:
		arr=line.split("\t")
		for i in range(9,len(arr),1):
			if "/" in arr[i].split(":")[0]:
				continue
			elif "|" in arr[i].split(":")[0]:
				#print(arr[i])
				gtarr=arr[i].split(":")
				tgt=gtarr[0].replace("|","/")
				gtarr[0]=tgt
				arr[i]=":".join(gtarr)
		print("\t".join(arr))
