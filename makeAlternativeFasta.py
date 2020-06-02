import sys,glob

ChrSNPpos = {}

def makeFastaSet(finp):
    FastaDict = {}
    FastaOrder= []
    with open(finp,'r') as fp:
        for line in fp:
            if line.startswith('>'):
                title = line.strip()[1:]
                FastaDict.setdefault(title,[])
                FastaOrder.append(title)
            else:
                seq = line.strip()
                FastaDict[title]+=seq
    return FastaDict,FastaOrder

def SwitchSeq(FastaDict,ListOfPos):
    with open(ListOfPos,'r') as PosList:
        for line in PosList:
            if line.strip() == "":continue
            title,pos,ref,alt = line.strip().split('\t')
            if FastaDict[title][int(pos)-1] !=ref:
                print ('reference seq is not matched its coordinate!')
            FastaDict[title][int(pos)-1] = alt
    return FastaDict

def NmaskSeq(FastaDict,ListOfPos):
    informed_dict = {}
    with open(ListOfPos,'r') as PosList:
        for line in PosList:
            if line.strip() == "":continue
            title,pos,ref,alt = line.strip().split('\t')
            FastaDict[title][int(pos)-1] = "N"
            informed_dict.setdefault(title,{})
            informed_dict[title][int(pos)] = ref+"/"+alt 
    return FastaDict,informed_dict

def main():
    ref_file = '/nv/hp10/xwu302/data2/moredata/apis/AM_plosP/apis_beebase_4.5/genome/Amel_4.5_scaffolds.fa' #reference genome
    flist = glob.glob("/nv/hp10/xwu302/data2/apis_PSM/RNA/snpsplit/snps/*informed_snps.txt") # can handle multiple snp files at once
    for informed_file in flist:
        sample_index = informed_file.split('/')[-1].split('_')[0]
        shared_file = informed_file.replace("_informed_","_shared_")
        conserved_share_file = shared_file
        conserved_informative_file = informed_file
        FastaDictionary, FastaOrderedList = makeFastaSet(ref_file)
        FastaDictionary = SwitchSeq(FastaDictionary,conserved_share_file)
        FastaDictionary,informed_pos_dict = NmaskSeq(FastaDictionary,conserved_informative_file)
        final_fasta = 'Nmasked_Amel_'+sample_index+'.fa'
        fSNP_split = open("SNP_"+sample_index+"_snpsplit.txt",'w')
        snp_index = 0
        with open(final_fasta,'w') as foutFas:
            for header in FastaOrderedList:
                if header in informed_pos_dict:
                    ordered_pos_list = sorted(informed_pos_dict[header].keys())
                    for snpidx in ordered_pos_list:
                        snp_index +=1
                        fSNP_split.write(str(snp_index)+'\t'+header+'\t'+str(snpidx)+'\t1\t'+informed_pos_dict[header][snpidx]+'\n')
                foutFas.write('>'+header+'\n')
                cnt = 0
                for i in FastaDictionary[header]:
                    cnt +=1
                    foutFas.write(i)
                    if cnt % 60 == 0:
                        foutFas.write('\n')
                foutFas.write('\n')
        fSNP_split.close()        

main()









